import numpy as np
import sys, obspy, os
import matplotlib.pyplot as plt
from obspy import read, read_inventory
from obspy.core.stream import Stream
import sys

'''
Test - Get inventory file
Code made by ICJG
'''
# Functions by seismolive
# collection of functions used in noise correlation processing

def normalize(tr, clip_factor=6, clip_weight=10, norm_win=None, norm_method="1bit"): 
    
    if norm_method == 'clipping':
        lim = clip_factor * np.std(tr.data)
        tr.data[tr.data > lim] = lim
        tr.data[tr.data < -lim] = -lim

    elif norm_method == "clipping_iter":
        lim = clip_factor * np.std(np.abs(tr.data))
        
        # as long as still values left above the waterlevel, clip_weight
        while tr.data[np.abs(tr.data) > lim] != []:
            tr.data[tr.data > lim] /= clip_weight
            tr.data[tr.data < -lim] /= clip_weight

    elif norm_method == 'ramn':
        lwin = tr.stats.sampling_rate * norm_win
        st = 0                                               # starting point
        N = lwin   
        N=int(N)                                          # ending point
        lwin=int(lwin)
        while N < tr.stats.npts:
            win = tr.data[st:N]
            w = np.mean(np.abs(win)) / (2. * lwin + 1)
            
            # weight center of window
            tr.data[st + lwin // 2] /= w

            # shift window
            st += 1
            N += 1

        # taper edges
        taper = get_window(tr.stats.npts)
        tr.data *= taper

    elif norm_method == "1bit":
        tr.data = np.sign(tr.data)
        tr.data = np.float32(tr.data)

    return tr


def get_window(N, alpha=0.2):

    window = np.ones(N)
    x = np.linspace(-1., 1., N)
    ind1 = (abs(x) > 1 - alpha) * (x < 0)
    ind2 = (abs(x) > 1 - alpha) * (x > 0)
    window[ind1] = 0.5 * (1 - np.cos(np.pi * (x[ind1] + 1) / alpha))
    window[ind2] = 0.5 * (1 - np.cos(np.pi * (x[ind2] - 1) / alpha))
    return window


def whiten(tr, freqmin, freqmax):
    
    nsamp = tr.stats.sampling_rate
    
    n = len(tr.data)
    if n == 1:
        return tr
    else: 
        frange = float(freqmax) - float(freqmin)
        nsmo = int(np.fix(min(0.01, 0.5 * (frange)) * float(n) / nsamp))
        f = np.arange(n) * nsamp / (n - 1.)
        JJ = ((f > float(freqmin)) & (f<float(freqmax))).nonzero()[0]
            
        # signal FFT
        FFTs = np.fft.fft(tr.data)
        FFTsW = np.zeros(n) + 1j * np.zeros(n)

        # Apodization to the left with cos^2 (to smooth the discontinuities)
        smo1 = (np.cos(np.linspace(np.pi / 2, np.pi, nsmo+1))**2)
        FFTsW[JJ[0]:JJ[0]+nsmo+1] = smo1 * np.exp(1j * np.angle(FFTs[JJ[0]:JJ[0]+nsmo+1]))

        # boxcar
        FFTsW[JJ[0]+nsmo+1:JJ[-1]-nsmo] = np.ones(len(JJ) - 2 * (nsmo+1))\
        * np.exp(1j * np.angle(FFTs[JJ[0]+nsmo+1:JJ[-1]-nsmo]))
        # Apodization to the right with cos^2 (to smooth the discontinuities)
        smo2 = (np.cos(np.linspace(0, np.pi/2, nsmo+1))**2)
        espo = np.exp(1j * np.angle(FFTs[JJ[-1]-nsmo:JJ[-1]+1]))
        FFTsW[JJ[-1]-nsmo:JJ[-1]+1] = smo2 * espo

        whitedata = 2. * np.fft.ifft(FFTsW).real
        
        tr.data = np.require(whitedata, dtype="float32")

        return tr


def correlateNoise(st, stations, corrwin):
    print ('correlating stations', (stations[0], stations[1]))

    # initialize sliding timewindow (length = corrwin) for correlation
    # start 1 corrwin after the start to account for different stream lengths
    timewin = st.select(station=stations[1])[0].stats.starttime + corrwin
    # loop over timewindows 
    # stop 1 corrwin before the end to account for different stream lengths
    while timewin < st.select(station=stations[0])[-1].stats.endtime - 2*corrwin:
        sig1 = st.select(station=stations[0]).slice(timewin, timewin+corrwin)
        sig1.split().merge(method=0, fill_value=0)
        sig2 = st.select(station=stations[1]).slice(timewin, timewin+corrwin)
        sig2.split().merge(method=0, fill_value=0)
        xcorr = np.correlate(sig1[0].data, sig2[0].data, 'same')
        try: 
            # build array with all correlations
            corr = np.vstack((corr, xcorr))
        except: 
            # if corr doesn't exist yet
            corr = xcorr
        # ncorr=corr  
        # print(ncorr)  
        # shift timewindow by one correlation window length
        timewin += corrwin

        # stack the correlations; normalize
    corr = np.nan_to_num(corr,0)
    stack = np.sum(corr, 0)
    stack = stack / float((np.abs(stack).max()))    
    print ("...done")
  
    return corr,stack


def plotStack(st, stack2, maxlag, sta1, sta2, figurename=None):

    # define the time vector for the correlation (length of corr = corrwin + 1)
    limit = (len(nstack) / 2.) * int(st[0].stats.delta * 1000.) / 1000
    #limit = (len(stack2) / 2.) * st[0].stats.delta
    timevec = np.arange(-limit, limit, st[0].stats.delta)
    fig = plt.figure(figsize=(18,4))
    fig = plt.plot(timevec, stack2, 'k', linewidth='0.5')
    #stations = list(set([_i.stats.station for _i in st]))
    fig = plt.title("Stacked correlation between %s and %s" % (sta1, sta2))
    fig = plt.xlim(-maxlag, maxlag)
    fig = plt.xlabel('time [s]')

    if figurename is not None:
        fig = plt.savefig(figurename, format="pdf")
    else:
        fig = plt.show()

        

if __name__ == '__main__':
    nst = Stream()
    jdays = range(152,160)     # Julian days from June-July two months didn't work
    # Reading files
    inv = read_inventory('/Users/home/juarezilma/PhD/Test/Python/ALL_inventory.xml')
    # stat1 = sys.argv[1] #GEONET
    # stat2 = sys.argv[2]	#COSA,DWARFS
    # print('Station1:',stat1)st
    # print('Station2:',stat2)
    for j in jdays:
        print('Reading Julian day:',j)
        nst += read('/Volumes/GeoPhysics_23/users-data/juarezilma/NoisePy/Raw_data/*' '*HHZ*'+ str(j) + '*' )
        nst += read('/Volumes/GeoPhysics_23/users-data/juarezilma/NoisePy/Raw_data/*' + str(j)  + '*HHZ*' )
    print('DONE READING!')
    print(len(nst))
    #nst = st.copy()
    nst.detrend('linear')
    nst.detrend('demean')
    nst.taper(max_percentage=0.05, type='cosine')                    # taper the edges
    #Remove instrumental response
    print('About to get rid of the response')
    pre_filt = (0.005, 0.006, 30.0, 35.0)
    st = Stream()
    for tr in nst:
        try:
            tr.remove_response(inventory=inv,output='VEL', pre_filt=pre_filt)
        except:
            print('Can NOT find response')
            print(tr.stats.station)
        tr.detrend('linear')
        tr.detrend('demean')
        tr.filter('bandpass', freqmin=0.04, freqmax=12, zerophase=True) # filter data of all traces in the streams
        print('DOWNSAMPLING')
        if tr.stats.sampling_rate == 200:
            tr.decimate(factor=10)
        elif tr.stats.sampling_rate == 100:
            tr.decimate(factor=5)
        print('Now lets normalize')
        tr = normalize(tr, norm_method='ramn',norm_win=12)
        tr = whiten(tr, 0.05, 2)
        st += tr       
    print ('done! - Time for cross-corelation')
    #Cross-correlation
    st.merge()
    for i in range(len(st)):
        for j in range(i+1,len(st)):
            station1 = st[i].stats.station
            station2 = st[j].stats.station
            ncorr, nstack = correlateNoise(st, [station1,station2], 3600)
            np.save('/Volumes/GeoPhysics_23/users-data/juarezilma/Test/V2/' + station1 + '_' + station2 + '_xcorr',ncorr.data)
            np.save('/Volumes/GeoPhysics_23/users-data/juarezilma/Test/V2/' + station1 + '_' + station2 + '_stack',nstack.data)
            plotStack(st,nstack,500,station1,station2,figurename='/Volumes/GeoPhysics_23/users-data/juarezilma/Test/V2/Stack_'+station1+'_'+station2+'.pdf')
    


