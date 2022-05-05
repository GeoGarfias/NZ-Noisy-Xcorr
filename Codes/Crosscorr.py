import numpy as np
import sys, obspy, os
import matplotlib.pyplot as plt
from obspy import read, read_inventory
from obspy.core.stream import Stream
from obspy import UTCDateTime
from obspy.clients.fdsn import Client as FDSN_Client

# Functions by seismolive
# https://krischer.github.io/seismo_live_build/html/Ambient%20Seismic%20Noise/NoiseCorrelation_wrapper.html
# collection of functions used in noise correlation processing

def normalize(tr, clip_factor=6, clip_weight=10, norm_win=None, norm_method="1bit"): 
    '''
    Temporal normalization of the traces, most after Bensen 2007. NB. before this treatment, traces must be
    demeaned, detrended and filtered. Description of argument:

    norm_method="clipping"
        signal is clipped to 'clip_factor' times the std
        clip_factor recommended: 1 (times std)

    norm_method="clipping_iter"
        the signal is clipped iteratively: values above 'clip_factor * std' 
        are divided by 'clip_weight'. until the whole signal is below 
        'clip_factor * std'
        clip_factor recommended: 6 (times std)


    norm_method="ramn"
        running absolute mean normalization: a sliding window runs along the 
        signal. The values within the window are used to calculate a 
        weighting factor, and the center of the window is scaled by this 
        factor. 
            weight factor: w = np.mean(np.abs(tr.data[win]))/(2. * norm_win + 1) 
        finally, the signal is tapered with a tukey window (alpha = 0.2).

        norm_win: running window length, in seconds.
          recommended: half the longest period

    norm_method="1bit"
        only the sign of the signal is conserved
    '''
    
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
    '''
    Return tukey window of length N
    N: length of window
    alpha: alpha parameter in case of tukey window.
    0 -> rectangular window
    1 -> cosine taper
    returns: window (np.array)
    '''

    window = np.ones(N)
    x = np.linspace(-1., 1., N)
    ind1 = (abs(x) > 1 - alpha) * (x < 0)
    ind2 = (abs(x) > 1 - alpha) * (x > 0)
    window[ind1] = 0.5 * (1 - np.cos(np.pi * (x[ind1] + 1) / alpha))
    window[ind2] = 0.5 * (1 - np.cos(np.pi * (x[ind2] - 1) / alpha))
    return window

def whiten(tr, freqmin, freqmax):
    '''
    spectral whitening of trace `tr` using a cosine tapered boxcar between `freqmin` and `freqmax`
    (courtesy Gaia Soldati & Licia Faenza, INGV)
    '''
    
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
    '''
    correlate two stations, using slices of 'corrwin' seconds at a time correlations are also stacked. 
    NB hardcoded: correlates 1st with 2nd station in the stream only signals are merged - any data gaps are
    filled with zeros.
    st : stream containing data from the two stations to correlate
    stations : list of stations
    corrwin : correlation window length
    returns 'corr' (all correlations) and 'stack' (averaged correlations)
    '''
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
    print ("DONE")
  
    return corr,stack

def plotStack(st, stack2, maxlag, sta1, sta2, figurename=None):

    '''
    plots stack of correlations with correct time axis
        st: stream containing noise (and station information)
        stack: array containing stack       
        maxlag: maximum length of correlation to plot (in seconds)
    '''

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
        fig = plt.savefig(figurename, format="png")
    else:
        fig = plt.show()

def geonetData (t1, t2, stats):
    '''
    This functions gets data from GeoNet.
    Input:
        t1(start time) 
        t2 (end time)
        stats (station list)
    Output: 
        st (Obspy stream)
    '''
    client = FDSN_Client("GEONET")
    print( ' Downloading data from GeoNet ')
    for stat in stats:
        nst = client.get_waveforms(network='NZ', station=stat,location='*',starttime=t1,endtime=t2,channel='HHZ',attach_response=True)
    nst.merge(method=0, fill_value=0)
    return nst

def removeIR (stIR, pre_fil, wl, d_path):
    '''
    This function remove the instrumental response of the files you just downloaded.
    Saves the data in 'Data_path'
    Input:
        stIR: Strems that coints the traces of data
        pre_fil: (Prefilters)
        wl: (water_level)
        d_path: the path of the data
    Output:
        str (Obspy stream with the instrumental response removed)
    '''
    stIR.detrend('linear')
    stIR.detrend('demean')
    print('Removing instrumental response')
    stIR.remove_response(output='VEL', pre_filt=pre_fil, water_level=wl)
    stIR.detrend('linear')
    stIR.detrend('demean')
    stIR.taper(max_percentage=0.05, type='cosine')
    # Save stream into one big file
    for tr in stIR:
        name = str(tr.stats.starttime)
        tr.write(d_path +'/'+ str(tr.stats.station) +'_'+ str(tr.stats.channel) +'_'+ name[0:4] +'_'+ name[5:7] +'_'+ name[8:10] +'T'+ name[11:13] +'.mseed',format='MSEED')
    print('Save data in the directory')
    return stIR

def mdir(m_path,d_dir,f_dir):
    data_path = os.path.join(m_path,d_dir)
    fig_path = os.path.join(m_path,f_dir)
    if os.path.isdir(data_path) == False:
        os.mkdir(data_path)
    if os.path.isdir(fig_path) == False:
        os.mkdir(fig_path)
    
    return data_path,fig_path

if __name__ == '__main__':
    '''
    Test - Perform cross-correlation with HHZ data
    Code made by ICJG
    Input:
        main_path: Add the a path where you want to keep the data and figures this code produce.
        data_dir: name of the data directory
        fig_dir: name of the figures directory
        t1, t2: The time you want to analise. Take into account that the code is slow. stat_list or stat_pair: Name of GeoNet stations you want to analise. You can choose a lot of stations (stat_list) or only a pair
        (stat_pair).
        pre_filt: A prefilter to remove the instrument response.
        waterl: A value for the water level filter to remove instrument response.

        You can also modified:

        Sample rate
        bandpass filter
        type of time normalization
        frequencys of spectral normalization (whitening)
        the crosscorrelation window lenght
        maximum length of correlation to plot

    Output:
        data files: Every trace is saved in your local computer following the 'data_path'. Files are MSEED format and they're save as: STAT_CHANNEL_YEAR_MONTH_DAYTHOUR.mseed
        stack figures: Stacks of pair-stations are saved as 'png' figures. 

    '''
    main_path = '/Volumes/GeoPhysics_23/users-data/juarezilma/Noise/'
    # Name of the data directory.
    data_dir = 'Data'
    fig_dir = 'Figures'
    data_path, fig_path = mdir(main_path,data_dir,fig_dir)
    # Choose the time you want to analise . Not more that a few days, this code isnt the fastest.
    t1 = UTCDateTime('2020-02-15T00:00:00')
    t2 = UTCDateTime('2020-02-16T00:00:00')
    # Choose the GeoNet stations you want to analise.
    #stat_list = ['MQZ,OPZ,NNZ,WEL,WKZ,JCZ,INZ,FOZ,MSZ,WHSZ,WVZ,GVZ,KHZ,CVZ,ODZ,TUZ']
    stat_pair = ['FOZ,MQZ']
    # This line downloads data from GeoNet
    st = geonetData(t1,t2,stat_pair)
    # Pre filters for the remove of instrumental response
    pre_filt = (0.00001, 0.0001, 40.0, 48.0)
    waterl = 20
    # This removes the instrumental response
    st1 = removeIR(st,pre_filt,waterl,data_path)
    # ------------------------------------------------------------------------------------------------------------------
    #st1 = read('')
    for tr in st1:
        # Downsampling to 20sps
        if tr.stats.sampling_rate == 200:
            tr.decimate(factor=10)
        elif tr.stats.sampling_rate == 100:
            tr.decimate(factor=5)
        tr.filter('bandpass', freqmin=0.05, freqmax=1, zerophase=True) # filter data of all traces in the streams
        tr = normalize(tr, norm_method='ramn',norm_win=10)
        tr = whiten(tr, 0.05, 2)
    # Time for the crosscorrelation
    print('Cross-correlation begins')
    for i in range(len(st1)):
        for j in range(i+1,len(st1)):
            station1 = st1[i].stats.station
            station2 = st1[j].stats.station
            ncorr, nstack = correlateNoise(st1, [station1,station2], 3600)
            #np.save('/Volumes/GeoPhysics_23/users-data/juarezilma/Test/V2/' + station1 + '_' + station2 + '_xcorr',ncorr.data)
            #np.save('/Volumes/GeoPhysics_23/users-data/juarezilma/Test/V2/' + station1 + '_' + station2 + '_stack',nstack.data)
            plotStack(st1,nstack,500,station1,station2,figurename=fig_path +'/stack_'+station1+'_'+station2+'.png')
    


