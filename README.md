# NZ-Noisy-Xcorr
So far there's only one code in this repository 'Crosscorr.py'
_______________________________________________________________________________________________________________________________________________________
1. Crosscorr.py:
   This codes contains function from Seismo-live** and from me.
   Variables you can modified from this code:
   - main_path: Add the a path where you want to keep the data and figures this code produce.
   - data_dir: name of the data directory
   - fig_dir: name of the figures directory
   - t1, t2: The time you want to analise. Take into account that the code is slow. 
              stat_list or stat_pair: Name of GeoNet stations you want to analise. You can choose a lot of stations (stat_list) or only  a pair       
              (stat_pair).
   - pre_filt: A prefilter to remove the instrument response.
   - waterl:  A value for the water level filter to remove instrument response.
  You can also modified:
   - Sample rate
   - bandpass filter
   - type of time normalization
   - frequencys of spectral normalization (whitening)
   - the crosscorrelation window lenght
   - maximum length of correlation to plot
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ** Seismo-live: https://krischer.github.io/seismo_live_build/html/Ambient%20Seismic%20Noise/NoiseCorrelation_wrapper.html
