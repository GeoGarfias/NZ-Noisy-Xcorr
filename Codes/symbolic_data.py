import os

''' 
Make symbolic links for GEONET data from 2019
'''

BroadGeonet_stats = ['DSZ', 'THZ','WEL', 'INZ', 'LTZ', 'WVZ', 'OXZ', 'FOZ', 'RPZ', 'JCZ', 'LBZ', 'MSZ', 'WKZ', 'MLZ','EAZ','QRZ','NNZ','KHZ','GVZ','MQZ','CVZ','ODZ','OPZ','TUZ','SYZ','WHZ','DCZ','PYZ']
#StrongGeonet_stats = ['MQZ','OPZ','NNZ,','WEL']
#ShortGeonet_stats = ['ARCZ',]
years = ['2019','2020','2021']
original_path = '/Volumes/QuakeArchive_01/users-data/chambeca/GeoNet_archive/'
destination_path = '/Volumes/GeoPhysics_23/users-data/juarezilma/Crosscorrelations_TEST/Data_1/'
for y in years:
    julian_list = os.listdir('/Volumes/QuakeArchive_01/users-data/chambeca/GeoNet_archive/' + y + '/')
    
    stat_list = os.listdir('/Volumes/QuakeArchive_01/users-data/chambeca/GeoNet_archive/' + y +'/')
    for stat in BroadGeonet_stats:
        stat_folder = stat + '.NZ/'
        if stat_folder in stat_list:





src = '/usr/bin/python'
dst = '/tmp/python'

# This creates a symbolic link on python in tmp directory
os.symlink(src, dst)