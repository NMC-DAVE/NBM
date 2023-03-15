from netCDF4 import Dataset
import numpy as np
from dateutils import daterange
import sys
import os
import os.path
from os import path
import numpy.ma as ma
import _pickle as cPickle

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# ---- commmand line inputs 

cstream = sys.argv[1]
clonlow = sys.argv[2]
clatlow = sys.argv[3]
clonhi = sys.argv[4]
clathi = sys.argv[5]
ctitle = sys.argv[6]
rlonlow = float(clonlow)
rlatlow = float(clatlow)
rlonhi = float(clonhi)
rlathi = float(clathi)


# --- read in sample lat/lon indices

infile = '/Volumes/Backup Plus/gefsv12/1999/bfg_2003123100_fhr00_control2.nc4'
nc = Dataset(infile)
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
nlons = len(lon)
nlats = len(lat)
nc.close()

# ---- determine the nearest grid index for box boundaries

imin = find_nearest(lon, rlonlow)
jmin = find_nearest(lat, rlatlow)
imax = find_nearest(lon, rlonhi)
jmax = find_nearest(lat, rlathi)

# ---- determine the date ranges to read in.

if cstream == '1999':
    date_list = daterange('2000010100','2003123100',24)
elif cstream == '2003':
    date_list = daterange('2004010100','2007123100',24)
elif cstream == '2007':
    date_list = daterange('2008010100','2011123100',24)
elif cstream == '2011':
    date_list = daterange('2012010100','2015123100',24)
elif cstream == '2015':
    date_list = daterange('2016010100','2019123100',24)
ndates = len(date_list)
swmean = np.zeros((ndates), dtype=np.float32)    
    
# ---- loop thru files and read in

for idate, date in enumerate(date_list):
    
    infile = '/Volumes/Backup Plus/gefsv12/'+cstream+'/bfg_'+date+'_fhr00_control2.nc4'
    print (infile)    
    does_soil_exist = path.exists(infile)
    if does_soil_exist == True:
        nc = Dataset(infile)
        sw = nc.variables['soilw10_40cmdow'][0,:,:]
        ls = nc.variables['landsfc'][0,:,:]
        nc.close()
        swmean[idate] = np.sum(sw[jmax:jmin,imin:imax]*ls[jmax:jmin,imin:imax]) / \
            np.sum(ls[jmax:jmin,imin:imax])
        print (idate, 'swmean = ', swmean[idate])
    else:  
        swmean[idate] = -99.99
        
        
swmean_masked = ma.masked_where(swmean <= 0.0, swmean)
date_list_vec = np.squeeze(np.asarray(date_list))

# ---- save to file.

outfile = 'gefsv12/'+ctitle+'_'+cstream+'_soilmoisture.dump'
print ('writing soil moisture time series to ',outfile)
swmean_masked.dump(outfile)

outfile = 'gefsv12/'+ctitle+'_'+cstream+'_datelist.dump'
print ('writing date_list time series to ',outfile)
date_list_vec.dump(outfile)
