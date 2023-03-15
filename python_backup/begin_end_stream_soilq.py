### begin_end_stream_soilq.py ###

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

cstreams = ['1999','2003','2007','2011','2015']

clonlow = sys.argv[1]
clatlow = sys.argv[2]
clonhi = sys.argv[3]
clathi = sys.argv[4]
ctitle = sys.argv[5]
rlonlow = float(clonlow)
rlatlow = float(clatlow)
rlonhi = float(clonhi)
rlathi = float(clathi)


# --- read in sample lat/lon indices

infile = '/Users/Tom/python/gefsv12/1999/bfg_2003123100_fhr00_control2.nc4'
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
nj = jmin - jmax 
ni = imax - imin 

# ---- loop thru streams

sw1_save = np.zeros((nj,ni,4), dtype=np.float32)
sw2_save = np.zeros((nj,ni,4), dtype=np.float32)
sw3_save = np.zeros((nj,ni,4), dtype=np.float32)

for istream, cstream in enumerate(cstreams[0:-1]):

    # ---- determine the dates

    cstream2 = cstreams[istream+1]
    print ('istream, cstream, cstream2 = ', istream, cstream, cstream2 )
    if cstream == '1999':
        date1 = '2003123100'
        date2 = '2004010100'
        date3 = '2004010200'
    elif cstream == '2003':
        date1 = '2007123100'
        date2 = '2008010100'
        date3 = '2008010200'
    elif cstream == '2007':
        date1 = '2011123100'
        date2 = '2012010100'
        date3 = '2012010200'
    elif cstream == '2011':
        date1 = '2015123100'
        date2 = '2016010100'
        date3 = '2016010200'
    
    infile1 = '/Users/Tom/python/gefsv12/'+cstream+'/bfg_'+date1+'_fhr00_control2.nc4'
    infile2 = '/Users/Tom/python/gefsv12/'+cstream2+'/bfg_'+date2+'_fhr00_control2.nc4'
    infile3 = '/Users/Tom/python/gefsv12/'+cstream2+'/bfg_'+date3+'_fhr00_control2.nc4'
   
    nc = Dataset(infile1)
    print (infile1)
    sw1 = nc.variables['soilw10_40cmdow'][0,jmax:jmin,imin:imax]
    ls = nc.variables['landsfc'][0,:,:]
    nc.close()

    nc = Dataset(infile2)
    print (infile2)
    sw2 = nc.variables['soilw10_40cmdow'][0,jmax:jmin,imin:imax]
    nc.close()
    
    nc = Dataset(infile3)
    print (infile3)
    sw3 = nc.variables['soilw10_40cmdow'][0,jmax:jmin,imin:imax]
    nc.close() 
    
    sw1_save[:,:,istream] = sw1[:,:] 
    sw2_save[:,:,istream] = sw2[:,:] 
    sw3_save[:,:,istream] = sw3[:,:]     

# ---- save to file.

outfile = 'gefsv12/'+ctitle+'_streamboundary_soilq.dump'
print ('writing to ', outfile)
ouf = open(outfile,'wb')
cPickle.dump(sw1_save, ouf)
cPickle.dump(sw2_save, ouf)
cPickle.dump(sw3_save, ouf)
ouf.close()

