import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import _pickle as cPickle

# --- open a sample forecast file to use for grib data template

filename = '../ecmwf/2t_2019123100_f120.grib2'
ffcst = pygrib.open(filename)
grb = ffcst.read(1)[0] 
latsf, lonsf = grb.latlons()
nyf, nxf = np.shape(latsf)
print ('nyf, nxf = ', nyf, nxf)
print ('sample forecast file: ', filename)
print ('latsf = ',latsf[:,0])
print ('lonsf = ',lonsf[0,:])

# ---- loop through all ERA5 records in grib file

#infile = 'era5_reanalyses.grib'
infile = '/Volumes/Backup Plus/ecmwf/Jan2020.grib'

print (infile)
fanal = pygrib.open(infile)
#fanal.seek(0)
for grb in fanal:
    #msg = grb.tostring()
    cdatadate = str(grb.dataDate)
    chour = str(grb.hour)
    if chour == '0': chour = '00'
    cdatadate = cdatadate + chour
    latsa, lonsa = grb.latlons()
    t2m = grb.values
    nya, nxa = np.shape(latsa)
    #print ('nya, nxa = ', nya, nxa)
    #cftime = str(grb.forecastTime)
    #print ('cdatadate, cftime = ', cdatadate, cftime)
    temp_quarterdegree = grb.values
    nyq, nxq = np.shape(temp_quarterdegree)
    #print ('quarter degree nyq, nxq = ', nyq, nxq)
    temp_halfdegree = temp_quarterdegree[0::2,0::2]
    latsa_halfdegree = latsa[0::2,0::2]
    lonsa_halfdegree = lonsa[0::2,0::2]
    nyh, nxh = np.shape(temp_quarterdegree)
    #print ('half degree nyh, nxh = ', nyh, nxh)
    #print ('temp_quarterdegree[:,nxq//2] = ',temp_quarterdegree[:,2])
    #print ('temp_halfdegree[:,nxh//2] = ',temp_halfdegree[:,1])
    #print ('latsa_halfdegree[:,0] = ', latsa_halfdegree[:,0])
    #print ('lonsa_halfdegree[0,:] = ', lonsa_halfdegree[0,:])
    #print ('t2m[:,0] = ', t2m[:,0])
    #sys.exit()
    
    outfilename = '/Volumes/Backup Plus/ecmwf/t2m_era5_halfdegree_'+cdatadate+'.cPick'
    print (outfilename)
    ouf = open(outfilename, 'wb')
    cPickle.dump(temp_halfdegree, ouf)
    cPickle.dump(latsa_halfdegree, ouf)
    cPickle.dump(lonsa_halfdegree, ouf)
    ouf.close()    



