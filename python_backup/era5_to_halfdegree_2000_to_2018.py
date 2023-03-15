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

for iyear in range(2000,2019):
    cyear = str(iyear)
    cyyyymmddhh_begin = cyear+'010100'
    cyyyymmddhh_end = cyear+'123100'
    date_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 24)
    for idate,date in enumerate(date_list):
        infile = '/Volumes/Backup Plus/ecmwf/'+cyear+'/2t_'+date+'_f0.grib2'
        print (infile)
        fanal = pygrib.open(infile)
        for grb in fanal:
            cdatadate = str(grb.dataDate)
            chour = str(grb.hour)
            if chour == '0': chour = '00'
            cdatadate = cdatadate + chour
            latsa, lonsa = grb.latlons()
            t2m = grb.values
            nya, nxa = np.shape(latsa)
            temp_quarterdegree = grb.values
            nyq, nxq = np.shape(temp_quarterdegree)
            temp_halfdegree = temp_quarterdegree[0::2,0::2]
            latsa_halfdegree = latsa[0::2,0::2]
            lonsa_halfdegree = lonsa[0::2,0::2]
            temp_halfdegree = np.flipud(temp_halfdegree)
            latsa_halfdegree = np.flipud(latsa_halfdegree)
            nyh, nxh = np.shape(temp_quarterdegree)
    
            outfilename = '/Volumes/Backup Plus/ecmwf/'+cyear+\
                '/t2m_era5_halfdegree_'+cdatadate+'.cPick'
            print (outfilename)
            ouf = open(outfilename, 'wb')
            cPickle.dump(temp_halfdegree, ouf)
            cPickle.dump(latsa_halfdegree, ouf)
            cPickle.dump(lonsa_halfdegree, ouf)
            ouf.close()    
            
            #sys.exit()



