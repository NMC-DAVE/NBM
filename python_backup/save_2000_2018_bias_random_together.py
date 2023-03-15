"""
save_2000_2018_bias_random.py

"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
import numpy as np
import numpy.ma as ma
import scipy.signal as signal
import scipy.stats as stats

# --------------------------------------------------------------   

# ---- various initialization

clead = sys.argv[1]
iskip = int(clead)//24
cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'

date_end = dateshift('2018123100',-2*int(clead))
date_list_anal = daterange('2000010100',date_end,24)

ndates = len(date_list_anal)
date_list_forecast = []
for i in range(ndates):
    date_list_forecast.append(dateshift(date_list_anal[i], int(clead)))
    
date_list_anal_ver = []
for i in range(ndates):
    date_list_anal_ver.append(dateshift(date_list_forecast[i], int(clead)))

# ---- loop through dates and process day by day

ndatektr = 0
ndatektr_yearly = np.zeros((18), dtype=np.int32)
for idate, date in enumerate(date_list_anal):
    
    # ---- read reanalysis appropriate at the time of this forecast for bias corr.
    
    rem = idate%30
    if rem == 0: print ('processing date = ', date)
    cyear = date[0:4]
    cmm = date[4:6]
    cmmdd = date[4:8]
    imm = int(cmm)
    iyear = int(cyear)-2000
    datef = date_list_forecast[idate]
    cyearf = datef[0:4]
    
    infile = cpath_era5 +cyearf+'/t2m_era5_halfdegree_'+datef+'.cPick'
    fexist1 = os.path.exists(infile)
    if fexist1 == True:
        inf = open(infile, 'rb')
        analysis = cPickle.load(inf) - 273.16
        analysis = np.flipud(analysis)
        if cmmdd == '0101': 
            lats = cPickle.load(inf)
            lons = cPickle.load(inf)
            nlats, nlons = np.shape(lats)
            lats = np.flipud(lats)
            differences_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            print (lats[:,0])
            print (analysis[:,0])
    else:
        print ('1. did not find file ', infile)
            
    # ---- read the forecast information for bias corr.
            
    cyear = date[0:4]    
    cpath_forecast = cpath_gefsv12+cyear+'/'
    infile = cpath_forecast + date + '_lead'+clead+'_conus_0.5deg_hour'+clead+'.cPick'
    fexist2 = os.path.exists(infile)
    if fexist2 == True:
        inf = open(infile,'rb')
        forecast = cPickle.load(inf)
        inf.close()
        if cmmdd == '0101':
            infile = '/Volumes/Backup Plus/gefsv12/t2m/gefsv12_latlon_subset.cPick'
            inf = open(infile,'rb')
            latsf = cPickle.load(inf)
            lonsf = cPickle.load(inf)
            inf.close()
            print (latsf[:,0])
            print (forecast[:,0])
            #sys.exit()
    else:
        print ('2. did not find file ', infile)
    
    # ---- compute the differences between the forecast and analysis, add to time series
    
    differences_3d[idate,:,:] = forecast[:,:] - analysis[:,:] 
    mae = np.mean(np.abs(forecast[:,:] - analysis[:,:]))
    if mae > 5.0: 
        print ('abnormally high mae ', mae,' for date ', date)
    #float_formatter = "{:.2f}".format
    #np.set_printoptions(precision=2, edgeitems=32, suppress=True)
    #if imm == 1: print (differences_3d[idate,:,:])
    

# ---- save the bias correction time series

bias_file = 'differences_2000_2018_lead'+clead+'.cPick'
ouf = open(bias_file, 'wb')
cPickle.dump(differences_3d, ouf)
cPickle.dump(date_list_anal, ouf)
ouf.close()



    
    

    
    
    