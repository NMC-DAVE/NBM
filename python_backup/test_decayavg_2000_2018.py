"""
test_decayavg_2000_2018.py

"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
import numpy as np
import numpy.ma as ma
import math
import scipy
from decayavg_biascorrection2 import decayavg_biascorrection2
import matplotlib.pyplot as plt
from verify_forecasts_2000_2018 import verify_forecasts_2000_2018

# -------------------------------------------------------------- 

def read_lsmask():

    infilename = 'lsmask_0p5.cPick'
    print (infilename)
    inf = open(infilename, 'rb')
    lsmask = cPickle.load(inf)
    inf.close()
    return lsmask
  
# --------------------------------------------------------------   


clead = sys.argv[1]
ilead = int(clead)
alphas = [0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20]

iskip = int(clead) // 24
if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
    iskip = iskip+1
lsmask = read_lsmask()

# ---- various initialization

cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'

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
    if rem == 0: print ('processing date = ', idate, date)

    datef = date_list_forecast[idate]
    cyearf = datef[0:4]
    iyear = int(cyearf)-2000
    cmm = datef[4:6]
    imm = int(cmm)
    idd = int(datef[6:8])
    iyear_full = int(cyearf)
    julday = dayofyear(iyear_full, imm, idd) - 1
    if julday > 364: julday = 364
    
    infile = cpath_era5 +cyearf+'/t2m_era5_halfdegree_'+datef+'.cPick'
    fexist1 = os.path.exists(infile)
    if fexist1 == True:
        inf = open(infile, 'rb')
        analysis = cPickle.load(inf) - 273.16
        analysis = np.flipud(analysis)
        if idate == 0: 
            lats = cPickle.load(inf)
            lons = cPickle.load(inf)
            nlats, nlons = np.shape(lats)
            lats = np.flipud(lats)
            differences_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            analyses_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            forecast_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            print ('arrays established')
    else:
        print ('1. did not find file ', infile)
            
    # ---- read the t2m forecast information for bias corr.
            
    cyear = date[0:4]    
    cpath_forecast = cpath_gefsv12+cyear+'/'
    infile = cpath_forecast + date + '_lead'+clead+'_conus_0.5deg_hour'+clead+'.cPick'
    fexist2 = os.path.exists(infile)
    if fexist2 == True:
        inf = open(infile,'rb')
        forecast = cPickle.load(inf)
        inf.close()
    else:
        print ('2. did not find file ', infile)       
    
    # ---- subtract off climatology, save info to 3D arrays.

    analyses_3d[idate,:,:] = analysis[:,:] 
    forecast_3d[idate,:,:] = forecast[:,:] 
    
    
for ialpha, alpha in enumerate(alphas):
    
    # ---- get 2000-2018 decaying average bias correction
    
    beta_3d = np.zeros((ndates, nlats, nlons), dtype=np.float64)
    npts = ndates
    beta_3d = decayavg_biascorrection2(npts, nlats, nlons, \
        forecast_3d, analyses_3d, beta_3d, lsmask, \
        date_list_forecast, alpha) 
        
    statsfile = cpath_errorstats + \
        'decayavg_2000_2018_alpha'+str(alpha)+'_lead'+clead+'.txt'
    rmse, bia, mae = verify_forecasts_2000_2018(ndates, \
        nlats, nlons, clead, analyses_3d, forecast_3d, \
        beta_3d, lsmask, lons, lats, iskip, \
        statsfile)
        
    print ('alpha, rmse = ', alpha, rmse)
    

    
