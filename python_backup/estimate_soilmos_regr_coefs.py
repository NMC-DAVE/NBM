"""
estimate_soilmos_regr_coefs.py

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
import statsmodels.api as sm
from astropy.convolution import convolve
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
import math
import scipy

# --------------------------------------------------------------   

clead = sys.argv[1]

infilename = 'lsmask_0p5.cPick'
inf = open(infilename, 'rb')
lsmask = cPickle.load(inf)
inf.close()
lsmask = lsmask.astype(int)

# ---- various initialization

cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_soilw = '/Volumes/Backup Plus/gefsv12/soilw/'
cpath_cloud = '/Volumes/Backup Plus/gefsv12/cloud/'

date_end = dateshift('2018123100',-2*int(clead))
date_list_anal = daterange('2000010100',date_end,24)
ndates = len(date_list_anal)

date_list_forecast = []
for i in range(ndates):
    date_list_forecast.append(dateshift(date_list_anal[i], int(clead)))
    
# ---- read estimated 2000-2018 climatology to file

if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
    infile = cpath_era5 + 'ERA5_temperature_climatology_12UTC.cPick'
else:
    infile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'

print ('reading from ', infile)
inf = open(infile,'rb')
climo_temps_estimated = cPickle.load(inf)
inf.close()

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
    #print (iyear_full, imm, idd)
    julday = dayofyear(iyear_full, imm, idd) - 1
    if julday > 364: julday = 364
    
    infile = cpath_era5 +cyearf+'/t2m_era5_halfdegree_'+datef+'.cPick'
    fexist1 = os.path.exists(infile)
    #print (idate, infile, fexist1)
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
            forecast_3d_cloud = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            forecast_3d_soilw = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            print ('arrays established')
    else:
        print ('1. did not find file ', infile)
            
    # ---- read the  t2m forecast information for MOS.
            
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
    
    # ---- read the soilw forecast information for MOS.
    
    cyear = date[0:4]    
    cpath_forecast = cpath_soilw+cyear+'/'
    infile = cpath_forecast + date + '_lead'+clead+\
        '_soilw_0.1_to_0.4_conus_0.5deg_hour'+clead+'.cPick'
    fexist2 = os.path.exists(infile)
    if fexist2 == True:
        inf = open(infile,'rb')
        forecast_soilw = cPickle.load(inf)
        inf.close()
    else:
        print ('2. did not find file ', infile)
        
    # ---- read the cloud forecast information for MOS.
    
    cyear = date[0:4]    
    cpath_forecast = cpath_cloud+cyear+'/'
    infile = cpath_forecast + date + '_lead'+clead+\
        '_cldcover_conus_0.5deg_hour'+clead+'.cPick'
    fexist2 = os.path.exists(infile)
    if fexist2 == True:
        inf = open(infile,'rb')
        forecast_cloud = cPickle.load(inf)
        inf.close()
    else:
        print ('2. did not find file ', infile)
        
    # ---- subtract off climatology.

    #analyses_3d[idate,:,:] = analysis[:,:] - climo_temps_estimated[julday,:,:]
    #forecast_3d[idate,:,:] = forecast[:,:] - climo_temps_estimated[julday,:,:]
    analyses_3d[idate,:,:] = analysis[:,:] 
    forecast_3d[idate,:,:] = forecast[:,:] 
    forecast_3d_cloud[idate,:,:] = forecast_cloud[:,:] 
    forecast_3d_soilw[idate,:,:] = forecast_soilw[:,:] 
    
# ---- loop thru months

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#slope = np.zeros((12,nlats,nlons), dtype=np.float64)
#intercept = np.zeros((12,nlats,nlons), dtype=np.float64)
regr_coefs = np.zeros((12,7,nlats,nlons), dtype=np.float64)

for imonth, cmonth in enumerate(cmonths):
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('**** PROCESSING month = ', cmonth, current_time)
    
    # ---- loop through dates and process day by day

    datevalid_indices = np.zeros(ndates,dtype=np.int32)
    for idate, date in enumerate(date_list_anal):
        thismonth = int(str(date)[4:6]) - 1

        if imonth == 0:
            if thismonth == imonth or thismonth == 11 or \
                thismonth == imonth+1: datevalid_indices[idate] = 1
        elif imonth == 11:
            if thismonth == imonth or thismonth == 0 or \
                thismonth == imonth-1: datevalid_indices[idate] = 1
        else:
            if thismonth == imonth or thismonth == imonth+1 or \
                thismonth == imonth-1: datevalid_indices[idate] = 1
                
    ndates_valid = np.sum(datevalid_indices)
    forecast_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)
    forecast_validdates_soilw = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)
    forecast_validdates_cloud = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)
    analysis_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)

    ktr = 0
    for idate, date in enumerate(date_list_anal):
        if datevalid_indices[idate] == 1:
            analysis_validdates[ktr,:,:] = analyses_3d[idate,:,:]
            forecast_validdates[ktr,:,:] = forecast_3d[idate,:,:]
            forecast_validdates_soilw[ktr,:,:] = forecast_3d_soilw[idate,:,:]
            forecast_validdates_cloud[ktr,:,:] = forecast_3d_cloud[idate,:,:]
            ktr = ktr+1
   
    # ---- perform linear regression analysis for each grid point
    
    for ix in range(nlons):
        for jy in range(nlats):
            if lsmask[jy,ix] == 1:
                X = np.zeros((ktr,7), dtype=np.float64)
                X[:,0] = 1.0
                X[:,1] = forecast_validdates[:,jy,ix]
                X[:,2] = forecast_validdates_soilw[:,jy,ix]
                X[:,3] = forecast_validdates_cloud[:,jy,ix]
                X[:,4] = X[:,1] * X[:,2]
                X[:,5] = X[:,1] * X[:,3]
                X[:,6] = X[:,2] * X[:,3]

                y = np.zeros((ktr), dtype=np.float64)
                y[:] = analysis_validdates[:,jy,ix]

                model = sm.OLS(y, X)
                results = model.fit()
                #tvalues = results.tvalues
                regr_coefs[imonth,:,jy,ix] = results.params
            
    # ---- write the slope and intercept for this month to file.
    
    outfile = cpath_gefsv12+'MOS_soil_cloud_'+cmonth+\
        '_lead='+clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(regr_coefs[imonth,:,:], ouf)
    ouf.close()
    print ('   done writing')
    

    
    
    
