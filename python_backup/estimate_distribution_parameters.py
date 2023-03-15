"""
estimate_distribution_parameters.py

estimate, for each month of the year, parameters of fitted
distributions for forecast and analyzed

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
from astropy.convolution import convolve
from det_params import det_params

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
    
    # ---- read reanalysis appropriate at the time.
    
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
        if idate == 0: 
            lats = cPickle.load(inf)
            lons = cPickle.load(inf)
            nlats, nlons = np.shape(lats)
            lats = np.flipud(lats)
            forecast_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            analysis_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
    else:
        print ('1. did not find file ', infile)
            
    # ---- read the forecast information .
            
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
    else:
        print ('2. did not find file ', infile)
    
    # ---- compute the differences between the forecast and analysis, add to time series
    
    forecast_3d[idate,:,:] = forecast[:,:] 
    analysis_3d[idate,:,:] = analysis[:,:] 
    
# ---- now estimate for each month and each grid point the distribution parameters.

power_analysis = np.zeros((12,nlats,nlons), dtype=np.float64)
mean_analysis = np.zeros((12,nlats,nlons), dtype=np.float64)
stddev_analysis = np.zeros((12,nlats,nlons), dtype=np.float64)
power_forecast = np.zeros((12,nlats,nlons), dtype=np.float64)
mean_forecast = np.zeros((12,nlats,nlons), dtype=np.float64)
stddev_forecast = np.zeros((12,nlats,nlons), dtype=np.float64)

for imonth in range(12):
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('**** PROCESSING month = ', imonth, current_time)
    
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
    forecast_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float32)
    analysis_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float32)

    ktr = 0
    for idate, date in enumerate(date_list_anal):
        if datevalid_indices[idate] == 1:
            forecast_validdates[ktr,:,:] = forecast_3d[idate,:,:]
            analysis_validdates[ktr,:,:] = analysis_3d[idate,:,:]
            ktr = ktr+1
    
    #print ('ktr, ndates_valid = ', ktr, ndates_valid)
    #print ('forecast_validdates[:,nlats//2,nlons//2] = ', forecast_validdates[:,nlats//2,nlons//2])
    #print ('analysis_validdates[:,nlats//2,nlons//2] = ', analysis_validdates[:,nlats//2,nlons//2])
    #sys.exit()
            
    # ---- call routine to test various power transformations, seeing which gives
    #      the best Lilliefors test.
    
    power_analysis[imonth,:,:], mean_analysis[imonth,:,:], \
        stddev_analysis[imonth,:,:] = det_params(ndates_valid, \
        nlats, nlons, analysis_validdates)
    power_forecast[imonth,:,:], mean_forecast[imonth,:,:], \
        stddev_forecast[imonth,:,:] = det_params(ndates_valid, \
        nlats, nlons, forecast_validdates)
    
# ---- save the optimal distributional data 
    
outfile = cpath_gefsv12+'GEFSv12_forecast_distribution_parameters_f'+clead+'.cPick'
ouf = open(outfile, 'wb')
cPickle.dump(power_forecast, ouf)
cPickle.dump(mean_forecast, ouf)
cPickle.dump(stddev_forecast, ouf)
ouf.close()

outfile = cpath_era5+'ERA5_analyzed_distribution_parameters_f'+clead+'.cPick'
ouf = open(outfile, 'wb')
cPickle.dump(power_analysis, ouf)
cPickle.dump(mean_analysis, ouf)
cPickle.dump(stddev_analysis, ouf)
ouf.close()
        





    
    

    
    
    