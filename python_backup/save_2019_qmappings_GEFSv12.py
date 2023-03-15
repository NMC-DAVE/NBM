"""
save_2019_qmappings_GEFSv12.py
"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
from read_reanalysis_timeseries import read_reanalysis_timeseries
from read_forecast_timeseries_GEFSv12 import read_forecast_timeseries_GEFSv12
from read_forecast_timeseries_GEFSv12_cloud import read_forecast_timeseries_GEFSv12_cloud
from read_forecast_timeseries_GEFSv12_soilw import read_forecast_timeseries_GEFSv12_soilw
from verify_forecasts import verify_forecasts
from verify_forecasts_dual import verify_forecasts_dual
from decayavg_biascorrection2 import decayavg_biascorrection2
import numpy as np
from quantile_mapping import quantile_mapping
from KFgain_GEFSv12_biascorr_2019 import \
    KFgain_GEFSv12_biascorr_2019
from MOS_forecast_2019 import MOS_forecast_2019
from MOS_multiple_regr_2019 import MOS_multiple_regr_2019
from MOS_multiple_regr_soilw_2019 import MOS_multiple_regr_soilw_2019
from read_cloudfcst_timeseries_GEFSv12 import read_cloudfcst_timeseries_GEFSv12
from analog_forecast_2019 import analog_forecast_2019
from lowess_forecast_2019 import lowess_forecast_2019

# --------------------------------------------------------------

def form_diagonal_matrix(npts, vary):
    B = vary*np.identity(npts, dtype=np.float64) 
    return B
# --------------------------------------------------------------

def set_alpha(clead):
    if clead == '12':
        alpha = 0.16
    elif clead == '24':
        alpha = 0.14  
    elif clead == '36':
        alpha = 0.08
    elif clead == '48':
        alpha = 0.06 
    elif clead == '60':
        alpha = 0.06 
    elif clead == '72':
        alpha = 0.06 
    elif clead == '84':
        alpha = 0.06
    elif clead == '96':
        alpha = 0.04 
    elif clead == '108':
        alpha = 0.04 
    elif clead == '120':
        alpha = 0.04  
    elif clead == '132':
        alpha = 0.04
    return alpha

# --------------------------------------------------------------

def initialize_date_lists(cyear, clead):

    start_date = cyear+'010100'
    end_date = cyear+'123100'
    date_list_anal = daterange(start_date, end_date, 24)
    ndates = len(date_list_anal)
    date_list_forecast = []
    for i in range(ndates):
        date_list_forecast.append(dateshift(date_list_anal[i], \
            int(clead)))
    return date_list_anal, date_list_forecast

# --------------------------------------------------------------   

def write_bias_corrections(outfile, beta_3d, lats, lons):

    print ('   writing beta estimates to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(beta_3d, ouf)  
    cPickle.dump(lats, ouf) 
    cPickle.dump(lons, ouf) 
    ouf.close()
    istat = 0
    return istat
    
# --------------------------------------------------------------   

def read_lsmask():

    infilename = 'lsmask_0p5.cPick'
    #print (infilename)
    inf = open(infilename, 'rb')
    lsmask = cPickle.load(inf)
    inf.close()
    lsmask = lsmask.astype(int)
    return lsmask

# --------------------------------------------------------------  

def read_era5_climatology(cpath_era5, clead):

    if clead == '12' or clead == '36' or clead == '60' \
    or clead == '84' or clead == '120':
        infile = cpath_era5 + 'ERA5_temperature_climatology_12UTC.cPick'
    else:
        infile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'
    print ('reading from ', infile)
    inf = open(infile,'rb')
    climo_temps_estimated = cPickle.load(inf)
    inf.close()
    return climo_temps_estimated

# --------------------------------------------------------------  

# ---- various initialization

cyear = '2019'
clead = sys.argv[1]
iskip = int(clead)//24
cvariable = '2t'

cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'
cpath_forecast = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_cloud = '/Volumes/Backup Plus/gefsv12/cloud/'
cpath_soilw = '/Volumes/Backup Plus/gefsv12/soilw/'
cpath_gain = '/Volumes/Backup Plus/gefsv12/t2m/'
    
# ---- start the processing.

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print (current_time,'   &&&&&&&&&& LEAD TIME = ',clead,'  &&&&&&&&&&')
    
lsmask = read_lsmask()
date_list_anal, date_list_forecast = initialize_date_lists(cyear, clead)
ndates = len(date_list_forecast)
alpha = set_alpha(clead)


climo_temps_estimated = read_era5_climatology(cpath_era5, clead)
    
# ---- read the reanalysis time series on the dates specified.  Note that
#      we pass in the dates of the forecast valid time, as we wish to 
#      upload the analyses to compare against the forecasts at this time.

analyses_3d, lats, lons = read_reanalysis_timeseries(cpath_era5, \
    date_list_forecast)
#print (date_list_forecast)
#print ('anal = ',analyses_3d[180:211,20,40])
nlats, nlons = np.shape(lons)
npts = nlats*nlons
    
# ---- read the forecast time series on the dates specified.   Pass in 
#      the lead time and the initial time of the analysis.

#print (date_list_anal)
forecast_3d, latsf, lonsf = read_forecast_timeseries_GEFSv12 ( \
    cpath_forecast, date_list_anal, clead)
#forecast_3d_cloud = read_forecast_timeseries_GEFSv12_cloud( \
#    nlats, nlons, cpath_cloud, date_list_anal, clead)
#forecast_3d_soilw = read_forecast_timeseries_GEFSv12_soilw( \
#    nlats, nlons, cpath_soilw, date_list_anal, clead)
#print ('fcst = ',forecast_3d[180:211,20,40])

# ---- verify the raw forecasts, i.e., with no bias correction. save stats

beta_3d = np.zeros((ndates, nlats, nlons), dtype=np.float64)
ftype = 'raw_forecast_errorstats_2019_'
statsfile = cpath_errorstats + ftype + clead+'h.txt'

rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, beta_3d, lsmask, lons, lats, \
    climo_temps_estimated, iskip, ftype, statsfile)
print ('raw rmse, bias, mae = ', rmse, bias, mae)
             
# ---- perform analog bias correction, write to file
        
#print ('performing analog bias correction')
#beta_analog_3d = np.zeros((ndates, nlats, nlons), dtype=np.float64)
#beta_analog_3d = analog_forecast_2019(npts, nlats, nlons, \
#    forecast_3d, analyses_3d, beta_analog_3d, lsmask, \
#    date_list_forecast, clead, cpath_forecast, \
#    cpath_era5)
#outfile = cpath_forecast + '2019_analog_lead'+clead+'.cPick'
#istat = write_bias_corrections(outfile, beta_analog_3d, lats, lons)

# ---- verify, save decay average stats to file   
    
#ftype = 'analog_forecast_errorstats_2019_'    
#statsfile = cpath_errorstats + ftype + clead+'h.txt'
#rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
#    analyses_3d, forecast_3d, beta_analog_3d, lsmask, lons, lats, \
#    climo_temps_estimated, iskip, ftype, statsfile)
#print ('analog rmse, bias, mae = ', rmse, bias, mae)        
#sys.exit()         
        
# ---- produce cloud-soilw MOS for 2019 

#print ('performing MOS')
#beta_3d[:,:,:] = 0.0
#beta_3d = MOS_multiple_regr_soilw_2019(npts, nlats, nlons, \
#    analyses_3d, forecast_3d, forecast_3d_soilw, forecast_3d_cloud, \
#    beta_3d, lsmask, date_list_forecast, clead, cpath_forecast)
    
# ---- verify the MOS bias correction forecasts.  

#ftype = 'MOS_soilw_GEFSv12_2019_'
#statsfile = cpath_errorstats + ftype + clead+'h.txt'
#rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
#    analyses_3d, forecast_3d, beta_3d, lsmask, lons, lats, \
#    climo_temps_estimated, iskip, ftype, statsfile)
#print ('MOS soilw/cloud:  rmse, bias, mae = ', rmse, bias, mae)
#outfile = cpath_forecast + '2019_MOS_soilw_lead'+clead+'.cPick'
#istat = write_bias_corrections(outfile, beta_3d, lats, lons)               
#sys.exit()        
        

# ---- perform decay avg bias correction, write to file
        
print ('performing decay avg type bias correction')
beta_3d = decayavg_biascorrection2(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, lsmask, \
    date_list_forecast, alpha)
beta_decay_3d = np.copy(beta_3d)
outfile = cpath_forecast + '2019_decayavg_lead'+clead+'.cPick'
istat = write_bias_corrections(outfile, beta_3d, lats, lons)

# ---- verify, save decay average stats to file   
    
ftype = 'decayavg_forecast_errorstats_2019_'    
statsfile = cpath_errorstats + ftype + clead+'h.txt'
rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, beta_3d, lsmask, lons, lats, \
    climo_temps_estimated, iskip, ftype, statsfile)
print ('decay avg rmse, bias, mae = ', rmse, bias, mae)   

# ---- perform quantile mapping to as alternate way to estimate bias
        
print ('performing quantile mapping')
beta_qmap_3d = np.copy(beta_3d)
beta_qmap_3d = quantile_mapping(nlats, nlons, \
    forecast_3d, analyses_3d, beta_qmap_3d, lsmask, \
    date_list_forecast, cpath_gefsv12, cpath_era5, clead)
outfile = cpath_forecast + '2019_qmap_lead'+clead+'.cPick'
istat = write_bias_corrections(outfile, beta_qmap_3d, lats, lons)

# ---- verify, save quantile mapping to file   
   
ftype = 'quantile_mapping_errorstats_2019_'    
statsfile = cpath_errorstats + ftype + clead+'h.txt'
rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, beta_qmap_3d, lsmask, lons, lats, \
    climo_temps_estimated, iskip, ftype, statsfile)
print ('quantile mapping rmse, bias, mae = ', rmse, bias, mae) 

# ---- produce MOS for 2019 

print ('performing MOS')
beta_3d[:,:,:] = 0.0
beta_3d = MOS_forecast_2019(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, lsmask, \
    date_list_forecast, clead, cpath_forecast, \
    cpath_era5)
outfile = cpath_forecast + '2019_MOS_lead'+clead+'.cPick'
istat = write_bias_corrections(outfile, beta_3d, lats, lons)
    
# ---- verify the MOS bias correction forecasts.  

ftype = 'MOS_GEFSv12_2019_'
statsfile = cpath_errorstats + ftype + clead+'h.txt'
rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, beta_3d, lsmask, lons, lats, \
    climo_temps_estimated, iskip, ftype, statsfile)
print ('MOS:  rmse, bias, mae = ', rmse, bias, mae)
outfile = cpath_forecast + '2019_MOS_lead'+clead+'.cPick'
istat = write_bias_corrections(outfile, beta_3d, lats, lons) 
    
# ---- produce multi-variate MOS for 2019 

print ('performing MOS multiple regression')
beta_3d = MOS_multiple_regr_2019(npts, nlats, nlons, \
    analyses_3d, forecast_3d, beta_decay_3d, beta_qmap_3d, \
    beta_3d, lsmask, date_list_forecast, clead, \
    cpath_forecast, cpath_era5)
                    
# ---- verify the multi-variate MOS forecasts.  

ftype = 'MOS_mvr_GEFSv12_2019_'
statsfile = cpath_errorstats + ftype + clead+'h.txt'
rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, beta_3d, lsmask, lons, lats, \
    climo_temps_estimated, iskip, ftype, statsfile)
print ('MOS multiple regression:  rmse, bias, mae = ', rmse, bias, mae)
outfile = cpath_forecast + '2019_MOS_mvr_lead'+clead+'.cPick'
istat = write_bias_corrections(outfile, beta_3d, lats, lons) 

# ---- produce Kalman filter bias correction for 2019 using 2018 gain estimates

#print ('performing Kalman filter')
#beta_3d = KFgain_GEFSv12_biascorr_2019(npts, nlats, nlons, \
#    forecast_3d, analyses_3d, beta_3d, lsmask, \
#    date_list_forecast, clead, cpath_gain)
    
# ---- verify the spatially varying Bx bias correction forecasts.  

#ftype = 'KF_corr_GEFSv12_2019_'
#statsfile = cpath_errorstats + ftype + clead+'h.txt'

#rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
#    analyses_3d, forecast_3d, beta_3d, lsmask, lons, lats, \
#    iskip, ftype, statsfile)
#print ('KF bias corr:  rmse, bias, mae = ', rmse, bias, mae)
#outfile = cpath_forecast + '2019_KF_lead'+clead+'.cPick'
#istat = write_bias_corrections(outfile, beta_3d, lats, lons)

# ---- verify linear combo of qmap and decay avg. 
    
#print ('performing avg of qmapping + decaying average')    
#statsfile = cpath_errorstats + 'dual_errorstats_2019_'+\
#    clead+'h.txt'
#ftype = 'dual'
#rmse, bias, mae = verify_forecasts_dual(ndates, nlats, \
#    nlons, clead, analyses_3d, forecast_3d, beta_3d_save, lsmask, \
#    beta_qmap, iskip, ftype, statsfile)
#print ('1/2 qmap, 1/2 decay rmse, bias, mae = ', rmse, bias, mae)    
    

    
    
    