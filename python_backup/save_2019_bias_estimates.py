"""
save_2019_bias_estimates.py

"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
from read_reanalysis_timeseries import read_reanalysis_timeseries
from read_forecast_timeseries import read_forecast_timeseries
#from forecast_error_covariance import forecast_error_covariance
from verify_forecasts import verify_forecasts
from decayavg_biascorrection2 import decayavg_biascorrection2
from Kalman_filter_biascorr_2019 import Kalman_filter_biascorr_2019
import numpy as np

# --------------------------------------------------------------
def form_diagonal_matrix(npts, vary):
    B = vary*np.identity(npts, dtype=np.float64) 
    return B
# --------------------------------------------------------------

def set_alpha(clead):
    if clead == '24':
        alpha = 0.16  # will yield approx 0.08 alpha coefficient if R=Bx=1.0
    elif clead == '48':
        alpha = 0.12  
    elif clead == '72':
        alpha = 0.08  
    elif clead == '96':
        alpha = 0.07  
    elif clead == '120':
        alpha = 0.06  

    return alpha

# --------------------------------------------------------------

def initialize_date_lists(cyear, clead):

    start_date = cyear+'010100'
    end_date = cyear+'123100'
    date_list_anal = daterange(start_date, end_date, 24)
    ndates = len(date_list_anal)
    date_list_forecast = []
    for i in range(ndates):
        date_list_forecast.append(dateshift(date_list_anal[i], int(clead)))
    return date_list_anal, date_list_forecast

# --------------------------------------------------------------   

# ---- various initialization

cyear = '2019'
clead = sys.argv[1]
iskip = int(clead)//24
cvariable = '2t'

cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_forecast = '/Volumes/Backup Plus/ecmwf/'
cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'
cpath_gain = '/Volumes/Backup Plus/python/KFgain/'
    
# ---- start the processing.

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print (current_time,'   &&&&&&&&&& LEAD TIME = ',clead,'  &&&&&&&&&&')
    
date_list_anal, date_list_forecast = initialize_date_lists(cyear, clead)
ndates = len(date_list_forecast)
alpha = set_alpha(clead)
    
# ---- read the reanalysis time series on the dates specified.  Note that
#      we pass in the dates of the forecast valid time, as we wish to 
#      upload the analyses to compare against the forecasts at this time.

analyses_3d, lats, lons = read_reanalysis_timeseries(cpath_era5, \
    date_list_forecast)
nlats, nlons = np.shape(lons)
npts = nlats*nlons
    
# ---- read the forecast time series on the dates specified.   Pass in 
#      the lead time and the initial time of the analysis.

forecast_3d, lats, lons = read_forecast_timeseries(cpath_forecast, \
    date_list_anal, clead, cvariable)
      
# ---- verify the raw forecasts, i.e., with no bias correction. save stats

beta_3d = np.zeros((ndates, nlats, nlons), dtype=np.float64)
beta_3d_save = np.copy(beta_3d)
statsfile = cpath_errorstats + 'raw_forecast_errorstats_2019_'+\
    clead+'h.txt'
rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
print ('raw rmse, bias, mae = ', rmse, bias, mae)
        
# ---- perform simplified Kalman filter of bias assuming diagonal Bias
#      covariance matrix
        
print ('performing decay avg type bias correction')
beta_3d = decayavg_biascorrection2(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, date_list_forecast, alpha)
    
# ---- verify, save decay average stats to file   
    
statsfile = cpath_errorstats + 'decayavg_forecast_errorstats_2019_'+\
    clead+'h_.txt'
rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
print ('decay avg rmse, bias, mae = ', rmse, bias, mae)    
  
# ---- produce a more careful estimate of the bias-corrected forecast
#      error covariance, using previously generated correlation model
#      based on land/sea, horiz and vert distances, and beta coeffient
#      standard deviations from previous decaying avg type approach.
                
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
beta_3d = np.copy(beta_3d_save) 

# ---- produce Kalman filter bias correction for 2019 using 2018 gain estimates
                    
beta_3d = Kalman_filter_biascorr_2019(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, \
    date_list_forecast, clead, cpath_gain)
    
# ---- verify the spatially varying Bx bias correction forecasts

statsfile = cpath_errorstats + 'KFbiascorrection_2019_'+\
    clead+'h.txt'
rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
print ('KF bias corr:  rmse, bias, mae = ',\
    rmse, bias, mae)
        
print ('Finished!')
    
    

    
    
    