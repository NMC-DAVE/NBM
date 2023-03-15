"""
save_2018_bias_estimates.py

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
from Bxmodel import Bxmodel
from Bbeta_model import Bbeta_model
from Kalman_filter_biascorrection import Kalman_filter_biascorrection
from Kalman_filter_biascorrection_savegain \
    import Kalman_filter_biascorrection_savegain
from numba import jit
import numpy as np

# --------------------------------------------------------------
def form_diagonal_matrix(npts, vary):
    B = vary*np.identity(npts, dtype=np.float64) 
    return B
# --------------------------------------------------------------

def set_error_variances(clead):
    if clead == '24':
        fcst_err_var = 1.0
        b_beta_var = 0.16  # will yield approx 0.08 alpha coefficient if R=Bx=1.0
    elif clead == '48':
        fcst_err_var = 2.0
        b_beta_var = 0.12  
    elif clead == '72':
        fcst_err_var = 3.0
        b_beta_var = 0.08  
    elif clead == '96':
        fcst_err_var = 4.0
        b_beta_var = 0.07  
    elif clead == '120':
        fcst_err_var = 5.0
        b_beta_var = 0.06  

    return fcst_err_var, b_beta_var

# --------------------------------------------------------------

def initialize_date_lists(warm_or_cold, cyear, clead):

    if warm_or_cold == 'warm' :
        start_date = cyear+'040100'
        end_date = cyear+'093000'
        date_list_anal = daterange(start_date, end_date, 24)
        ndates = len(date_list_anal)
        date_list_forecast = []
        for i in range(ndates):
            date_list_forecast.append(dateshift(date_list_anal[i], int(clead)))
    else:
        
        start_date = cyear+'010100'
        end_date = cyear+'033100'
        date_list_anal1 = daterange(start_date, end_date, 24)
        ndates = len(date_list_anal1)
        date_list_forecast1 = []
        for i in range(ndates):
            date_list_forecast1.append(dateshift(date_list_anal1[i], int(clead)))
            
        start_date = cyear+'100100'
        end_date = cyear+'123100'
        date_list_anal2 = daterange(start_date, end_date, 24)
        ndates = len(date_list_anal2)
        date_list_forecast2 = []
        for i in range(ndates):
            date_list_forecast2.append(dateshift(date_list_anal2[i], int(clead)))
            
        date_list_anal = date_list_anal2 + date_list_anal1
        date_list_forecast = date_list_forecast2 + date_list_forecast1

    return date_list_anal, date_list_forecast

# --------------------------------------------------------------   

# ---- various initialization

clead = sys.argv[1]
warm_or_cold = sys.argv[2]
cefold_random = sys.argv[3]
cefold_bias = sys.argv[4]
cyear = '2018'
iskip = int(clead)//24
cvariable = '2t'
efold_random = float(cefold_random)
efold_bias = float(cefold_bias)

cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_forecast = '/Volumes/Backup Plus/ecmwf/'
cpath_Bx = '/Volumes/Backup Plus/ecmwf/Bx/'
cpath_Bbeta = '/Volumes/Backup Plus/ecmwf/Bbeta/'
cpath_decay = '/Volumes/Backup Plus/ecmwf/biascorr/'
cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'
cpath_gain = '/Volumes/Backup Plus/python/KFgain/'
cpath_bias_est = '/Volumes/Backup Plus/python/bias_est/'

anal_err_var = 1.0
exponenty = 2.0
already = False
already_decay = False
save_beta = True
    
date_list_anal, date_list_forecast = \
    initialize_date_lists(warm_or_cold, cyear, clead)
ndates = len(date_list_forecast)
fcst_err_var, b_beta_var =  set_error_variances(clead)
    
# ---- read the reanalysis time series on the dates specified.  Note that
#      we pass in the dates of the forecast valid time, as we wish to 
#      upload the analyses to compare against the forecasts at this time.

print ('reading analyses')
analyses_3d, lats, lons = read_reanalysis_timeseries(cpath_era5, \
    date_list_forecast)
nlats, nlons = np.shape(lons)
npts = nlats*nlons
    
# ---- read the forecast time series on the dates specified.   Pass in 
#      the lead time and the initial time of the analysis.

print ('reading forecasts')
forecast_3d, lats, lons = read_forecast_timeseries(cpath_forecast, \
    date_list_anal, clead, cvariable)
      
# ---- verify the raw forecasts, i.e., with no bias correction

print ('decay avg bias correction')
beta_3d = np.zeros((ndates, nlats, nlons), dtype=np.float64)
beta_3d_save = np.copy(beta_3d)
Bx = form_diagonal_matrix(npts, fcst_err_var)
Bbeta = form_diagonal_matrix(npts, b_beta_var)
R = form_diagonal_matrix(npts, anal_err_var)
        
savefile = cpath_decay+'2018_'+warm_or_cold+'_beta_3d.cPick'
beta_3d = Kalman_filter_biascorrection( \
    npts, nlats, nlons, forecast_3d, analyses_3d, beta_3d, \
    date_list_forecast, R, Bx, Bbeta, savefile, already_decay)

# ---- produce a more careful estimate of the bias-corrected forecast
#      error covariance
                
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('processing random = ',efold_random,' bias = ',efold_bias,current_time)
strr = str(efold_random)
strb = str(efold_bias)
beta_3d_save = np.copy(beta_3d)
Bx_localized = Bxmodel(nlats, nlons, ndates, lats, lons, cyear, clead, \
    warm_or_cold, cpath_Bx, efold_random, exponenty, forecast_3d, \
    analyses_3d, beta_3d, already)
 
# ---- formulate and localize model for covariance of bias-correction
#      estimates
 
print ('calling Bbeta_model')
Bbeta_localized = Bbeta_model(nlats, nlons, ndates, npts, \
    cyear, clead, warm_or_cold, cpath_Bbeta, efold_bias, \
    exponenty, beta_3d, already)
 
# ---- apply the Kalman filter configured to use a more careful
#      estimate of Bx, but still diagonal Bbeta.  
    
print ('producing Kalman filter bias correction with improved Bx')
gain_outfile = cpath_gain + '2018_KFgain_flocal'+strr+\
    '_blocal'+strb+'_'+cyear+'_'+warm_or_cold+'_lead'+\
    clead+'.cPick'
beta_3d = Kalman_filter_biascorrection_savegain(\
    npts, nlats, nlons, forecast_3d, analyses_3d, beta_3d, \
    date_list_forecast, R, Bx_localized, Bbeta_localized, \
    gain_outfile, already)
print ('max, min beta_3d = ', np.max(beta_3d), np.min(beta_3d))
                    
if save_beta == True:
    beta_outfile = cpath_bias_est + '2018_bias_est_flocal'+strr+\
        '_blocal'+strb+'_'+cyear+'_'+warm_or_cold+'_lead'+\
        clead+'.cPick'
    print (beta_outfile)
    ouf = open(beta_outfile,'wb')
    cPickle.dump(beta_3d, ouf)
    cPickle.dump(lats,ouf)
    cPickle.dump(lons,ouf)
    ouf.close()
    
# ---- verify the spatially varying Bx bias correction forecasts

statsfile = cpath_errorstats + '2018_KF_forecast_errorstats_'+\
    clead+'h_flocal'+strr+'_blocal'+strb+'.txt'
rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
    analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
print ('   Localization, random and bias = ', efold_random, efold_bias,\
    ' KF with Bx_localized, Bbeta_localized: rmse, bias, mae = ',\
    rmse, bias, mae)
        
print ('Finished!')
    
    

    
    
    