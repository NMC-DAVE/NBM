"""
control_KF_biascorr_v2.py

controls sequential execution of a process for the potential
refinement of location- and date-dependent estimates of 2-m
temperature forecast bias using a more formal Kalman-filter
extension of the decaying-average bias correction.

The data used in this process are deterministic forecast data 
from the ECMWF prediction system in a grid encompassing the 
CONUS (-125 to -60 W, 20 to 50 N, 1/2-degree grid spacing).
The 2-m temperature analysis produced in the ERA5 analysis
system on the same grid is used as the verification and 
training data.

The sequential process follows these steps:

1.  Apply a simple decaying-average bias correction independently, 
grid point by grid point.   Actually, the formulation of the 
decaying average bias correction is expressed as a Kalman-filter 
bias correction, strongly inspired by Dick Dee's article,

https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1256/qj.05.137

and in particular eqs. (37) - (39) therein.  However, the B_x and
B_beta matrices therein are diagonal, and thus effectively 
become the decaying-average bias correction.  Validate these 
forecasts.

2.  With the bias estimates produced, correct the time series 
of forecasts, and use this corrected time series of forecasts
to produce a covariance matrix of the bias-corrected forecast 
errors.   From data assimilation experience, e.g., 
https://psl.noaa.gov/people/tom.hamill/covlocal_mwr.pdf
the localization of the matrix improves its accuracy, so 
apply localization.   Return the localized forecast error
covariance model.

3. Recompute the Kalman filter bias correction with a
spatially correlated forecast-error covariance matrix from 
step 2, but still independent estimates of the 
bias-correction coefficients. Validate these forecasts.

4.  Estimate the error covariance matrix of the bias
correction coefficients from the time series of bias estimates
at every grid point.   As with forecast-error covariances,
apply localization and return the localized bias-correction
error covariance.

5.  Apply a full Kalman-filter bias correction with the
more carefully validated forecast- and bias-correction error
covariances. Validate these forecasts.

"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
from read_reanalysis_timeseries import read_reanalysis_timeseries
from read_forecast_timeseries import read_forecast_timeseries
from forecast_error_covariance import forecast_error_covariance
from verify_forecasts import verify_forecasts
from Bxmodel import Bxmodel
from Bbeta_model import Bbeta_model
from Kalman_filter_biascorrection import Kalman_filter_biascorrection
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

cyear = '2018'
clead = '24'
iskip = int(clead)//24
cvariable = '2t'

cpath_era5 = '/Users/Tom/python/ecmwf/'
cpath_forecast = '/Volumes/Backup Plus/ecmwf/'
cpath_Bx = '/Volumes/Backup Plus/ecmwf/Bx/'
cpath_Bbeta = '/Volumes/Backup Plus/ecmwf/Bbeta/'
cpath_errorstats = '/Users/Tom/python/fcst_stats/'

anal_err_var = 1.0
exponenty = 2.0
already = False
    
for clead in ['24','48','72','96','120']:
    
    for warm_or_cold in ['warm', 'cold']:
        
        print ('&&&&&&&&&&&&&&&&&&&&&&&& LEAD TIME = ',clead,' SEASON = ',\
             warm_or_cold,'  &&&&&&&&&&&&&&&&&&&&&&&&&')
    
        date_list_anal, date_list_forecast = \
            initialize_date_lists(warm_or_cold, cyear, clead)
        ndates = len(date_list_forecast)
        fcst_err_var, b_beta_var =  set_error_variances(clead)
    
        # ---- read the reanalysis time series on the dates specified.  Note that
        #      we pass in the dates of the forecast valid time, as we wish to 
        #      upload the analyses to compare against the forecasts at this time.

        analyses_3d, lons, lats = read_reanalysis_timeseries(cpath_era5, \
            date_list_forecast)
        nlats, nlons = np.shape(lons)
        npts = nlats*nlons
    
        # ---- form a simple diagonal estimate for the Bx and Bbeta error
        #      covariance matrices.

        npts = nlons*nlats
        Bx = form_diagonal_matrix(npts, fcst_err_var)
        Bbeta = form_diagonal_matrix(npts, b_beta_var)
        R = form_diagonal_matrix(npts, anal_err_var)
    
        # ---- read the forecast time series on the dates specified.   Pass in 
        #      the lead time and the initial time of the analysis.

        forecast_3d, lons, lats = read_forecast_timeseries(cpath_forecast, \
            date_list_anal, clead, cvariable)
      
        # ---- verify the raw forecasts, i.e., with no bias correction

        beta_3d = np.zeros((ndates, nlats, nlons), dtype=np.float64)
        statsfile = cpath_errorstats + 'raw_forecast_errorstats_'+\
            clead+'h_'+warm_or_cold+'.txt'
        rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
            analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
        print ('   raw forecast rmse, bias, mae = ', rmse, bias, mae)

        # ---- apply the Kalman filter configured as a simple decaying-average
        #      bias correction.

        beta_3d = Kalman_filter_biascorrection(\
            npts, nlats, nlons, forecast_3d, analyses_3d, beta_3d, \
            date_list_forecast, R, Bx, Bbeta)
    
        # ---- verify, save to file   
    
        statsfile = cpath_errorstats + 'decayavg_forecast_errorstats_'+\
            clead+'h_'+warm_or_cold+'.txt'
        rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
            analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
        print ('   decay avg rmse, bias, mae = ', rmse, bias, mae)    
    
        beta_3d_save = np.copy(beta_3d)
        
        #for efold_random in [50.0,100.0,150.0,200.0]:
        #    for efold_bias in [400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0]:
        for efold_random in [250.0,300.0]:
            for efold_bias in [1700.0, 2000.0, 2500.0]:
    
                # ---- produce a more careful estimate of the bias-corrected forecast
                #      error covariance
    
                beta_3d = np.copy(beta_3d_save)
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
                beta_3d = Kalman_filter_biascorrection(\
                    npts, nlats, nlons, forecast_3d, analyses_3d, beta_3d, \
                    date_list_forecast, R, Bx_localized, Bbeta_localized)
    
                # ---- verify the spatially varying Bx bias correction forecasts

                strr = str(efold_random)
                strb = str(efold_bias)
                statsfile = cpath_errorstats + 'KF_forecast_errorstats_'+\
                    clead+'h_flocal'+strr+'_blocal'+strb+'.txt'
                rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
                    analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
                print ('   Localization, random and bias = ', efold_random, efold_bias,\
                    ' KF with Bx_localized, Bbeta_localized: rmse, bias, mae = ',\
                    rmse, bias, mae)
        
print ('Finished!')
    
    

    
    
    