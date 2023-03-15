"""
save_2018_KFgain_v2.py

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
from Bbeta_model_4d import Bbeta_model_4d
from Kalman_filter_biascorrection2 import Kalman_filter_biascorrection2
from weighted_avg_biascorrection \
    import weighted_avg_biascorrection
from numba import jit
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

cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_forecast = '/Volumes/Backup Plus/ecmwf/'
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
    
for clead in ['24','48','72','96','120']:
    
    for warm_or_cold in ['warm', 'cold']:
        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print (current_time,'   &&&&&&&&&& LEAD TIME = ',clead,' SEASON = ',\
             warm_or_cold,'  &&&&&&&&&&')
    
        date_list_anal, date_list_forecast = \
            initialize_date_lists(warm_or_cold, cyear, clead)
        ndates = len(date_list_forecast)
        alpha =  set_alpha(clead)
    
        # ---- read the reanalysis time series on the dates specified.  Note that
        #      we pass in the dates of the forecast valid time, as we wish to 
        #      upload the analyses to compare against the forecasts at this time.

        analyses_3d, lats, lons = read_reanalysis_timeseries(cpath_era5, \
            date_list_forecast)
        nlats, nlons = np.shape(lons)
        npts = nlats*nlons
        R = form_diagonal_matrix(npts, anal_err_var)
    
        # ---- read the forecast time series on the dates specified.   Pass in 
        #      the lead time and the initial time of the analysis.

        forecast_3d, lats, lons = read_forecast_timeseries(cpath_forecast, \
            date_list_anal, clead, cvariable)
      
        # ---- verify the raw forecasts, i.e., with no bias correction. save stats

        beta_3d = np.zeros((ndates, nlats, nlons), dtype=np.float64)
        statsfile = cpath_errorstats + 'raw_forecast_errorstats_'+\
            clead+'h_'+warm_or_cold+'.txt'
        rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
            analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
        print ('   raw rmse, bias, mae = ', rmse, bias, mae)
        
        # ---- perform simplified Kalman filter of bias assuming diagonal Bias
        #      covariance matrix
        
        print ('   performing decay avg type bias correction')
        beta_3d = decay_avg_biascorrection2( npts, nlats, nlons, \
            forecast_3d, analyses_3d, beta_3d, date_list_forecast, \
            alpha)
    
        # ---- verify, save stats to file   
    
        statsfile = cpath_errorstats + 'decayavg_forecast_errorstats_'+\
            clead+'h_'+warm_or_cold+'.txt'
        rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
            analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
        print ('   decay avg rmse, bias, mae = ', rmse, bias, mae)    
        beta_3d_save = np.copy(beta_3d)
        
        # ---- produce a more careful estimate of the bias-corrected forecast
        #      error covariance, using previously generated correlation model
        #      based on land/sea, horiz and vert distances, and beta coeffient
        #      standard deviations from previous decaying avg type approach.
                
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        beta_3d = np.copy(beta_3d_save) 
        print ('   calling Bbeta_model ')
        Bbeta = Bbeta_model_4d(nlats, nlons, ndates, npts, \
            cyear, clead, warm_or_cold, cpath_Bbeta, beta_3d, already)
 
        # ---- apply the Kalman filter configured to use a more careful
        #      estimate of Bx, but still diagonal Bbeta.  
    
        print ('   producing Kalman filter type bias correction with improved Bbeta')
        beta_3d = weighted_avg_biascorrection(\
            npts, nlats, nlons, forecast_3d, analyses_3d, beta_3d, \
            date_list_forecast, Bbeta, alpha) 
                    
        if save_beta == True:
            beta_outfile = cpath_bias_est + '2018_bias_weightedsum_'+\
                cyear+'_'+warm_or_cold+'_lead'+clead+'.cPick'
            print ('   writing ', beta_outfile)
            ouf = open(beta_outfile,'wb')
            cPickle.dump(beta_3d, ouf)
            cPickle.dump(lats,ouf)
            cPickle.dump(lons,ouf)
            ouf.close()
    
        # ---- verify the spatially varying Bx bias correction forecasts

        statsfile = cpath_errorstats + '2018_weightsum_forecast_errorstats_'+\
            clead+'h.txt'
        rmse, bias, mae = verify_forecasts(ndates, nlats, nlons, clead, \
            analyses_3d, forecast_3d, beta_3d, iskip, statsfile)
        print ('   weighted sum type bias corr:  rmse, bias, mae = ',\
            rmse, bias, mae)
        
print ('Finished!')
    
    

    
    
    