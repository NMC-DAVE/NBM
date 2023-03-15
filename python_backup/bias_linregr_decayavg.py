"""
using ECMWF forecast and ERA5 verification data, implement Kalman-filter
type regression bias correction approach to estimate forecast bias.   Grid covers
CONUS, approximately.  Coded by Tom Hamill, May-Jun 2020.
tom.hamill@noaa.gov
"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
import statsmodels.api as sm
from numba import jit
from os import path
    
# =====================================================================

#def Kalman_filter_regression(R, sigma2_intercept, sigma2_slope, \
#    regr_slope, regr_intercept, analysis, forecast):   
#
#    # --- Kalman-filter update our estimates of the regression coefficients with the most
#    #     recent forecast and analysis deviations from climatology
#
#    regr_corr_forecast = regr_intercept + regr_slope*forecast
#    Kalman_gain_intercept = sigma2_intercept / \
#        (sigma2_intercept + sigma2_slope*regr_corr_forecast + R )   # **** regression corrected or not?
#    Kalman_gain_slope = sigma2_slope*forecast / \
#        (sigma2_intercept + sigma2_slope*regr_corr_forecast + R )
#    regr_slope = regr_slope +  Kalman_gain_slope*(analysis-regr_corr_forecast)
#    regr_intercept = regr_intercept +  Kalman_gain_intercept*(analysis-regr_corr_forecast)
#    return regr_slope, regr_intercept


@jit
def Kalman_filter_regression(R, var_const, var_slope, cov_const_slope, \
    regr_slope, regr_intercept, analysis, forecast):   

    # --- Kalman-filter update our estimates of the regression coefficients with the most
    #     recent forecast and analysis deviations from climatology

    regr_corr_forecast = regr_intercept + regr_slope*forecast
    ny, nx = np.shape(var_const)
    P_beta = np.zeros((2, 2), dtype=np.float32)
    P_beta_inv = np.zeros((2, 2), dtype=np.float32)
    for ix in range(nx):
        for jy in range(ny):
            H = np.array([1.0, regr_corr_forecast[jy,ix]])
            P_beta[0,0] = var_const[jy,ix]
            P_beta[1,1] = var_slope[jy,ix]
            P_beta[1,0] = cov_const_slope[jy,ix]
            P_beta[0,1] = P_beta[1,0]
            #P_beta_inv = np.linalg.inv(P_beta)
            P_beta_HT = np.matmul(P_beta, np.transpose(H))
            H_Pbeta_HT = np.matmul(H, P_beta_HT)
            H_Pbeta_HT_plusR_inv = 1.0 / (H_Pbeta_HT + R)
            K_beta = P_beta_HT * H_Pbeta_HT_plusR_inv
            if jy == 27 and ix == 13:
                print ('prior intercept, slope = ', regr_intercept[jy,ix] , regr_slope[jy,ix] )
            regr_intercept[jy,ix] = regr_intercept[jy,ix] + \
                K_beta[0]*(analysis[jy,ix]-regr_corr_forecast[jy,ix])
            regr_slope[jy,ix] = regr_slope[jy,ix] + \
                K_beta[1]*(analysis[jy,ix]-regr_corr_forecast[jy,ix])
            if jy == 27 and ix == 13:
                print ('K_beta = ', K_beta)
                print ('posterior intercept, slope = ', regr_intercept[jy,ix] , regr_slope[jy,ix] )
                print ('analysis[27,13], regr_corr_forecast[jy,ix], forecast[jy,ix] = ', \
                    analysis[27,13], regr_corr_forecast[jy,ix], forecast[jy,ix])
            
    return regr_slope, regr_intercept

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
ilead = int(clead)
ndstart = ilead // 24
R = 1.0
datadir_reanl = '/Users/Tom/python/ecmwf/'
datadir = '/Users/Tom/python/ecmwf/'
cvariable = '2t'
tally_statistics = False
#datestart = dateshift('2018110100',ilead)
dateend = dateshift('2019123100',-ilead)
date_list_anal = daterange('2018110100',dateend,24) # initial time of the current forecast
ndates = len(date_list_anal)
date_list_fcst = []
for idate in range(ndates):
    date_list_fcst.append(dateshift(date_list_anal[idate],ilead)) # initial times of fcst

# ---- read in the time-dependent ERA5 climatology of t2m

infilename = 'ecmwf/t2m_climo_daily_era5_halfdegree.cPick'
print (infilename)
inf = open(infilename, 'rb')
climo_yearly = cPickle.load(inf)
print ('shape climo_yearly = ', np.shape(climo_yearly))
latsa_halfdegree = cPickle.load(inf)
lonsa_halfdegree = cPickle.load(inf)
nlats, nlons = np.shape(latsa_halfdegree)
inf.close()

infile = 'covariance_regression_coefficients_lead='+clead+'.cPick'
inf = open(infile, 'rb')
var_const = cPickle.load(inf)
var_slope = cPickle.load(inf)
cov_const_slope = cPickle.load(inf)

#var_const = var_const*2.0
#var_slope = var_slope*2.0
#cov_const_slope[:,:] = 0.0

inf.close()
mfact = 1.
#var_const = var_const * mfact
#var_slope = var_slope * mfact
#cov_const_slope = cov_const_slope*mfact

regr_intercept = np.zeros((nlats, nlons), dtype=np.float32)
regr_slope = np.ones((nlats, nlons), dtype=np.float32)
bias_estimate = np.zeros((nlats, nlons), dtype=np.float32)

#mfact = 4.0
#if clead == '24':
#    #sigma2_intercept = mfact*0.0054
#    #sigma2_slope = mfact*0.00027
#    sigma2_intercept = mfact*7.77
#    sigma2_slope = mfact*9.0316564e-05
#elif clead == '48':
#    sigma2_intercept = mfact*0.0073
#    sigma2_slope = mfact*0.00034
#elif clead == '72':
#    sigma2_intercept = mfact*0.0093
#    sigma2_slope = mfact*0.00043
#elif clead == '96':
#    sigma2_intercept = mfact*0.0122
#    sigma2_slope = mfact*0.00054

# ---- loop over dates and update bias estimates

for idate, datea in enumerate(date_list_anal):
    
    # ---- determine the julian day of the year
    
    yyyy,mm,dd,hh = splitdate(datea)
    doy = dayofyear(yyyy,mm,dd)
    datef = date_list_fcst[idate]
    print ('------ processing analysis, forecast dates = ', datea, datef)
    if int(datea) >= 2019010100: tally_statistics = True

    # ---- read the ECMWF ERA5 reanalysis at valid at the forecast date.
    #      Convert this to deviation from climatology
    
    infile = datadir_reanl + 't2m_era5_halfdegree_'+datef+'.cPick'
    #print (infile)
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats)
        print (nlats, nlons)
        analyzed_yearly = np.zeros((ndates,nlats,nlons), dtype=np.float32) 
        forecast_yearly = np.zeros((ndates,nlats,nlons), dtype=np.float32)    
    inf.close()
    analysis[:,:] = analysis[:,:] - climo_yearly[doy,:,:]
    #print ('max, min analysis = ', np.max(analysis), np.min(analysis))
    
    # ---- read the control forecast at this lead time and initial date.
    #      Convert this to deviation from climatology
 
    infile = datadir + cvariable+'_'+datea+'_f'+clead+'.grib2'   
    fexist = path.exists(infile)
    #print (infile, fexist)
    if fexist:
        grbfile = pygrib.open(infile) 
        grb = grbfile.select()[0] 
        forecast = grb.values
        grbfile.close()
        forecast[:,:] = forecast[:,:] - climo_yearly[doy,:,:]
        #print ('max, min forecast = ', np.max(forecast), np.min(forecast))
        #print ('max, min F-A = ', np.max(forecast-analysis), np.min(forecast-analysis))
    
        # ---- update the estimates of the regression coefficients with Kalman-filter
        #      type machinery.
    
        #regr_slope, regr_intercept = Kalman_filter_regression(R, \
        #    sigma2_intercept, sigma2_slope, \
        #    regr_slope, regr_intercept, analysis, forecast)
        
        regr_slope, regr_intercept = Kalman_filter_regression(R, \
            var_const, var_slope, cov_const_slope, \
            regr_slope, regr_intercept, analysis, forecast) 
        #print ('       max, min regr_slope = ', ma.max(regr_slope), ma.min(regr_slope)) 
        #print ('       max, min regr_intercept = ', ma.max(regr_intercept), ma.min(regr_intercept))  
        #print ('       mean slope, intercept = ', np.mean(regr_slope), np.mean(regr_intercept))
        #print ('       std slope, intercept = ', np.std(regr_slope), np.std(regr_intercept))
        #print ('       mean O-F  = ', np.mean(analysis)-np.mean(forecast))  
        #ind = np.unravel_index(np.argmin(regr_intercept, axis=None), regr_intercept.shape)
    
    # ---- write bias estimates to file if in 2019.

    if tally_statistics == True:
        outfilename = datadir + 'bias_KFregression_'+datef+'_f'+clead+'.cPick'
        #print ('       writing bias estimates to ', outfilename)

        ouf = open(outfilename, 'wb')
        cPickle.dump(regr_slope, ouf)
        cPickle.dump(regr_intercept, ouf)
        ouf.close()