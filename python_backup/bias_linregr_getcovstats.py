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
from os import path
import statsmodels.api as sm
    
# =====================================================================

def KF_linear_regres(Bbeta, B, R, KF_betahat, \
    obs, fcst, bias_estimate):

    # ---- estimate the Kalman gain for the bias correction.

    L = np.array([1.0, ])
    BbetaLT = np.matmul(Bbeta[:,:], np.transpose(L))
    LBbetaLT = np.matmul(L,BbetaLT)
    LBbetaLT_plus_B_plus_R = LBbetaLT + B + R
    LBbetaLT_plus_B_plus_R_inv = 1.0 / LBbetaLT_plus_B_plus_R
    Kfgain_beta = BbetaLT * LBbetaLT_plus_B_plus_R_inv

    # ---- update bias estimate with new data

    for i in range(3):
        KF_betahat[i,:,:] = KF_betahat[i,:,:] - \
            Kfgain_beta[i]*(obs[:,:] - (forecast[:,:] - bias_estimate[:,:]))
    bias_estimate = L[0]*KF_betahat[0,:,:] + \
        L[1]*KF_betahat[1,:,:] + L[2]*KF_betahat[2,:,:] + \
        L[3]*KF_betahat[3,:,:] + L[4]*KF_betahat[4,:,:]
    return KF_betahat, bias_estimate


# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
ilead = int(clead)
ndstart = ilead // 24
datadir_reanl = '/Users/Tom/python/ecmwf/'
datadir = '/Users/Tom/python/ecmwf/'

cvariable = '2t'
tally_statistics = False
#datestart = dateshift('2018110100',ilead)
dateend = dateshift('2018123100',-ilead)
date_list_anal = daterange('2018010100',dateend,24) # initial time of the current forecast
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
inf.close()

# ---- loop over dates and update bias estimates

for idate, datea in enumerate(date_list_anal):
    
    datef = date_list_fcst[idate]
    #print ('------ processing analysis, forecast dates = ', datea, datef)
    if int(datea) >= 2019010100: tally_statistics = True

    # ---- read the ECMWF ERA5 reanalysis at valid at the forecast date.
    
    infile = datadir_reanl + 't2m_era5_halfdegree_'+datef+'.cPick'
    #print (infile)
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats)
        print (nlats, nlons)
        analyzed_yearly = ma.zeros((ndates,nlats,nlons), dtype=np.float32) 
        forecast_yearly = ma.zeros((ndates,nlats,nlons), dtype=np.float32)    
    inf.close()
    analyzed_yearly[idate,:,:] = analysis[:,:] 
    
    # ---- read the control forecast at this lead time and initial date

 
    infile = datadir + cvariable+'_'+datea+'_f'+clead+'.grib2'  
    fexist = path.exists(infile)
    print (infile, fexist)
    if fexist == True:
        #print (infile)
        grbfile = pygrib.open(infile) 
        grb = grbfile.select()[0] 
        forecast = grb.values
        grbfile.close()
        forecast_yearly[idate,:,:] = forecast[:,:] 
    else:
        forecast_yearly[idate,:,:] = ma.masked
    
# ---- for each grid point, first produce a linear regression in order 
#      to get a sense of the standard error of the regression coefficients.
    
forecast_deviation = forecast_yearly - climo_yearly[ndstart:ndstart+ndates,:,:]
analyzed_deviation = analyzed_yearly - climo_yearly[ndstart:ndstart+ndates,:,:]

var_const = np.zeros((nlats,nlons), dtype = np.float32)
var_slope = np.zeros((nlats,nlons), dtype = np.float32)
cov_const_slope = np.zeros((nlats,nlons), dtype = np.float32)
regr_const = np.zeros((nlats,nlons), dtype = np.float32)
regr_slope = np.zeros((nlats,nlons), dtype = np.float32)

for ix in range(nlons):
    for jy in range(nlats):
        a = analyzed_deviation[:,jy,ix]
        f = forecast_deviation[:,jy,ix]
        x = ma.zeros((ndates,2), dtype=np.float32)
        x[:,0] = 1.0
        x[:,1] = f[:]
        model = sm.OLS(a, x, hasconst=True)
        results = model.fit()
        params = results.params
        cov = results.cov_params()
        var_const[jy,ix] = cov[0,0]
        var_slope[jy,ix] = cov[1,1]
        cov_const_slope[jy,ix] = cov[1,0]
        regr_const[jy,ix] = params[0]
        regr_slope[jy,ix] = params[1]


print ('mean slope, intercept = ', np.mean(regr_slope), np.mean(regr_const))
print ('std slope, intercept = ', np.std(regr_slope), np.std(regr_const))
ind = np.unravel_index(np.argmin(regr_const, axis=None), regr_const.shape)
print ('ind = ', ind)
print ('regr_const at min = ', regr_const[ind])

        

outfile = 'covariance_regression_coefficients_lead='+clead+'.cPick'
ouf = open(outfile, 'wb')
cPickle.dump(var_const, ouf)
cPickle.dump(var_slope, ouf)
cPickle.dump(cov_const_slope, ouf)
ouf.close()

print ('mean const, slope = ', np.mean(regr_const), np.mean(regr_slope))
print ('variance of const, slope = ', np.mean(var_const), np.mean(var_slope))
        
#KF_betahat, bias_seasonalKF = seasonalKFbias(cosfac, sinfac, \
#    cos2fac, sin2fac, Bbeta, B, R, KF_betahat, obs, \
#    forecast, bias_seasonalKF)
