"""
using ECMWF forecast and ERA5 verification data, implement Kalman-filter
type bias correction approach to estimate forecast bias.   Grid covers
CONUS, approximately.  Coded by Tom Hamill, May-Jun 2020.
tom.hamill@noaa.gov
"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import _pickle as cPickle

# =====================================================================

def initialize(ilead, alpha, nlats, nlons):
    
    # ---- initialize various
    
    R = 1.**2 # 1.0**2 # observation-error variance
    B = (0.75 + (ilead/48.))**2 # estimated forecast random-error variance, which grows with lead time
    Bbeta = np.zeros((5,5), dtype=np.float32) # simplified error covariance for seas'ly depdt coefficients
    Bbeta[0,0] = 0.12 
    Bbeta[1,1] = 0.12
    Bbeta[2,2] = 0.12
    Bbeta[3,3] = 0.12 
    Bbeta[4,4] = 0.12
    KF_betahat = np.zeros((5, nlats, nlons), dtype=np.float32)
    decay_betahat = np.zeros((nlats, nlons), dtype=np.float32) 
    return R, B, Bbeta, KF_betahat, decay_betahat

# =====================================================================

def cosfac_sinfac (date):
    
    # ---- compute cos, sin of Julian day/365.
    
    yyyy,mm,dd,hh = splitdate(date) 
    doy = dayofyear(yyyy,mm,dd)
    fac = 2.*3.14159*(np.real(doy)/365.)
    cosfac = np.cos(fac)
    sinfac = np.sin(fac)
    cos2fac = np.cos(2.*fac)
    sin2fac = np.sin(2.*fac)
    return cosfac, sinfac, cos2fac, sin2fac
    
# =====================================================================

def decayavg_bias(alpha, obs, forecast, bias_decayavg):
    
    # ---- compute the bog-standard decaying average bias correction estimate
       
    bias_decayavg = (1-alpha)*bias_decayavg[:,:] + alpha*(forecast[:,:]-obs[:,:])
    return bias_decayavg

# =====================================================================

def seasonalKFbias(cosfac, sinfac, cos2fac, sin2fac, Bbeta, B, R, KF_betahat, \
    obs, fcst, bias_estimate):

    # ---- estimate the Kalman gain for the bias correction.
    
    L = np.array([1.0, sinfac, cosfac, sin2fac, cos2fac])
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
calpha = sys.argv[2]  # alpha, specifying weighting of new vs. old data in 
    # decaying average bias correction and Kalman filter.
alpha = float(calpha)
ilead = int(clead)
datadir = '/Users/Tom/python/ecmwf/'
cvariable = '2t'
tally_statistics = False
datestart = dateshift('2018110100',ilead)
date_list_anal = daterange(datestart,'2019123100',24)
ndates = len(date_list_anal)
date_list_fcst = []
for idate in range(ndates):
    date_list_fcst.append(dateshift(date_list_anal[idate],-ilead)) # initial times of fcst

# ---- loop over dates and update bias estimates

for idate, datea in enumerate(date_list_anal):
    
    datef = date_list_fcst[idate]
    
    #print ('------ processing analysis, forecast dates = ', datea, datef)
    if int(datea) >= 2019010100: tally_statistics = True

    # ---- read the ECMWF ERA5 reanalysis at this analysis date.
    
    infile = datadir + 't2m_era5_halfdegree_'+datea+'.cPick'
    #print (infile)
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats)
        bias_decayavg = np.zeros((nlats, nlons), dtype=np.float32)
        bias_seasonalKF = np.zeros((nlats, nlons), dtype=np.float32)
        
        # ---- call initialization routine

        R, B, Bbeta, KF_betahat, decay_betahat = \
            initialize(ilead, alpha, nlats, nlons)
            
    inf.close()
    
    # ---- read the ECMWF control forecast at this lead time and initial date
 
    infile = datadir + cvariable+'_'+datef+'_f'+clead+'.grib2'  
    #print (infile)
    grbfile = pygrib.open(infile) 
    grb = grbfile.select()[0] 
    forecast = grb.values
    grbfile.close()
    
    # ---- read the ERA5 analysis valid at this date.
    
    infilename = datadir+'t2m_era5_halfdegree_'+datea+'.cPick'
    #print (infilename)
    inf = open(infilename, 'rb')
    obs = cPickle.load(inf)
    inf.close()    
    
    # ---- define cos, sin of day of year.
    
    cosfac, sinfac, cos2fac, sin2fac = cosfac_sinfac(datea)
    
    # ---- produce estimate of standard decaying-average bias correction
    
    bias_decayavg = decayavg_bias(alpha, obs, forecast, bias_decayavg)
    
    # ---- produce estimate of Kalman filter bias correction with seasonal variability.
    
    KF_betahat, bias_seasonalKF = seasonalKFbias(cosfac, sinfac, \
        cos2fac, sin2fac, Bbeta, B, R, KF_betahat, obs, \
        forecast, bias_seasonalKF)
        
    bias_decay_mean = np.mean(bias_decayavg)
    bias_KF_mean = np.mean(bias_seasonalKF)
    #print ('       mean decay, KF bias corrs= ', bias_decay_mean, bias_KF_mean)
        
    # ---- write bias estimates to file if in 2019.
    
    if tally_statistics == True:
        outfilename = datadir + 'bias_est_'+datef+'_f'+clead+'.cPick'
        #print ('       writing bias estimates to ', outfilename)
        ouf = open(outfilename, 'wb')
        cPickle.dump(bias_decayavg, ouf)
        cPickle.dump(bias_seasonalKF, ouf)
        ouf.close()
