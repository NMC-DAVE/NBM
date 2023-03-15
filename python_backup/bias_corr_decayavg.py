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
import numpy.ma as ma
from os import path
import _pickle as cPickle
    
# =====================================================================

def decayavg_bias(alpha, obs, forecast, bias_decayavg):
    
    # ---- compute the bog-standard decaying average bias correction estimate
       
    bias_decayavg = (1-alpha)*bias_decayavg[:,:] + alpha*(forecast[:,:]-obs[:,:])
    return bias_decayavg

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
calpha = sys.argv[2]  # alpha, specifying weighting of new vs. old data in 
    # decaying average bias correction and Kalman filter.

# --- other initialization stuff    
    
alpha = float(calpha)
ilead = int(clead)
datadir_reanl = '/Users/Tom/python/ecmwf/'
datadir = '/Users/Tom/python/ncep/'
cvariable = '2t'
tally_statistics = False
#datestart = dateshift('2018110100',ilead)
dateend = dateshift('2019123100',-ilead)
date_list_anal = daterange('2018110100',dateend,24) # initial time of the current forecast
ndates = len(date_list_anal)
date_list_fcst = []
for idate in range(ndates):
    date_list_fcst.append(dateshift(date_list_anal[idate],ilead)) # initial times of fcst

# ---- loop over dates and update bias estimates

for idate, datea in enumerate(date_list_anal):
    
    datef = date_list_fcst[idate]
    #print ('------ processing analysis, forecast dates = ', datea, datef)
    if int(datea) >= 2019010100: tally_statistics = True

    # ---- read the ECMWF ERA5 reanalysis at valid at the forecast date.
    
    infile = datadir_reanl + 't2m_era5_halfdegree_'+datef+'.cPick'
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats)
        bias_decayavg = ma.zeros((nlats, nlons), dtype=np.float32)       
    inf.close()
    
    infile = datadir + cvariable+'_'+datea+'_f'+clead+'.grib2'  
    fexist = path.exists(infile)
    #print (infile, fexist)
    if fexist:
        
        # ---- read the control forecast at this lead time and initial date

        #print (infile)
        grbfile = pygrib.open(infile) 
        grb = grbfile.select()[0] 
        forecast = grb.values
        grbfile.close()
    
        # ---- produce estimate of standard decaying-average bias correction
    
        bias_decayavg = decayavg_bias(alpha, analysis, forecast, bias_decayavg)
        
    # ---- write bias estimates to file if in 2019.
    
    if tally_statistics == True:
        outfilename = datadir + 'bias_decayavg_alpha'+calpha+'_'+datef+'_f'+clead+'.cPick'
        print ('       writing bias estimates to ', outfilename)
        ouf = open(outfilename, 'wb')
        cPickle.dump(bias_decayavg, ouf)
        ouf.close()
