"""

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.signal as signal
import scipy.stats as stats
from astropy.convolution import convolve


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
ilead = int(clead)
datadir = '/Users/Tom/python/ecmwf/'
cvariable = '2t'
datestart = dateshift('2018010100',ilead)
date_list_anal = daterange(datestart,'2018123100',24)
ndates = len(date_list_anal)
date_list_fcst = []
for idate in range(ndates):
    date_list_fcst.append(dateshift(date_list_anal[idate],-ilead)) # initial times of fcst
    
# ---- call initialization routine

for idate, datea in enumerate(date_list_anal):
    
    datef = date_list_fcst[idate]
    if datea == '2018010100': dstart = idate
    print ('------ processing analysis, forecast dates = ', datea, datef)

    # ---- read the ECMWF ERA5 reanalysis at this analysis date.
    
    infile = datadir + 't2m_era5_halfdegree_'+datea+'.cPick'
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats) 
        complete_fcst = np.zeros((ndates,nlats,nlons), dtype=np.float32)
        complete_anal = np.zeros((ndates,nlats,nlons), dtype=np.float32)
    inf.close()
    complete_anal[idate,:,:] = analysis[:,:]
    
    # ---- read the ECMWF control forecast at this lead time and initial date
 
    infile = datadir + cvariable+'_'+datef+'_f'+clead+'.grib2'  
    grbfile = pygrib.open(infile) 
    grb = grbfile.select()[0] 
    forecast = grb.values
    grbfile.close()
    complete_fcst[idate,:,:] = forecast[:,:]
 
# ---- estimate the slowly time varying bias.
  
window_std = 5. # 10
window = signal.gaussian(201,window_std)  
complete_difference = complete_fcst - complete_anal
bias_estimate = np.copy(complete_difference)    
random_error_estimate = np.copy(complete_difference) 
for ix in range(nlons):
    for jy in range(nlats):
        timeseries = complete_difference[:,jy,ix]
        bias_estimate[:,jy,ix]= convolve(timeseries, window)
        random_error_estimate[:,jy,ix] = timeseries[:] - bias_estimate[:,jy,ix]
    
# --- generate statistics 

fmean = np.mean(complete_fcst,axis=0)
amean = np.mean(complete_anal,axis=0)
bias = fmean-amean

# --- compute covariance statistics

print ('computing covariance statistics')
print ('nlons, nlats = ', nlons, nlats)
cov_bias = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float32)
cov_random = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float32)
var_bias = np.zeros((nlats, nlons), dtype=np.float32)
var_random = np.zeros((nlats, nlons), dtype=np.float32)
for j1 in range(nlats):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('processing j1 = ',j1,' of ',nlats,'. Current time = ',current_time)
    for i1 in range(nlons): 
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('processing i1 = ',i1,' of ',nlons,'. Current time = ',current_time)
        a1 = bias_estimate[:,j1,i1]
        a2 = random_error_estimate[:,j1,i1]
        for j2 in range(nlats):
            for i2 in range(nlons):
                b1 = bias_estimate[:,j2,i2]
                b2 = random_error_estimate[:,j2,i2]
                cov = np.cov(a1,b1) * 364./2.
                var_bias[j1,i1] = cov[0,0]* 364./2.
                cov_bias[j1,i1,j2,i2] = cov[1,0] * 364./2.
                cov = np.cov(a2,b2)
                var_random[j1,i1] = cov[0,0]
                cov_random[j1,i1,j2,i2] = cov[1,0]
                
# ---- write to cPickle file

outfile = 'covstats_bias_random_ecmwf2018_lead='+clead+'.cPick'
ouf = open(outfile,'wb')
print ('writing to ',outfile)
cPickle.dump(cov_bias, ouf)
cPickle.dump(cov_random , ouf)
cPickle.dump(var_bias, ouf)
cPickle.dump(var_random, ouf)
cPickle.dump(lats, ouf)
cPickle.dump(lons, ouf)
ouf.close()              


