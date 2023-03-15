"""
compare_errors_decay_conuskf.py

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

datadir = '/Users/Tom/python/ecmwf/'
cvariable = '2t'

rmse_kf_mean = np.zeros(4,dtype=np.float32)
rmse_decayavg_mean = np.zeros(4,dtype=np.float32)
rmse_raw_mean = np.zeros(4,dtype=np.float32)
calpha_optimal = ['0.10', '0.06', '0.04', '0.04']

for ilead in range(24,97,24):
    clead = str(ilead)
    datestart = dateshift('2019010100',ilead)
    dateend = dateshift('2019123100',-ilead)
    date_list_initial = daterange(datestart,dateend,24)
    print (date_list_initial)
    #date_list_anal = daterange(datestart,'2019013100',24)
    ndates = len(date_list_initial)
    rmse_kf = np.zeros((ndates), dtype=np.float32)
    rmse_decayavg = np.zeros((ndates), dtype=np.float32)
    rmse_raw = np.zeros((ndates), dtype=np.float32)
    
    ilead_idx = ilead//24 -1
    calpha = calpha_optimal[ilead_idx]
    date_list_valid = []
    for idate in range(ndates):
        date_list_valid.append(dateshift(date_list_initial[idate],ilead)) # valid times
        
    for idate, datea in enumerate(date_list_initial):
    
        datev = date_list_valid[idate]
        if datev == '2019010100': dstart = idate

        # ---- read the ECMWF ERA5 reanalysis at this analysis date.
    
        infile = datadir + 't2m_era5_halfdegree_'+datev+'.cPick'
        print (infile)
        inf = open(infile, 'rb')
        analysis = cPickle.load(inf)
        if idate == 0:
            lats = cPickle.load(inf)
            lons = cPickle.load(inf)
            nlats, nlons = np.shape(lats)
            npts = nlats*nlons 
        inf.close()
    
        # ---- read the ECMWF control forecast at this lead time and initial date
 
        infile = datadir + cvariable+'_'+datea+'_f'+clead+'.grib2'  
        grbfile = pygrib.open(infile) 
        grb = grbfile.select()[0] 
        forecast = grb.values
        grbfile.close()
    
        # ---- read the decaying average and Kalman filter bias corrections estimates
        #      for this date.
    
        infilename = datadir + 'bias_decayavg_alpha'+calpha+'_'+datea+'_f'+clead+'.cPick'
        inf = open(infilename, 'rb')
        bias_decayavg = cPickle.load(inf)
        inf.close()
        
        infilename = datadir + 'bias_est_conusKF'+datea+'_f'+clead+'.cPick'
        inf = open(infilename, 'rb')
        bias_estimate = cPickle.load(inf)
        inf.close()
    
        #frac2019[idate] = fracyear = doy/365.
    
        rmse_raw[idate] = np.sqrt(np.sum((forecast-analysis)**2)/(npts-1.))
        rmse_decayavg[idate] = np.sqrt(np.sum(((forecast-bias_decayavg)-analysis)**2)/(npts-1.))
        rmse_kf[idate] = np.sqrt(np.sum(((forecast-bias_estimate)-analysis)**2)/(npts-1.))

    rmse_raw_mean[ilead_idx] = np.mean(rmse_raw)
    rmse_decayavg_mean[ilead_idx]  = np.mean(rmse_decayavg)
    rmse_kf_mean[ilead_idx]  = np.mean(rmse_kf)

# ---- plot errors
    
f = plt.figure(figsize=(6.5,4))

ax = f.add_axes([.13,.13,.83,.79])
plt.title(r'RMSE of 2019 ECMWF T$_{2m}$ forecasts over CONUS',fontsize=14)
ax.plot([1,2,3,4], rmse_raw_mean, 'o-', color='Black',\
    lw=1.5, markersize=0.6, label='Raw', markerfacecolor='Red')
ax.plot([1,2,3,4], rmse_decayavg_mean, 'o-', color='RoyalBlue',\
    lw=1.5, markersize=0.6, label='Decaying average', markerfacecolor='RoyalBlue')
ax.plot([1,2,3,4], rmse_kf_mean, 'o-', color='Red',\
    lw=1.5, markersize=0.6, label='CONUS Kalman filter', markerfacecolor='Red')
ax.set_ylim(0,2.5)
ax.set_ylabel('RMSE (deg C)',fontsize=13)
ax.legend(loc=0)
plt.grid(True, lw=0.25)
ax.set_xlim(0,5)
ax.set_xticks([0,1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=13)

imagefile = 'rmse.pdf'
plt.savefig(imagefile)
print ('Plot done', imagefile)
    
            

