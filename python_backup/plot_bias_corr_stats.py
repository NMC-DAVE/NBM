#!/usr/local/bin/python2.7
import numpy as np
import numpy.ma as ma
import dateutils
from netCDF4 import Dataset
from datetime import datetime, date, time, timedelta
import matplotlib as mpl
mpl.use('Agg') #for web 
from dateutils import datetohrs
import matplotlib.pyplot as plt
import sys
import _pickle as cPickle

clead = sys.argv[1]


rmse_raw = np.array([1.931, 2.062, 2.224, 2.436, 2.698])
rmse_decay = np.array([1.417, 1.641, 1.862, 2.120, 2.417])
rmse_kf = np.array([1.410, 1.628, 1.864, 2.125, 2.417])
rmse_qmap = np.array([1.485, 1.657, 1.861, 2.109, 2.387])
rmse_both = np.array([1.400, 1.605, 1.828, 2.088, 2.384])

bia_raw = np.array([0.690, 0.750, 0.755, 0.757, 0.769])
bia_decay = np.array([0.010, 0.015, 0.021, 0.025, 0.038])
bia_kf = np.array([0.058, 0.054, 0.048, 0.079, 0.108])
bia_qmap = np.array([0.114, 0.160, 0.160, 0.154, 0.151])
bia_both = np.array([0.062, 0.088, 0.092, 0.091, 0.097])


infile = '/Volumes/Backup Plus/python/fcst_stats/raw_forecast_errorstats_2019_'+clead+'h.cPick'
inf = open(infile,'rb')
rmse_raw_timeseries = cPickle.load(inf)
bia_raw_timeseries = cPickle.load(inf)
inf.close()

infile = '/Volumes/Backup Plus/python/fcst_stats/decayavg_forecast_errorstats_2019_'+clead+'h.cPick'
inf = open(infile,'rb')
rmse_decay_timeseries = cPickle.load(inf)
bia_decay_timeseries = cPickle.load(inf)
inf.close()

infile = '/Volumes/Backup Plus/python/fcst_stats/KF_corr_GEFSv12_2019_'+clead+'h.cPick'
inf = open(infile,'rb')
rmse_KF_timeseries = cPickle.load(inf)
bia_KF_timeseries = cPickle.load(inf)
inf.close()

infile = '/Volumes/Backup Plus/python/fcst_stats/quantile_mapping_errorstats_2019_'+clead+'h.cPick'
inf = open(infile,'rb')
rmse_qmap_timeseries = cPickle.load(inf)
bia_qmap_timeseries = cPickle.load(inf)
inf.close()

infile = '/Volumes/Backup Plus/python/fcst_stats/dual_errorstats_2019_'+clead+'h.cPick'
inf = open(infile,'rb')
rmse_both_timeseries = cPickle.load(inf)
bia_both_timeseries = cPickle.load(inf)
inf.close()



# ---- plot data
    
f = plt.figure(figsize=(9.,4.5))

ax = f.add_axes([.07,.12,.23,.78])
ax.set_title(r'(a) RMSEs',fontsize=15)
ax.plot([1,2,3,4,5], rmse_raw, 'o-', color='Red',\
    lw=2, markersize=4, label='Raw', markerfacecolor='Red')
ax.plot([1,2,3,4,5], rmse_decay,'o-', color='Blue',\
    lw=2, markersize=4, label='Decaying average', markerfacecolor='Blue')
ax.set_ylabel('RMSE (deg C)',fontsize=13)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0.5,5.5)
ax.set_ylim(0.0, 4.0)
ax.set_xticks([1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=13)
    
ax = f.add_axes([.43,.12,.23,.78])
ax.set_title(r'(b) RMSE differences',fontsize=15)
ax.plot([1,2,3,4,5], rmse_kf - rmse_decay,'o-', color='Red',\
    lw=2, markersize=4, label='Kalman filter -\ndecaying average', markerfacecolor='Red')
ax.plot([1,2,3,4,5], rmse_qmap - rmse_decay, 'o-', color='Blue',\
    lw=2, markersize=4, label='Quantile mapping -\ndecaying average', markerfacecolor='Blue')
ax.plot([1,2,3,4,5], rmse_both - rmse_decay, 'o-', color='Green',\
    lw=2, markersize=4, label='Mixture -\ndecaying average', markerfacecolor='Green')
ax.plot([0.5,5.5],[0.0, 0.0],lw=2,color='Black')
ax.set_ylabel('RMSE difference (deg C)',fontsize=13)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0.5,5.5)
ax.set_ylim(-0.1, 0.1)
ax.set_xticks([1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=13)


ax = f.add_axes([.76,.12,.23,.78])
ax.set_title(r'(c) Biases',fontsize=15)
ax.plot([1,2,3,4,5], bia_raw, 'o-', color='Black',\
    lw=2, markersize=4, label='Raw', markerfacecolor='Black')
ax.plot([1,2,3,4,5], bia_decay, 'o-', color='Gray',\
    lw=2, markersize=4, label='Decaying average', markerfacecolor='Gray')
ax.plot([1,2,3,4,5], bia_kf, 'o-', color='Red',\
    lw=2, markersize=4, label='Kalman filter', markerfacecolor='Red')
ax.plot([1,2,3,4,5], bia_qmap, 'o-', color='Blue',\
    lw=2, markersize=4, label='Quantile mapping', markerfacecolor='Blue')
ax.plot([1,2,3,4,5], bia_both,'o-', color='Green',\
    lw=2, markersize=4, label='Mixture', markerfacecolor='Green')
ax.set_ylabel('Bias (deg C)',fontsize=13)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0.5,5.5)
ax.set_ylim(0.0, 1.0)
ax.set_xticks([1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=13)

imagefile = 'rmse_bias_biascorrs.pdf' 
plt.savefig(imagefile, dpi=400)
print ('Plot done', imagefile)
plt.close()



# ---- plot time series of errors data

ndates = len(rmse_raw_timeseries)
ibracket = int(clead) // 24
    
f = plt.figure(figsize=(9.,4.5))

ax = f.add_axes([.1,.15,.85,.78])
ax.set_title('Daily RMSEs, lead = '+clead+' h',fontsize=15)
ax.plot(range(ibracket,ndates-ibracket), rmse_KF_timeseries[ibracket:ndates-ibracket], 'o-', color='Red',\
    lw=2, markersize=2.5, label='Kalman filter', markerfacecolor='Red')
ax.plot(range(ibracket,ndates-ibracket), rmse_decay_timeseries[ibracket:ndates-ibracket],'o-', color='Blue',\
    lw=2, markersize=2.5, label='Decaying average', markerfacecolor='Blue')
ax.set_ylabel('RMSE (deg C)',fontsize=13)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0,366)
ax.set_ylim(0.0, 6.0)
ax.set_xticks([0,30,60,90,120,150,180,210,240,270,300,330,360])
ax.set_xlabel('Julian day of the year', fontsize=13)
    
imagefile = 'rmse_KF_decay_timeseries_f'+clead+'.pdf' 
plt.savefig(imagefile, dpi=400)
print ('Plot done', imagefile)
plt.close()


# ---- plot time series of errors data

ndates = len(rmse_raw_timeseries)
    
f = plt.figure(figsize=(9.,4.5))

ax = f.add_axes([.1,.15,.85,.78])
ax.set_title('Daily RMSEs, lead = '+clead+' h',fontsize=15)
ax.plot(range(ibracket,ndates-ibracket), rmse_qmap_timeseries[ibracket:ndates-ibracket], 'o-', color='Red',\
    lw=2, markersize=2.5, label='Quantile mapping', markerfacecolor='Red')
ax.plot(range(ibracket,ndates-ibracket), rmse_decay_timeseries[ibracket:ndates-ibracket],'o-', color='Blue',\
    lw=2, markersize=2.5, label='Decaying average', markerfacecolor='Blue')
ax.set_ylabel('RMSE (deg C)',fontsize=13)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0,366)
ax.set_ylim(0.0, 6.0)
ax.set_xticks([0,30,60,90,120,150,180,210,240,270,300,330,360])
ax.set_xlabel('Julian day of the year', fontsize=13)
    
imagefile = 'rmse_qmap_decay_timeseries_f'+clead+'.pdf' 
plt.savefig(imagefile, dpi=400)
print ('Plot done', imagefile)
plt.close()

