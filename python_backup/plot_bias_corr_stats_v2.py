# plot_bias_corr_stats_v2.py 

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
import scipy.stats as stats
from paired_bootstrap import paired_bootstrap
from matplotlib import rcParams
rcParams['legend.fontsize']='small'

def read_bia_mae_rmse(cpath_errorstats, ftype, clead):
    statsfile = cpath_errorstats+ftype+clead+'h.txt'
    inf = open(statsfile, 'r')
    txtline = inf.read() 
    x = txtline.split()
    bia = float(x[0])
    mae = float(x[1])
    rmse = float(x[2])
    inf.close()
    return bia, mae, rmse
    
    
def read_bia_mae_rmse_bydate(cpath_errorstats, ftype, clead):
    statsfile = cpath_errorstats+ftype+clead+'h.cPick'
    inf = open(statsfile, 'rb')
    rmse_bydate = cPickle.load(inf)
    bia_bydate = cPickle.load(inf)
    mae_bydate = cPickle.load(inf)
    inf.close()
    return bia_bydate, mae_bydate, rmse_bydate

cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'
cleads = ['12','24','36','48','60','72','84','96','108','120']
nleads = len(cleads)
ndates = 365

rmse_raw = np.zeros((nleads), dtype=np.float32)
bia_raw = np.zeros((nleads), dtype=np.float32)
mae_raw = np.zeros((nleads), dtype=np.float32)

rmse_decay = np.zeros((nleads), dtype=np.float32)
bia_decay = np.zeros((nleads), dtype=np.float32)
mae_decay = np.zeros((nleads), dtype=np.float32)

rmse_qmap = np.zeros((nleads), dtype=np.float32)
bia_qmap = np.zeros((nleads), dtype=np.float32)
mae_qmap = np.zeros((nleads), dtype=np.float32)

rmse_KF = np.zeros((nleads), dtype=np.float32)
bia_KF = np.zeros((nleads), dtype=np.float32)
mae_KF = np.zeros((nleads), dtype=np.float32)

rmse_analog = np.zeros((nleads), dtype=np.float32)
bia_analog = np.zeros((nleads), dtype=np.float32)
mae_analog = np.zeros((nleads), dtype=np.float32)

rmse_MOS = np.zeros((nleads), dtype=np.float32)
bia_MOS = np.zeros((nleads), dtype=np.float32)
mae_MOS = np.zeros((nleads), dtype=np.float32)

rmse_mvMOS = np.zeros((nleads), dtype=np.float32)
bia_mvMOS = np.zeros((nleads), dtype=np.float32)
mae_mvMOS = np.zeros((nleads), dtype=np.float32)


rmse_raw_bydate = np.zeros((nleads, ndates), dtype=np.float32)
bia_raw_bydate = np.zeros((nleads, ndates), dtype=np.float32)
mae_raw_bydate = np.zeros((nleads, ndates), dtype=np.float32)

rmse_decay_bydate = np.zeros((nleads, ndates), dtype=np.float32)
bia_decay_bydate = np.zeros((nleads, ndates), dtype=np.float32)
mae_decay_bydate = np.zeros((nleads, ndates), dtype=np.float32)

rmse_qmap_bydate = np.zeros((nleads, ndates), dtype=np.float32)
bia_qmap_bydate = np.zeros((nleads, ndates), dtype=np.float32)
mae_qmap_bydate = np.zeros((nleads, ndates), dtype=np.float32)

rmse_KF_bydate = np.zeros((nleads, ndates), dtype=np.float32)
bia_KF_bydate = np.zeros((nleads, ndates), dtype=np.float32)
mae_KF_bydate = np.zeros((nleads, ndates), dtype=np.float32)

rmse_analog_bydate = np.zeros((nleads, ndates), dtype=np.float32)
bia_analog_bydate = np.zeros((nleads, ndates), dtype=np.float32)
mae_analog_bydate = np.zeros((nleads, ndates), dtype=np.float32)

rmse_MOS_bydate = np.zeros((nleads, ndates), dtype=np.float32)
bia_MOS_bydate = np.zeros((nleads, ndates), dtype=np.float32)
mae_MOS_bydate = np.zeros((nleads, ndates), dtype=np.float32)

rmse_mvMOS_bydate = np.zeros((nleads, ndates), dtype=np.float32)
bia_mvMOS_bydate = np.zeros((nleads, ndates), dtype=np.float32)
mae_mvMOS_bydate = np.zeros((nleads, ndates), dtype=np.float32)

q05_decay = np.zeros((nleads), dtype=np.float32)
q05_qmap = np.zeros((nleads), dtype=np.float32)
q05_KF = np.zeros((nleads), dtype=np.float32)
q05_analog = np.zeros((nleads), dtype=np.float32)
q05_MOS = np.zeros((nleads), dtype=np.float32)
q05_mvMOS = np.zeros((nleads), dtype=np.float32)

q95_decay = np.zeros((nleads), dtype=np.float32)
q95_qmap = np.zeros((nleads), dtype=np.float32)
q95_KF = np.zeros((nleads), dtype=np.float32)
q95_analog = np.zeros((nleads), dtype=np.float32)
q95_MOS = np.zeros((nleads), dtype=np.float32)
q95_mvMOS = np.zeros((nleads), dtype=np.float32)

for ilead, clead in enumerate(cleads):

    ftype = 'raw_forecast_errorstats_2019_'
    bia, mae, rmse = read_bia_mae_rmse(cpath_errorstats, ftype, clead)
    bia_raw[ilead] = bia
    mae_raw[ilead] = mae
    rmse_raw[ilead] = rmse
    
    bia_bydate0, mae_bydate0, rmse_bydate0 = read_bia_mae_rmse_bydate(cpath_errorstats, ftype, clead)
    bia_raw_bydate[ilead,:] = bia_bydate0[:]
    mae_raw_bydate[ilead,:] = mae_bydate0[:]
    rmse_raw_bydate[ilead,:] = rmse_bydate0[:]
    
    ftype = 'decayavg_forecast_errorstats_2019_'
    bia, mae, rmse = read_bia_mae_rmse(cpath_errorstats, ftype, clead)
    bia_decay[ilead] = bia
    mae_decay[ilead] = mae
    rmse_decay[ilead] = rmse
    
    bia_bydate, mae_bydate, rmse_bydate = read_bia_mae_rmse_bydate(cpath_errorstats, ftype, clead)
    bia_decay_bydate[ilead,:] = bia_bydate[:]
    mae_decay_bydate[ilead,:] = mae_bydate[:]
    rmse_decay_bydate[ilead,:] = rmse_bydate[:]
    
    ftype = 'quantile_mapping_errorstats_2019_'
    bia, mae, rmse = read_bia_mae_rmse(cpath_errorstats, ftype, clead)
    bia_qmap[ilead] = bia
    mae_qmap[ilead] = mae
    rmse_qmap[ilead] = rmse
    
    bia_bydate2, mae_bydate2, rmse_bydate2 = read_bia_mae_rmse_bydate(cpath_errorstats, ftype, clead)
    bia_qmap_bydate[ilead,:] = bia_bydate2[:]
    mae_qmap_bydate[ilead,:] = mae_bydate2[:]
    rmse_qmap_bydate[ilead,:] = rmse_bydate2[:]
    
    q05_qmap[ilead], q95_qmap[ilead] = paired_bootstrap(rmse_bydate, rmse_bydate2)
    
    ftype = 'KF_corr_GEFSv12_2019_'
    bia, mae, rmse = read_bia_mae_rmse(cpath_errorstats, ftype, clead)
    bia_KF[ilead] = bia
    mae_KF[ilead] = mae
    rmse_KF[ilead] = rmse
    
    bia_bydate2, mae_bydate2, rmse_bydate2 = read_bia_mae_rmse_bydate(cpath_errorstats, ftype, clead)
    bia_KF_bydate[ilead,:] = bia_bydate2[:]
    mae_KF_bydate[ilead,:] = mae_bydate2[:]
    rmse_KF_bydate[ilead,:] = rmse_bydate2[:]
    
    ftype = 'analog_forecast_errorstats_2019_'
    bia, mae, rmse = read_bia_mae_rmse(cpath_errorstats, ftype, clead)
    bia_analog[ilead] = bia
    mae_analog[ilead] = mae
    rmse_analog[ilead] = rmse
    
    bia_bydate2, mae_bydate2, rmse_bydate2 = read_bia_mae_rmse_bydate(cpath_errorstats, ftype, clead)
    bia_analog_bydate[ilead,:] = bia_bydate2[:]
    mae_analog_bydate[ilead,:] = mae_bydate2[:]
    rmse_analog_bydate[ilead,:] = rmse_bydate2[:]    
    q05_analog[ilead], q95_analog[ilead] = paired_bootstrap(rmse_bydate, rmse_bydate2)
    
    ftype = 'MOS_GEFSv12_2019_'
    bia, mae, rmse = read_bia_mae_rmse(cpath_errorstats, ftype, clead)
    bia_MOS[ilead] = bia
    mae_MOS[ilead] = mae
    rmse_MOS[ilead] = rmse
    
    bia_bydate2, mae_bydate2, rmse_bydate2 = read_bia_mae_rmse_bydate(cpath_errorstats, ftype, clead)
    bia_MOS_bydate[ilead,:] = bia_bydate2[:]
    mae_MOS_bydate[ilead,:] = mae_bydate2[:]
    rmse_MOS_bydate[ilead,:] = rmse_bydate2[:]
    
    q05_MOS[ilead], q95_MOS[ilead] = paired_bootstrap(rmse_bydate, rmse_bydate2)
    
    ftype = 'MOS_mvr_GEFSv12_2019_'
    bia, mae, rmse = read_bia_mae_rmse(cpath_errorstats, ftype, clead)
    bia_mvMOS[ilead] = bia
    mae_mvMOS[ilead] = mae
    rmse_mvMOS[ilead] = rmse
    
    bia_bydate2, mae_bydate2, rmse_bydate2 = read_bia_mae_rmse_bydate(cpath_errorstats, ftype, clead)
    bia_mvMOS_bydate[ilead,:] = bia_bydate2[:]
    mae_mvMOS_bydate[ilead,:] = mae_bydate2[:]
    rmse_mvMOS_bydate[ilead,:] = rmse_bydate2[:]
    
    q05_mvMOS[ilead], q95_mvMOS[ilead] = paired_bootstrap(rmse_bydate, rmse_bydate2)

    q05_decay[ilead], q95_decay[ilead] = paired_bootstrap(rmse_bydate, rmse_bydate0)

# ---- plot data
    
f = plt.figure(figsize=(9.,4.5))

leads = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
ax = f.add_axes([.07,.12,.23,.78])
ax.set_title(r'(a) RMSEs',fontsize=13)

yup = rmse_decay + q95_decay
ydn = rmse_decay + q05_decay
print ('yup = ', yup)
print ('ydn = ', ydn)
ax.fill_between(leads, yup, ydn, color='NavajoWhite',zorder=0,alpha=0.6)

ax.plot(leads, rmse_raw, 'o-', color='Red',\
    lw=2, markersize=4, label='Raw', markerfacecolor='Red')
ax.plot(leads, rmse_decay,'o', linestyle='--',color='Peru',\
    lw=2, markersize=4, label='DAV', markerfacecolor='Peru')
ax.set_ylabel('RMSE (deg C)',fontsize=12)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0.4,5.1)
ax.set_ylim(0.0, 4.0)
ax.set_xticks([1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=12)
    
ax = f.add_axes([.4,.12,.23,.78])
ax.set_title(r'(b) RMSE differences',fontsize=13)

yup = (rmse_qmap - rmse_decay) + q05_qmap
ydn = (rmse_qmap - rmse_decay) + q95_qmap
ax.fill_between(leads, yup, ydn, color='PowderBlue',zorder=0,alpha=0.6)

yup = (rmse_MOS - rmse_decay) + q05_MOS
ydn = (rmse_MOS - rmse_decay) + q95_MOS
ax.fill_between(leads, yup, ydn, color='PaleGreen',zorder=1,alpha=0.6)

yup = (rmse_mvMOS - rmse_decay) + q05_mvMOS
ydn = (rmse_mvMOS - rmse_decay) + q95_mvMOS
ax.fill_between(leads, yup, ydn, color='Gainsboro',zorder=2,alpha=0.6)

yup = (rmse_analog - rmse_decay) + q05_analog
ydn = (rmse_analog - rmse_decay) + q95_analog
ax.fill_between(leads, yup, ydn, color='Plum',zorder=2,alpha=0.6)

ax.plot(leads, rmse_qmap - rmse_decay, 'o-', color='Blue',zorder=4,\
    lw=2, markersize=4, label='QM - DAV', markerfacecolor='Blue')
ax.plot(leads, rmse_MOS - rmse_decay, 'o-', color='Green',zorder=5,\
    lw=2, markersize=4, label='uMOS - DAV', markerfacecolor='Green')
ax.plot(leads, rmse_mvMOS - rmse_decay, 'o', linestyle='--', color='Gray',\
    lw=2, markersize=4, label='mvMOS - DAV', \
    zorder=6,markerfacecolor='Gray')
ax.plot(leads, rmse_analog - rmse_decay, 'o', linestyle='--', color='DarkMagenta',\
    lw=2, markersize=4, label='Analog - DAV', \
    zorder=7,markerfacecolor='DarkMagenta')
    
ax.plot([0.4,5.1],[0.0, 0.0],lw=2,color='Black',zorder=0)
ax.set_ylabel('RMSE difference (deg C)',fontsize=12)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0.4,5.1)
ax.set_ylim(-0.4, 0.15)
ax.set_xticks([1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=12)

ax = f.add_axes([.73,.12,.23,.78])
ax.set_title(r'(c) Unconditional biases',fontsize=13)
ax.plot(leads, bia_raw, 'o-', color='Red',\
    lw=2, markersize=4, label='Raw', markerfacecolor='Red')
ax.plot(leads, bia_decay, 'o-', color='Peru',\
    lw=2, markersize=4, label='DAV', markerfacecolor='Peru')
ax.plot(leads, bia_analog, 'o', linestyle='--', color='DarkMagenta',\
    lw=2, markersize=4, label='Analog', markerfacecolor='DarkMagenta')
ax.plot(leads, bia_qmap, 'o-', color='Blue',\
    lw=2, markersize=4, label='QM', markerfacecolor='Blue')
ax.plot(leads, bia_MOS,'o-', color='Green',\
    lw=2, markersize=4, label='uMOS', markerfacecolor='Green')
ax.plot(leads, bia_mvMOS,'o', linestyle='--', color='Gray',\
    lw=2, markersize=4, label='mvMOS', markerfacecolor='Gray')
ax.plot([0.4,5.1],[0,0],color='Black',lw=2,zorder=0)
ax.set_ylabel('Bias (deg C)',fontsize=12)
ax.legend(loc=9)
ax.grid(True, lw=0.25)
ax.set_xlim(0.4,5.1)
ax.set_ylim(-0.2, 1.5)
ax.set_xticks([1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=12)

imagefile = 'rmse_bias_biascorrs.png' 
plt.savefig(imagefile, dpi=400)
print ('Plot done', imagefile)
plt.close()








# ---- plot data
    
f = plt.figure(figsize=(9.,4.5))

leads = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
ax = f.add_axes([.07,.12,.23,.78])
ax.set_title(r'(a) MAEs',fontsize=13)

#yup = rmse_decay + q95_decay
#ydn = rmse_decay + q05_decay
#print ('yup = ', yup)
#print ('ydn = ', ydn)
#ax.fill_between(leads, yup, ydn, color='NavajoWhite',zorder=0,alpha=0.6)

ax.plot(leads, mae_raw, 'o-', color='Red',\
    lw=2, markersize=4, label='Raw', markerfacecolor='Red')
ax.plot(leads, mae_decay,'o', linestyle='--',color='Peru',\
    lw=2, markersize=4, label='Decaying average', markerfacecolor='Peru')
ax.set_ylabel('MAE (deg C)',fontsize=12)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0.4,5.1)
ax.set_ylim(0.0, 4.0)
ax.set_xticks([1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=12)
    
ax = f.add_axes([.4,.12,.23,.78])
ax.set_title(r'(b) MAE differences',fontsize=13)

#yup = (rmse_qmap - rmse_decay) + q05_qmap
#ydn = (rmse_qmap - rmse_decay) + q95_qmap
#ax.fill_between(leads, yup, ydn, color='PowderBlue',zorder=0,alpha=0.6)

#yup = (rmse_MOS - rmse_decay) + q05_MOS
#ydn = (rmse_MOS - rmse_decay) + q95_MOS
#ax.fill_between(leads, yup, ydn, color='PaleGreen',zorder=1,alpha=0.6)

#yup = (rmse_mvMOS - rmse_decay) + q05_mvMOS
#ydn = (rmse_mvMOS - rmse_decay) + q95_mvMOS
#ax.fill_between(leads, yup, ydn, color='Gainsboro',zorder=2,alpha=0.6)

#yup = (rmse_KF - rmse_decay) + q05_KF
#ydn = (rmse_KF - rmse_decay) + q95_KF
#ax.fill_between(leads, yup, ydn, color='Plum',zorder=2,alpha=0.6)

ax.plot(leads, mae_qmap - mae_decay, 'o-', color='Blue',zorder=4,\
    lw=2, markersize=4, label='Quantile mapping -\ndecaying average', markerfacecolor='Blue')
ax.plot(leads, mae_MOS - mae_decay, 'o-', color='Green',zorder=5,\
    lw=2, markersize=4, label='MOS -\ndecaying average', markerfacecolor='Green')
ax.plot(leads, mae_mvMOS - mae_decay, 'o', linestyle='--', color='Gray',\
    lw=2, markersize=4, label='Multi-variate MOS -\ndecaying average', \
    zorder=6,markerfacecolor='Gray')
#ax.plot(leads, mae_KF - mae_decay, 'o', linestyle='--', color='DarkMagenta',\
#    lw=2, markersize=4, label='Kalman filter -\ndecaying average', \
#    zorder=7,markerfacecolor='DarkMagenta')
    
ax.plot([0.4,5.1],[0.0, 0.0],lw=2,color='Black',zorder=0)
ax.set_ylabel('RMSE difference (deg C)',fontsize=12)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0.4,5.1)
ax.set_ylim(-0.4, 0.15)
ax.set_xticks([1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=12)

ax = f.add_axes([.73,.12,.23,.78])
ax.set_title(r'(c) Unconditional biases',fontsize=13)
ax.plot(leads, bia_raw, 'o-', color='Red',\
    lw=2, markersize=4, label='Raw', markerfacecolor='Red')
ax.plot(leads, bia_decay, 'o-', color='Peru',\
    lw=2, markersize=4, label='Decaying average', markerfacecolor='Peru')
#ax.plot(leads, bia_KF, 'o', linestyle='--', color='DarkMagenta',\
#    lw=2, markersize=4, label='Kalman filter', markerfacecolor='DarkMagenta')
ax.plot(leads, bia_qmap, 'o-', color='Blue',\
    lw=2, markersize=4, label='Quantile mapping', markerfacecolor='Blue')
ax.plot(leads, bia_MOS,'o-', color='Green',\
    lw=2, markersize=4, label='MOS', markerfacecolor='Green')
ax.plot(leads, bia_mvMOS,'o', linestyle='--', color='Gray',\
    lw=2, markersize=4, label='Multi-variate MOS', markerfacecolor='Gray')
ax.plot([0.4,5.1],[0,0],color='Black',lw=2,zorder=0)
ax.set_ylabel('Bias (deg C)',fontsize=12)
ax.legend(loc=9)
ax.grid(True, lw=0.25)
ax.set_xlim(0.4,5.1)
ax.set_ylim(-0.2, 1.5)
ax.set_xticks([1,2,3,4,5])
ax.set_xlabel('Forecast lead (days)', fontsize=12)

imagefile = 'mae_bias_biascorrs.png' 
plt.savefig(imagefile, dpi=400)
print ('Plot done', imagefile)
plt.close()



