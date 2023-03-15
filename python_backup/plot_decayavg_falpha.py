# plot_decayavg_falpha.py

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
rcParams['legend.fontsize']='x-small'

def read_bia_mae_rmse(statsfile):
    inf = open(statsfile, 'r')
    txtline = inf.read() 
    x = txtline.split()
    bia = float(x[0])
    mae = float(x[1])
    rmse = float(x[2])
    inf.close()
    return bia, mae, rmse
    

def set_optimal_alpha(clead):
    if clead == '12':
        oalpha = 0.16
        ialpha = 7
    elif clead == '24':
        oalpha = 0.14  
        ialpha = 6
    elif clead == '36':
        oalpha = 0.08
        ialpha = 3
    elif clead == '48':
        oalpha = 0.06 
        ialpha = 2
    elif clead == '60':
        oalpha = 0.06 
        ialpha = 2
    elif clead == '72':
        oalpha = 0.06 
        ialpha = 2
    elif clead == '84':
        oalpha = 0.06
        ialpha = 2
    elif clead == '96':
        oalpha = 0.04 
        ialpha = 1
    elif clead == '108':
        oalpha = 0.04 
        ialpha = 1
    elif clead == '120':
        oalpha = 0.04 
        ialpha = 1 
    elif clead == '132':
        oalpha = 0.04
        ialpha = 1
    return oalpha, ialpha    
    


cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'
cleads = ['12','24','36','48','60','72','84','96','108','120']
cleads_r = ['120','108','96','84','72','60','48','36','24','12']
calphas = ['0.02','0.04','0.06','0.08','0.1','0.12','0.14','0.16','0.18','0.2']
alphas = np.arange(0.02,0.201,0.02)
nalpha = len(calphas)
nleads = len(cleads)

rmse_decay = np.zeros((nleads,nalpha), dtype=np.float32)
bia_decay = np.zeros((nleads,nalpha), dtype=np.float32)
mae_decay = np.zeros((nleads,nalpha), dtype=np.float32)

for ilead, clead in enumerate(cleads):
    for ialpha, calpha in enumerate(calphas):

        statsfile = cpath_errorstats + \
            'decayavg_2000_2018_alpha'+calpha+'_lead'+clead+'.txt'
        
        bia, mae, rmse = read_bia_mae_rmse(statsfile)
        bia_decay[ilead,ialpha] = bia
        mae_decay[ilead,ialpha] = mae
        rmse_decay[ilead,ialpha] = rmse
    

# ---- plot data
    
f = plt.figure(figsize=(3.1,4.5))

leads = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
ileads = [0,1,2,3,4,5,6,7,8,9]
ileads_r = ileads.reverse()
ax = f.add_axes([.2,.12,.77,.78])
ax.set_title(r'Decaying-average RMSE',fontsize=13)

for clead in cleads_r:

    if clead == '12':
        color = 'Black'
        markerfacecolor = 'Black'
        linestyle = '--'
        label = 'Day +0.5'
        ilead = 0
    elif clead == '24':
        color = 'Black'
        markerfacecolor = 'Black'
        linestyle = '-'
        label = 'Day +1'
        ilead = 1
    elif clead == '36':
        color = 'Gray'
        markerfacecolor = 'Gray'
        linestyle = '--'
        label = 'Day +1.5'
        ilead = 2
    elif clead == '48':
        color = 'Gray'
        markerfacecolor = 'Gray'
        linestyle = '-'
        label = 'Day +2'
        ilead = 3
    elif clead == '60':
        color = 'Red'
        markerfacecolor = 'Red'
        linestyle = '--'
        label = 'Day +2.5'
        ilead = 4
    elif clead == '72':
        color = 'Red'
        markerfacecolor = 'Red'
        linestyle = '-'
        label = 'Day +3'
        ilead = 5
    elif clead == '84':
        color = 'RoyalBlue'
        markerfacecolor = 'RoyalBlue'
        linestyle = '--'
        label = 'Day +3.5'
        ilead = 6
    elif clead == '96':
        color = 'RoyalBlue'
        markerfacecolor = 'RoyalBlue'
        linestyle = '-'
        label = 'Day +4'
        ilead = 7
    elif clead == '108':
        color = 'Peru'
        markerfacecolor = 'Peru'
        linestyle = '--'
        label = 'Day +4.5'
        ilead = 8
    elif clead == '120':
        color = 'Peru'
        markerfacecolor = 'Peru'
        linestyle = '-'
        label = 'Day +5'
        ilead = 9

    ax.plot(alphas, rmse_decay[ilead,:], 'o-', color=color, linestyle = linestyle, \
        lw=2, markersize=4, label=label, markerfacecolor=markerfacecolor)
        
        
    oalpha,ialpha = set_optimal_alpha(clead)
    ax.plot([oalpha,oalpha], [rmse_decay[ilead,ialpha],rmse_decay[ilead,ialpha]], 'o', color=color, \
        markersize=6, markerfacecolor=markerfacecolor)

ax.set_ylabel('RMSE (deg C)',fontsize=12)
ax.legend(loc=0)
ax.grid(True, lw=0.25)
ax.set_xlim(0.01,0.21)
ax.set_ylim(0.0, 3.5)
ax.set_xticks([0.04,0.08,0.12,0.16,0.20])
ax.set_xlabel(r'$\alpha$', fontsize=12)
    
imagefile = 'rmse_decayavg_falpha.png' 
plt.savefig(imagefile, dpi=400)
print ('Plot done', imagefile)
plt.close()






