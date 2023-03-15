"""
plot_gain_2018.py

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
from read_reanalysis_timeseries import read_reanalysis_timeseries
from reformat_2d_to_4d_f90 import reformat_2d_to_4d_f90

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
cwarm = sys.argv[2] # warm or cold
cefold_random = sys.argv[3]
cefold_bias = sys.argv[4]
cday = sys.argv[5] # day past start of "season"
iday = int(cday) - 1
ilead = int(clead)

cpath_bias_est = '/Volumes/Backup Plus/python/bias_est/'
cpath_decay = '/Volumes/Backup Plus/ecmwf/biascorr/'

# --- read in the decaying average bias correction estimates.

savefile = cpath_decay+'2018_'+cwarm+'_beta_3d.cPick'
inf = open(savefile, 'rb')
beta_3d_decay = cPickle.load(inf)
inf.close()

# --- read in the Kalman filter bias estimates

beta_infile = cpath_bias_est + '2018_bias_est_flocal'+cefold_random+\
    '_blocal'+cefold_bias+'_2018_'+cwarm+'_lead'+\
    clead+'.cPick'
print (beta_infile)
inf = open(beta_infile,'rb')
beta_3d_KF = cPickle.load(inf)
print ('max, min beta_3d_KF = ', np.max(beta_3d_KF), np.min(beta_3d_KF))
lats = cPickle.load(inf)
lons =  cPickle.load(inf)
inf.close()

# --- now plot the localized bias error correlation for the selected point
   
#clevs = [-0.99,-0.9, -0.7, -0.5, -0.3, -0.1,0.1,0.3,0.5,0.7,0.9,0.99]
clevs = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000'] 
    
    
fig = plt.figure(figsize=(6.5,9.0))

axloc = [0.02,0.57,0.96,0.36]
ax1 = fig.add_axes(axloc)
title = r'(a) Decaying-avg. bias estimate for day '+cday+' of\n'+cwarm+\
    '-season 2018, '+clead+'-hour forecast'
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,beta_3d_decay[iday,:,:],clevs,\
    cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')


axloc = [0.02,0.12,0.96,0.36]
ax1 = fig.add_axes(axloc)
title = r'(b) Kalman-filter bias estimate for day '+cday+' of\n'+cwarm+\
    '-season 2018, '+clead+'-hour forecast'
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,beta_3d_KF[iday,:,:],clevs,\
    cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.02,0.07,0.96,0.02])
cbar = plt.colorbar(CS2, orientation='horizontal',\
    cax=cax, extend='both', ticks=clevs, format='%g')
cbar.ax.tick_params(labelsize=9)
cbar.set_label('Bias (deg C)')

# ---- set plot title

plot_title = 't2m_bias_estimate_day'+cday+'_'+cwarm+'2018_f'+\
    clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


