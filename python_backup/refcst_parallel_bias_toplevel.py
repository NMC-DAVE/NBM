""" refcst_parallel_bias.py

"""


import csv
import sys
import _pickle as cPickle
import numpy as np
from datetime import datetime
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams 
from netCDF4 import Dataset
rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'

# ---- read data from Sherrie Fredrick's netCDF file.

ncfile = 'lat_time_all.nc'
ncdata  = Dataset(ncfile, 'r')
plevs = ncdata['plev'][:]/100.
forecast_hours = ncdata['forecast_hours'][:]
lats = ncdata['lat'][:]
RH = ncdata['RH'][:,:]
T = ncdata['T'][:,:]
U = ncdata['U'][:,:]
V = ncdata['V'][:,:]
ncdata.close()
print (plevs)
print (np.shape(plevs), np.shape(lats), np.shape(T))
print (np.max(T), np.min(T))

#print (forecast_hours)
#sys.exit()

T_bias_top = T[:,0,:]  # fhr, lat
T_bias_250 = T[:,4,:]  # fhr, lat

# ----- make plots

fig = plt.figure(figsize=(6.5,7.))

a1 = fig.add_axes([.16,0.17,.8,0.78])
a1.set_title('Differences, 10 hPa T, GEFSv12 pre-production parallel minus reforecast',fontsize=11)
clevels = [-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
print (np.shape(forecast_hours),np.shape(lats), np.shape(T_bias_top))
cf = a1.contourf(forecast_hours/24,lats, np.transpose(T_bias_top),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours/24,lats, np.transpose(T_bias_top),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(-90,90)
a1.set_ylabel('Latitude (degrees)', fontsize=11)
a1.set_xlim(0,10.)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=11)

cax = fig.add_axes([0.16,0.06,0.8,0.03])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('Bias (deg C)',fontsize=8)
cb.ax.tick_params(labelsize=8)

plot_title = 't-bias-top.png'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title,dpi=600)
print ('Plot done')


# ----- make plots

fig = plt.figure(figsize=(6.5,7.))

a1 = fig.add_axes([.16,0.17,.8,0.78])
a1.set_title('Differences, 250 hPa T, GEFSv12 pre-production parallel minus reforecast',fontsize=10)
clevels = [-1.5,-1.0,-0.7,-0.5,-0.2,0.2,0.5,0.7,1.0,1.5]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
print (np.shape(forecast_hours),np.shape(lats), np.shape(T_bias_250))
cf = a1.contourf(forecast_hours/24,lats, np.transpose(T_bias_250),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours/24,lats, np.transpose(T_bias_250),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(-90,90)
a1.set_ylabel('Latitude (degrees)', fontsize=11)
a1.set_xlim(0,10.)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=11)

cax = fig.add_axes([0.16,0.06,0.8,0.03])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('Bias (deg C)',fontsize=8)
cb.ax.tick_params(labelsize=8)

plot_title = 't-bias-250.png'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title,dpi=600)
print ('Plot done')

        
        
