""" refcst_parallel_rmse_difference.py

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

ncfile = 'zonal_UVTRH.nc'
ncdata  = Dataset(ncfile, 'r')
plevs = ncdata['plev'][:]/100.
forecast_hours = ncdata['forecast_hours'][:]
RH = ncdata['RH'][:,:]
T = ncdata['T'][:,:]
U = ncdata['U'][:,:]
V = ncdata['V'][:,:]
ncdata.close()
print (plevs)
print (np.shape(plevs), np.shape(forecast_hours), np.shape(T))
print (np.max(RH))
print (np.min(U), np.max(U))
print (np.min(V), np.max(V))
#print (forecast_hours)
#sys.exit()

# ----- make plots

fig = plt.figure(figsize=(6.5,7.))

# ----- Bias CFSR

a1 = fig.add_axes([.09,.61,.4,.3])
a1.set_title('RMSE temperature differences, GEFSv12\npre-production parallel minus reforecast',fontsize=11)
clevels = [-1.0,-0.7,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.7,1.0]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
cf = a1.contourf(forecast_hours/24,plevs,np.transpose(T),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours/24,plevs,np.transpose(T),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.09,0.54,0.4,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('Bias (deg C)',fontsize=8)
cb.ax.tick_params(labelsize=8)



plot_title = 'retro_refcst_bias_3panel.png'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title,dpi=600)
print ('Plot done')


        
        
