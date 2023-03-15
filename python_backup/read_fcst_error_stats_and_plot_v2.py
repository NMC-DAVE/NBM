""" read_fcst_error_stats_and_plot_v2.py

now re-done with Sherrie's verification against ERA5

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

ncfile = 'Temp_bias_rmse.nc'
ncdata  = Dataset(ncfile, 'r')
plevs = ncdata['plev'][:]
forecast_hours = ncdata['forecast_hours'][:]
regions = ncdata['regions'][:]
rmse_diff = ncdata['rmse_diff'][:,:,:]
print (np.shape(rmse_diff))
print (np.shape(rmse_diff[0,:,:]))
bias_diff = ncdata['bias_diff'][:,:,:]
GEFSv12_bias = ncdata['GEFSv12_bias'][:,:,:]
GEFSv12_rmse = ncdata['GEFSv12_rmse'][:,:,:]
CFSR_bias = ncdata['CFSR_bias'][:,:,:]
CFSR_rmse = ncdata['CFSR_rmse'][:,:,:]
ncdata.close()

print (plevs)
print (forecast_hours)
print (rmse_diff[0,:,:])

print ()
print (np.flipud(rmse_diff[0,:,:]))
#sys.exit()


# ----- make plots

fig = plt.figure(figsize=(6.5,7.))

fig.suptitle('2000-2019 temperature forecast RMSE difference and biases',fontsize=13)


hours = [6,24,48,72,96,120]
plevels = [1000, 925,850,800,750,700,600,500,400,300,250,200,150,100,50,10]
#clevels = [-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.01,1.5,2.0,2.5]
clevels = [-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000'] 
    
a1 = fig.add_axes([.10,.72,.24,.20])
a1.set_title('(a) NH RMSE diff (GEFSv12-CFSR)',fontsize=7)
cf = a1.contourf(hours,plevels, np.flipud(rmse_diff[0,:,:]),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,np.flipud(rmse_diff[0,:,:]),\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=7)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

a1 = fig.add_axes([.42,.72,.24,.20])
a1.set_title('(b) TR RMSE diff(GEFSv12-CFSR)',fontsize=7)
a1.contourf(hours,plevels,np.flipud(rmse_diff[1,:,:]),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,np.flipud(rmse_diff[1,:,:]),\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

a1 = fig.add_axes([.74,.72,.24,.20])
a1.set_title('(c) SH RMSE diff (GEFSv12-CFSR)',fontsize=7)
a1.contourf(hours,plevels,np.flipud(rmse_diff[2,:,:]),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,np.flipud(rmse_diff[2,:,:]),\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

cax = fig.add_axes([0.1,0.67,0.8,0.013])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('RMSE, GEFSv12-CFSR (deg C)',fontsize=6)
cb.ax.tick_params(labelsize=6)

# ----- GEFS bias

clevels = [-2,-1.5,-1,-0.5,-0.2,0.2,0.5,1,1.5,2]

a1 = fig.add_axes([.10,.38,.24,.2])
a1.set_title('(d) NH bias, GEFSv12 init',fontsize=7.5)
cf = a1.contourf(hours,plevels,np.flipud(GEFSv12_bias[0,:,:]),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,np.flipud(GEFSv12_bias[0,:,:]),\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=7)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

a1 = fig.add_axes([.42,.38,.24,.20])
a1.set_title('(e) TR bias, GEFSv12 initial',fontsize=7.5)
a1.contourf(hours,plevels,np.flipud(GEFSv12_bias[1,:,:]),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,np.flipud(GEFSv12_bias[1,:,:]),\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

a1 = fig.add_axes([.74,.38,.24,.20])
a1.set_title('(f) SH bias, GEFSv12 initial',fontsize=7.5)
a1.contourf(hours,plevels,np.flipud(GEFSv12_bias[2,:,:]),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,np.flipud(GEFSv12_bias[2,:,:]),\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])



# ----- Bias CFSR

a1 = fig.add_axes([.10,.12,.24,.20])
a1.set_title('(g) NH bias, CFSR initial',fontsize=7.5)
cf = a1.contourf(hours,plevels,np.flipud(CFSR_bias[0,:,:]),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,np.flipud(CFSR_bias[0,:,:]),\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=7)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])
a1.set_xlabel('Forecast lead time (h)', fontsize=7)

a1 = fig.add_axes([.42,.12,.24,.2])
a1.set_title('(h) TR bias, CFSR initial',fontsize=7.5)
a1.contourf(hours,plevels,np.flipud(CFSR_bias[1,:,:]),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,np.flipud(CFSR_bias[1,:,:]),\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])
a1.set_xlabel('Forecast lead time (h)', fontsize=7)

a1 = fig.add_axes([.74,.12,.24,.2])
a1.set_title('(i) SH bias, CFSR initial',fontsize=7.5)
a1.contourf(hours,plevels,np.flipud(CFSR_bias[2,:,:]),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,np.flipud(CFSR_bias[2,:,:]),\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])
a1.set_xlabel('Forecast lead time (h)', fontsize=7)

cax = fig.add_axes([0.1,0.05,0.8,0.013])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('Bias (deg C)',fontsize=6)
cb.ax.tick_params(labelsize=6)


plot_title = 'temperature_rmse_difference_2000-2019.png'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title,dpi=600)
print ('Plot done')



        
        
