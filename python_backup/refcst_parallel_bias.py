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
rcParams['legend.fontsize']='xx-small'

# ---- read data from Sherrie Fredrick's netCDF file 
#      the differences between retro and reforecast

ncfile = 'zonal_UVTRH.nc'
ncdata  = Dataset(ncfile, 'r')
plevs = ncdata['plev'][:]/100.
#print (plevs)
#sys.exit()
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
print (forecast_hours)

# ----- get the temp, U, and RH RMSEs from retro and reforecast

ncfile = 'Temp_RETRO_rmse.nc'
ncdata  = Dataset(ncfile, 'r')
T_retro_rmse = ncdata['Trmse_pt'][:,:] # forcast hours, plev
forecast_hours2 = np.arange(0.25,10.01,0.25) #ncdata['forecast_hours'][:]
ncdata.close()

ncfile = 'Temp_Reforecast_rmse.nc'
ncdata  = Dataset(ncfile, 'r')
T_reforecast_rmse = ncdata['Trmse_pt'][:,:] # forcast hours, plev
ncdata.close()

ncfile = 'U_RETRO_rmse.nc'
ncdata  = Dataset(ncfile, 'r')
U_retro_rmse = ncdata['Urmse_pt'][:,:] # forcast hours, plev
print ('np.shape(U_retro_rmse)', np.shape(U_retro_rmse))
ncdata.close()

ncfile = 'U_Reforecast_rmse.nc'
ncdata  = Dataset(ncfile, 'r')
U_reforecast_rmse = ncdata['Urmse_pt'][:,:] # forcast hours, plev
ncdata.close()

ncfile = 'RH_RETRO_rmse.nc'
ncdata  = Dataset(ncfile, 'r')
RH_retro_rmse = ncdata['RHrmse_pt'][:,:] # forcast hours, plev
print ('np.shape(RH_retro_rmse)', np.shape(RH_retro_rmse))
ncdata.close()

ncfile = 'RH_Reforecast_rmse.nc'
ncdata  = Dataset(ncfile, 'r')
RH_reforecast_rmse = ncdata['RHrmse_pt'][:,:] # forcast hours, plev
ncdata.close()


print ('lead (h)   [U250 RMSE retro]    [U250 RMSE reforecast]')
for i in range(len(forecast_hours2)):
    print (forecast_hours2[i],U_retro_rmse[i,4], U_reforecast_rmse[i,4])
    
# ----- make plots

fig = plt.figure(figsize=(6.5,9.))

plt.suptitle('GEFSv12 retro and reforecast differences and RMSE', fontsize=13)

# ----- Differences

a1 = fig.add_axes([.1,.75,.22,.18])
a1.set_title('(a) Mean T difference',fontsize=11)
clevels = [-1.0,-0.7,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.7,1.0]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
cf = a1.contourf(forecast_hours/24,plevs,-np.transpose(T),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours/24,plevs,-np.transpose(T),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.08,0.69,0.26,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('Difference, reforecast - retro (deg C)',fontsize=7)
cb.ax.tick_params(labelsize=6)


a1 = fig.add_axes([.425,.75,.22,.18])
a1.set_title('(b) Mean RH difference',fontsize=11)
clevels = [-10,-7,-5,-3,-1,1,3,5,7,10]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
cf = a1.contourf(forecast_hours/24,plevs,-np.transpose(RH),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours/24,plevs,-np.transpose(RH),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.405,0.69,0.26,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('Difference, reforecast - retro (%)',fontsize=7)
cb.ax.tick_params(labelsize=6)


a1 = fig.add_axes([.74,.75,.22,.18])
a1.set_title('(c) Mean U difference',fontsize=11)
clevels = [-0.7,-0.5,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.5,0.7]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
cf = a1.contourf(forecast_hours/24,plevs,-np.transpose(U),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours/24,plevs,-np.transpose(U),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.72,0.69,0.26,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label(r'Difference, reforecast - retro (m $s^{-1}$)',fontsize=7)
cb.ax.tick_params(labelsize=6)


# ---- Retro RMSE


a1 = fig.add_axes([.1,.43,.22,.18])
a1.set_title('(d) T retro RMSE',fontsize=11)
clevels = [0,.3,.6,1,1.3,1.6,2,2.5,3.0,3.5]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
cf = a1.contourf(forecast_hours2,plevs,np.transpose(T_retro_rmse),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours2,plevs,np.transpose(T_retro_rmse),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.08,0.37,0.26,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label(r'RMSE (deg C)',fontsize=7)
cb.ax.tick_params(labelsize=6)


a1 = fig.add_axes([.425,.43,.21,.18])
a1.set_title('(e) RH retro RMSE',fontsize=11)
clevels = [2, 5, 7, 10, 15, 20, 25, 30, 35, 40]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
cf = a1.contourf(forecast_hours2,plevs,np.transpose(RH_retro_rmse),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours2,plevs,np.transpose(RH_retro_rmse),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.405,0.37,0.26,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label(r'RMSE (%)',fontsize=7)
cb.ax.tick_params(labelsize=6)


a1 = fig.add_axes([.74,.43,.22,.18])
a1.set_title('(f) U retro RMSE',fontsize=11)
clevels = [0,1,2,3,4,5,6,7,8,9,10,12,14]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
cf = a1.contourf(forecast_hours2,plevs,np.transpose(U_retro_rmse),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours2,plevs,np.transpose(U_retro_rmse),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.72,0.37,0.26,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label(r'RMSE (m $s^{-1}$)',fontsize=7)
cb.ax.tick_params(labelsize=6)


# ---- RMSE differences


a1 = fig.add_axes([.1,.11,.22,.18])
a1.set_title('(g) T RMSE difference',fontsize=11)
clevels = [-1.0,-0.7,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.7,1.0]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
cf = a1.contourf(forecast_hours2,plevs,\
    np.transpose(T_reforecast_rmse - T_retro_rmse),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours2,plevs,\
    np.transpose(T_reforecast_rmse - T_retro_rmse),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.08,0.04,0.26,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') 
cb.set_label('RMSE diff., reforecast - retro (deg C)',fontsize=6)
cb.ax.tick_params(labelsize=5)


a1 = fig.add_axes([.425,.11,.22,.18])
a1.set_title('(h) RH RMSE difference',fontsize=11)
clevels = [-10,-7,-5,-3,-1,1,3,5,7,10]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
cf = a1.contourf(forecast_hours2,plevs,\
    np.transpose(RH_reforecast_rmse - RH_retro_rmse),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours2,plevs,\
    np.transpose(RH_reforecast_rmse - RH_retro_rmse),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.405,0.04,0.26,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('RMSE diff., reforecast - retro(%)',fontsize=6)
cb.ax.tick_params(labelsize=5)


a1 = fig.add_axes([.74,.11,.22,.18])
a1.set_title('(i) U RMSE difference',fontsize=11)
clevels = [-0.7,-0.5,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.5,0.7]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
cf = a1.contourf(forecast_hours2,plevs,\
    np.transpose(U_reforecast_rmse - U_retro_rmse),\
    levels = clevels, colors=colorst,extend='both')
a1.contour(forecast_hours2,plevs,\
    np.transpose(U_reforecast_rmse - U_retro_rmse),\
    levels = clevels,colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=9)
a1.set_xlim(0,10.25)
a1.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
a1.set_xlabel('Forecast lead time (days)', fontsize=9)

cax = fig.add_axes([0.72,0.04,0.26,0.01])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label(r'RMSE diff., reforecast - retro (m $s^{-1}$)',fontsize=6)
cb.ax.tick_params(labelsize=5)


plot_title = 'retro_refcst_bias_and_rmse.png'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title,dpi=600)
print ('Plot done')



