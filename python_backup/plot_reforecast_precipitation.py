"""
plot_reforecast_precipitation.py

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
from netCDF4 import Dataset

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# ---- inputs from command line

cmonth = sys.argv[1]
cyear = sys.argv[2]
clead = sys.argv[3]
cyyyymmddhh_init = sys.argv[4]
iyyyymmddhh_init = int(cyyyymmddhh_init)

master_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'
ncfile = master_directory + cmonth + cyear + '_h' + clead + '.nc'
print (ncfile)
nc = Dataset(ncfile)
yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
idx = np.where(yyyymmddhh_init==iyyyymmddhh_init)[0]

print ('matching index at ', idx)
apcp_ens = np.squeeze(nc.variables['apcp_fcst'][idx,:,:,:])
lons_1d = nc.variables['lons_fcst'][:]
lats_1d = nc.variables['lats_fcst'][:]
lats_1d = np.flipud(lats_1d)
print (lons_1d)
print (lats_1d)
nc.close()
lons_2d, lats_2d = np.meshgrid(lons_1d,lats_1d)

# ---- calculate ensemble mean 

print (np.shape(apcp_ens))
print ('max, min apcp_ens = ', np.max(apcp_ens), np.min(apcp_ens))
apcp_mean = np.flipud(np.mean(apcp_ens, axis=0))
print (np.shape(apcp_mean))
print ('max, min apcp_mean = ', np.max(apcp_mean), np.min(apcp_mean))

# --- now plot the precipitation

clevs = [0.0,0.1,0.5,0.7,1.0,2.0,3.0,4.0,5.0,\
    6.0,8.0,10.0,12.0,15.0,20.0,25.0,30.0,50.0,100.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','Gray','Black']

fig = plt.figure(figsize=(9.,6.5))
axloc = [0.02,0.11,0.96,0.8]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 ensemble-mean precipitation, IC =  '+cyyyymmddhh_init+\
    ', '+clead+'-hour forecast'
ax1.set_title(title, fontsize=14,color='Black')
m = Basemap(llcrnrlon=lons_2d[0,0],llcrnrlat=lats_2d[0,0],\
    urcrnrlon=lons_2d[-1,-1],urcrnrlat=lats_2d[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons_2d, lats_2d)
CS2 = m.contourf(x,y,apcp_mean, clevs,alpha=0.6, \
    cmap=None, colors=colorst, extend='both') #
m.drawcoastlines(linewidth=0.5,color='Gray',zorder=32)
m.drawcountries(linewidth=0.5,color='Gray',zorder=32)
m.drawstates(linewidth=0.5,color='Gray',zorder=32)

# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.02,0.08,0.96,0.03])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,format='%g') # ticks=clevs,
cb.ax.tick_params(labelsize=7)
cb.ax.set_xticks(clevs)
#cb.ax.set_xticklabels(['0.0','0.1','0.5','0.7','1.0','2.0','3.0','4.0','5.0',\
#    '6.0','8.0','10.0','12.0','15.0','20.0','25.0','30.0','50.0','100.0'])  # horizontal colorbar
cb.set_label('Ensemble-mean precipitation amount (mm)')

# ---- set plot title

plot_title = 'precip_ens_mean_IC'+cyyyymmddhh_init+'_h'+\
    clead+'.png'
fig.savefig(plot_title, dpi=300,fontsize=9)

print ('saving plot to file = ',plot_title)
print ('Done!')