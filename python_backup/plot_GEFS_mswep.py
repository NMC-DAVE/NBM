"""
plot_GEFS_mswep.py cmonth clead jy ix

"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import _pickle as cPickle
import scipy.stats as stats
rcParams['legend.fontsize']='medium'

# =====================================================================
    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clead = sys.argv[2] # 03, 06, 12, etc.
cyyyymmddhh_fcst = sys.argv[3]
cyear = cyyyymmddhh_fcst[0:4]
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
imonth = cmonths.index(cmonth)+1
if imonth < 10:
    cmonthnum = '0'+str(imonth)
else:
    cmonthnum = str(imonth)

master_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'

# ---- read in the previously generated netCDF file with precipitation
#      forecast for this month and lead time.  All members, dates for this 
#      month have been smushed into one leading index, dimension
#      nsamps, since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution fitting.

jmin = 93
jmax = 246
imin = 368
imax = 686
ncfile = master_directory + cmonth + '_apcp' '_h' + clead + '.nc'
print (ncfile)
nc = Dataset(ncfile)
lons_1d = nc.variables['lons_fcst'][imin:imax]
lats_1d = nc.variables['lats_fcst'][jmin:jmax]
yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
yyyymmddhh_fcst = nc.variables['yyyymmddhh_fcst'][:]
indices = np.where(yyyymmddhh_fcst == int(cyyyymmddhh_fcst))	
precip_fcst = nc.variables['apcp_fcst'][indices[0],jmin:jmax,imin:imax]
precip_mean = np.mean(precip_fcst, axis=0)
nlats_gefs, nlon_gefs = np.shape(precip_mean)
nc.close()
lons_fcst_2d, lats_fcst_2d = np.meshgrid(lons_1d,lats_1d)

# ---- read in the mswep precipitation analyses for this time

mswep_directory = '/Volumes/Backup Plus/mswep/'
ncfile = mswep_directory + cyear + cmonthnum + '_on_ndfd_grid_6hourly.nc'
print (ncfile)
nc = Dataset(ncfile)
lons_mswep = nc.variables['lons'][:,:]
lats_mswep = nc.variables['lats'][:,:]
print ('min, max lon = ', lons_mswep[0,0],lons_mswep[-1,-1])
print ('min, max lat = ', lats_mswep[0,0],lats_mswep[-1,-1])
nlats_mswep, nlons_mswep = np.shape(lons_mswep)
yyyymmddhh_begin = nc.variables['yyyymmddhh_begin'][:]
yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
idx = np.where(yyyymmddhh_end == int(cyyyymmddhh_fcst))[0]	
precip_mswep = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
nlats_mswep, nlons_mswep = np.shape(precip_mswep)
nc.close()

# ---- first plot GEFS control forecast amount

m = Basemap(llcrnrlon=lons_fcst_2d[0,0],llcrnrlat=lats_fcst_2d[-1,-1],\
    urcrnrlon=lons_fcst_2d[-1,-1],urcrnrlat=lats_fcst_2d[0,0],\
    resolution='l', projection='mill')
x, y = m(lons_fcst_2d, lats_fcst_2d)

# ---- make plots of fraction positive precip

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

clevs = [0.0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20]
fig = plt.figure(figsize=(9,6.))
axloc = [0.02,0.1,0.96,0.8]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Ensemble-mean precipitation amount (mm) for '+clead+\
    '-h forecast valid '+cyyyymmddhh_fcst
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, precip_mean, clevs,\
    cmap=None, colors=colorst, extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.02,0.07,0.96,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Mean precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'forecast_precip_'+clead+'_h_valid_'+cyyyymmddhh_fcst+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')




# --- now plot the MSWEP data interpolated to the NDFD grid
# map parameters from 
# grib_list /Volumes/Backup Plus/blend_domains/blend.t00z.qmd.f001.co.grib2


m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_mswep, lats_mswep)

fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = '6-h MSWEP precipitation amount (mm) for valid '+cyyyymmddhh_fcst
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, precip_mswep, clevs,\
    cmap=None, colors=colorst, extend='both')

#for ilon in range(0,nlons_mswep,10): 
#    for ilat in range(0,nlats_mswep,10): 
#        xdot, ydot = m(lons_mswep[ilat,ilon], lats_mswep[ilat,ilon])
#        m.plot(xdot,ydot,marker='.',markersize=0.2,color='Black')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'mswep_precip_6h_accumulation_valid_'+cyyyymmddhh_fcst+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')
