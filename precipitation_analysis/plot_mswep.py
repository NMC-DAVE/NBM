"""
plot_mswep.py cmonth 

Plot MSWEP data in a domain surrounding the CONUS that has 
previously been interpolated to the NDFD 3-km grid and
stored in a netCDF file.  

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

# ==========================================================
    
# ---- inputs from command line

cyyyymmddhh = sys.argv[1]
cyear = cyyyymmddhh[0:4]
cmonthnum = cyyyymmddhh[4:6]
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
imonth = int(cmonthnum)-1

# ---- read in the mswep precipitation analyses for this time.
#      These were downloaded from gloh2o.org, v 260 data, then 
#      interpolated to the NDFD grid and saved in netcdf format
#      by the python routine mswep_to_netcdf.py

mswep_directory = '/Volumes/Backup Plus/mswep/'
ncfile = mswep_directory + cyear + \
    cmonthnum + '_on_ndfd_grid_6hourly.nc'
print (ncfile)
nc = Dataset(ncfile)
lons_mswep = nc.variables['lons'][:,:]
lats_mswep = nc.variables['lats'][:,:]
print ('min, max lon = ', lons_mswep[0,0],lons_mswep[-1,-1])
print ('min, max lat = ', lats_mswep[0,0],lats_mswep[-1,-1])
nlats_mswep, nlons_mswep = np.shape(lons_mswep)
yyyymmddhh_begin = nc.variables['yyyymmddhh_begin'][:]
yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
idx = np.where(yyyymmddhh_end == int(cyyyymmddhh))[0]	
precip_mswep = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
nlats_mswep, nlons_mswep = np.shape(precip_mswep)
nc.close()

# --- now plot the MSWEP data interpolated to the NDFD grid
# map parameters from 
# grib_list /Volumes/Backup Plus/blend_domains/blend.t00z.qmd.f001.co.grib2

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
clevs = [0,0.1,0.3,0.5,1.0,2.0,3.0,5.0,7.0,10.0,15.0,20.0,25.0,35.0,50.0]
m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
    urcrnrlon = 300.95782, urcrnrlat = 54.37,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_mswep, lats_mswep)

fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
title = '6-h MSWEP precipitation amount (mm) for valid '+cyyyymmddhh
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, precip_mswep, clevs,\
    cmap=None, colors=colorst, extend='both')

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

plot_title = 'mswep_precip_6h_accumulation_valid_'+\
    cyyyymmddhh+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')
