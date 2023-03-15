"""

python plot_ccpa_netcdf.py cyyyymmddhh 

Plot the 24-h accumulated precipitation in the merged CCPA/MSWEP
netCDF database on the NDFD 3-km grid.  Save output to .png file.

If you want something other than the 24-h period, lines 50-55 can
be modified to accommodate this.

Uses Jeff Whitaker's dateutils routine.

"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.stats as stats
import pygrib
from netCDF4 import Dataset
from dateutils import hrs_since_day1CE_todate, \
    dateto_hrs_since_day1CE, hrstodate, datetohrs, dateshift
from mpl_toolkits.basemap import Basemap, interp

# ---- get the month and end time from the commmand line

cyyyymmddhh = sys.argv[1] # 01 etc
iyyyymmddhh = int(cyyyymmddhh)
cyear = cyyyymmddhh[0:4]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cday = cyyyymmddhh[6:8]
chour = cyyyymmddhh[8:10]
cyyyymmdd = cyyyymmddhh[0:8]
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']

# ---- get the lat/lons of the output NDFD CONUS grid.   These are
#      oriented S to N as interp requires

infile = '/Volumes/NBM/conus_panal/'+cyear+\
    cmonth+'_ccpa_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile,'r')
lons = nc.variables['lons'][:,:]
lats = nc.variables['lats'][:,:]
conusmask = nc.variables['conusmask'][:,:]
yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
idx = int(np.where(yyyymmddhh_end == iyyyymmddhh)[0])
apcp_anal_0 = nc.variables['apcp_anal'][idx,:,:]
apcp_anal_1 = nc.variables['apcp_anal'][idx-1,:,:]
apcp_anal_2 = nc.variables['apcp_anal'][idx-2,:,:]
apcp_anal_3 = nc.variables['apcp_anal'][idx-3,:,:]
nc.close()
apcp_anal = apcp_anal_0 + apcp_anal_1 + apcp_anal_2 + apcp_anal_3
        
# ===========================================================

m = Basemap(llcrnrlon=-125.,llcrnrlat=32.5,
    urcrnrlon = -113.5, urcrnrlat = 42.,\
    projection='lcc',lat_1=35.,lat_2=35.,lon_0=-120.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons, lats)  
title = 'Analyzed CCPA/MSWEP precipitation,\n24 h ending '+\
    chour+' UTC '+cday+' '+cmonths[imonth]+' '+cyear
    
clevs = [0.0,0.01,0.1, 0.25, 0.5, 0.75, 1.0, \
    1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0,10.0,15.0]
colorst = ['LightGray','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','LightGray']

xdim = 5.
ydim = 6.3
fig = plt.figure(figsize=(xdim,ydim))
axloc = [0.02,0.1,0.96,0.79]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, apcp_anal/25.4, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=5)
cb.set_label('Precipitation (in)',fontsize=9)

# ---- set plot title

plot_title = 'CCPA_daily_precipitation_'+cyyyymmdd+'_'+chour+'UTC.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')   

