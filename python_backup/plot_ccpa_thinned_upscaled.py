"""

Plot the netCDF file of full field, thinned, or upscaled precipitation
for the date chosen

python plot_ccpa_thinned_upscaled.py cyyyymmddhh ctype

Tom Hamill

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

cyyyymmddhh_end = sys.argv[1] # 01 etc
ctype = sys.argv[2] # 'fullfield','thinned','upscaled'
print ('ctype = ', ctype)
cyyyymm_anal = cyyyymmddhh_end[0:6]

# ---- read the precipitation analyses for the chosen date

master_directory = '/Volumes/NBM/conus_panal/'
if ctype != 'fullfield':
    infile = master_directory + cyyyymm_anal + \
        '_ccpa_on_ndfd_grid_6hourly_'+ctype+'.nc'
else:
    infile = master_directory + cyyyymm_anal + \
        '_ccpa_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile)
yyyymmddhh_end_in = nc.variables['yyyymmddhh_end'][:]
idx = np.where(yyyymmddhh_end_in == int(cyyyymmddhh_end))[0]
conusmask_in = nc.variables['conusmask'][:,:]
precip_anal = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
ny_anal, nx_anal = np.shape(precip_anal)
ones = np.ones((ny_anal, nx_anal), dtype=np.float64)
zeros = np.zeros((ny_anal, nx_anal), dtype=np.float64)
lons_in = nc.variables['lons'][:,:]
lats_in = nc.variables['lats'][:,:]
weight = np.where(conusmask_in == 1.0, \
    ones*np.cos(lats_in*3.1415926/180.), zeros)

nc.close()    
        
# ===========================================================

m = Basemap(llcrnrlon=233.7234-360.,llcrnrlat=19.229,
    urcrnrlon = 300.95782-360., urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_in, lats_in)  
title = ctype + ', merged CCPA/MSWEP precipitation, 6 h ending '+cyyyymmddhh_end 
    
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
clevs = [0.0,0.25,0.5,1,2,3,5,7,10,15,20,25,35,50,75,100]
cclevs = ['0.0','0.25\n(0.01)','0.5\n(0.02)','1\n(0.04)','2\n(0.08)',\
    '3\n(0.12)','5\n(0.2)','7\n(0.28)','10\n(0.4)','15\n(0.6)',\
    '20\n(0.8)','25\n(1.0)','35\n(1.4)','50\n(2.0)','75\n(3.0)',\
    '100\n(4.0)']

fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.11,0.96,0.84]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, precip_anal, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.08,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.ax.set_xticklabels(cclevs) 
cb.set_label('Precipitation in mm and (in)',fontsize=9)

# ---- set plot title

plot_title = 'CCPA_'+ctype+'_'+cyyyymmddhh_end+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')   
  
