"""

Plot the ccpa data from its original grib file 

python plot_ccpa_timeslice.py cyyyymmddhh

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

# ---- get the month and end time from the commmand line.
#      Set input directory

cyyyymmddhh = sys.argv[1] # 01-12 
cyyyymm = cyyyymmddhh[0:6]
chour = cyyyymmddhh[8:10]
cyyyymmdd = cyyyymmddhh[0:8]
input_directory = '/Volumes/Backup Plus/ccpa/'

# ---- get the lat/lons of the NDFD CONUS grid.   These should be
#      oriented S to N as interp requires

infile = input_directory + 'ccpa.'+cyyyymmdd+\
    '/00/ccpa.t00z.06h.ndgd2p5.conus.gb2'
print (infile)
flatlon = pygrib.open(infile)
fcst = flatlon.select()[0]
lats_ndfd, lons_ndfd = fcst.latlons()
ny, nx = np.shape(lats_ndfd)
dy = 111.*(lats_ndfd[ny//2,nx//2]- lats_ndfd[ny//2-1,nx//2])
coslat = np.cos(3.14159*lats_ndfd[ny//2,nx//2]/180.)
dx = dy*coslat
dxy = np.sqrt(dx**2 + dy**2)

if lats_ndfd[0,0] > lats_ndfd[-1,0]: 
    flipud = True
else:
    flipud = False
if flipud == True:
    lats_ndfd = np.flipud(lats_ndfd)
    lons_ndfd = np.flipud(lons_ndfd)
nlats_ndfd, nlons_ndfd = np.shape(lons_ndfd)
flatlon.close()
print ('min, max lons_ndfd = ', np.min(lons_ndfd), np.max(lons_ndfd))

# ---- read in the CONUS mask. 

master_directory = '/Volumes/NBM/conus_panal/'
infile = master_directory + cyyyymm + \
    '_ccpa_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile)
conusmask_in = nc.variables['conusmask'][:,:]
nc.close()

# ---- Read in the precipitation analysis data

try:
    infile = input_directory + 'ccpa.'+cyyyymmdd + \
        '/'+chour+'/ccpa.t'+chour+'z.06h.ndgd2p5.conus.gb2'
    print (infile)
    grb = pygrib.open(infile)
    print (grb)
    panal = grb.select()[0]
    precip_ccpa = panal.values * conusmask_in
    print ('max, min precip_ccpa = ', np.max(precip_ccpa), np.min(precip_ccpa))
    if flipud == True:
        precip_ccpa = np.flipud(precip_ccpa)
    grb.close()
    zeros = np.zeros((nlats_ndfd, nlons_ndfd), dtype=np.float64)
    ones = np.ones((nlats_ndfd, nlons_ndfd), dtype=np.float64)
    print (precip_ccpa[nlats_ndfd//2,0:nlons_ndfd//4:5])
    apcp_mask = ma.getmask(precip_ccpa)
    apcp_mask_ones = ma.where(apcp_mask == True, zeros, ones)                
except:
    print ('whoops!   some problem with ', infile)
        
# ===========================================================


# ---- make plot of precipitation amount.

m = Basemap(llcrnrlon=233.7234-360.,llcrnrlat=19.229,
    urcrnrlon = 300.95782-360., urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ndfd, lats_ndfd)  
title = 'Masked CCPA precipitation, 6 h ending '+chour+' UTC '+cyyyymmdd
    
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
clevs = [0.0,0.1,0.3,0.6,1,2,3,5,7,10,15,20,25]
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, conusmask_in*precip_ccpa, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=5)
cb.set_label('Precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'CCPA_daily_precipitation_'+cyyyymmdd+'_'+chour+'UTC.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')   






# ------ plot the mask

title = 'CCPA implicit CONUS mask'  
clevs = [-100,-0.00001,0.00001,0.05,0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,5,7,10]
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, conusmask_in, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=5)
cb.set_label('mask',fontsize=9)

# ---- set plot title

plot_title = 'CCPA_mask.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')    
