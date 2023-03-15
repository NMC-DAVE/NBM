"""

upscale_ccpa_multiple_timeslice.py cyyyymmddhh_begin cyyyymmddhh_end  

"""

import os, sys
from datetime import datetime
from dateutils import daterange, dateshift
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
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import _pickle as cPickle
import scipy.stats as stats
import cartopy.io.shapereader as shpreader
from PIL import Image

# ---- get the month and end time from the commmand line

cyyyymmddhh_begin = sys.argv[1] # 01 etc
cyyyymmddhh_end = sys.argv[2] # 01 etc


cyyyymmddhh_realbegin = dateshift(cyyyymmddhh_begin, -6)

# ---- read Cartopy shapes 

reader = shpreader.Reader('countyl010g_shp_nt00964/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cf.ShapelyFeature(counties, ccrs.PlateCarree())

reader = shpreader.Reader('statesl010g/statesp010g.shp')
states = list(reader.geometries())
STATES = cf.ShapelyFeature(states, ccrs.PlateCarree())

# ---- loop over analysis times to read in, and sum.

cyyyymmddhh_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 6)

for itime, cyyyymmddhh in enumerate(cyyyymmddhh_list):
    cyyyymmdd = cyyyymmddhh[0:8]
    chh = cyyyymmddhh[8:10]

    # ---- get the lat/lons of the output NDFD CONUS grid.   These are
    #      oriented S to N as interp requires

    infile = '/Volumes/NBM/ccpa/ccpa.'+cyyyymmdd+'/'+\
        chh+'/ccpa.t'+chh+'z.06h.ndgd5p0.conus.gb2'
    print (infile)
    flatlon = pygrib.open(infile)
    anal = flatlon.select(shortName='tp')[0]
    lats_ndfd, lons_ndfd = anal.latlons()
    anal_values = anal.values
    ny, nx = np.shape(lats_ndfd)
    dy = 111.*(lats_ndfd[ny//2,nx//2]- lats_ndfd[ny//2-1,nx//2])
    coslat = np.cos(3.14159*lats_ndfd[ny//2,nx//2]/180.)
    dx = dy*coslat
    dxy = np.sqrt(dx**2 + dy**2)
    if itime == 0:
        accum_precip = np.copy(anal_values)
    else:
        accum_precip = accum_precip + anal_values
        
        
a = np.where(accum_precip > 500.)
accum_precip[a] = 0.0

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

lats_upscaled = np.zeros((ny//5, nx//5), dtype=np.float32)
lons_upscaled = np.zeros((ny//5, nx//5), dtype=np.float32)
accum_precip_upscaled = np.zeros((ny//5, nx//5), dtype=np.float32)
for jy in range(5,ny-6):
    jyout = jy//5
    for ix in range(5,nx-6):
        ixout = ix//5
        lats_upscaled[jyout,ixout] = lats_ndfd[jy,ix]
        lons_upscaled[jyout,ixout] = lons_ndfd[jy,ix]
        accum_precip_upscaled[jyout,ixout] = np.mean(accum_precip[jy-3:jy+3,ix-3:ix+3])

#print ('performing upscaling of lat/lon arrays, initializing work arrays.')
#im = Image.fromarray(lons_ndfd)
#imu = im.resize(size=(ny//10, nx//10),resample=Image.BOX)
#lons_upscaled = np.transpose(np.asarray(imu))
#im = Image.fromarray(lats_ndfd)
#imu = im.resize(size=(ny//10, nx//10),resample=Image.BOX)
#lats_upscaled = np.transpose(np.asarray(imu))
#ny_upscaled, nx_upscaled = np.shape(lons_upscaled)


#im = Image.fromarray(accum_precip)
#imu = im.resize(size=(ny//10, nx//10),resample=Image.BOX)
#ccpa_upscaled = np.transpose(np.asarray(imu))

# ---- read in the CONUS mask.  Not sure about accuracy.

infile = '/Volumes/Backup Plus/ccpa/supplemental_locations_ndfd2p5_Jan.nc'
nc = Dataset(infile)
conusmask_in = nc.variables['conusmask'][:,:]
conusmask = ma.masked_equal(conusmask_in, 0)
nc.close()
        
# ===========================================================

#latb = 36.7 # Colorado domain
#late = 41.2 # Colorado domain 4.5
#lonb = -110. # Colorado domain
#lone = -100.9 # Colorado domain 9.1

latb = 39. # NY
late = 44. # NY
lonb = -72. # NY
lone = -79. # NY

xdim = 6.
ydim = 5.8
drawcoasts = False

proj = ccrs.LambertConformal(\
    central_latitude = (latb+late)/2.,
    central_longitude = (lonb+lone)/2,
    standard_parallels = (latb, late))
 
title = 'CCPA precipitation, '+cyyyymmddhh_realbegin+' to '+cyyyymmddhh_end
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
#clevs = [-100,-0.00001,0.00001,0.1,0.3,0.6,1,2,3,5,7,10,15,20,25]
clevs = [0,5,10,20,30,40,50,60,70,85,100,120,140,160,180,200]

fig = plt.figure(figsize=(xdim, ydim))
axloc = [0.02,0.14,0.96,0.78]
ax = plt.axes(axloc,projection = proj)
ax.set_extent([lonb,lone,latb,late])
ax.coastlines(resolution='50m',lw=0.5)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
ax.add_feature(cf.BORDERS,lw=0.5)
ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
ax.set_title(title, fontsize=13,color='Black')
CS = ax.contourf(lons_ndfd, lats_ndfd, \
    accum_precip, clevs, cmap=None, colors=colorst, \
    extend='both', transform=ccrs.PlateCarree())

cax = fig.add_axes([0.02,0.11,0.96,0.02])
cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=9)
cb.set_label('precipitation amount (mm)',\
    fontsize=11)

plot_title = 'CCPA_NY.png'
fig.savefig(plot_title, dpi=400)
plt.close()
print ('saving plot to file = ',plot_title)



title = 'Upscaled CCPA precipitation, '+cyyyymmddhh_realbegin+' to '+cyyyymmddhh_end
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
#clevs = [-100,-0.00001,0.00001,0.1,0.3,0.6,1,2,3,5,7,10,15,20,25]
clevs = [0,5,10,20,30,40,50,60,70,85,100,120,140,160,180,200]

fig = plt.figure(figsize=(xdim, ydim))
axloc = [0.02,0.14,0.96,0.78]
ax = plt.axes(axloc,projection = proj)
ax.set_extent([lonb,lone,latb,late])
ax.coastlines(resolution='50m',lw=0.5)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
ax.add_feature(cf.BORDERS,lw=0.5)
ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
ax.set_title(title, fontsize=13,color='Black')
CS = ax.contourf(lons_upscaled, lats_upscaled, \
    accum_precip_upscaled, clevs, cmap=None, colors=colorst, \
    extend='both', transform=ccrs.PlateCarree())

cax = fig.add_axes([0.02,0.11,0.96,0.02])
cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=9)
cb.set_label('precipitation amount (mm)',\
    fontsize=11)

plot_title = 'CCPA_upscaled_NY.png'
fig.savefig(plot_title, dpi=400)
plt.close()
print ('saving plot to file = ',plot_title)
