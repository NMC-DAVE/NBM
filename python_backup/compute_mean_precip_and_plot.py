"""
compute_mean_precip_and_plot.py cmonth cleadb cleade

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
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import _pickle as cPickle
import scipy.stats as stats
import cartopy.io.shapereader as shpreader

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'


cyyyymmddhh = sys.argv[1] # 
cleadb = sys.argv[2] # 018 not 18
cleade = sys.argv[3] # 018 not 18
ileadb = int(cleadb)
ileade = int(cleade)
cyear = cyyyymmddhh[0:4]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
master_directory = '/Volumes/NBM/conus_gefsv12/qmapped/'

reader = shpreader.Reader('countyl010g_shp_nt00964/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cf.ShapelyFeature(counties, ccrs.PlateCarree())

reader = shpreader.Reader('statesl010g/statesp010g.shp')
states = list(reader.geometries())
STATES = cf.ShapelyFeature(states, ccrs.PlateCarree())


for ilead in range(ileadb, ileade+1, 6):

    if ilead < 10: 
        clead = '00'+str(ilead)
    elif ilead >= 10 and ilead < 100:
        clead = '0'+str(ilead)
    else:
        clead = str(ilead)

    # ---- read ensemble data

    infile = master_directory + cyyyymmddhh+'_use99_lead='+clead+'.cPick'
    print ('reading raw and qmapped ens from ', infile) 
    inf = open(infile,'rb')
    precip_ens_raw_ndfd = cPickle.load(inf)
    precip_ens_qmapped_ndfd = cPickle.load(inf)
    lons_ndfd = cPickle.load(inf)
    lats_ndfd = cPickle.load(inf)
    inf.close()
    
    if ilead == ileadb:
        precip_ens_raw = precip_ens_raw_ndfd
        precip_ens_qmapped = precip_ens_qmapped_ndfd
    else:
        precip_ens_raw = precip_ens_raw + precip_ens_raw_ndfd
        precip_ens_qmapped =  precip_ens_qmapped + precip_ens_qmapped_ndfd


precip_raw_mean = np.mean(precip_ens_raw, axis=0)
precip_qmapped_mean = np.mean(precip_ens_qmapped, axis=0)

# ---- plot the stamp map

clevs = [0, 1, 3, 5, 10, 15, 20, 25, 30, 40, 50, 65, 80, 100.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

# ---- plot the storm total ens mean over the Rockies

#latb = 33.
#late = 47.
#lonb = -117.
#lone = -101.9

latb = 36.7 # Colorado domain
late = 41.2 # Colorado domain 4.5
lonb = -110. # Colorado domain
lone = -100.9 # Colorado domain 9.1

xdim = 6.
ydim = 4.8
drawcoasts = False

proj = ccrs.LambertConformal(\
    central_latitude = (latb+late)/2.,
    central_longitude = (lonb+lone)/2,
    standard_parallels = (latb, late))

fig = plt.figure(figsize=(xdim, ydim))
axloc = [0.02,0.14,0.96,0.78]
ax = plt.axes(axloc,projection = proj)
ax.set_extent([lonb,lone,latb,late])
ax.coastlines(resolution='50m',lw=0.5)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
ax.add_feature(cf.BORDERS,lw=0.5)
ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)

title = 'Raw ensemble-mean forecast storm-total meltwater'
ax.set_title(title, fontsize=15,color='Black')
CS = ax.contourf(lons_ndfd, lats_ndfd, \
    precip_raw_mean, clevs, cmap=None, colors=colorst, \
    extend='both', transform=ccrs.PlateCarree())

cax = fig.add_axes([0.02,0.11,0.96,0.02])
cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=9)
cb.set_label('precipitation amount (mm)',\
    fontsize=11)

plot_title = 'storm_total_raw_Colorado.png'
fig.savefig(plot_title, dpi=400)
plt.close()
print ('saving plot to file = ',plot_title)



fig = plt.figure(figsize=(xdim, ydim))
axloc = [0.02,0.14,0.96,0.78]
ax = plt.axes(axloc,projection = proj)
ax.set_extent([lonb,lone,latb,late])
ax.coastlines(resolution='50m',lw=0.5)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
ax.add_feature(cf.BORDERS,lw=0.5)
ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)

title = 'Statistically modified forecast storm-total meltwater'
ax.set_title(title, fontsize=15,color='Black')
CS = ax.contourf(lons_ndfd, lats_ndfd, \
    precip_qmapped_mean, clevs, cmap=None, colors=colorst, \
    extend='both', transform=ccrs.PlateCarree())

cax = fig.add_axes([0.02,0.11,0.96,0.02])
cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=9)
cb.set_label('precipitation amount (mm)',\
    fontsize=11)

plot_title = 'storm_total_qmapped_Colorado.png'
fig.savefig(plot_title, dpi=400)
plt.close()
print ('saving plot to file = ',plot_title)

              
                










