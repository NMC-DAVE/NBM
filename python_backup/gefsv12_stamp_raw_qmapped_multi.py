"""
gefsv12_stamp_raw_qmapped_multi.py cyyyymmddhh cleadb cleade
"""

import os, sys
from datetime import datetime
from dateutils import dateshift
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
from pyproj import Proj
import cartopy.io.shapereader as shpreader

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'


reader = shpreader.Reader('countyl010g_shp_nt00964/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cf.ShapelyFeature(counties, ccrs.PlateCarree())

reader = shpreader.Reader('statesl010g/statesp010g.shp')
states = list(reader.geometries())
STATES = cf.ShapelyFeature(states, ccrs.PlateCarree())

# =====================================================================
# =====================================================================

# ---- inputs from command line

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


for ilead in range(ileadb, ileade+1, 6):

    if ilead < 10: 
        clead = '00'+str(ilead)
    elif ilead >= 10 and ilead < 100:
        clead = '0'+str(ilead)
    else:
        clead = str(ilead)

    # ---- read ensemble data

    infile = master_directory + cyyyymmddhh+'_lead='+clead+'.cPick'
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
        
    print(ilead, clead, precip_ens_qmapped_ndfd[19,809,923])

# ---- plot the stamp map

#clevs = [0.0,0.254, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0]
clevs = [0,5,10,20,30,40,50,60,70,85,100,120,140,160,180,200]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
#clevs = [0,5,10,20,30,50,75,100,125,150,175,200]

# ----- plot raw data.

fig = plt.figure(figsize=(9.,6.5))

ileadbb = int(cleadb)-6
if ileadbb < 10: 
    cleadbb = '00'+str(ileadbb)
elif ileadbb >= 10 and ileadbb < 100:
    cleadbb = '0'+str(ileadbb)
else:
    cleadbb = str(ileadbb)

fig.suptitle('Raw GEFSv12 total precipitation, IC = '+cyyyymmddhh+\
    ', period '+cleadbb+' to '+cleade+' h',\
    fontsize=12,color='Black')

xbegin = [0.025,0.185,0.345,0.505,0.665,0.825]
xlen = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15]
ybegin = [.76, .59, .42, .25, .08]
ylen = [0.17, 0.17,0.17, 0.17,0.17]

latb = 36.7 # Colorado domain
late = 41.2 # Colorado domain
lonb = -110. # Colorado domain
lone = -100.9 # Colorado domain

proj = ccrs.LambertConformal(\
    central_latitude = (latb+late)/2.,
    central_longitude = (lonb+lone)/2,
    standard_parallels = (latb, late))

process_raw = True
if process_raw == True:
    
    fig.suptitle('Raw GEFSv12 total precipitation, IC = '+cyyyymmddhh+\
        ', period '+cleadbb+' to '+cleade+' h',\
        fontsize=12,color='Black')
    
    for imem in range(30):
        print ('processing ',imem)
        irow = imem//6
        icol = imem - irow*6
        axloc = [xbegin[icol],ybegin[irow],xlen[icol],ylen[irow]]
        ax = plt.axes(axloc, projection = proj)
        if imem == 0:
            ax.set_title('Control', fontsize=8,color='Black')
        else:
            ax.set_title('Perturbed member '+str(imem), fontsize=8,color='Black')
        ax.set_extent([lonb,lone,latb,late])
        #ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        CS = ax.contourf(lons_ndfd, lats_ndfd, precip_ens_raw[imem,:,:], clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.01,0.065,0.98,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Precipitation amount (mm)',fontsize=9)

    plot_title = 'raw_stamp_'+\
        cyyyymmddhh+'_lead'+cleadbb+'_to_'+cleade+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)

# ----- plot quantile-mapped data.

process_qmapped = True
if process_qmapped == True:
    fig = plt.figure(figsize=(9.,6.5))
    fig.suptitle('Quantile-mapped GEFSv12 total precipitation, IC = '+cyyyymmddhh+\
        ', period '+cleadbb+' to '+cleade+' h',\
        fontsize=12,color='Black')

    proj = ccrs.LambertConformal(\
        central_latitude = (latb+late)/2.,
        central_longitude = (lonb+lone)/2,
        standard_parallels = (latb, late))

    for imem in range(30):
        print ('processing ',imem)
        irow = imem//6
        icol = imem - irow*6
        axloc = [xbegin[icol],ybegin[irow],xlen[icol],ylen[irow]]
        ax = plt.axes(axloc, projection = proj)
        if imem == 0:
            ax.set_title('Control', fontsize=8,color='Black')
        else:
            ax.set_title('Perturbed member '+str(imem), fontsize=8,color='Black')
        ax.set_extent([lonb,lone,latb,late])
        #ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        CS = ax.contourf(lons_ndfd, lats_ndfd, precip_ens_qmapped[imem,:,:], clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.01,0.065,0.98,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Precipitation amount (mm)',fontsize=9)

    plot_title = 'qmapped_stamp_'+\
        cyyyymmddhh+'_lead'+cleadbb+'_to_'+cleade+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)
    
    



