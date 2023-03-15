"""
gefsv12_ensmean_cartopy_v2.py cyyyymmddhh cleadb cleade
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

#states_provinces = cf.NaturalEarthFeature( category='cultural', \
#    name='admin_1_states_provinces_lines', scale='50m', facecolor='none') 


# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
cleadb = sys.argv[2] # 018 not 18
cleade = sys.argv[3] # 018 not 18

ileadb = int(cleadb)
ileade = int(cleade)
clead_start = '0'+str(ileadb - 6)
for ilead in range(ileadb, ileade+1, 6):
    if ilead < 10:
        clead = '00'+str(ilead)
    elif ilead > 10 and ilead < 100:
        clead = '0'+str(ilead)
    else:
        clead = str(ilead)
        
    master_directory = '/Volumes/NBM/conus_gefsv12/qmapped/'
    infile = master_directory + \
        cyyyymmddhh+'_mean_lead='+clead+'.cPick'
    print ('reading raw and qmapped mean from ', infile)
    inf = open(infile,'rb')
    ensmean_raw_ndfd = cPickle.load(inf)
    ensmean_qmapped_ndfd = cPickle.load(inf) 
    print ('max(ensmean_qmapped_ndfd) = ', np.max(ensmean_qmapped_ndfd))
    lons_ndfd = cPickle.load(inf)    
    lats_ndfd = cPickle.load(inf)   
    inf.close() 
    
    if ilead == ileadb:
        ensmean_raw_sum = np.copy(ensmean_raw_ndfd)
        ensmean_qmapped_sum = np.copy(ensmean_qmapped_ndfd)
    else:
        ensmean_raw_sum = ensmean_raw_sum + ensmean_raw_ndfd
        ensmean_qmapped_sum = ensmean_qmapped_sum + ensmean_qmapped_ndfd

print ('max(ensmean_qmapped_sum) = ', np.max(ensmean_qmapped_sum))
#sys.exit()

# ---- plot the raw ensemble mean by sector

clevs = [0.0,0.254, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
sectors = ['CONUS','California','Northwest_US','Central_Rockies',\
        'Northeast_US','Southeast_US','Great_Plains']
nsectors = len(sectors)
plot_CONUS = False
plot_California = False
plot_northwest = False
plot_central_rockies = True
plot_northeast = False
plot_southeast = False
plot_plains = False

for isector,sector in enumerate(sectors):
    print (isector, sector, plot_central_rockies)
    plotit = False
    if isector == 0 and plot_CONUS == True:
        plotit = True
        latb = 23.
        late = 50.
        lonb = -125.
        lone = -71.
        xdim = 9.
        ydim = 6.5
        drawcoasts = True
    if isector == 1 and plot_California == True:
        plotit = True
        latb = 31.
        late = 43.
        lonb = -125.
        lone = -112.
        xdim = 6.
        ydim = 8.
        drawcoasts = True
    if isector == 2 and plot_northwest == True:
        plotit = True
        latb = 40.
        late = 51.
        lonb = -126.
        lone = -105.
        xdim = 8.
        ydim = 6.
        drawcoasts = True
    if isector == 3 and plot_central_rockies == True:
        print ('plotting central Rockies')
        plotit = True
        latb = 33.
        late = 47.
        lonb = -117.
        lone = -101.9
        xdim = 6.
        ydim = 8.
        drawcoasts = False
    if isector == 4 and plot_northeast == True:
        plotit = True
        latb = 36.
        late = 49.
        lonb = -90.
        lone = -65.
        xdim = 9.
        ydim = 6.5
        drawcoasts = True
    if isector == 5 and plot_southeast == True:
        plotit = True
        latb = 24.
        late = 38.
        lonb = -95.
        lone = -75.
        xdim = 6.5
        ydim = 6.
        drawcoasts = True
    if isector == 6 and plot_plains == True:
        plotit = Tue
        latb = 25.
        late = 50.
        lonb = -108.
        lone = -88.
        xdim = 6.5
        ydim = 9.
        drawcoasts = True

    if plotit == True:
        
        proj = ccrs.LambertConformal(\
            central_latitude = (latb+late)/2., 
            central_longitude = (lonb+lone)/2, 
            standard_parallels = (latb, late))

        # ---- plot the raw ensemble mean

        fig = plt.figure(figsize=(xdim,ydim))
        axloc = [0.02,0.1,0.96,0.81]
        ax = plt.axes(axloc,projection = proj)
        ax.set_extent([lonb,lone,latb,late])
        if drawcoasts == True: ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(cf.BORDERS,lw=0.5)
        #ax.add_feature(cf.STATES,lw=0.5)
        #ax.states(cf.STATES,'50m',lw=0.5)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        #ax.add_feature(states_provinces.with_scale('50m'), edgecolor='black') 
        #ax.add_feature(cf.BORDERS.with_scale('50m'))
        #ax.add_feature(cfeature.STATES.with_scale('50m'))
        
        cleadb = str(int(clead)-6)
        title = 'Raw GEFSv12 ensemble mean,\nIC = '+\
            cyyyymmddhh+' leads '+clead_start+' to '+cleade+' h'
        ax.set_title(title, fontsize=14,color='Black')
        CS = ax.contourf(lons_ndfd, lats_ndfd, ensmean_raw_sum, clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

        cax = fig.add_axes([0.08,0.07,0.84,0.02])
        cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Precipitation amount (mm)',fontsize=9)

        plot_title = 'raw_ensmean_'+sector+'_'+\
            cyyyymmddhh+'_lead'+clead+'.png'
        fig.savefig(plot_title, dpi=900)
        plt.close()
        print ('saving plot to file = ',plot_title)
        
        
        
        # ---- plot the quantile-mapped ensemble mean

        fig = plt.figure(figsize=(xdim,ydim))
        axloc = [0.02,0.1,0.96,0.81]
        ax = plt.axes(axloc,projection = proj)
        ax.set_extent([lonb,lone,latb,late])
        if drawcoasts == True: ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(cf.BORDERS,lw=0.5)
        #ax.add_feature(cf.STATES,lw=0.5)
        #ax.states(cf.STATES,'50m',lw=0.5)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        #ax.add_feature(states_provinces.with_scale('50m'), edgecolor='black') 
        #ax.add_feature(cf.BORDERS.with_scale('50m'))
        #ax.add_feature(cfeature.STATES.with_scale('50m'))
        
        cleadb = str(int(clead)-6)
        title = 'Quantile-mapped GEFSv12 ensemble mean,\nIC = '+\
            cyyyymmddhh+' leads '+clead_start+' to '+cleade+' h'
        ax.set_title(title, fontsize=14,color='Black')
        CS = ax.contourf(lons_ndfd, lats_ndfd, ensmean_qmapped_sum, clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

        cax = fig.add_axes([0.08,0.07,0.84,0.02])
        cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Precipitation amount (mm)',fontsize=9)

        plot_title = 'qmapped_ensmean_'+sector+'_'+\
            cyyyymmddhh+'_lead'+clead+'.png'
        fig.savefig(plot_title, dpi=900)
        plt.close()
        print ('saving plot to file = ',plot_title)


