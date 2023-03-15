"""
gefsv12_ensmean_multilead_cartopy.py cyyyymmddhh cleadb cleade
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
cleadbb = '0'+str(int(cleadb)-6)
cleade = sys.argv[3] # 018 not 18
cyear = cyyyymmddhh[0:4]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'

ileadb = int(cleadb)
ileade = int(cleade) +1
for ilead in range(ileadb, ileade, 6):
    
    if ilead < 10:
        clead = '00'+str(ilead)
    elif ilead >= 10 and ilead < 100:
        clead = '0'+str(ilead)
    else:
        clead = str(ilead)

    # ---- read the stored netCDF ensemble mean, raw 
    #      and quantile mapped

    infile = master_directory_fullfield_output +\
        ccmonth+cyear+'_use99_lead'+clead+\
        '_probabilities_fullfield.nc'
    print ('reading from ', infile)
    nc = Dataset(infile,'r')
    yyyymmddhh = nc.variables['yyyymmddhh_init'][:]
    #print (yyyymmddhh)
    thresholds = nc.variables['thresholdv'][:]
    #print (thresholds)
    idd = int(np.where(yyyymmddhh == int(cyyyymmddhh))[0])
    lons_ndfd = nc.variables['lons'][:,:]
    lats_ndfd = nc.variables['lats'][:,:]
    ensmean_raw_fullfield = nc.variables['ensmean_raw'][idd,:,:]
    ensmean_qmapped_fullfield = nc.variables['ensmean_qmapped'][idd,:,:]
    nc.close()
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
    lons_ndfd = np.flipud(lons_ndfd)
    lats_ndfd = np.flipud(lats_ndfd)
    ensmean_raw_fullfield = np.flipud(ensmean_raw_fullfield)
    ensmean_qmapped_fullfield = np.flipud(ensmean_qmapped_fullfield)

    if ilead == ileadb:
        ensmean_raw = ensmean_raw_fullfield
        ensmean_qmapped = ensmean_qmapped_fullfield
    else:
        ensmean_raw = ensmean_raw + ensmean_raw_fullfield
        ensmean_qmapped = ensmean_qmapped + ensmean_qmapped_fullfield

# ---- plot the raw ensemble mean by sector

#clevs = [0.0,0.254, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0]
#clevs = [0,5,10,20,30,50,75,100,125,150,175,200]

clevs = [0,2,5,10,15,20,25,30,40,50,60,75,100,125]

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
sectors = ['CONUS','California','Northwest_US','Central_Rockies',\
    'Northeast_US','Southeast_US','Great_Plains']
nsectors = len(sectors)
plot_CONUS = False
plot_California = False
plot_northwest = True
plot_central_rockies = False
plot_northeast = False
plot_southeast = False
plot_plains = False

for isector,sector in enumerate(sectors):

    plotit = False
    print ('isector, sector, plot_northwest = ', isector, sector, plot_northwest)
    if isector == 0 and plot_CONUS == True:
        plotit = True
        latb = 23.
        late = 50.
        lonb = -125.
        lone = -71.
        xdim = 9.
        ydim = 6.
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
            cyyyymmddhh+', '+cleadbb+' to '+cleade+' h'
        ax.set_title(title, fontsize=14,color='Black')
        CS = ax.contourf(lons_ndfd, lats_ndfd, ensmean_raw, clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

        cax = fig.add_axes([0.08,0.07,0.84,0.02])
        cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Precipitation amount (mm)',fontsize=9)

        plot_title = 'raw_ensmean_'+sector+'_'+\
            cyyyymmddhh+'_lead'+cleadb+'_to_'+cleade+'.png'
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
            cyyyymmddhh+', '+cleadbb+' to '+cleade+' h'
        ax.set_title(title, fontsize=14,color='Black')
        CS = ax.contourf(lons_ndfd, lats_ndfd, ensmean_qmapped, clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

        cax = fig.add_axes([0.08,0.07,0.84,0.02])
        cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Precipitation amount (mm)',fontsize=9)

        plot_title = 'qmapped_ensmean_'+sector+'_'+\
            cyyyymmddhh+'_lead'+cleadb+'_to_'+cleade+'.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)


