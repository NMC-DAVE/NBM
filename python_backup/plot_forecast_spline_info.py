"""
plot_forecast_spline_info.py cmonth clead 

plot the forecast spline parameters for the chosen month and lead

"""

import os, sys
from datetime import datetime # Jeff Whitaker's datetime library
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import scipy.stats as stats
from mpl_toolkits.basemap import Basemap, interp
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
from pyproj import Proj
import cartopy.io.shapereader as shpreader
from matplotlib import rcParams

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# ===========================================================

def read_forecast_spline_info(cmonth, clead, \
    master_directory_forecast_spline):

    """ load spline information for forecast netCDF file.  
        For initial times other than 00 UTC, we need to find the
        appropriate 00 UTC spline file to read in.   Assuming that
        the CDF characteristics are primarily a function of the 
        diurnal cycle, make this adjustment."""
    
    infile = master_directory_forecast_spline + \
        cmonth + '_conus_GEFSv12_spline_info_h' + clead + '.nc' 
    print ('reading forecast spline information from ', infile)
    nc = Dataset(infile)    
    lons_spline_gefsv12_1d = nc.variables['lons'][:]
    print ('min, max lons_spline_gefsv12_1d = ', \
        np.min(lons_spline_gefsv12_1d), \
        np.max(lons_spline_gefsv12_1d))
    lats_spline_gefsv12_1d = nc.variables['lats'][:]
    spline_info_gefsv12 = nc.variables['spline_info'][:,:,:,:]
    fraction_zero_gefsv12 = nc.variables['fzero'][:,:]
    usegamma_gefsv12 = nc.variables['usegamma'][:,:]
    quantile_98 = nc.variables['quantile_98'][:,:]
    print ('max, min quantile_98 = ', np.min(quantile_98), \
        np.max(quantile_98))
    ny_gefsv12, nx_gefsv12 = np.shape(usegamma_gefsv12)
    nc.close()
    lons_fcst_2d, lats_fcst_2d = \
        np.meshgrid(lons_spline_gefsv12_1d,lats_spline_gefsv12_1d)
        
    return lons_spline_gefsv12_1d, lats_spline_gefsv12_1d, \
        spline_info_gefsv12, fraction_zero_gefsv12, \
        usegamma_gefsv12, quantile_98, ny_gefsv12, nx_gefsv12, \
        lons_fcst_2d, lats_fcst_2d
        
# ===========================================================

# ---- inputs from command line

cmonth = sys.argv[1] # Jan, etc.
clead = sys.argv[2]  # lead time as 3-digit number, e.g., 
                     # forecast ending at +18h is 018.
                     
master_directory_forecast_spline = '/Volumes/NBM/conus_gefsv12/CDF_spline/'

lons_1d_realtime, lats_1d_realtime, \
    spline_info_gefsv12, fraction_zero_gefsv12, \
    usegamma_gefsv12, quantile_98, ny_gefsv12, nx_gefsv12, \
    lons_fcst_2d, lats_fcst_2d = \
    read_forecast_spline_info(cmonth, clead, \
    master_directory_forecast_spline)

# ---- cartopy stuff

reader = shpreader.Reader('countyl010g_shp_nt00964/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cf.ShapelyFeature(counties, ccrs.PlateCarree())

reader = shpreader.Reader('statesl010g/statesp010g.shp')
states = list(reader.geometries())
STATES = cf.ShapelyFeature(states, ccrs.PlateCarree())

# ---- plot the q98

clevs = [0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 25, 30, 35, 40]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
sectors = ['CONUS','California','Northwest_US','Central_Rockies',\
        'Northeast_US','Southeast_US','Great_Plains']
nsectors = len(sectors)
plot_CONUS = True
plot_California = False
plot_northwest = False
plot_central_rockies = True
plot_northeast = False
plot_southeast = False
plot_plains = False

print ('min_max lons_fcst_2d = ', np.min(lons_fcst_2d), np.max(lons_fcst_2d))
print ('min_max lats_fcst_2d = ', np.min(lats_fcst_2d), np.max(lats_fcst_2d))
print ('max, min quantile_98 = ', np.min(quantile_98), \
    np.max(quantile_98))
for isector,sector in enumerate(sectors):

    plotit = False
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

    proj = ccrs.LambertConformal(\
        central_latitude = (latb+late)/2.,
        central_longitude = (lonb+lone)/2,
        standard_parallels = (latb, late))

    if plotit == True:

        # ---- plot the 98th percentile

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

        title = '98th percentile of '+cmonth+'\n6-h GEFSv12 rainfall ending '+clead+' h'
        ax.set_title(title, fontsize=16,color='Black')
        CS = ax.contourf(lons_fcst_2d, lats_fcst_2d, \
            quantile_98, clevs, cmap=None, colors=colorst, \
            extend='both', transform=ccrs.PlateCarree())

        cax = fig.add_axes([0.02,0.08,0.96,0.02])
        cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=9)
        cb.set_label('98th percentile precipitation amount (mm)',fontsize=11)

        plot_title = 'q98_'+cmonth+'_'+sector+'_lead'+clead+'h.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)




