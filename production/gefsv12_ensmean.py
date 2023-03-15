"""
gefsv12_ensmean.py cyyyymmddhh clead 
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

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2] # 018 not 18
cyear = cyyyymmddhh[0:4]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'

# ---- read the stored netCDF ensemble mean, raw 
#      and quantile mapped

infile = master_directory_fullfield_output +\
    ccmonth+cyear+'_use99_lead'+clead+\
    '_probabilities_fullfield.nc'
print ('reading from ', infile)
nc = Dataset(infile,'r')
yyyymmddhh = nc.variables['yyyymmddhh_init'][:]
print (yyyymmddhh)
thresholds = nc.variables['thresholdv'][:]
print (thresholds)
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

ensmean_raw = ensmean_raw_fullfield
ensmean_qmapped = ensmean_qmapped_fullfield
print ('np.shape(ensmean_raw) = ', np.shape(ensmean_raw) )
print ('np.shape(ensmean_qmapped) = ', np.shape(ensmean_qmapped) )
print ('np.shape(lons_ndfd) = ', np.shape(lons_ndfd))
print ('np.shape(lats_ndfd) = ', np.shape(lats_ndfd))

# ---- plot the raw ensemble mean by sector


m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ndfd, lats_ndfd)

clevs = [0.0,0.254, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 35.0, 50.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
sectors = ['CONUS','California','Northwest_US','Central_Rockies',\
        'Northeast_US','Southeast_US','Great_Plains']
nsectors = len(sectors)
plot_CONUS = True
plot_California = False
plot_northwest = False
plot_central_rockies = False
plot_northeast = False
plot_southeast = False
plot_plains = True

for isector,sector in enumerate(sectors):

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
        latb = 38.
        late = 51.
        lonb = 233.7234 - 360.
        lone = -100.
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
        plotit = True
        latb = 25.
        late = 50.
        lonb = -108.
        lone = -88.
        xdim = 6.5
        ydim = 9.
        drawcoasts = True


    if plotit == True:

        print ('lonb, latb, lone, late = ',lonb, latb, lone, late )
        print ('np.max(lons_ndfd), np.min(lons_ndfd) = ', \
            np.max(lons_ndfd), np.min(lons_ndfd))
        m = Basemap(llcrnrlon=lonb,llcrnrlat=latb,
            urcrnrlon = lone, urcrnrlat = late,\
            projection='mill',resolution ='i')
        x, y = m(lons_ndfd, lats_ndfd)

        # ---- plot the raw ensemble mean

        fig = plt.figure(figsize=(xdim,ydim))
        axloc = [0.02,0.1,0.96,0.81]
        ax1 = fig.add_axes(axloc)
        cleadb = str(int(clead)-6)
        title = 'Raw GEFSv12 ensemble mean,\nIC = '+\
            cyyyymmddhh+' lead = '+clead+' h'
        ax1.set_title(title, fontsize=14,color='Black')
        CS2 = m.contourf(x, y, ensmean_raw[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')

        if drawcoasts == True: m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
        m.drawcounties(linewidth=0.1,color='Gray')

        cax = fig.add_axes([0.08,0.07,0.84,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Precipitation amount (mm)',fontsize=9)

        plot_title = 'raw_ensmean_'+sector+'_'+\
            cyyyymmddhh+'_lead'+clead+'.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)


        # ---- plot the quantile-mapped ensemble mean

        fig = plt.figure(figsize=(xdim,ydim))
        axloc = [0.02,0.1,0.96,0.81]
        ax1 = fig.add_axes(axloc)
        cleadb = str(int(clead)-6)
        title = 'Quantile-mapped GEFSv12 ensemble mean,\nIC = '+\
            cyyyymmddhh+' lead = '+clead+' h'
        ax1.set_title(title, fontsize=14,color='Black')
        CS2 = m.contourf(x, y, ensmean_qmapped[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')

        if drawcoasts == True: m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
        m.drawcounties(linewidth=0.1,color='Gray',zorder=20)

        cax = fig.add_axes([0.08,0.07,0.84,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Precipitation amount (mm)',fontsize=9)

        plot_title = 'qmapped_ensmean_'+sector+'_'+\
            cyyyymmddhh+'_lead'+clead+'.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)
        #print ('Done!')



    
