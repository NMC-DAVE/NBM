"""
gefsv12_ensemble_probability_v2.py cyyyymmddhh clead 
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
import pygrib
#from quantile_mapping_gamma_mixture_v2_f90 import \
#    quantile_mapping_gamma_mixture_v2_f90
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
import _pickle as cPickle
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
iyyyymmddhh = int(cyyyymmddhh)
clead = sys.argv[2] # 018 not 18
cyyyy = cyyyymmddhh[0:4]
use98 = True
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
#master_directory_probability_output = \
#    '/Volumes/NBM/conus_gefsv12/probabilities/'
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
nthresh = len(thresholds)

# ---- read the stored netCDF probabilities, raw 
#      and quantile mapped

if use98 == False:
    infile = master_directory_fullfield_output + ccmonth + cyyyy + \
        '_lead'+clead+'_probabilities_fullfield.nc'
else:
    infile = master_directory_fullfield_output + ccmonth + cyyyy + \
        '_use98_lead'+clead+'_probabilities_fullfield.nc'
    
print ('reading from ', infile)
nc = Dataset(infile,'r')
lons_ndfd = nc.variables['lons'][:,:]
lats_ndfd = nc.variables['lats'][:,:]
yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
idx = np.where(yyyymmddhh_init == iyyyymmddhh)[0]
prob_raw = np.squeeze(nc.variables['probability_raw'][idx,:,:,:])
prob_qmapped = np.squeeze(nc.variables['probability_qmapped'][idx,:,:,:])

thresholds_in = nc.variables['thresholdv'][:]
nthresh = len(thresholds_in)
nc.close()
ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
print (np.shape(prob_raw))
print (np.shape(prob_qmapped))
print (thresholds_in)

zeros = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
ones = np.ones((ny_ndfd, nx_ndfd), dtype=np.float64)


clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.93,0.96,1.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

thresholds_in = [0.254, 1.0, 5.0]
for ithresh, thresh in enumerate(thresholds_in):
    
    
    m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
        urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
        projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
        resolution ='i',area_thresh=1000.)
    x, y = m(lons_ndfd, lats_ndfd)
    
    # ---- plot the raw probability
    
    fig = plt.figure(figsize=(9,6.5))
    axloc = [0.06,0.15,0.92,0.8]
    ax1 = fig.add_axes(axloc)
    cleadb = str(int(clead)-6)
    title = 'Raw GEFSv12 P(obs > '+str(thresh)+' mm), IC = '+\
        cyyyymmddhh+', lead = '+clead+' h'
    ax1.set_title(title, fontsize=14,color='Black')
    CS2 = m.contourf(x, y, prob_raw[ithresh,:,:], clevs,\
        cmap=None, colors=colorst, extend='both')
    
    parallels = np.arange(20.,60.01,5.)
    meridians = np.arange(-140.,-49.95,5.)
    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    
    m.drawparallels(parallels,labels=[1,0,0,0],color='Gray',fontsize=8,linewidth=0.2)
    meridians = np.arange(220.,360.,5.)
    m.drawmeridians(meridians,labels=[0,0,0,1],color='Gray',fontsize=8,linewidth=0.2)
    
    cax = fig.add_axes([0.06,0.07,0.88,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Probability',fontsize=9)

    plot_title = 'raw_probability_'+str(thresh)+'mm_'+\
        cyyyymmddhh+'_lead'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)
    #print ('Done!')

    # ---- plot the quantile-mapped probability
    
    fig = plt.figure(figsize=(9,6.5))
    axloc = [0.06,0.15,0.92,0.8]
    ax1 = fig.add_axes(axloc)
    cleadb = str(int(clead)-6)
    title = 'Quantile-mapped GEFSv12 P(obs > '+str(thresh)+' mm), IC = '+\
        cyyyymmddhh+', lead = '+clead+' h'
    ax1.set_title(title, fontsize=14,color='Black')
    CS2 = m.contourf(x, y, prob_qmapped[ithresh,:,:], clevs,\
        cmap=None, colors=colorst, extend='both')
    
    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    
    m.drawparallels(parallels,labels=[1,0,0,0],color='Gray',fontsize=8,linewidth=0.2)
    meridians = np.arange(220.,360.,5.)
    m.drawmeridians(meridians,labels=[0,0,0,1],color='Gray',fontsize=8,linewidth=0.2)
    
    cax = fig.add_axes([0.06,0.07,0.88,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Probability',fontsize=9)

    plot_title = 'qmapped_probability_'+str(thresh)+'mm_'+\
        cyyyymmddhh+'_lead'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)

    
    sectors = ['California','Northwest_US','Central_Rockies',\
        'Northeast_US','Southeast_US','Great_Plains']
    nsectors = len(sectors)
    plot_California = False
    plot_northwest = False
    plot_central_rockies = True
    plot_northeast = False
    plot_southeast = False
    plot_plains = False

    for isector,sector in enumerate(sectors):
        plotit = False 
        if isector == 0 and plot_California == True:
            plotit = True
            latb = 31.
            late = 43.
            lonb = -125.
            lone = -112.
            xdim = 6.
            ydim = 8.
            drawcoasts = True
            area_thresh = 1000
        if isector == 1 and plot_northwest == True:
            plotit = True
            latb = 40.
            late = 51.
            lonb = -126.
            lone = -105.
            xdim = 8.
            ydim = 6.
            drawcoasts = True
            area_thresh = 1000
        if isector == 2 and plot_central_rockies == True:
            plotit = True
            latb = 33.
            late = 47.
            lonb = -117.
            lone = -101.9
            xdim = 6.
            ydim = 8.
            drawcoasts = False
            area_thresh = 0
        if isector == 3 and plot_northeast == True:
            plotit = True
            latb = 36.
            late = 49.
            lonb = -90.
            lone = -65.
            xdim = 9.
            ydim = 6.5
            drawcoasts = True
            area_thresh = 1000
        if isector == 4 and plot_southeast == True:
            plotit = True
            latb = 24.
            late = 38.
            lonb = -95.
            lone = -75.
            xdim = 6.5
            ydim = 6.
            drawcoasts = True
            area_thresh = 1000
        if isector == 5 and plot_plains == True:
            plotit = True
            latb = 25.
            late = 50.
            lonb = -108.
            lone = -88.
            xdim = 6.5
            ydim = 9.
            drawcoasts = True
            area_thresh = 1000

            
        if plotit == True:
            print (lonb, latb, lone, late)
            m = Basemap(llcrnrlon=lonb,llcrnrlat=latb,
                urcrnrlon = lone, urcrnrlat = late,\
                projection='mill',resolution ='l')
            x, y = m(lons_ndfd, lats_ndfd)
    
            # ---- plot the raw probability
    
            fig = plt.figure(figsize=(xdim,ydim))
            axloc = [0.02,0.1,0.96,0.81]
            ax1 = fig.add_axes(axloc)
            cleadb = str(int(clead)-6)
            title = 'Raw GEFSv12 P(obs > '+str(thresh)+' mm),\nIC = '+\
                cyyyymmddhh+' lead = '+clead+' h'
            ax1.set_title(title, fontsize=14,color='Black')
            print ('max prob_raw = ', np.max(prob_raw[ithresh,:,:]))
            print ('prob_raw[ithresh,ny_ndfd//2, 0:nx_ndfd:10] = ',prob_raw[ithresh,ny_ndfd//2, 0:nx_ndfd:10])
            CS2 = m.contourf(x, y, prob_raw[ithresh,:,:], clevs,\
                cmap=None, colors=colorst, extend='both')
    
            if drawcoasts == True: m.drawcoastlines(linewidth=0.8,color='Gray')
            m.drawcountries(linewidth=0.8,color='Gray')
            m.drawstates(linewidth=0.8,color='Gray')
            m.drawcounties(linewidth=0.1,color='Gray')
    
            cax = fig.add_axes([0.08,0.07,0.84,0.02])
            cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
                drawedges=True,ticks=clevs,format='%g')
            cb.ax.tick_params(labelsize=7)
            cb.set_label('Probability',fontsize=9)

            plot_title = 'raw_probability_'+sector+'_'+str(thresh)+'mm_'+\
                cyyyymmddhh+'_lead'+clead+'.png'
            fig.savefig(plot_title, dpi=300)
            plt.close()
            print ('saving plot to file = ',plot_title)
            
            
            # ---- plot the quantile-mapped probability
    
            fig = plt.figure(figsize=(xdim,ydim))
            axloc = [0.02,0.1,0.96,0.81]
            ax1 = fig.add_axes(axloc)
            cleadb = str(int(clead)-6)
            title = 'Quantile-mapped GEFSv12 P(obs > '+str(thresh)+' mm),\nIC = '+\
                cyyyymmddhh+' lead = '+clead+' h'
            ax1.set_title(title, fontsize=14,color='Black')
            print ('ithresh = ', ithresh)
            print ('max prob_qmapped = ', np.max(prob_raw[ithresh,:,:]))
            print ('prob_qmapped[ithresh,ny_ndfd//2, 0:nx_ndfd:10] = ',prob_qmapped[ithresh,ny_ndfd//2, 0:nx_ndfd:10])
            CS2 = m.contourf(x, y, prob_qmapped[ithresh,:,:], clevs,\
                cmap=None, colors=colorst, extend='both')
    
            if drawcoasts == True: m.drawcoastlines(linewidth=0.8,color='Gray')
            m.drawcountries(linewidth=0.8,color='Gray')
            m.drawstates(linewidth=0.8,color='Gray')
            m.drawcounties(linewidth=0.1,color='Gray',zorder=20)
    
            cax = fig.add_axes([0.08,0.07,0.84,0.02])
            cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
                drawedges=True,ticks=clevs,format='%g')
            cb.ax.tick_params(labelsize=7)
            cb.set_label('Probability',fontsize=9)

            if use98 == True:
                plot_title = 'qmapped_probability_use98_'+sector+'_'+str(thresh)+'mm_'+\
                    cyyyymmddhh+'_lead'+clead+'.png'
            else:
                plot_title = 'qmapped_probability_'+sector+'_'+str(thresh)+'mm_'+\
                    cyyyymmddhh+'_lead'+clead+'.png'    
            fig.savefig(plot_title, dpi=300)
            plt.close()
            print ('saving plot to file = ',plot_title)
            #print ('Done!')



