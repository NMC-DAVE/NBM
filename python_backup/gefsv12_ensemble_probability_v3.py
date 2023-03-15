"""
gefsv12_ensemble_probability_v3.py cyyyymmddhh clead cthresh
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
cthresh = sys.argv[3]  #0.254, 1.0,  5.0, 10.0, 25.0

cyear = cyyyymmddhh[0:4]
rthresh = float(cthresh)
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'
master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/thinned/'
master_directory_upscaled_output = '/Volumes/NBM/conus_gefsv12/upscaled/'
plot_fullfield = False
plot_thinned = True
plot_upscaled = False

# ---- read the stored netCDF full-field probabilities, raw 
#      and quantile mapped

if plot_fullfield == True:
    infile = master_directory_fullfield_output +\
        ccmonth+cyear+'_use99_lead'+clead+\
        '_probabilities_fullfield.nc'
    print ('reading from ', infile)
    nc = Dataset(infile,'r')
    yyyymmddhh = nc.variables['yyyymmddhh_init'][:]
    print (yyyymmddhh)
    thresholds = nc.variables['thresholdv'][:]
    print (thresholds)
    #idd = int(np.where(yyyymmddhh == int(cyyyymmddhh))[0])
    idd = np.where(yyyymmddhh == int(cyyyymmddhh))[0]
    print ('idd = ', idd)
    if rthresh == 0.254:
        itt = 0
    else:
        itt = int(np.where(thresholds == rthresh)[0])
    print ('idd,itt = ',idd,itt)
    lons_ndfd = nc.variables['lons'][:,:]
    lats_ndfd = nc.variables['lats'][:,:]
    prob_raw_fullfield = nc.variables['probability_raw'][idd,itt,:,:]
    prob_qmapped_fullfield = nc.variables['probability_qmapped'][idd,itt,:,:]
    nc.close()
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)

# ---- read the stored netCDF thinned-field probabilities, raw 
#      and quantile mapped

if plot_thinned == True:
    infile = master_directory_thinned_output+\
        ccmonth+cyear+'_use99_lead'+clead+'_probabilities_thinned.nc'
    print ('reading from ', infile)
    nc = Dataset(infile,'r')
    yyyymmddhh = nc.variables['yyyymmddhh_init'][:]
    thresholds = nc.variables['thresholdv'][:]
    idd = np.where(yyyymmddhh == int(cyyyymmddhh))[0]
    print ('idd = ', idd)
    if rthresh == 0.254:
        itt = 0
    else:
        itt = int(np.where(thresholds == rthresh)[0])
    print ('itt =',itt)
    lons_thinned = nc.variables['lons'][:,:]
    lats_thinned = nc.variables['lats'][:,:]
    prob_raw_thinned = nc.variables['probability_raw'][idd,itt,:,:]
    prob_qmapped_thinned = nc.variables['probability_qmapped'][idd,itt,:,:]
    nc.close()
    ny_thinned, nx_thinned = np.shape(lons_thinned)

# ---- read the stored netCDF thinned-field probabilities, raw 
#      and quantile mapped

if plot_upscaled == True:
    infile = master_directory_upscaled_output+ccmonth+\
        cyear+'_use99_lead'+clead+'_probabilities_upscaled.nc'
    print ('reading from ', infile)
    nc = Dataset(infile,'r')
    yyyymmddhh = nc.variables['yyyymmddhh_init'][:]
    thresholds = nc.variables['thresholdv'][:]
    idd = np.where(yyyymmddhh == int(cyyyymmddhh))[0]
    if rthresh == 0.254:
        itt = 0
    else:
        itt = int(np.where(thresholds == rthresh)[0])
    print ('itt =',itt)
    lons_upscaled = nc.variables['lons'][:,:]
    lats_upscaled = nc.variables['lats'][:,:]
    prob_raw_upscaled = nc.variables['probability_raw'][idd,itt,:,:]
    prob_qmapped_upscaled = nc.variables['probability_qmapped'][idd,itt,:,:]
    nc.close()
    ny_upscaled, nx_upscaled = np.shape(lons_upscaled)

# CONUS
#m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
#    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
#    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
#    resolution ='l',area_thresh=1000.)
   
# Western US 
m = Basemap(llcrnrlon=235.,llcrnrlat=30.,\
    urcrnrlon = 260., urcrnrlat = 45.,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
    
# Pac Northwest
#m = Basemap(llcrnrlon=233.7234,llcrnrlat=37.,\
#    urcrnrlon = 260., urcrnrlat = 54.37279,\
#    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
#    resolution ='l',area_thresh=1000.)
    
#m = Basemap(llcrnrlon=245.,llcrnrlat=35.,
#    urcrnrlon = 260., urcrnrlat = 43.,\
#    projection='lcc',lat_1=25.,lat_2=25.,lon_0=255.,\
#    resolution ='l',area_thresh=1000.)
    
# Cali
#m = Basemap(llcrnrlon=235.,llcrnrlat=31.,
#    urcrnrlon = 245., urcrnrlat = 42.,\
#    projection='lcc',lat_1=25.,lat_2=25.,lon_0=245.,\
#    resolution ='l',area_thresh=1000.)
    
#m = Basemap(llcrnrlon=230.,llcrnrlat=22.,
#    urcrnrlon = 290., urcrnrlat = 51.,\
#    projection='lcc',lat_1=25.,lat_2=25.,lon_0=255.,\
#    resolution ='l',area_thresh=1000.)

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.93,0.97,1.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
parallels = np.arange(20.,60.01,5.)
meridians = np.arange(-140.,-49.95,5.)
    

for itype in range(3):
    makeplot = False    
    if itype == 0 and plot_fullfield == True :
        lons = lons_ndfd
        lats = lats_ndfd
        prob_raw = np.squeeze(prob_raw_fullfield)
        prob_qmapped = np.squeeze(prob_qmapped_fullfield)
        ctype_out = 'Fullfield_'
        #ctype = 'Full field '
        ctype = ''
        makeplot = True
    elif itype == 1 and plot_thinned == True:
        lons = lons_thinned
        lats = lats_thinned
        prob_raw = np.squeeze(prob_raw_thinned)
        prob_qmapped = np.squeeze(prob_qmapped_thinned)
        ctype_out = 'Thinned_'
        ctype = 'Thinned '
        makeplot = True
    elif itype == 2 and plot_upscaled == True:
        lons = lons_upscaled
        lats = lats_upscaled
        prob_raw = np.squeeze(prob_raw_upscaled)
        prob_qmapped = np.squeeze(prob_qmapped_upscaled)
        ctype_out = 'Upscaled_' 
        ctype = 'Upscaled ' 
        makeplot = True
   
    if makeplot == True:
        x, y = m(lons, lats)
    
        # ---- plot the raw probability
    
        #fig = plt.figure(figsize=(9,6.5))
        fig = plt.figure(figsize=(9,6.5))
        #fig = plt.figure(figsize=(6.,7.7)) # cali
        #axloc = [0.06,0.15,0.92,0.8]
        axloc = [0.06,0.11,0.92,0.8]
        ax1 = fig.add_axes(axloc)
        cleadb = str(int(clead)-6)
        title = ctype+' GEFSv12 P(obs > '+cthresh+' mm), IC = '+\
            cyyyymmddhh+', lead = '+clead+' h'
        ax1.set_title(title, fontsize=15,color='Black')
        print ('np.shape(prob_raw) = ', np.shape(prob_raw))
        CS2 = m.contourf(x, y, prob_raw[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')

        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
        m.drawcounties(linewidth=0.3,color='LightGray')
    
        #m.drawparallels(parallels,labels=[1,0,0,0],\
        #    color='Gray',fontsize=8,linewidth=0.2)  # uncomment for lat/lon lines
        #meridians = np.arange(220.,360.,5.)
        #m.drawmeridians(meridians,labels=[0,0,0,1],\
        #    color='Gray',fontsize=8,linewidth=0.2)
    
        cax = fig.add_axes([0.06,0.06,0.92,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Probability',fontsize=9)

        plot_title = ctype_out+'raw_probability_'+cthresh+'mm_'+\
            cyyyymmddhh+'_lead'+clead+'.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)

        # ---- plot the quantile-mapped probability
    
        #fig = plt.figure(figsize=(9,6.5))
        fig = plt.figure(figsize=(9.,6.5))
        #fig = plt.figure(figsize=(6.,7.7))
        #axloc = [0.06,0.15,0.92,0.8]
        axloc = [0.06,0.11,0.92,0.8]
        ax1 = fig.add_axes(axloc)
        cleadb = str(int(clead)-6)
        #title = ctype+' Quantile-mapped GEFSv12 P(obs > '+\
        #    cthresh+' mm), IC = '+cyyyymmddhh+', lead = '+clead+' h'
        title = ctype+' Quantile-mapped GEFSv12 P(obs > '+\
            cthresh+' mm),\nIC = '+cyyyymmddhh+', lead = '+clead+' h'
        ax1.set_title(title, fontsize=15,color='Black')
        CS2 = m.contourf(x, y, prob_qmapped[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')
    
        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
        m.drawcounties(linewidth=0.3,color='LightGray')
    
        #m.drawparallels(parallels,labels=[1,0,0,0],\
        #    color='Gray',fontsize=8,linewidth=0.2)
        #meridians = np.arange(220.,360.,5.)
        #m.drawmeridians(meridians,labels=[0,0,0,1],\
        #    color='Gray',fontsize=8,linewidth=0.2)
    
        cax = fig.add_axes([0.06,0.06,0.88,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Probability',fontsize=9)

        plot_title = ctype_out+'qmapped_probability_'+cthresh+'mm_'+\
            cyyyymmddhh+'_lead'+clead+'.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)

    
