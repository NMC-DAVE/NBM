"""
python gefsv12_ensemble_probability_weighted.py cyyyymmddhh clead cthresh
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
jnd = 1350
ind = 1110

cyear = cyyyymmddhh[0:4]
rthresh = float(cthresh)
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
#master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'
master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/thinned/'
dressed = False
plot_fullfield = True
weighted = True
makeplot = True

# ---- read the stored netCDF full-field probabilities, quantile mapped

#if weighted == True:
#    if dressed == True:
#        infile = master_directory_fullfield_output +\
#            ccmonth+cyear+'_use99_lead'+clead+\
#            '_probabilities_weighted_dressed_fullfield.nc'
#    else:
#        infile = master_directory_fullfield_output +\
#            ccmonth+cyear+'_use99_lead'+clead+\
#            '_probabilities_weighted_fullfield.nc'
#else:
#    infile = master_directory_fullfield_output +\
#        ccmonth+cyear+'_use99_lead'+clead+\
#        '_probabilities_fullfield.nc'
        
if weighted == True:
    if dressed == True:
        infile = master_directory_thinned_output +\
            ccmonth+cyear+'_use99_lead'+clead+\
            '_probabilities_weighted_dressed_thinned.nc'
    else:
        infile = master_directory_thinned_output +\
            ccmonth+cyear+'_use99_lead'+clead+\
            '_probabilities_weighted_thinned.nc'
else:
    infile = master_directory_thinned_output +\
        ccmonth+cyear+'_use99_lead'+clead+\
        '_probabilities_thinned.nc'
        
print ('reading from ', infile)
nc = Dataset(infile,'r')
yyyymmddhh = nc.variables['yyyymmddhh_init'][:]
idd = np.where(yyyymmddhh == int(cyyyymmddhh))[0]
print ('yyyymmddhh = ', yyyymmddhh)
thresholds = nc.variables['thresholdv'][:]
if rthresh == 0.254:
    itt = 0
else:
    itt = int(np.where(thresholds == rthresh)[0])
print ('itt = ', itt)
lons_ndfd = nc.variables['lons'][:,:]
lats_ndfd = nc.variables['lats'][:,:]
prob_qmapped_fullfield = \
    np.squeeze(nc.variables['probability_qmapped'][idd,itt,:,:])
nc.close()
ny_ndfd, nx_ndfd = np.shape(lons_ndfd)



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
lons = lons_ndfd
lats = lats_ndfd
prob_qmapped = np.squeeze(prob_qmapped_fullfield)
ctype = ''


if makeplot == True:
    x, y = m(lons, lats)

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
    if weighted == True:
        if dressed == True:
            title = 'Weighted, dressed quantile-mapped GEFSv12 P(obs > '+\
                cthresh+' mm),\nIC = '+cyyyymmddhh+', lead = '+clead+' h'
            plot_title = 'qmapped_chist_dressed_probability_'+cthresh+'mm_'+\
                cyyyymmddhh+'_lead'+clead+'.png'
        else:
            title = 'Weighted, quantile-mapped GEFSv12 P(obs > '+\
                cthresh+' mm),\nIC = '+cyyyymmddhh+', lead = '+clead+' h'
            plot_title = 'qmapped_chist_probability_'+cthresh+'mm_'+\
                cyyyymmddhh+'_lead'+clead+'.png' 
    else:
        title = 'Quantile-mapped GEFSv12 P(obs > '+\
            cthresh+' mm),\nIC = '+cyyyymmddhh+', lead = '+clead+' h'
        plot_title = 'qmapped_probability_'+cthresh+'mm_'+\
            cyyyymmddhh+'_lead'+clead+'.png'

        
    ax1.set_title(title, fontsize=14,color='Black')
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

    cax = fig.add_axes([0.06,0.06,0.92,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Probability',fontsize=9)


    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)




#