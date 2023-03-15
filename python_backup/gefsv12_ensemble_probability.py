"""
gefsv12_ensemble_probability.py cyyyymmddhh clead 

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
from quantile_mapping_gamma_mixture_v2_f90 import \
    quantile_mapping_gamma_mixture_v2_f90
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import LSQUnivariateSpline, splrep, splev

import _pickle as cPickle
import scipy.stats as stats


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'



def set_domain_boundaries(cdomain):

    """ used grib file of 2.5-km blend output grid to determine bounding
        lat and lon, and from that, the domain bounding indices for the
        0.25 GEFSv12 reforecast data that will encompass the domain.
    """
    if cdomain == 'conus':
        jmin = 93
        jmax = 246
        imin = 368
        imax = 686
    elif cdomain == 'pr':
        jmin = 243
        jmax = 256
        imin = 649
        imax = 667
    elif cdomain == 'ak':
        jmin = 19
        jmax = 161
        imin = 201
        imax = 967
    else:
        print ('invalid domain.  Exiting.')
        sys.exit()

    return jmin, jmax, imin, imax


# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2] # 18 not 18
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
jmin, jmax, imin, imax = set_domain_boundaries(cdomain)

cmembers = ['c00','p01', 'p02','p03','p04','p05','p06','p07','p08','p09','p10',\
    'p11', 'p12','p13','p14','p15','p16','p17','p18','p19','p20',\
    'p21', 'p22','p23','p24','p25','p26','p27','p28','p29','p30']
nmembers = len(cmembers)
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
nthresh = len(thresholds)

# ---- read in the previously generated netCDF file with precipitation
#      to get lat/lon of NDFD grid

ccpa_directory = '/Volumes/Backup Plus/ccpa/'
infile = ccpa_directory + '200201_ccpa_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile)
lons_ndfd = nc.variables['lons'][:,:]
lons_ndfd = lons_ndfd
lats_ndfd = nc.variables['lats'][:,:]
ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
nc.close()

zeros = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
ones = np.ones((ny_ndfd, nx_ndfd), dtype=np.float64)

sum_binary_raw = np.zeros((nthresh, ny_ndfd, nx_ndfd), dtype=np.float64)
sum_binary_qmap = np.zeros((nthresh, ny_ndfd, nx_ndfd), dtype=np.float64)
probability_raw = np.zeros((nthresh, ny_ndfd, nx_ndfd), dtype=np.float64)
probability_qmap = np.zeros((nthresh, ny_ndfd, nx_ndfd), dtype=np.float64)

for cmem in cmembers: 
    
    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print ('processing member = ',cmem,time)
    
    # ---- read the quantile mapped array from file
    
    input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
    infile = input_directory + cyyyymmddhh + \
        '_qmap_'+cmem+'.t00z.pgrb2s.0p25.f' + clead+'.cPick'
    print (infile)
    inf = open(infile,'rb')
    precip_realtime = cPickle.load(inf)
    qmapped_precip = cPickle.load(inf)
    inf.close()
  
    for ithresh, thresh in enumerate(thresholds):
        binary_map_raw = np.where(precip_realtime > thresh, ones, zeros)
        binary_map_qmap = np.where(qmapped_precip > thresh, ones, zeros)
        sum_binary_raw[ithresh,:,:] = sum_binary_raw[ithresh,:,:] + \
            binary_map_raw[:,:]
        sum_binary_qmap[ithresh,:,:] = sum_binary_qmap[ithresh,:,:] + \
            binary_map_qmap[:,:]
            
# ---- finalize probabilities

for ithresh, thresh in enumerate(thresholds):
    probability_raw[ithresh,:,:] = \
        sum_binary_raw[ithresh,:,:] / float(nmembers)
    probability_qmap[ithresh,:,:] = \
        sum_binary_qmap[ithresh,:,:] / float(nmembers)
    

    m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
        urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
        projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
        resolution ='l',area_thresh=1000.)
    x, y = m(lons_ndfd, lats_ndfd)

    clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.0]
    colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
        '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
        '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
    # ---- plot the raw probability
    
    fig = plt.figure(figsize=(9,6.5))
    axloc = [0.02,0.1,0.96,0.84]
    ax1 = fig.add_axes(axloc)
    cleadb = str(int(clead)-6)
    title = 'Raw GEFSv12 P(obs > '+str(thresh)+' mm), IC = '+\
        cyyyymmddhh+', lead = '+clead+' h'
    ax1.set_title(title, fontsize=14,color='Black')
    CS2 = m.contourf(x, y, probability_raw[ithresh,:,:], clevs,\
        cmap=None, colors=colorst, extend='both')
    
    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    
    cax = fig.add_axes([0.06,0.07,0.88,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Probability',fontsize=9)

    plot_title = 'raw_probability_'+str(thresh)+'mm_'+\
        cyyyymmddhh+'_lead'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    print ('Done!')

    # ---- plot the quantile-mapped probability
    
    fig = plt.figure(figsize=(9,6.5))
    axloc = [0.02,0.1,0.96,0.84]
    ax1 = fig.add_axes(axloc)
    cleadb = str(int(clead)-6)
    title = 'Quantile-mapped GEFSv12 P(obs > '+str(thresh)+' mm), IC = '+\
        cyyyymmddhh+', lead = '+clead+' h'
    ax1.set_title(title, fontsize=14,color='Black')
    CS2 = m.contourf(x, y, probability_qmap[ithresh,:,:], clevs,\
        cmap=None, colors=colorst, extend='both')
    
    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    
    cax = fig.add_axes([0.06,0.07,0.88,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Probability',fontsize=9)

    plot_title = 'qmapped_probability_'+str(thresh)+'mm_'+\
        cyyyymmddhh+'_lead'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    print ('Done!')
    
    
    
    # ---- colorado
    
    #m = Basemap(llcrnrlon=-114.,llcrnrlat=36.0,
    #    urcrnrlon = -103., urcrnrlat = 43.0,\
    #    projection='mill',resolution ='l')
        
    # ---- california
    m = Basemap(llcrnrlon=-127.,llcrnrlat=32.0,
        urcrnrlon = -110., urcrnrlat = 49.0,\
        projection='mill',resolution ='l')
    x, y = m(lons_ndfd, lats_ndfd)

    clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.0]
    colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
        '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
        '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
    # ---- plot the raw probability
    
    #fig = plt.figure(figsize=(9,6.5))
    fig = plt.figure(figsize=(6,6.5))
    axloc = [0.02,0.1,0.96,0.81]
    ax1 = fig.add_axes(axloc)
    cleadb = str(int(clead)-6)
    title = 'Raw GEFSv12 P(obs > '+str(thresh)+' mm),\nIC = '+\
        cyyyymmddhh+' lead = '+clead+' h'
    ax1.set_title(title, fontsize=14,color='Black')
    CS2 = m.contourf(x, y, probability_raw[ithresh,:,:], clevs,\
        cmap=None, colors=colorst, extend='both')
    
    m.drawcoastlines(linewidth=0.8,color='Gray')
    #m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    #m.drawcounties(linewidth=0.8,color='Gray')
    
    cax = fig.add_axes([0.08,0.07,0.84,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Probability',fontsize=9)

    #plot_title = 'raw_probability_mountains_'+str(thresh)+'mm_'+\
    #    cyyyymmddhh+'_lead'+clead+'.png'
    plot_title = 'raw_probability_westcoast_'+str(thresh)+'mm_'+\
        cyyyymmddhh+'_lead'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    print ('Done!')

    # ---- plot the quantile-mapped probability
    
    fig = plt.figure(figsize=(6,6.5))
    axloc = [0.02,0.1,0.96,0.81]
    ax1 = fig.add_axes(axloc)
    cleadb = str(int(clead)-6)
    title = 'Quantile-mapped GEFSv12 P(obs > '+str(thresh)+' mm),\nIC = '+\
        cyyyymmddhh+' lead = '+clead+' h'
    ax1.set_title(title, fontsize=14,color='Black')
    CS2 = m.contourf(x, y, probability_qmap[ithresh,:,:], clevs,\
        cmap=None, colors=colorst, extend='both')
    
    m.drawcoastlines(linewidth=0.8,color='Gray')
    #m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    #m.drawcounties(linewidth=0.8,color='Gray')
    
    cax = fig.add_axes([0.08,0.07,0.84,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Probability',fontsize=9)

    plot_title = 'qmapped_probability_westcoast_'+str(thresh)+'mm_'+\
        cyyyymmddhh+'_lead'+clead+'.png'
    #plot_title = 'qmapped_probability_mountains_'+str(thresh)+'mm_'+\
    #    cyyyymmddhh+'_lead'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    print ('Done!')



