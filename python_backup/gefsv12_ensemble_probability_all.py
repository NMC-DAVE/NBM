"""
python gefsv12_ensemble_probability_all.py cyyyymmddhh clead cthresh
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
master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/all_thinned/'
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'

# ---- read the stored netCDF thinned probabilities

#infile = master_directory_thinned_output +\
#    ccmonth + cyear+'_lead'+clead+\
#            '_probabilities_all_thinned.nc'
            
infile = master_directory_fullfield_output +\
    ccmonth + cyear+'_lead'+clead+\
            '_probabilities_all_fullfield.nc'
            
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
probability_raw = np.squeeze(nc.variables['probability_raw'][idd,itt,:,:])
probability_qmapped = np.squeeze(nc.variables[\
    'probability_qmapped'][idd,itt,:,:])
probability_qmapped_weighted = np.squeeze(\
    nc.variables['probability_qmapped_weighted'][idd,itt,:,:])
probability_qmapped_weighted_dressed = np.squeeze(\
    nc.variables['probability_qmapped_weighted_dressed'][idd,itt,:,:])
nc.close()
ny_ndfd, nx_ndfd = np.shape(lons_ndfd)

diff = probability_qmapped_weighted_dressed - probability_qmapped_weighted
print ('max diff = ', np.max(diff))
print ('min diff = ', np.min(diff))

# CONUS
#m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
#    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
#    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
#    resolution ='l',area_thresh=1000.)
    
# Western US
m = Basemap(llcrnrlon=235.,llcrnrlat=25.,\
    urcrnrlon = 260., urcrnrlat = 50.,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=250.,\
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
lons = lons_ndfd
lats = lats_ndfd
    
x, y = m(lons, lats)

for itype in range(4):
    
    if itype == 0:
        prob = probability_raw
        cftype = '(a) Raw '
        cptype = 'raw_'
    elif itype == 1:
        prob = probability_qmapped
        cftype = '(b) Quantile-mapped '
        cptype = 'qmapped_'
    elif itype == 2:
        prob = probability_qmapped_weighted
        cftype = '(c) Weighted, quantile-mapped\n'
        cptype = 'qmapped_weighted_'
    elif itype == 3:
        prob = probability_qmapped_weighted_dressed
        cftype = '(d) Weighted, quantile-mapped,\ndressed '
        cptype = 'qmapped_weighted_dressed_'

    # ---- plot the quantile-mapped probability

    #fig = plt.figure(figsize=(9.,6.5))
    fig = plt.figure(figsize=(5.,7.4))
    #axloc = [0.06,0.15,0.92,0.8]
    axloc = [0.04,0.09,0.92,0.8]
    ax1 = fig.add_axes(axloc)
    cleadb = str(int(clead)-6)
    title = cftype + 'GEFSv12 P(obs > '+\
        cthresh+' mm),\nIC = '+cyyyymmddhh+', lead = '+clead+' h'
    plot_title = cptype+cthresh+'mm_'+cyyyymmddhh+'_lead'+clead+'.png'
        
    ax1.set_title(title, fontsize=14,color='Black')
    CS2 = m.contourf(x, y, prob[:,:], clevs,\
        cmap=None, colors=colorst, extend='both')

    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    m.drawcounties(linewidth=0.3,color='LightGray')

    cax = fig.add_axes([0.04,0.06,0.92,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Probability',fontsize=9)

    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)




#