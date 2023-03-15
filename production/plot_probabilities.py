"""
python plot_probabilities.py cyyyymmddhh clead cthresh ctype
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
ctype = sys.argv[4] # fullfield, thinned
fullfield = False
thinned = False
if ctype == 'fullfield': fullfield = True
if ctype == 'thinned' : thinned = True
jnd = 1350
ind = 1110
ileadb = int(clead)-6
if ileadb < 10:
    cleadb = '00'+str(ileadb)
elif ileadb >= 10 and ileadb < 100:
    cleadb = '0'+str(ileadb)
else:
    cleadb = str(ileadb)
    
cyear = cyyyymmddhh[0:4]
rthresh = float(cthresh)
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'
    # where the full field probabilities are stored, if netcdf_fullfield == True
master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/all_thinned/'
    # where the thinned field probabilities are stored, if netcdf_thinned == True

# ---- read the stored netCDF full-field probabilities, quantile mapped

if fullfield == True:
    infile = master_directory_fullfield_output +\
        ccmonth+cyear+'_lead'+clead+\
            '_probabilities_all_fullfield.nc'
if thinned == True:
    infile = master_directory_thinned_output +\
        ccmonth+cyear+'_lead'+clead+\
            '_probabilities_all_thinned.nc'
print ('reading from ', infile)
nc = Dataset(infile,'r')
yyyymmddhh = nc.variables['yyyymmddhh_init'][:]
idd = np.where(yyyymmddhh == int(cyyyymmddhh))[0]
print ('yyyymmddhh = ', yyyymmddhh, int(cyyyymmddhh))
thresholds = nc.variables['thresholdv'][:]
if rthresh == 0.254:
    itt = 0
else:
    itt = int(np.where(thresholds == rthresh)[0])
print ('itt = ', itt)
print ('idd = ', idd)
lons_ndfd = nc.variables['lons'][:,:]
lats_ndfd = nc.variables['lats'][:,:]
probability_raw = nc.variables['probability_raw'][idd,itt,:,:][0]
probability_qmapped = nc.variables['probability_qmapped'][idd,itt,:,:]
probability_qmapped_weighted = nc.variables['probability_qmapped_weighted'][idd,itt,:,:][0]
probability_qmapped_weighted_dressed = nc.variables['probability_qmapped_weighted_dressed'][idd,itt,:,:][0]
nc.close()
ny_ndfd, nx_ndfd = np.shape(lons_ndfd)

# CONUS
llcrnrlon=235.
llcrnrlat=20.
urcrnrlon = 295.
urcrnrlat = 51.
drawtheborders = True
makeplot = True
    
# Western US
#llcrnrlon=235.
#llcrnrlat=30.
#urcrnrlon = 260.
#urcrnrlat = 45.
#drawtheborders = True

# Pac Northwest
#llcrnrlon=233.7234
#llcrnrlat=37.
#urcrnrlon = 260.
#urcrnrlat =  54.37279
#drawtheborders = True

#Colorado?
#llcrnrlon=245
#llcrnrlat=35.
#urcrnrlon = 260.
#urcrnrlat =  43.
#drawtheborders = False

# Cali
#llcrnrlon=235.
#llcrnrlat=31.
#urcrnrlon = 245.
#urcrnrlat =  42.
# drawtheborders = False

m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,
    urcrnrlon = urcrnrlon, urcrnrlat = urcrnrlat,\
    projection='lcc',lat_1=25.,lat_2=25., lon_0=(llcrnrlon+urcrnrlon)/2.,\
    resolution ='l',area_thresh=1000.)

plot_ratio = (urcrnrlon -llcrnrlon) / (urcrnrlat - llcrnrlat)

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.93,0.97,1.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
parallels = np.arange(20.,60.01,5.)
meridians = np.arange(-140.,-49.95,5.)

if makeplot == True:
    x, y = m(lons_ndfd, lats_ndfd)

    # ---- plot the quantile-mapped probability

    for itype in range(4):
        if itype == 0:
            plot_type = 'raw'
            title_header = 'Raw ensemble '
            prob = np.squeeze(probability_raw)
        elif itype == 1:
            plot_type = 'qmapped_unweighted'
            title_header = 'Quantile-mapped, unweighted '
            prob = np.squeeze(probability_qmapped)
        if itype == 2:
            plot_type = 'qmapped_weighted'
            title_header = 'Quantile-mapped, weighted '
            prob = np.squeeze(probability_qmapped_weighted)
        if itype == 3:
            plot_type = 'qmapped_weighted_dressed'
            title_header = 'Quantile-mapped, weighted, and dressed '
            prob = np.squeeze(probability_qmapped_weighted_dressed)
            
        plot_title = plot_type + '_linux_probability_'+cthresh+'mm_'+\
            cyyyymmddhh+'_lead'+clead+'.png'
        title = title_header + 'probability of > '+\
            cthresh+' mm,\nIC = '+cyyyymmddhh+\
            ', lead = '+cleadb+' to '+clead+' h'    
        print(np.shape(prob))
            
        fig = plt.figure(figsize=(7.,1.3+7./plot_ratio))
        axloc = [0.04,0.12,0.92,0.79]
        ax1 = fig.add_axes(axloc)
        cleadb = str(int(clead)-6)
        
        ax1.set_title(title, fontsize=22/plot_ratio,color='Black')
        CS2 = m.contourf(x, y, prob, clevs,\
        cmap=None, colors=colorst, extend='both')

        if drawtheborders == True:
            m.drawcoastlines(linewidth=0.5,color='Gray')
            m.drawcountries(linewidth=0.5,color='Gray')
        m.drawstates(linewidth=0.3,color='Gray')
        m.drawcounties(linewidth=0.15,color='LightGray')

        cax = fig.add_axes([0.1,0.075,0.8,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=10/plot_ratio)
        cb.set_label('Probability',fontsize=13/plot_ratio)

        fig.savefig(plot_title, dpi=400)
        plt.close()
        print ('saving plot to file = ',plot_title)




#