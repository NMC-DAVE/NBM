# plot_landwater_ndfd_conus.py

import numpy as np
import _pickle as cPickle
import sys
from netCDF4 import Dataset
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.stats as stats
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
import pygrib

rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'


for iyear in range(2002,2020):
    cyear = str(iyear)
    for imonth in range(1,13):
        if imonth < 10:
            cmonth = '0'+str(imonth)
        else:
            cmonth = str(imonth)
            
        infile = '/Volumes/Backup Plus/ccpa/'+cyear+cmonth+'_ccpa_sum_on_ndfd_grid_6hourly.nc'
        print (infile)        
        nc = Dataset(infile,'r')
        apcp_sum = nc.variables['apcp_sum'][:,:]
        apcp_count = nc.variables['apcp_count'][:,:]
        if iyear == 2002 and imonth == 1:
            lons = nc.variables['lons'][:,:]
            print (np.min(lons), np.max(lons))
            lats = nc.variables['lats'][:,:]
            ny, nx = np.shape(lons)
            zeros = np.zeros((ny,nx), dtype=np.float64)
            ones = np.ones((ny,nx), dtype=np.float64)
            ktrarr = np.zeros((ny,nx), dtype=np.float64)
            apcp_sum_month_year = np.zeros((ny,nx), dtype=np.float64)
            apcp_sum_count = np.zeros((ny,nx), dtype=np.float32)
        nc.close()
    
        apcp_sum_month_year = apcp_sum_month_year + apcp_sum
        apcp_sum_count = apcp_sum_count + apcp_count

print (apcp_sum_month_year[ny//2, 0:nx:10])
print (apcp_sum_count[ny//2, 0:nx:10])
apcp_mean = apcp_sum_month_year / apcp_sum_count
apcp_mask = ma.getmask(apcp_mean)
apcp_mean_nomask = np.where(apcp_mask == True, -99.99*ones, apcp_mean)
apcp_mask_ones = np.where(apcp_mask == True, zeros, ones)


# ===========================================================

m = Basemap(llcrnrlon=233.7234-360.,llcrnrlat=19.229,
    urcrnrlon = 300.95782-360., urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
colorst = ['White','Black','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
x, y = m(lons, lats)  

title = 'CCPA precipitation daily mean'     
clevs = [-100,-0.00001,0.00001,0.05,0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,5,7,10]
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, apcp_mean_nomask, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=5)
cb.set_label('Mean precipitation / 6 h (mm)',fontsize=9)

# ---- set plot title

plot_title = 'CCPA_daily_precipitation_mean.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')   




# ------ plot the mask

title = 'CCPA implicit mask '     
clevs = [-100,-0.00001,0.00001,0.05,0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,5,7,10]
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, apcp_mask_ones, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=5)
cb.set_label('mask',fontsize=9)

# ---- set plot title

plot_title = 'CCPA_mask.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')    
    
