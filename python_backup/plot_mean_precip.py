"""
plot_fitted_precip_cdfs.py cmonth clead jy ix

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
rcParams['legend.fontsize']='medium'

# =====================================================================

def fraczero_possamps(nsamps, precip):
    """
    from the vector input sample precip_ens, define the fraction of
    samples with zero precipitation.   For the positive samples, add
    a small random number to deal with the fact that the data was 
    discretized to 0.1 mm, so that when later creating CDFs we don't 
    have values with lots of tied amounts.   Sort the nonzero amounts 
    and return.
    """
    number_zeros = 0
    precip_nonzero = np.delete(precip, \
        np.where(precip <= 0.0))  # censor at 0.1 mm
    nz = len(precip_nonzero)
    # data discretized, so add random component of this magnitude
    precip_nonzero = precip_nonzero + \
        np.random.uniform(low=-0.005,high=0.005,size=nz) 
    precip_nonzero = np.sort(precip_nonzero)  
    #print (precip_ens_nonzero[0:10]) 
    ntotal = len(precip)
    nzero = ntotal - len(precip_nonzero)
    fraction_zero = float(nzero) / float(ntotal)
    return fraction_zero, precip_nonzero, nz


# =====================================================================
    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clead = sys.argv[2] # 03, 06, 12, etc.


master_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time.  All members, dates for this 
#      month have been smushed into one leading index, dimension
#      nsamps, since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution fitting.

jmin = 93
jmax = 246
imin = 368
imax = 686
        
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
ncfile = master_directory + cmonth + '_apcp' '_h' + clead + '.nc'
print (ncfile)
nc = Dataset(ncfile)
lons_1d = nc.variables['lons_fcst'][imin:imax]
lats_1d = nc.variables['lats_fcst'][jmin:jmax]
precip = nc.variables['apcp_fcst'][:,jmin:jmax,imin:imax]
precip_mean = np.mean(precip, axis=0)
nsamps, nlats, nlons = np.shape(precip)
print ('nlats, lats[0], lats[-1] = ', nlats, lats_1d[0], lats_1d[-1])
print ('nlons, lons[0], lons[-1] = ', nlons, lons_1d[0], lons_1d[-1])
print ('nlats, nlons = ', nlats,nlons)
nc.close()
lons_2d, lats_2d = np.meshgrid(lons_1d,lats_1d)

fraction_zero = np.zeros((nlats,nlons), dtype=np.float32)
precip_mean = np.zeros((nlats,nlons), dtype=np.float32)
for ix in range(nlons):
    for jy in range(nlats):
        precip_1d = precip[:,jy,ix]
        fz, precip_nonzero, nz = fraczero_possamps(nsamps, precip_1d)
        fraction_zero[jy,ix] = fz
        precip_mean[jy,ix] = np.mean(precip_nonzero)

m = Basemap(llcrnrlon=lons_2d[0,0],llcrnrlat=lats_2d[-1,-1],\
    urcrnrlon=lons_2d[-1,-1],urcrnrlat=lats_2d[0,0],\
    resolution='l', projection='mill')
x, y = m(lons_2d, lats_2d)

# ---- make plots of fraction positive precip

clevs = [0.0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

fig = plt.figure(figsize=(9,6.))
axloc = [0.08,0.1,0.89,0.8]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Climatological probability of non-zero precipitation, '+\
    cmonth+' '+cleadb+' to '+clead+' h forecast'
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, 1.0-fraction_zero, clevs,\
    cmap=None, colors=colorst, extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
parallels = np.arange(-20.,80,10.)
m.drawparallels(parallels,labels=[1,0,0,0],color='LightGray')
meridians = np.arange(0.,360.,20.)
m.drawmeridians(meridians,labels=[0,0,0,1],color='LightGray')
    
# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.28)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('1 - fraction of samples with zero precipitation',fontsize=9)

# ---- set plot title

plot_title = 'fraction_zero_'+cmonth+'_h'+clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')




# --- plot the ensemble mean

clevs = [0.0,0.3,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10]
fig = plt.figure(figsize=(9,6.))
axloc = [0.08,0.1,0.89,0.8]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Mean of non-zero precipitation amounts, '+\
    cmonth+' '+cleadb+' to '+clead+' h forecast'
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, precip_mean, clevs,\
    cmap=None, colors=colorst, extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
parallels = np.arange(-20.,80,10.)
m.drawparallels(parallels,labels=[1,0,0,0],color='LightGray')
meridians = np.arange(0.,360.,20.)
m.drawmeridians(meridians,labels=[0,0,0,1],color='LightGray')
    
# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.28)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('mean of non-zero precipitation amounts (mm)',fontsize=9)

# ---- set plot title

plot_title = 'mean_precip_'+cmonth+'_h'+clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


