"""
plot_GEFS_beta_cov.py

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
cseason = sys.argv[2] # warm or cold
clon = sys.argv[3]
clat = sys.argv[4]
rlon = float(clon)
rlat = float(clat)
ilead = int(clead)
efold = 800.0

cpath_beta = '/Volumes/Backup Plus/gefsv12/t2m/beta/'

# --- read in the decaying average bias correction estimates.

infile = cpath_beta+'Localized_Bbeta_'+cseason+\
    '_lead='+clead+'_'+str(efold)+'.cPick'
inf = open(infile,'rb')
Bbeta_localized_4D = cPickle.load(inf)
inf.close()

infile = '/Volumes/Backup Plus/gefsv12/t2m/gefsv12_latlon_subset.cPick'
inf = open(infile,'rb')
lats = cPickle.load(inf)
lons = cPickle.load(inf)
nlats, nlons = np.shape(lats)
npts = nlats*nlons
inf.close()
lats = np.flipud(lats)
print ('min, max lons = ',np.min(lons), np.max(lons))

lats_1d = lats[:,0]
lons_1d = lons[0,:]
ilon = find_nearest(lons_1d, rlon)
ilat = find_nearest(lats_1d, rlat)

# --- now plot the localized bias error correlation for the selected point
   
clevs = [-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000'] 
    
bias_error_cov_map = Bbeta_localized_4D[ilat,ilon,:,:]  
#bias_error_cov_map = np.flipud(bias_error_cov_map)

# --- now plot the localized bias error covariance for the selected point

#clevs = [-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
#colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
#    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']
clevs = [-0.3,0.01,0.03,0.05,0.07,0.1,0.12,0.14,0.17,0.2,0.25,0.3,0.5,0.75,1.0,1.5,2.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

print ('min, max covariance = ', np.min(bias_error_cov_map), np.max(bias_error_cov_map))
fig = plt.figure(figsize=(6.,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 T$_{2m}$ bias covariance map, '+cseason+\
    ' '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,bias_error_cov_map, clevs, \
    cmap=None, colors=colorst, extend='both') #
xdot, ydot = m(rlon,rlat)
m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.1)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,format='%g') # ticks=clevs,
cb.ax.tick_params(labelsize=7)
cb.set_label(r'Temperature bias covariance (deg C$^2$)')

# ---- set plot title

plot_title = 't2m_bias_error_covariance_'+cseason+'_GEFSv12_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300,fontsize=9)

print ('saving plot to file = ',plot_title)
print ('Done!')




