"""
display_kalman_gain.py 

"""
import pygrib
import os, sys
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
clon = sys.argv[2]
clat = sys.argv[3]
rlon = float(clon)
rlat = float(clat)
ilead = int(clead)

# --- get lat/lon from a sample file

datadir = '/Users/Tom/python/ecmwf/'
infile = datadir + 't2m_era5_halfdegree_2019010100.cPick'
inf = open(infile, 'rb')
analysis = cPickle.load(inf)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
nlats, nlons = np.shape(lats)
inf.close()
lats_1d = lats[:,0]
lons_1d = lons[0,:]

ilon = find_nearest(lons_1d, rlon)
print ('ilon = ', ilon)
ilat = find_nearest(lats_1d, rlat)
print ('ilat = ', ilat)

# --- read in the Kalman gain generated by invert_localized_covs_2019.py

infile = 'Kalman_gain_4D_lead='+clead+'.cPick'
print (infile)
inf = open(infile, 'rb')
Kalman_gain_4D = cPickle.load(inf)
inf.close()        
gain_2d = Kalman_gain_4D[ilat,ilon,:,:]
            
# --- now plot the bias correlation map for the selected point

#clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1.0] 
clevs = [0.0,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.085,0.09,0.095,0.10] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.07,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF Kalman gain, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,gain_2d,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 'Kalman_gain_f'+clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')

