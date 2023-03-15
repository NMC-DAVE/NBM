# plot_various_masks.py

import numpy as np
import sys
import pygrib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset

rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'

# ===========================================================

infile = '/Volumes/Backup Plus/ccpa/various_nbm_plus_mask.nc'
nc = Dataset(infile)
lats = nc.variables['latitude'][:,:]
lons = nc.variables['longitude'][:,:]
landmask = nc.variables['landmask'][:,:]
landmask_nolakes = nc.variables['landmask_nolakes'][:,:]
validmask = nc.variables['validmask'][:,:]
validmask_ccpa = nc.variables['validmask_ccpa'][:,:]
nlats, nlons = np.shape(lats)

if lats[0,0] > lats[-1,0]: 
    landmask = np.flipud(landmask)
    landmask_nolakes = np.flipud(landmask_nolakes)
    validmask = np.flipud(validmask)
    validmask_ccpa = np.flipud(validmask_ccpa)
    lats = np.flipud(lats)
    lons = np.flipud(lons)
nc.close()

for imask in range(4):
    if imask == 0:
        mask = landmask
        title = 'landmask'
        plot_title = 'landmask.png'
    elif imask == 1:
        mask = landmask_nolakes
        title = 'landmask_nolakes'
        plot_title = 'landmask_nolakes.png'
    elif imask == 2:
        mask = validmask
        title = 'validmask'
        plot_title = 'validmask.png'
    else:
        mask = validmask_ccpa
        title = 'validmask_ccpa'
        plot_title = 'validmask_ccpa.png'
    print ('processing ',title)
    
    m = Basemap(llcrnrlon=233.7234-360.,llcrnrlat=19.229,
        urcrnrlon = 300.95782-360., urcrnrlat = 54.37279,\
        projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
        resolution ='l',area_thresh=1000.)
    x, y = m(lons, lats)    

    fig = plt.figure(figsize=(8,6.5))
    axloc = [0.08,0.1,0.87,0.88]
    ax1 = fig.add_axes(axloc)
    ax1.set_title(title, fontsize=16,color='Black')

    for jy in range(0,nlats,5):
        for ix in range(0,nlons,5):

            if mask[jy,ix] > 0:
                xdot, ydot = m(lons[jy,ix],lats[jy,ix])
                m.plot(xdot,ydot,marker='s',markersize=0.1,color='Black')

    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    
    # draw parallels and meridians.
    # label on left and bottom of map.
    parallels = np.arange(20.,60.01,5.)
    m.drawparallels(parallels,labels=[1,0,0,0],color='Black',fontsize=8,linewidth=0.3)
    meridians = np.arange(220.,360.,5.)
    m.drawmeridians(meridians,labels=[0,0,0,1],color='Black',fontsize=8,linewidth=0.3)
   
    # ---- set plot title

    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    
print ('Done!')






