# plot_ccpa_landmask.py

import numpy as np
import sys
import pygrib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'

# ===========================================================

infile = 'ndfd_terrain_landwater.grib2.gb2'
print (infile)
flatlon = pygrib.open(infile)
fcst = flatlon.select()[0]
landmask = fcst.values
lats_ndfd, lons_ndfd = fcst.latlons()
if lats_ndfd[0,0] > lats_ndfd[-1,0]: 
    flipud = True
else:
    flipud = False
if flipud == True:
    lats_ndfd = np.flipud(lats_ndfd)
    lons_ndfd = np.flipud(lons_ndfd)
nlats_ndfd, nlons_ndfd = np.shape(lons_ndfd)
flatlon.close()
print (landmask[nlats_ndfd//2,0:nlons_ndfd:10])


title = 'CCPA on NDFD - land mask' 
    
m = Basemap(llcrnrlon=233.7234-360.,llcrnrlat=19.229,
    urcrnrlon = 300.95782-360., urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ndfd, lats_ndfd)    

fig = plt.figure(figsize=(8,6.5))
axloc = [0.08,0.1,0.87,0.88]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=16,color='Black')

for jy in range(0,nlats_ndfd,2):
    for ix in range(0,nlons_ndfd,2):

        if landmask[jy,ix] > 0:
            xdot, ydot = m(lons_ndfd[jy,ix],lats_ndfd[jy,ix])
            #print ('xdot, ydot = ',xdot, ydot)
            m.plot(xdot,ydot,marker='.',markersize=0.1,color='Black')

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

plot_title = 'ccpa_landmask.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')






