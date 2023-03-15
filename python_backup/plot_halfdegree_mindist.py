# plot_halfdegree_mindist.py

import numpy as np
import sys, os
import pygrib
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'

# ===========================================================


infile = 'nearest_NDFD_to_halfdegree.cPick'
fexist = os.path.exists(infile)
if fexist == True:
    inf = open(infile,'rb')
    jnear = cPickle.load(inf)
    inear = cPickle.load(inf)
    dist = cPickle.load(inf)
    ny, nx = np.shape(jnear)
    inf.close()

print ('max, min dist = ', np.max(dist), np.min(dist))
print ('jnear[0:ny,nx//2] = ',jnear[0:ny,nx//2])
print ('dist[0:ny,nx//2] = ',dist[0:ny,nx//2])
print ('inear[ny//2,0:nx] = ',inear[ny//2,0:nx])
print ('dist[ny//2,0:nx] = ',dist[ny//2,0:nx])
lats_halfdegree = np.arange(57.5, 18.99, -0.5)
lons_halfdegree = np.arange(-128.5, -57.5, 0.5)
lons_2d, lats_2d = np.meshgrid(lons_halfdegree,lats_halfdegree)

title = 'distance' 
    
m = Basemap(llcrnrlon=-128.5,llcrnrlat=19.,urcrnrlon =-57.5, urcrnrlat = 57.5,\
    projection='mill',resolution ='l')
x, y = m(lons_2d, lats_2d)    

fig = plt.figure(figsize=(8,6.5))
axloc = [0.08,0.1,0.87,0.88]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=16,color='Black')
clevs = [0,1,2,3,4,5]
cs = m.contour(x,y,dist,levels=clevs)
plt.clabel(cs,levels=clevs)
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# draw parallels and meridians.
# label on left and bottom of map.

parallels = np.arange(20.,55.01,5.)
m.drawparallels(parallels,labels=[1,0,0,0],color='Black',fontsize=8,linewidth=0.3)
meridians = np.arange(-125,-55,5.)
m.drawmeridians(meridians,labels=[0,0,0,1],color='Black',fontsize=8,linewidth=0.3)
   
# ---- set plot title

plot_title = 'halfdegree_distance.png'
fig.savefig(plot_title, dpi=400)
print ('saving plot to file = ',plot_title)
print ('Done!')






