import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, interp



infilename = '/Volumes/Backup Plus/ecmwf/lsmask_era5_0p5.cPick'
print (infilename)
inf = open(infilename, 'rb')
lsmask_ecmwf = cPickle.load(inf)
lats_ecmwf = cPickle.load(inf)
lons_ecmwf = cPickle.load(inf)
inf.close()    
#lsmask_ecmwf = np.flipud(lsmask_ecmwf)
#lats_ecmwf = np.flipud(lats_ecmwf)

infilename = '/Volumes/Backup Plus/gefsv12/lsmask_gefsv12_0p5.cPick'
print (infilename)
inf = open(infilename, 'rb')
lsmask_gefsv12 = cPickle.load(inf)
lats_gefsv12 = cPickle.load(inf)
lons_gefsv12 = cPickle.load(inf)
nlats, nlons = np.shape(lats_gefsv12)
inf.close()    

lsmask = lsmask_ecmwf*lsmask_gefsv12

            

outfilename = 'lsmask_0p5.cPick'
print (outfilename)
ouf = open(outfilename, 'wb')
cPickle.dump(lsmask, ouf)
cPickle.dump(lats_gefsv12, ouf)
cPickle.dump(lons_gefsv12, ouf)
ouf.close()
            
        
# ---- plot the final land-sea mask

f = plt.figure(figsize=(6,4))

ax = f.add_axes([.02,.02,.96,.93])
ax.set_title('1/2-degree land grid points',fontsize=13)

m = Basemap(llcrnrlon=lons_gefsv12[0,0],llcrnrlat=lats_gefsv12[-1,-1],\
    urcrnrlon=lons_gefsv12[-1,-1],urcrnrlat=lats_gefsv12[0,0],\
    resolution='l', projection='mill')

x, y = m(lons_gefsv12, lats_gefsv12)
for i in range(nlons):
    for j in range(nlats):
        if lsmask[j,i] == 1:
            xdot, ydot = m(lons_gefsv12[j,i],lats_gefsv12[j,i])
            m.plot(xdot,ydot,marker='.',markersize=1.5,color='Gray')
            
m.drawstates(linewidth=0.4,color='Black')
m.drawcountries(linewidth=1,color='Black')
m.drawcoastlines(linewidth=1,color='Black')

figname = 'postproc/Fig1.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()

