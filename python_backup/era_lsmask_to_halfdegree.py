import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, interp

# --- open a sample forecast file to use for grib data template

filename = '../ecmwf/lsmask_0p25.grib'
ffcst = pygrib.open(filename)
grb = ffcst.read(1)[0] 
lsmask = grb.values
lats, lons = grb.latlons()
ny, nx = np.shape(lats)
print ('ny, nx = ', ny, nx)
print ('sample forecast file: ', filename)
print ('lats = ',lats[:,0])
print ('lons = ',lons[0,:])

# ---- loop through all ERA5 records in grib file

lats_halfdegree = lats[0::2,0::2]
lons_halfdegree = lons[0::2,0::2]
lfract_halfdegree = lsmask[0::2,0::2]
nlats, nlons = np.shape(lats_halfdegree)
ones = np.ones((nlats,nlons), dtype=np.float32)
zeros = np.zeros((nlats,nlons), dtype=np.float32)
lsmask_halfdegree = np.where(lfract_halfdegree > 0.5, ones, zeros)

#print (lsmask_halfdegree[0,:])

outfilename = '../ecmwf/lsmask_era5_0p5.cPick'
print (outfilename)
ouf = open(outfilename, 'wb')
cPickle.dump(lsmask_halfdegree, ouf)
cPickle.dump(lats_halfdegree, ouf)
cPickle.dump(lons_halfdegree, ouf)
ouf.close()    
            
           
# ---- Q-Q plot of forecast data

f = plt.figure(figsize=(6,4))

ax = f.add_axes([.02,.02,.96,.93])
ax.set_title('ERA-5 Land-sea mask',fontsize=13)

m = Basemap(llcrnrlon=lons_halfdegree[0,0],llcrnrlat=lats_halfdegree[-1,-1],\
    urcrnrlon=lons_halfdegree[-1,-1],urcrnrlat=lats_halfdegree[0,0],\
    resolution='l', projection='mill')

x, y = m(lons_halfdegree, lats_halfdegree)
for i in range(nlons):
    for j in range(nlats):
        if lsmask_halfdegree[j,i] == 1:
            xdot, ydot = m(lons_halfdegree[j,i],lats_halfdegree[j,i])
            m.plot(xdot,ydot,marker='.',markersize=1,color='LightGray')
            
m.drawstates(linewidth=0.2,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawcoastlines(linewidth=0.8,color='Gray')

figname = 'lsmask_era5.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()

