"""
save_lsmask.py -- GEFS land-sea mask
"""
from dateutils import daterange # utility of Jeff Whitaker's
import os, sys
import pygrib
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, interp

infile = 'hgt_sfc.grib2'

# read in grib file

fexist = os.path.exists(infile)
if fexist == True:
    grbfile = pygrib.open(infile)
    grb = grbfile.select()[0] # forecastTime
    height = grb.values 
    lats, lons = grb.latlons()
    nlats, nlons = np.shape(lats)
    height_subset = height[160:281:2,940:1202:2]
    lats_subset = lats[160:281:2,940:1202:2]
    lons_subset = lons[160:281:2,940:1202:2] - 360.
    nlats, nlons = np.shape(lats_subset)
    for i in range(131):
        print (height_subset[0,i], height_subset[-1,i])
             
land_sea_mask = np.ones((nlats,nlons), dtype=np.int32)
zeros = np.zeros((nlats,nlons), dtype=np.int32)
land_sea_mask = np.where(np.logical_and(height_subset > 0.0, height_subset < 0.03), zeros, land_sea_mask)
        

outfilename = '../gefsv12/lsmask_gefsv12_0p5.cPick'
print (outfilename)
ouf = open(outfilename, 'wb')
cPickle.dump(land_sea_mask, ouf)
cPickle.dump(lats_subset, ouf)
cPickle.dump(lons_subset, ouf)
ouf.close()

# ---- plot mask

f = plt.figure(figsize=(6,4))

ax = f.add_axes([.02,.02,.96,.93])
ax.set_title('GEFSv12 Land-sea mask',fontsize=13)

m = Basemap(llcrnrlon=lons_subset[0,0],llcrnrlat=lats_subset[-1,-1],\
    urcrnrlon=lons_subset[-1,-1],urcrnrlat=lats_subset[0,0],\
    resolution='l', projection='mill')

x, y = m(lons, lats)
for i in range(nlons):
    for j in range(nlats):
        if land_sea_mask[j,i] == 1:
            xdot, ydot = m(lons_subset[j,i],lats_subset[j,i])
            m.plot(xdot,ydot,marker='.',markersize=1,color='LightGray')
            
m.drawstates(linewidth=0.2,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawcoastlines(linewidth=0.8,color='Gray')

figname = 'lsmask.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()



# ==== write cPickle file

#outfile = 'land-sea-mask.cPick'
#print (outfile)
#ouf = open(outfile, 'wb')
#cPickle.dump(height_subset, ouf)
#ouf.close()


