# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import datetime
import time
#import pygrib
import sys
#append the path to shapefile.py in quotes, or just have it in your working directory
#sys.path.append("")
import shapefile
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon

plt.ion()



sf = shapefile.Reader("ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
recs    = sf.records()
shapes  = sf.shapes()
records = []
shp_vertices   = []
for nshp in range(0,len(shapes)): 
    pts     = np.array(shapes[nshp].points,dtype=np.float32)
    prt     = shapes[nshp].parts
    par     = list(prt) + [len(pts)]
    records.append(sf.record)
    for pij in range(0,len(prt)):
        shp_vertices.append(pts[par[pij]:par[pij+1]])
shapes2 = []
for v in range(0,len(shp_vertices)):
    shapes2.append(shp_vertices[v])

shapes2East = []
for v in range(0,len(shapes2)):
    if v == 718 or v == 2802:
        pass
    else:
        currarray = np.array(shapes2[v])
        currarray[currarray[:,0]<0,0] = currarray[currarray[:,0]<0,0]+360
        #currarray[:,0] += 180
        fnd1 = np.nonzero(currarray[:,0]>180)[0]
        fnd2 = np.nonzero(currarray[:,0]<180)[0]
        currarray2 = np.array(currarray)
        if len(fnd1) > 0 and len(fnd2) > 0:
            currarray[fnd2,0] += 360
            currarray2[fnd1,0] -= 360        
            shapes2East.append(currarray2)
        shapes2East.append(currarray)

#################################################################
##### Uncomment the next section to obtain the NCEP/NCAR reanalysis data used in the example
#################################################################
#import ftplib
#ftp = ftplib.FTP("ftp.cdc.noaa.gov")
#ftp.login()
#ftp.cwd("Datasets/ncep.reanalysis.derived/surface")
#ftp.retrbinary("RETR slp.mon.mean.nc", open("slp_mon_mean.nc", 'wb').write)
#ftp.quit()
#################################################################

f1 = netCDF4.Dataset("slp.mon.mean.nc")
slp = np.array(f1.variables['slp'])
f1.close()
slp = (slp[:,0:-1,:]+slp[:,1:,:])/2


###Plot 1
###This one uses lines for the Countries

fig1 = plt.figure(figsize=(10,5.))
ax1 = fig1.add_subplot(1,1,1)
ax1.set_position([.08,.1,.85,.85])
linecoll1 = LineCollection(shapes2East, edgecolor=[.2,.2,.2],linewidths=.6,zorder=1,alpha=1.)  
ax1.add_collection(linecoll1)
cs1 = ax1.imshow(np.mean(slp[852:,:,:],axis=0),cmap=plt.cm.jet,vmin=1000,vmax=1025,extent=[0,360,-90,90])
#potentially use this to show statistical signficance 
fndGreaterThan = np.nonzero(np.mean(slp[852:,:,:],axis=0)>1020)
lon,lat = np.meshgrid(np.arange(1.25,360,2.5),np.arange(88.75,-90,-2.5))
ax1.scatter(lon[fndGreaterThan],lat[fndGreaterThan],marker='s',s=18,lw=.5,facecolor='None',edgecolor=[0,0,0])
ax1.set_xticks(np.arange(0,361,30))
ax1.set_yticks(np.arange(-60,91,30))
ax1.axis([0,360,-60,90])
ax1.set_title("Mean Sea Level Pressure Jan-Sep 2019 (hPa)",fontsize=16)
cax = fig1.add_axes([0.08,0.07,0.85,0.02])
cbar = fig1.colorbar(cs1, orientation='horizontal',cax=cax,extend='both')###,drawedges=True,ticks=np,format='%g')extend='both'
cbar.set_ticks(np.arange(-1000,1026,5))
cbar.ax.tick_params(labelsize=9)


###Plot 2
###This one uses polygons for the Countries

fig1 = plt.figure(figsize=(10,5.))
ax1 = fig1.add_subplot(1,1,1)
ax1.set_position([.08,.1,.85,.85])
patches = []
for v in range(0,len(shapes2East)):
    patches.append(Polygon(shapes2East[v]))
patchColl = PatchCollection(patches, edgecolor=[.2,.2,.2],facecolor=[.5,.5,.5],linewidths=.4,alpha=.5,zorder=2)
ax1.add_collection(patchColl)
lon,lat = np.meshgrid(np.arange(0,361,2.5),np.arange(90,-91,-2.5))
#Use pcolormesh if you want to have the image or mesh conform to a curved surface, like if you are using some specified projection
cs1 = ax1.pcolormesh(lon,lat,np.mean(slp[852:,:,:],axis=0),cmap=plt.cm.jet,vmin=1000,vmax=1025,alpha=1.,zorder=1)
ax1.set_xticks(np.arange(0,361,30))
ax1.set_yticks(np.arange(-60,91,30))
ax1.set_aspect(1.0)
ax1.axis([0,360,-60,90])
ax1.set_title("Mean Sea Level Pressure Jan-Sep 2019 (hPa)",fontsize=16)
cax = fig1.add_axes([0.08,0.07,0.85,0.02])
cbar = fig1.colorbar(cs1, orientation='horizontal',cax=cax,extend='both')###,drawedges=True,ticks=np,format='%g')extend='both'
cbar.set_ticks(np.arange(-1000,1026,5))
cbar.ax.tick_params(labelsize=9)



