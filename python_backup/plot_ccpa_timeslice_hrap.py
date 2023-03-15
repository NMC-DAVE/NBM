"""

plot_ccpa_timeslice_hrap.py cyyyymmddhh chour

"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.stats as stats
import pygrib
from netCDF4 import Dataset
from dateutils import hrs_since_day1CE_todate, \
    dateto_hrs_since_day1CE, hrstodate, datetohrs, dateshift
from mpl_toolkits.basemap import Basemap, interp

# ---- get the month and end time from the commmand line

cyyyymmddhh = sys.argv[1] # 01 etc
chour = sys.argv[2]
cyyyymmdd = cyyyymmddhh[0:8]

try:
    infile = '/Volumes/NBM/ccpa/ccpa.'+cyyyymmdd+\
        '/'+chour+'/ccpa.t'+chour+'z.06h.hrap.conus.gb2'
    print (infile)
    grb = pygrib.open(infile)
    print ('got grb')
    panal = grb.select()[0]
    precip_anal = panal.values
    print ('got panal')
    lats, lons = panal.latlons()
    lons = lons+360.
    print ('min, max lons = ', np.min(lons), np.max(lons))
    grb.close()           
except:
    print ('whoops!   some problem with ', infile)
        
# ===========================================================

m = Basemap(llcrnrlon=235.,llcrnrlat=31.,
    urcrnrlon = 245., urcrnrlat = 42.,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=245.,\
    resolution ='l',area_thresh=1000.) 
       
x, y = m(lons, lats)  
title = 'CCPA precipitation, 6 h ending '+chour+' UTC '+cyyyymmddhh
    
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
clevs = [0.0,0.1,0.3,0.6,1,2,3,5,7,10,15,20,25,30,40,50]

fig = plt.figure(figsize=(6.,7.7))
axloc = [0.06,0.11,0.92,0.8]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, precip_anal, clevs,\
    cmap=None, colors=colorst, extend='both')
CS3 = m.contour(x, y, precip_anal, levels=[25.0],\
    colors='Black', linewidths=1.5)
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
m.drawcounties(linewidth=0.3,color='LightGray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.06,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=5)
cb.set_label('Precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'CCPA_6hourly_precipitation_'+cyyyymmdd+'_'+chour+'UTC.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')   
