"""
python gefsv12_ensmean_Dec2021.py cyyyymmddhh
"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import _pickle as cPickle
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 

cyear = cyyyymmddhh[0:4]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'

# ---- read the 66-h forecast

infile = '/Volumes/NBM/conus_gefsv12/qmapped/2021121100_mean_lead=066.cPick'
inf = open(infile,'rb')
ensmean_raw_ndfd_066 = cPickle.load(inf)
ensmean_qmapped_ndfd_066 = cPickle.load(inf) 
lons_ndfd = cPickle.load(inf)     
lats_ndfd = cPickle.load(inf)   
inf.close()      

# ---- read the 72-h forecast

infile = '/Volumes/NBM/conus_gefsv12/qmapped/2021121100_mean_lead=072.cPick'
inf = open(infile,'rb')
ensmean_raw_ndfd_072 = cPickle.load(inf)
ensmean_qmapped_ndfd_072 = cPickle.load(inf) 
inf.close()      

# ---- read the 78-h forecast

infile = '/Volumes/NBM/conus_gefsv12/qmapped/2021121100_mean_lead=078.cPick'
inf = open(infile,'rb')
ensmean_raw_ndfd_078 = cPickle.load(inf)
ensmean_qmapped_ndfd_078 = cPickle.load(inf) 
inf.close() 

# ---- read the 84-h forecast

infile = '/Volumes/NBM/conus_gefsv12/qmapped/2021121100_mean_lead=084.cPick'
inf = open(infile,'rb')
ensmean_raw_ndfd_084 = cPickle.load(inf)
ensmean_qmapped_ndfd_084 = cPickle.load(inf) 
inf.close()     

ensmean_raw = ensmean_raw_ndfd_066 + ensmean_raw_ndfd_072 + \
    ensmean_raw_ndfd_078 + ensmean_raw_ndfd_084
    
ensmean_qmapped = ensmean_qmapped_ndfd_066 + ensmean_qmapped_ndfd_072 + \
    ensmean_qmapped_ndfd_078 + ensmean_qmapped_ndfd_084

# ---- plot the raw ensemble mean by sector


m = Basemap(llcrnrlon=-125.,llcrnrlat=32.5,
    urcrnrlon = -113.5, urcrnrlat = 42.,\
    projection='lcc',lat_1=35.,lat_2=35.,lon_0=-120.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ndfd, lats_ndfd)

#clevs = [0.0,0.254, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 35.0, 50.0]
clevs = [0.0,0.01,0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0,10.0,15.0]
colorst = ['LightGray','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','LightGray']

drawcoasts = True

# ---- plot the raw ensemble mean

xdim = 5.
ydim = 6.3
fig = plt.figure(figsize=(xdim,ydim))
axloc = [0.02,0.1,0.96,0.79]
ax1 = fig.add_axes(axloc)
title = 'Raw GEFSv12 ensemble mean,\n'+\
    '60 to 84 h accumulations for the forecast\n'+\
    'starting at 00 UTC 11 December 2021'
ax1.set_title(title, fontsize=11,color='Black')
CS2 = m.contourf(x, y, ensmean_raw[:,:]/25.4, clevs,\
    cmap=None, colors=colorst, extend='max')

if drawcoasts == True: m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
m.drawcounties(linewidth=0.1,color='Gray')

cax = fig.add_axes([0.04,0.07,0.92,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Precipitation amount (in)',fontsize=9)

plot_title = 'raw_ensmean_cali_'+\
    cyyyymmddhh+'_60_to_84h.png'
fig.savefig(plot_title, dpi=300)
plt.close()
print ('saving plot to file = ',plot_title)


# ---- plot the quantile-mapped ensemble mean

fig = plt.figure(figsize=(xdim,ydim))
axloc = [0.02,0.1,0.96,0.79]
ax1 = fig.add_axes(axloc)
title = 'Statistically adjusted GEFSv12 ensemble mean,\n'+\
    '60 to 84 h accumulations for the forecast\n'+\
    'starting at 00 UTC 11 December 2021'
ax1.set_title(title, fontsize=11,color='Black')
CS2 = m.contourf(x, y, ensmean_qmapped[:,:]/25.4, clevs,\
    cmap=None, colors=colorst, extend='max')

if drawcoasts == True: m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
m.drawcounties(linewidth=0.1,color='Gray',zorder=20)

cax = fig.add_axes([0.04,0.07,0.92,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Precipitation amount (in)',fontsize=9)

plot_title = 'qmapped_ensmean_cali_'+\
    cyyyymmddhh+'_60_to_84h.png'
fig.savefig(plot_title, dpi=300)
plt.close()
print ('saving plot to file = ',plot_title)
#print ('Done!')



    
