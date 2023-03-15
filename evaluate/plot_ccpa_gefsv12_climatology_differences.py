"""
plot_ccpa_gefsv12_climatology_differences.py cmonth cend_hour clead

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
rcParams['legend.fontsize']='medium'

# =====================================================================
    
    
# ---- inputs from command line

cmonth = sys.argv[1] # '01', '02' etc.
cend_hour = sys.argv[2] # 06, 12, 18, 00 -- end hour of 6-h period
clead = sys.argv[3] # 072
imonth = int(cmonth) - 1
nstride = 1

# ---- set parameters

pflag = False # for print statements
master_directory = '/Volumes/Backup Plus/ccpa/'
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']


   
   
   
# --- now read CCPA gridded.

infile = 'quantiles_ccpa_'+cmonths[imonth]+'_'+cend_hour+'UTC.cPick'
print (infile)
inf = open(infile, 'rb')
precip_mean_ccpa = cPickle.load(inf)
precip_q95_ccpa = cPickle.load(inf)
precip_q98_ccpa = cPickle.load(inf)
lons_ccpa = cPickle.load(inf)
lats_ccpa = cPickle.load(inf)
inf.close()
   

infile = 'quantiles_gefsv12_'+cmonths[imonth]+'_lead'+clead+'.cPick'
print (infile)
inf = open(infile, 'rb')
precip_mean_gefsv12 = cPickle.load(inf)
precip_q95_gefsv12 = cPickle.load(inf)
precip_q98_gefsv12 = cPickle.load(inf)
lons_1d = cPickle.load(inf)
lats_1d = cPickle.load(inf)
lons_gefsv12= cPickle.load(inf)
lats_gefsv12 = cPickle.load(inf)
inf.close()





print ('interpolating GEFSv12 to NDFD grid. ')
precip_mean_gefsv12 = np.flipud(precip_mean_gefsv12)
precip_q95_gefsv12 = np.flipud(precip_q95_gefsv12)
precip_q98_gefsv12 = np.flipud(precip_q98_gefsv12)
lats_1d_flipud = np.flipud(lats_1d)
precip_mean_gefsv12_on_ndfd = interp(precip_mean_gefsv12, \
    lons_1d, lats_1d_flipud, lons_ccpa, lats_ccpa, \
    checkbounds=False, masked=False, order=1)
    
precip_q95_gefsv12_on_ndfd = interp(precip_q95_gefsv12, \
    lons_1d, lats_1d_flipud, lons_ccpa, lats_ccpa, \
    checkbounds=False, masked=False, order=1)
    
precip_q98_gefsv12_on_ndfd = interp(precip_q98_gefsv12, \
    lons_1d, lats_1d_flipud, lons_ccpa, lats_ccpa, \
    checkbounds=False, masked=False, order=1)

precip_mean_diff = precip_mean_gefsv12_on_ndfd - precip_mean_ccpa
precip_q95_diff = precip_q95_gefsv12_on_ndfd - precip_q95_ccpa
precip_q98_diff = precip_q98_gefsv12_on_ndfd - precip_q98_ccpa


# ---- plot mean differences

#colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
#    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
colorst = ['#b6773a','Peru','SandyBrown', 'BurlyWood','Bisque','White','#edfced', \
    '#d3f8d3','LightGreen','#48f648','#00b824']
clevels = [-1,-0.7,-0.5,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.5,0.7,1.0]
m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
    urcrnrlon = 300.95782, urcrnrlat = 54.37,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ccpa, lats_ccpa)

fig = plt.figure(figsize=(8.,6.5))
axloc = [0.02,0.1,0.96,0.81]
ax1 = fig.add_axes(axloc)
title = cmonths[imonth]+'GEFSv12 - CCPA mean differences, '+clead+'-h forecast'
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, precip_mean_diff, clevels,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevels,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Mean 6-hourly precipitation difference (mm)',fontsize=9)

# ---- set plot title

plot_title = 'gefsv12_minus_ccpa_mean_difference_precip_'+cmonths[imonth]+\
    '_'+clead+'UTC.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# ---- plot q95 differences

#colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
#    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
colorst = ['#b6773a','Peru','SandyBrown', 'BurlyWood','Bisque','White','#edfced', \
    '#d3f8d3','LightGreen','#48f648','#00b824']
clevels = [-9,-7,-4,-2,-1,-0.5,0.5,1,2,4,7,9]
m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
    urcrnrlon = 300.95782, urcrnrlat = 54.37,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)


fig = plt.figure(figsize=(8.,6.5))
axloc = [0.02,0.1,0.96,0.81]
ax1 = fig.add_axes(axloc)
title = cmonths[imonth]+' GEFSv12 - CCPA 95th percentile differences, '+clead+'-h forecast'
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, precip_q95_diff, clevels,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevels,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('95th percentile 6-hourly precipitation difference (mm)',fontsize=9)

# ---- set plot title

plot_title = 'gefsv12_minus_ccpa_q95_difference_precip_'+cmonths[imonth]+\
    '_'+clead+'h.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



# ---- plot q98 differences

#colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
#    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
    
colorst = ['#b6773a','Peru','SandyBrown', 'BurlyWood','Bisque','White','#edfced', \
    '#d3f8d3','LightGreen','#48f648','#00b824']   
    
m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
    urcrnrlon = 300.95782, urcrnrlat = 54.37,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)


fig = plt.figure(figsize=(8.,6.5))
axloc = [0.02,0.1,0.96,0.81]
ax1 = fig.add_axes(axloc)
title = cmonths[imonth]+' GEFSv12 - CCPA 98th percentile differences, '+clead+'-h forecast'
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, precip_q98_diff, clevels,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevels,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('98th percentile 6-hourly precipitation difference (mm)',fontsize=9)

# ---- set plot title

plot_title = 'gefsv12_minus_ccpa_q98_difference_precip_'+cmonths[imonth]+\
    '_'+clead+'h.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')






# ---- plot 3-panel q98 in western US

#colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
#    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
    

    
m = Basemap(llcrnrlon=237,llcrnrlat=30,\
    urcrnrlon = 260., urcrnrlat = 51,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=250.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ccpa, lats_ccpa)
    
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
clevs = [0,0.1,0.3,0.5,1.0,1.3,1.6,2.0,3.0,4.0,6.0,8.0,10.0,15.0,20.0, 25.0]
fig = plt.figure(figsize=(9.,4.5))

axloc = [0.02,0.1,0.3,0.81]
ax1 = fig.add_axes(axloc)
title = '(a) '+cmonths[imonth]+' GEFSv12 98th %ile,\n'+clead+' h forecast'
ax1.set_title(title, fontsize=11,color='Black')
CS2 = m.contourf(x, y, precip_q98_gefsv12_on_ndfd, clevs,\
    cmap=None, colors=colorst, extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

cax = fig.add_axes([0.02,0.09,0.3,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Precipitation (mm)',fontsize=9)


axloc = [0.35,0.1,0.3,0.81]
ax1 = fig.add_axes(axloc)
title = '(b) '+cmonths[imonth]+' CCPA 98th %ile,\n6-h ending '+cend_hour+' UTC'
ax1.set_title(title, fontsize=11,color='Black')
CS2 = m.contourf(x, y, precip_q98_ccpa, clevs,\
    cmap=None, colors=colorst, extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

cax = fig.add_axes([0.35,0.09,0.3,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Precipitation (mm)',fontsize=9)


axloc = [0.68,0.1,0.3,0.81]
ax1 = fig.add_axes(axloc)
colorst = ['#b6773a','Peru','SandyBrown', \
    'BurlyWood','Bisque','White','#edfced', \
    '#d3f8d3','LightGreen','#48f648','#00b824']   
clevels = [-9,-7,-4,-2,-1,-0.5,0.5,1,2,4,7,9]
title = '(c) '+cmonths[imonth]+' GEFSv12 - CCPA diff.'
ax1.set_title(title, fontsize=11,color='Black')
CS2 = m.contourf(x, y, precip_q98_diff, clevels,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.68,0.09,0.3,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevels,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Precipitation difference (mm)',fontsize=9)

# ---- set plot title

plot_title = 'gefsv12_minus_ccpa_q98_difference_precip_westUS'+cmonths[imonth]+\
    '_'+clead+'h.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



