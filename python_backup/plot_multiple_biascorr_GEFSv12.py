"""
plot_multiple_biascorr_GEFSv12.py
"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
from read_reanalysis_timeseries import read_reanalysis_timeseries
from read_forecast_timeseries_GEFSv12 import read_forecast_timeseries_GEFSv12
from decayavg_biascorrection2 import decayavg_biascorrection2
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable


# --------------------------------------------------------------

def initialize_date_lists(cyear, clead):

    start_date = cyear+'010100'
    end_date = cyear+'123100'
    date_list_anal = daterange(start_date, end_date, 24)
    ndates = len(date_list_anal)
    date_list_forecast = []
    for i in range(ndates):
        date_list_forecast.append(dateshift(date_list_anal[i], int(clead)))
    return date_list_anal, date_list_forecast

# --------------------------------------------------------------   

def read_bias_corrections(infile):

    print ('   reading beta estimates from ', infile)
    inf = open(infile,'rb')
    beta = cPickle.load(inf)
    lats = cPickle.load(inf)
    lons = cPickle.load(inf)
    inf.close()
    istat = 0
    return beta, lats, lons

# --------------------------------------------------------------  

# ---- various initialization

cyyyymmddhh = sys.argv[1]
clead = sys.argv[2]
cvariable = '2t'

cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'
cpath_forecast = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_gain = '/Volumes/Backup Plus/gefsv12/t2m/'
    
# ---- start the processing.

cyear = '2019'
date_list_anal, date_list_forecast = initialize_date_lists(cyear, clead)
ndates = len(date_list_forecast)
print (date_list_forecast)

#timeindex = np.where(date_list_forecast == cyyyymmddhh)[0]
timeindex = date_list_forecast.index(cyyyymmddhh)
print ('timeindex = ', timeindex)
    
# ---- read the reanalysis time series on the dates specified.  Note that
#      we pass in the dates of the forecast valid time, as we wish to 
#      upload the analyses to compare against the forecasts at this time.

analyses_3d, lats, lons = read_reanalysis_timeseries(cpath_era5, \
    date_list_forecast)

nlats, nlons = np.shape(lons)
npts = nlats*nlons

# ---- read the decaying average bias correction.

forecast_3d, latsf, lonsf = read_forecast_timeseries_GEFSv12 ( \
    cpath_forecast, date_list_anal, clead)     
        
infile = cpath_forecast + '2019_decayavg_lead'+clead+'.cPick'        
beta_decayavg, lats, lons = read_bias_corrections(infile)

# ---- read quantile mapping bias correction.

infile = cpath_forecast + '2019_qmap_lead'+clead+'.cPick'
beta_qmap, lats, lons = read_bias_corrections(infile)

# ---- read the univariate MOS correction

infile = cpath_forecast + '2019_MOS_lead'+clead+'.cPick'
beta_mos, lats, lons = read_bias_corrections(infile)

# ---- read the multi-variate MOS correction

infile = cpath_forecast + '2019_MOS_mvr_lead'+clead+'.cPick'
beta_mvmos, lats, lons = read_bias_corrections(infile)
 

forecast_2d = forecast_3d[timeindex,:,:]
analyzed_2d = analyses_3d[timeindex,:,:]
beta_decayavg_2d = beta_decayavg[timeindex,:,:]
beta_qmap_2d = beta_qmap[timeindex,:,:]
beta_mos_2d = beta_mos[timeindex,:,:]
beta_mvmos_2d = beta_mvmos[timeindex,:,:]

# ---- make 2-panel plot, forecast and analyzed


fig = plt.figure(figsize=(9,3.5))
plt.suptitle (r'T$_{2m}$ forecast and analyzed, IC = '+\
    cyyyymmddhh+', '+clead+'-hour forecast',fontsize=15)
    
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],\
    resolution='l', projection='mill')

# (a) Forecast

#colorst = ['#8FB3FF','#C4E8FF','#E4FFFF','White','#D8F9D8',\
#    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
#    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
colorst = ['DarkBlue','Blue', 'RoyalBlue', 'LightSkyBlue','PaleGreen',\
    'YellowGreen','Yellow','Orange','OrangeRed','Red','DarkRed']
    
lmin = (int(np.min(forecast_2d))//10)*10 - 5
lmax = (int(np.max(forecast_2d))//10)*10 + 5
clevs = range(lmin, lmax+1, 5)
clevs = range(-25, 26, 5)



axloc = [0.02,0.16,0.46,0.74]
ax1 = fig.add_axes(axloc)
title = '(a) Forecast'
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,forecast_2d,clevs,\
    cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='White')
m.drawcountries(linewidth=0.8,color='White')
m.drawstates(linewidth=0.8,color='White')

divider = make_axes_locatable(ax1)
cax = fig.add_axes([0.02,0.13,0.46,0.04])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Forecast temperature (deg C)',fontsize=9)

# (b) Decay avg

axloc = [0.52,0.16,0.46,0.74]
ax1 = fig.add_axes(axloc)
title = '(b) ERA-5 analyzed'
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,analyzed_2d,clevs,\
    cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='White')
m.drawcountries(linewidth=0.8,color='White')
m.drawstates(linewidth=0.8,color='White')

divider = make_axes_locatable(ax1)
cax = fig.add_axes([0.52,0.13,0.46,0.04])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Analyzed temperature (deg C)',fontsize=9)

# ---- set plot title

plot_title = 'forecast_analyzed_'+cyyyymmddhh+'_f'+\
    clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')





# ---- make 4-panel plot


fig = plt.figure(figsize=(9,6.5))
plt.suptitle (r'T$_{2m}$ bias corrections, IC = '+\
    cyyyymmddhh+', '+clead+'-hour forecast',fontsize=15)
    
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],\
    resolution='l', projection='mill')

# (a) Forecast

clevs = [-5, -4,-3,-2,-1,-0.5,0.5,1,2,3,4,5]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']
lmin = (int(np.min(forecast_2d))//10)*10 - 5
lmax = (int(np.max(forecast_2d))//10)*10 + 5

axloc = [0.02,0.56,0.46,0.37]
ax1 = fig.add_axes(axloc)
title = '(a) Decaying average'
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,beta_decayavg_2d,clevs,\
    cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

divider = make_axes_locatable(ax1)
cax = fig.add_axes([0.02,0.55,0.46,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Bias correction (deg C)',fontsize=9)

# (b) Decay avg

clevs = [-5, -4,-3,-2,-1,-0.5,0.5,1,2,3,4,5]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']
axloc = [0.52,0.56,0.46,0.37]
ax1 = fig.add_axes(axloc)
title = '(b) Quantile mapping'
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,beta_qmap_2d,clevs,\
    cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

divider = make_axes_locatable(ax1)
cax = fig.add_axes([0.52,0.55,0.46,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Bias correction (deg C)',fontsize=9)

# (c) Univariate MOS

clevs = [-5, -4,-3,-2,-1,-0.5,0.5,1,2,3,4,5]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']
axloc = [0.02,0.09,0.46,0.37]
ax1 = fig.add_axes(axloc)
title = '(c) Univariate MOS'
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,beta_mos_2d,clevs,\
    cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

divider = make_axes_locatable(ax1)
cax = fig.add_axes([0.02,0.08,0.46,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Bias correction (deg C)',fontsize=9)


# (d) Multi-variate MOS

clevs = [-5, -4,-3,-2,-1,-0.5,0.5,1,2,3,4,5]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']
axloc = [0.52,0.09,0.46,0.37]
ax1 = fig.add_axes(axloc)
title = '(d) Multi-variate MOS'
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,beta_mvmos_2d,clevs,\
    cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

divider = make_axes_locatable(ax1)
cax = fig.add_axes([0.52,0.08,0.46,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Bias correction (deg C)',fontsize=9)

# ---- set plot title

plot_title = 'bias_correction_4panel_'+cyyyymmddhh+'_f'+\
    clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



    
    
    