"""
display_bias_corrs_2019.py

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.signal as signal
import scipy.stats as stats
from astropy.convolution import convolve


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
clon = sys.argv[2]
clat = sys.argv[3]
rlon = float(clon)
rlat = float(clat)
ilead = int(clead)
datadir = '/Users/Tom/python/ecmwf/'
cvariable = '2t'
datestart = dateshift('2019010100',ilead)
date_list_anal = daterange(datestart,'2019123100',24)
ndates = len(date_list_anal)
date_list_fcst = []
for idate in range(ndates):
    date_list_fcst.append(dateshift(date_list_anal[idate],-ilead)) # initial times of fcst
    
# ---- call initialization routine

for idate, datea in enumerate(date_list_anal):
    
    datef = date_list_fcst[idate]
    if datea == '2019010100': dstart = idate
    print ('------ processing analysis, forecast dates = ', datea, datef)

    # ---- read the ECMWF ERA5 reanalysis at this analysis date.
    
    infile = datadir + 't2m_era5_halfdegree_'+datea+'.cPick'
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lats_1d = lats[:,0]
        lons = cPickle.load(inf)
        lons_1d = lons[0,:]
        nlats, nlons = np.shape(lats) 
        complete_fcst = np.zeros((ndates,nlats,nlons), dtype=np.float32)
        complete_anal = np.zeros((ndates,nlats,nlons), dtype=np.float32)
        ilon = find_nearest(lons_1d, rlon)
        print ('ilon = ', ilon)
        ilat = find_nearest(lats_1d, rlat)
        print ('ilat = ', ilat)
    inf.close()

    # ---- read the ECMWF control forecast at this lead time and initial date
 
    infile = datadir + cvariable+'_'+datef+'_f'+clead+'.grib2'  
    grbfile = pygrib.open(infile) 
    grb = grbfile.select()[0] 
    forecast = grb.values
    grbfile.close()
    complete_fcst[idate,:,:] = forecast[:,:]
    
    # ---- read the ERA5 analysis valid at this date.
    
    infilename = datadir+'t2m_era5_halfdegree_'+datea+'.cPick'
    inf = open(infilename, 'rb')
    obs = cPickle.load(inf)
    inf.close()
    complete_anal[idate,:,:] = obs[:,:]
 
# ---- estimate the slowly time varying bias.
  
window_std = 10
window = signal.gaussian(201,window_std)  
complete_difference = complete_fcst - complete_anal
bias_estimate = np.copy(complete_difference)    
random_error_estimate = np.copy(complete_difference) 
for ix in range(nlons):
    for jy in range(nlats):
        timeseries = complete_difference[:,jy,ix]
        bias_estimate[:,jy,ix]= convolve(timeseries, window)
        random_error_estimate[:,jy,ix] = timeseries[:] - bias_estimate[:,jy,ix]
    
# --- generate statistics 

fmean = np.mean(complete_fcst,axis=0)
amean = np.mean(complete_anal,axis=0)
bias = fmean-amean
                
# ---- read 2019 from cPickle file

infile = 'covstats_bias_random_ecmwf2019_lead='+clead+'.cPick'
inf = open(infile,'rb')
cov_bias_2019 = cPickle.load(inf)
cov_random_2019 = cPickle.load(inf)
var_bias_2019 = cPickle.load(inf)
var_random_2019 = cPickle.load(inf)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
inf.close()     

# ---- read localized 2019 from cPickle file

infile = 'covstats_bias_random_ecmwf2019_localized_lead='+clead+'.cPick'
inf = open(infile,'rb')
cov_bias_localized_2019 = cPickle.load(inf)
cov_random_localized_2019 = cPickle.load(inf)
inf.close()   


# ---- read from cPickle file

infile = 'covstats_bias_random_ecmwf2018_lead='+clead+'.cPick'
inf = open(infile,'rb')
cov_bias_2018 = cPickle.load(inf)
cov_random_2018 = cPickle.load(inf)
var_bias_2018 = cPickle.load(inf)
var_random_2018 = cPickle.load(inf)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
inf.close()     

# ---- read localized from cPickle file

infile = 'covstats_bias_random_ecmwf2018_localized_lead='+clead+'.cPick'
inf = open(infile,'rb')
cov_bias_localized_2018 = cPickle.load(inf)
cov_random_localized_2018 = cPickle.load(inf)
inf.close()   

            

rmse = np.sqrt(np.sum((complete_fcst - complete_anal)**2,  axis=0)/ndates)
analvar = np.ones((nlats, nlons), dtype=np.float32)
random_error = np.sqrt(rmse**2 - analvar - bias**2)
    
correlation_bias_2019 = np.zeros((nlats, nlons), dtype=np.float32)
correlation_random_error_2019 = np.zeros((nlats, nlons), dtype=np.float32)
correlation_bias_local_2019 = np.zeros((nlats, nlons), dtype=np.float32)
correlation_random_error_local_2019 = np.zeros((nlats, nlons), dtype=np.float32)
correlation_bias_2018 = np.zeros((nlats, nlons), dtype=np.float32)
correlation_random_error_2018 = np.zeros((nlats, nlons), dtype=np.float32)
correlation_bias_local_2018 = np.zeros((nlats, nlons), dtype=np.float32)
correlation_random_error_local_2018 = np.zeros((nlats, nlons), dtype=np.float32)
for ilat2 in range(nlats):
    for ilon2 in range(nlons):
        correlation_bias_2018[ilat2,ilon2] = cov_bias_2018[ilat,ilon,ilat2,ilon2] / \
            (np.sqrt(var_bias_2018[ilat,ilon])*np.sqrt(var_bias_2018[ilat2,ilon2]) )
        correlation_random_error_2018[ilat2,ilon2] = cov_random_2018[ilat,ilon,ilat2,ilon2] / \
            (np.sqrt(var_random_2018[ilat,ilon])*np.sqrt(var_random_2018[ilat2,ilon2]) )    
        correlation_bias_local_2018[ilat2,ilon2] = cov_bias_localized_2018[ilat,ilon,ilat2,ilon2] / \
            (np.sqrt(var_bias_2018[ilat,ilon])*np.sqrt(var_bias_2018[ilat2,ilon2]) )
        correlation_random_error_local_2018[ilat2,ilon2] = cov_random_localized_2018[ilat,ilon,ilat2,ilon2] / \
            (np.sqrt(var_random_2018[ilat,ilon])*np.sqrt(var_random_2018[ilat2,ilon2]) ) 
        correlation_bias_2019[ilat2,ilon2] = cov_bias_2019[ilat,ilon,ilat2,ilon2] / \
            (np.sqrt(var_bias_2019[ilat,ilon])*np.sqrt(var_bias_2019[ilat2,ilon2]) )
        correlation_random_error_2019[ilat2,ilon2] = cov_random_2019[ilat,ilon,ilat2,ilon2] / \
            (np.sqrt(var_random_2019[ilat,ilon])*np.sqrt(var_random_2019[ilat2,ilon2]) )    
        correlation_bias_local_2019[ilat2,ilon2] = cov_bias_localized_2019[ilat,ilon,ilat2,ilon2] / \
            (np.sqrt(var_bias_2019[ilat,ilon])*np.sqrt(var_bias_2019[ilat2,ilon2]) )
        correlation_random_error_local_2019[ilat2,ilon2] = cov_random_localized_2019[ilat,ilon,ilat2,ilon2] / \
            (np.sqrt(var_random_2019[ilat,ilon])*np.sqrt(var_random_2019[ilat2,ilon2]) ) 
            
# --- now plot the bias correlation map for the selected point

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.07,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ bias correlation map, Jan-Dec 2019, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,correlation_bias_2019,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 't2m_bias_correlation_JanDec2019_f'+clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now plot the random error correlation map for the selected point

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.07,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ random error correlation map, Jan-Dec 2019, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,correlation_random_error_2019,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 't2m_random_correlation_JanDec2019_f'+clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')




# ---- localized!



# --- now plot the bias correlation map for the selected point

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.07,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'Localized ECMWF T$_{2m}$ bias correlation map, Jan-Dec 2019, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,correlation_bias_local_2019,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 't2m_bias_correlation_local_JanDec2019_f'+clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now plot the random error correlation map for the selected point

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.07,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'Localized ECMWF T$_{2m}$ random error correlation map, Jan-Dec 2019, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,correlation_random_error_local_2019,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 't2m_random_correlation_local_JanDec2019_f'+clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')





# --- now plot the bias correlation map for the selected point

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.07,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ bias correlation map, Jan-Dec 2018, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,correlation_bias_2018,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 't2m_bias_correlation_JanDec2018_f'+clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now plot the random error correlation map for the selected point

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.07,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ random error correlation map, Jan-Dec 2018, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,correlation_random_error_2018,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 't2m_random_correlation_JanDec2018_f'+clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')




# ---- localized!



# --- now plot the bias correlation map for the selected point

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.07,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'Localized ECMWF T$_{2m}$ bias correlation map, Jan-Dec 2018, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,correlation_bias_local_2018,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 't2m_bias_correlation_local_JanDec2018_f'+clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now plot the random error correlation map for the selected point

clevs = [0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,1.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.07,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'Localized ECMWF T$_{2m}$ random error correlation map, Jan-Dec 2018, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,correlation_random_error_local_2018,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 't2m_random_correlation_local_JanDec2018_f'+clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')







