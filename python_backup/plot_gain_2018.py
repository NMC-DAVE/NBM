"""
plot_gain_2018.py

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
from read_reanalysis_timeseries import read_reanalysis_timeseries
from reformat_2d_to_4d_f90 import reformat_2d_to_4d_f90

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
cwarm = sys.argv[2]
clon = sys.argv[3]
clat = sys.argv[4]
rlon = float(clon)
rlat = float(clat)
ilead = int(clead)


if clead == '24':
    efold_random = '400.0'
    efold_bias = '1400.0'
elif clead == '120':
    if cwarm == 'warm':
        efold_random = '400.0'
        efold_bias = '1400.0'
    else:
        efold_random = '800.0' # '800.0'
        efold_bias = '1400.0' # '1400.0'
else:
    print ('invalid lead time. Stopping.')
    sys.exit()
        

cpath_Bx = '/Volumes/Backup Plus/ecmwf/Bx/'
cpath_Bbeta = '/Volumes/Backup Plus/ecmwf/Bbeta/'
cpath_gain = '/Volumes/Backup Plus/python/KFgain/'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
date_list_forecast = ['2018010100']

# --- read in a single reanalysis in order to get lats and lons.  Find the 
#     index of the nearest grid point.

analyses_3d, lats, lons = read_reanalysis_timeseries(cpath_era5, \
    date_list_forecast)
nlats, nlons = np.shape(lons)
print (nlats,nlons)
print ('min, max lons ', np.min(lons), np.max(lons))
print ('min, max lats ', np.min(lats), np.max(lats))
lons = lons - 360.
lats_1d = lats[:,0]
lons_1d = lons[0,:]
ilon = find_nearest(lons_1d, rlon)
ilat = find_nearest(lats_1d, rlat)

# ---- read the localized Bx forecast error covariance

infile = cpath_Bx+'Localized_Bx_'+cwarm+\
    'season_year2018_lead='+clead+'_efold'+\
    str(efold_random)+'.cPick'
print ('reading localized forecast error covariances from ', infile)
inf = open(infile,'rb')
Bx_localized_2d = cPickle.load(inf)
inf.close()
Bx_localized = reformat_2d_to_4d_f90(Bx_localized_2d, nlats, nlons)

# ---- read the bias correction error covariance localized

infile = cpath_Bbeta+'Localized_Bbeta_'+cwarm+\
    'season_year2018_lead='+clead+'_'+str(efold_bias)+'.cPick'
inf = open(infile,'rb')
Bbeta_localized_2d = cPickle.load(inf)
inf.close()
Bbeta_localized = reformat_2d_to_4d_f90(Bbeta_localized_2d, nlats, nlons)

# ---- extract forecast error variance and bias correction error variance

Bx_variance = np.zeros((nlats,nlons), dtype=np.float32)
Bbeta_variance = np.zeros((nlats,nlons), dtype=np.float32)

for j in range(nlats):
    for i in range(nlons):
        Bx_variance[j,i] = Bx_localized[j,i,j,i]
        Bbeta_variance[j,i] = Bbeta_localized[j,i,j,i]
        

# ---- compute 1-point correlation maps

#bias_error_corr_map = Bbeta_localized[:,:,ilat,ilon] / \
#    ( np.sqrt(Bbeta_variance[ilat,ilon]) * np.sqrt(Bbeta_variance[:,:]) ) 
bias_error_corr_map = np.zeros((nlats, nlons), dtype=np.float32)
bias_error_cov_map = np.zeros((nlats, nlons), dtype=np.float32)
for i in range(nlons):
    for j in range(nlats):
        bias_error_corr_map[j,i] = Bbeta_localized[ilat,ilon,j,i] / \
            ( np.sqrt(Bbeta_variance[ilat,ilon]) * np.sqrt(Bbeta_variance[j,i]) )
        bias_error_cov_map[j,i] = Bbeta_localized[ilat,ilon,j,i]
Bmax = np.max(Bbeta_variance)
        
forecast_error_corr_map = Bx_localized[:,:,ilat,ilon] / \
    (np.sqrt(Bx_variance[ilat,ilon])*np.sqrt(Bx_variance[:,:]))

# --- read the Kalman filter gain for bias

gain_infile = cpath_gain + '2018_KFgain_flocal'+str(efold_random)+\
    '_blocal'+str(efold_bias)+'_2018_'+cwarm+'_lead'+\
    clead+'.cPick'
inf = open(gain_infile, 'rb')
Kalman_gain_beta_4d = cPickle.load(inf)
Kalman_gain_beta_map = Kalman_gain_beta_4d[:,:, ilat,ilon]
inf.close()
print ('min, max Kalman_gain_beta_map = ',\
    np.min(Kalman_gain_beta_map), np.max(Kalman_gain_beta_map) )
Kmax = np.max(np.abs(Kalman_gain_beta_map))  
                    
# --- now plot the localized bias error correlation for the selected point 

clevs = [-0.99,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000'] 
    
fig = plt.figure(figsize=(6.,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ bias correlation map, '+cwarm+\
    '-season 2018, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,bias_error_corr_map,clevs,\
    cmap=None,colors=colorst,extend='both')
xdot, ydot = m(rlon,rlat)
m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.1)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias correlation')

# ---- set plot title

plot_title = 't2m_bias_error_correlation_'+cwarm+'2018_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300,fontsize=9)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now plot the localized bias error covariance for the selected point 

#clevs = [-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
#colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
#    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000'] 
clevs = [-0.3,0.01,0.03,0.05,0.07,0.1,0.12,0.14,0.17,0.2,0.25,0.3,0.5,0.75,1.0,1.5,2.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

print ('min, max covariance = ', np.min(bias_error_cov_map), np.max(bias_error_cov_map))    
fig = plt.figure(figsize=(6.,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ bias covariance map, '+cwarm+\
    '-season 2018, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,bias_error_cov_map, clevs, \
    cmap=None, colors=colorst, extend='both') # 
xdot, ydot = m(rlon,rlat)
m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.1)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,format='%g') # ticks=clevs,
cb.ax.tick_params(labelsize=7) 
cb.set_label('Temperature bias covariance')

# ---- set plot title

plot_title = 't2m_bias_error_covariance_'+cwarm+'2018_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300,fontsize=9)

print ('saving plot to file = ',plot_title)
print ('Done!')




      
                    
# --- now plot the Kalman gain map for the selected point

clevs = [-0.99,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']    
    
fig = plt.figure(figsize=(6,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ Kalman gain map, '+cwarm+\
    '-season 2018, '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,Kalman_gain_beta_map/Kmax,clevs,\
    cmap=None,colors=colorst,extend='both')
xdot, ydot = m(rlon,rlat)
m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.1)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Kalman gain / max(abs(Kalman gain))',fontsize=9)

# ---- set plot title

plot_title = 't2m_Kalman_gain_'+cwarm+'2018_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now plot the forecast-error variance

bmax = np.max(np.sqrt(Bx_variance))
if bmax > 100:
    bmaxplot = 100.*100*int(bmax/100.)
elif bmax > 10 and bmax <=100.:
    bmaxplot = 10.*10*int(bmax/10.)
else:
    bmaxplot = 1 + int(bmax)
clevs = [0., bmaxplot/10., 2*bmaxplot/10., 3*bmaxplot/10.,\
    4*bmaxplot/10.,5*bmaxplot/10., 6*bmaxplot/10.,\
    7*bmaxplot/10.,8*bmaxplot/10., 9*bmaxplot/10.,\
    bmaxplot]
    
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
fig = plt.figure(figsize=(6,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ forecast-error standard deviation, '+cwarm+\
    '-season 2018,\n'+clead+'-hour forecast, ' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,np.sqrt(Bx_variance),clevs,\
    cmap=None,colors=colorst,extend='both')
xdot, ydot = m(rlon,rlat)
m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.1)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Forecast background error standard deviation',fontsize=9)

# ---- set plot title

plot_title = 't2m_forecast_spread_'+cwarm+'2018_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



# --- now plot the bias variance

bmax = np.max(np.sqrt(Bbeta_variance))
if bmax < 0.1:
    bmaxplot = int(100. + bmax*100.) / 100.
elif bmax >=0.1 and bmax < 1:
    bmaxplot = int(10. + bmax*10.) / 10.
else:
    bmaxplot = 1 + int(bmax)
    
clevs = [-0.3,0.01,0.03,0.05,0.07,0.1,0.12,0.14,0.17,0.2,0.25,0.3,0.5,0.75,1.0,1.5,2.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
fig = plt.figure(figsize=(6,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ bias background variance, '+cwarm+\
    '-season 2018,\n'+clead+'-hour forecast, ' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x, y, Bbeta_variance, clevs,\
    cmap=None, colors=colorst, extend='both')
xdot, ydot = m(rlon,rlat)
m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.1)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7) 
cb.set_label('Bias background error variance',fontsize=9)

# ---- set plot title

plot_title = 't2m_beta_variance_'+cwarm+'2018_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')





# --- now plot the bias spread

bmax = np.max(np.sqrt(Bbeta_variance))
if bmax < 0.1:
    bmaxplot = int(100. + bmax*100.) / 100.
elif bmax >=0.1 and bmax < 1:
    bmaxplot = int(10. + bmax*10.) / 10.
else:
    bmaxplot = 1 + int(bmax)
    
clevs = [-0.3,0.01,0.03,0.05,0.07,0.1,0.12,0.14,0.17,0.2,0.25,0.3,0.5,0.75,1.0,1.5,2.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
fig = plt.figure(figsize=(6,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'ECMWF T$_{2m}$ bias background spread, '+cwarm+\
    '-season 2018,\n'+clead+'-hour forecast, ' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x, y, np.sqrt(Bbeta_variance), clevs,\
    cmap=None, colors=colorst, extend='both')
xdot, ydot = m(rlon,rlat)
m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.1)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7) 
cb.set_label('Bias background error spread',fontsize=9)

# ---- set plot title

plot_title = 't2m_beta_spread_'+cwarm+'2018_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')




# --- make plots of covariance along cross sections

BX_covariance_SW_NE = np.zeros((41,41),dtype=np.float32)
BX_covariance_NW_SE = np.zeros((41,41),dtype=np.float32)        
BX_covariance_N_S = np.zeros((41,41),dtype=np.float32)         
BX_covariance_W_E = np.zeros((41,41),dtype=np.float32)

for ktr1 in range(-20,21):  
    for ktr2 in range(-20,21): 
        BX_covariance_SW_NE[ktr1,ktr2] = Bbeta_localized[ilat+ktr1, ilon+ktr1    ,ilat+ktr2, ilon+ktr2]
        BX_covariance_NW_SE[ktr1,ktr2] = Bbeta_localized[ilat-ktr1, ilon+ktr1    ,ilat-ktr2, ilon+ktr2]
        BX_covariance_N_S[ktr1,ktr2]   = Bbeta_localized[ilat-ktr1, ilon         ,ilat-ktr2, ilon]
        BX_covariance_W_E[ktr1,ktr2]   = Bbeta_localized[ilat     , ilon+ktr1    ,ilat     , ilon+ktr2]



for isection in range(4):
    if isection == 0:
        xtitle = 'SW to NE cross section through point'
        xptit = 'SW_to_NE_xsect_'
        xs = BX_covariance_SW_NE
    elif isection == 1:
        xtitle = 'NW to SE cross section through point'
        xptit = 'NW_to_SE_xsect_'
        xs = BX_covariance_NW_SE
    elif isection == 2:
        xtitle = 'North to south cross section through point'
        xptit = 'N_to_S_xsect_'
        xs = BX_covariance_N_S
    else:
        xtitle = 'West to East cross section through point'
        xptit = 'W_to_E_xsect_'
        xs = BX_covariance_W_E
    
    #clevs = [-0.99,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
    #colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    #    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']
                     
    colorst = ['#8FB3FF','#C4E8FF','#E4FFFF','White','#D8F9D8',\
        '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
        '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    clevs = [-0.16,-0.09,-0.04,-0.01,0.01,0.02,0.04,0.09,0.16,0.25,0.36,0.49,0.64,0.81,1.0]
        
    fig = plt.figure(figsize=(6,6.5))
    axloc = [0.12,0.17,0.84,0.75]
    ax1 = fig.add_axes(axloc)
    title = r'Bias error covariance '+xtitle+',\n'+cwarm+\
        '-season 2018, '+clead+'-hour forecast, ' + \
        ' lon = '+clon+', lat = '+clat
    ax1.set_title(title, fontsize=11,color='Black')
    
    CS2 = ax1.contourf(range(-20,21), range(-20,21),\
        xs,clevs,\
        cmap=None,colors=colorst,extend='both')
    xdot, ydot = m(rlon,rlat)
    ax1.plot([0,0],[0,0],marker='.',markersize=5,color='Black')
    ax1.set_xlabel('Grid point 1 relative coordinate')
    ax1.set_ylabel('Grid point 2 relative coordinate')
    
    # ---- use axes_grid toolkit to make colorbar axes.

    rcParams['legend.fontsize']='xx-small'
    divider = make_axes_locatable(ax1)
    #cax = divider.append_axes("bottom", size="3%", pad=0.01)
    cax = fig.add_axes([0.05,0.07,0.9,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7) 
    cb.set_label(r'Error covariance (deg C$^2$)')

    # ---- set plot title

    plot_title = xptit+'bias_cov_'+cwarm+'2018_f'+\
        clead+'_lon='+clon+'_lat='+clat+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    print ('Done!')
    
    
    


        
        


