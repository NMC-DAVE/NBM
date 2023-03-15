"""
plot_GEFSv12_gain_corr.py

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
cmonth = sys.argv[2]
clon = sys.argv[3]
clat = sys.argv[4]
rlon = float(clon)
rlat = float(clat)
ilead = int(clead)
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_beta = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_random = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_gain = '/Volumes/Backup Plus/gefsv12/t2m/'


# ---- read the Bbeta_localized 2D to pickle file.

infile = cpath_beta+'Localized_Bbeta_4D_'+cmonth+\
    '_lead='+clead+'.cPick'
print ('   reading from ', infile)
inf = open(infile,'rb')
Bbeta_localized_4D = cPickle.load(inf)
inf.close()
    
# ---- write the Brandom_localized 4D to pickle file.

infile = cpath_random+'Localized_Brandom_4D_'+cmonth+\
    '_lead='+clead+'.cPick'
print ('   reading from ', infile)
inf = open(infile,'rb')
Brandom_localized_4D = cPickle.load(inf)
inf.close()
        
gain_infile = cpath_gain + 'GEFSv12_KFgain_'+cmonth+'_lead'+clead+'.cPick'
print ('   reading Kalman gain from ', gain_infile)
inf = open(gain_infile, 'rb')
Kalman_gain_beta_4D = cPickle.load(inf)
inf.close()


# --- read in a single reanalysis in order to get lats and lons.  Find the 
#     index of the nearest grid point.

date_list_forecast = ['2018010100']
analyses_3d, lats, lons = read_reanalysis_timeseries(cpath_era5, \
    date_list_forecast)
nlats, nlons = np.shape(lons)
print (nlats,nlons)
print ('min, max lons ', np.min(lons), np.max(lons))
print ('min, max lats ', np.min(lats), np.max(lats))
lons = lons - 360.
lats_1d = lats[:,0]
print (lats_1d)
lons_1d = lons[0,:]
ilon = find_nearest(lons_1d, rlon)
ilat = find_nearest(lats_1d, rlat)
print (ilat,ilon)

Kalman_gain_beta_map = Kalman_gain_beta_4D[:,:,ilat,ilon]
Kalman_gain_beta_map2 = Kalman_gain_beta_4D[ilat,ilon,:,:]
print ('min, max Kalman_gain_beta_max = ',\
    np.min(Kalman_gain_beta_map), np.max(Kalman_gain_beta_map) )
    
Kalman_gain_beta_map2_1deg = np.zeros((30,65), dtype=np.float32)
for j in range(0,nlats-2,2):
    for i in range(0,nlats-2,2):
        Kalman_gain_beta_map2_1deg[j//2,i//2] = \
           np.sum(Kalman_gain_beta_map2[j:j+2,i:i+2])/4
lats_1deg = lats[0:-1:2,0:-1:2]
lons_1deg = lons[0:-1:2,0:-1:2] 
    
Kmax = np.max(np.abs(Kalman_gain_beta_map)) 
Kmax2 = np.max(np.abs(Kalman_gain_beta_map2_1deg))  

# ---- extract forecast error variance and bias correction error variance

Bx_variance = np.zeros((nlats,nlons), dtype=np.float32)
Bbeta_variance = np.zeros((nlats,nlons), dtype=np.float32)

for j in range(nlats):
    for i in range(nlons):
        Bx_variance[j,i] = Brandom_localized_4D[j,i,j,i]
        Bbeta_variance[j,i] = Bbeta_localized_4D[j,i,j,i]
        

# ---- compute 1-point correlation maps

bias_error_corr_map = np.zeros((nlats, nlons), dtype=np.float32)
bias_error_corr_map2 = np.zeros((nlats, nlons), dtype=np.float32)
bias_error_cov_map = np.zeros((nlats, nlons), dtype=np.float32)
for i in range(nlons):
    for j in range(nlats):
        bias_error_corr_map[j,i] = Bbeta_localized_4D[ilat,ilon,j,i] / \
            ( np.sqrt(Bbeta_variance[ilat,ilon]) * np.sqrt(Bbeta_variance[j,i]) )
        bias_error_corr_map2[j,i] = Bbeta_localized_4D[j,i,ilat,ilon] / \
            ( np.sqrt(Bbeta_variance[ilat,ilon]) * np.sqrt(Bbeta_variance[j,i]) )
        bias_error_cov_map[j,i] = Bbeta_localized_4D[ilat,ilon,j,i]
Bmax = np.max(Bbeta_variance)
        
forecast_error_corr_map = Brandom_localized_4D[:,:,ilat,ilon] / \
    (np.sqrt(Bx_variance[ilat,ilon])*np.sqrt(Bx_variance[:,:]))


# --- now plot the localized bias error correlation for the selected point 

clevs = [-0.99,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000'] 
    
fig = plt.figure(figsize=(6.,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 T$_{2m}$ random-error correlation map, '+\
    cmonth+', '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,forecast_error_corr_map,clevs,\
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
cb.set_label('Temperature random error correlation')

# ---- set plot title

plot_title = 't2m_GEFSv12_forecast_error_correlation_'+cmonth+'_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300,fontsize=9)
print ('saving plot to file = ',plot_title)
print ('Done!')



                    
# --- now plot the localized bias error correlation for the selected point 

clevs = [-0.99,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000'] 
    
fig = plt.figure(figsize=(6.,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 T$_{2m}$ bias correlation map, '+\
    cmonth+', '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
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

plot_title = 't2m_GEFSv12_bias_error_correlation_'+cmonth+'_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300,fontsize=9)
print ('saving plot to file = ',plot_title)
print ('Done!')
print ('bias_err_correlation at location = ',bias_error_corr_map[ilat-2:ilat+2,ilon])



# --- now plot the localized bias error correlation for the selected point 

clevs = [-0.99,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000'] 
    
fig = plt.figure(figsize=(6.,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 T$_{2m}$ bias correlation map (flipped), '+\
    cmonth+', '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,bias_error_corr_map2,clevs,\
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
cb.set_label('Temperature bias correlation (flipped)')

# ---- set plot title

plot_title = 't2m_GEFSv12_bias_error_correlation_flipped_'+cmonth+'_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300,fontsize=9)
print ('saving plot to file = ',plot_title)
print ('Done!')
print ('bias_err_correlation at location = ',bias_error_corr_map2[ilat-2:ilat+2,ilon])





# --- now plot the localized bias error covariance for the selected point 

clevs = [-0.3,0.01,0.03,0.05,0.07,0.1,0.12,0.14,0.17,0.2,0.25,0.3,0.5,0.75,1.0,1.5,2.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

print ('min, max covariance = ', np.min(bias_error_cov_map), \
    np.max(bias_error_cov_map))    
fig = plt.figure(figsize=(6.,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 T$_{2m}$ bias covariance map, '+cmonth+\
    ', '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
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

plot_title = 't2m_GEFSv12_bias_error_covariance_'+cmonth+'_f'+\
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
title = r'GEFSv12 T$_{2m}$ Kalman gain map, '+cmonth+\
    ', '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
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

plot_title = 't2m_GEFSv12_Kalman_gain_'+cmonth+'_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')

# --- now plot the Kalman gain map2 for the selected point

clevs = [-0.99,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
    'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']    
    
fig = plt.figure(figsize=(6,4.2))
axloc = [0.02,0.09,0.96,0.82]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 T$_{2m}$ data-impact on center point, '+cmonth+\
    ', '+clead+'-hour forecast,\n' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons_1deg, lats_1deg)
CS2 = m.contourf(x,y,Kalman_gain_beta_map2_1deg/Kmax2,clevs,\
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
cb.set_label('Data impact / max(abs(data impact))',fontsize=9)

# ---- set plot title

plot_title = 't2m_GEFSv12_Kalman_gain2_'+cmonth+'_f'+\
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
title = r'GEFSv12 T$_{2m}$ forecast-error standard deviation, '+cmonth+\
    ',\n'+clead+'-hour forecast, ' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,np.sqrt(Bx_variance),clevs,\
    cmap=None,colors=colorst,extend='both')
#xdot, ydot = m(rlon,rlat)
#m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
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

plot_title = 't2m_GEFSv12_forecast_spread_'+cmonth+'_f'+\
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
title = r'GEFSv12 T$_{2m}$ bias background variance, '+cmonth+\
    ',\n'+clead+'-hour forecast, ' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
x, y = m(lons, lats)
CS2 = m.contourf(x, y, Bbeta_variance, clevs,\
    cmap=None, colors=colorst, extend='both')
#xdot, ydot = m(rlon,rlat)
#m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
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

plot_title = 't2m_GEFSv12_beta_variance_'+cmonth+'2018_f'+\
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
title = r'GEFSv12 T$_{2m}$ bias background spread, '+cmonth+\
    ',\n'+clead+'-hour forecast, ' + \
    ' lon = '+clon+', lat = '+clat
ax1.set_title(title, fontsize=11,color='Black')
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

plot_title = 't2m_GEFSv12_beta_spread_'+cmonth+'2018_f'+\
    clead+'_lon='+clon+'_lat='+clat+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')





    
    
    


        
        


