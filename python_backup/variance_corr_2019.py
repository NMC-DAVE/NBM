"""

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

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
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
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats) 
        complete_fcst = np.zeros((ndates,nlats,nlons), dtype=np.float32)
        complete_anal = np.zeros((ndates,nlats,nlons), dtype=np.float32)
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

# --- compute covariance statistics

print ('computing covariance statistics')
print ('nlons, nlats = ', nlons, nlats)
cov_bias = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float32)
cov_random = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float32)
var_bias = np.zeros((nlats, nlons), dtype=np.float32)
var_random = np.zeros((nlats, nlons), dtype=np.float32)
for j1 in range(nlats):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('processing j1 = ',j1,' of ',nlats,'. Current time = ',current_time)
    for i1 in range(nlons): 
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('processing i1 = ',i1,' of ',nlons,'. Current time = ',current_time)
        a1 = bias_estimate[:,j1,i1]
        a2 = random_error_estimate[:,j1,i1]
        for j2 in range(nlats):
            for i2 in range(nlons):
                b1 = bias_estimate[:,j2,i2]
                b2 = random_error_estimate[:,j2,i2]
                cov = np.cov(a1,b1)
                var_bias[j1,i1] = cov[0,0]
                cov_bias[j1,i1,j2,i2] = cov[1,0]
                cov = np.cov(a2,b2)
                var_random[j1,i1] = cov[0,0]
                cov_random[j1,i1,j2,i2] = cov[1,0]
                
# ---- write to cPickle file

outfile = 'covstats_bias_random_ecmwf2019_lead='+clead+'.cPick'
ouf = open(outfile,'wb')
cPickle.dump(cov_bias, ouf)
cPickle.dump(cov_random , ouf)
cPickle.dump(var_bias, ouf)
cPickle.dump(var_random, ouf)
cPickle.dump(lats, ouf)
cPickle.dump(lons, ouf)
ouf.close()              

sys.exit()

rmse = np.sqrt(np.sum((complete_fcst - complete_anal)**2,  axis=0)/ndates)
analvar = np.ones((nlats, nlons), dtype=np.float32)
random_error = np.sqrt(rmse**2 - analvar - bias**2)
    
# --- make plot    

# --- now plot the bias

clevs = [-3,-2,-1,-0.5,0.5,1,2,3] 
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.1,0.96,0.82]
ax1 = fig.add_axes(axloc)
ax1.set_title('ECMWF forecast bias, Nov-Dec 2018, '+clead+'-hour forecast', fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,bias,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature bias (deg C)')

# ---- set plot title

plot_title = 't2m_bias_NovDec2018_f'+clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now plot the rmse

clevs = [0.0,0.5,1.0,1.5,2.0,2.5,3,3.5,4.0,5.0,6.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.1,0.96,0.82]
ax1 = fig.add_axes(axloc)
ax1.set_title('ECMWF forecast RMSE, Nov-Dec 2018, '+clead+'-hour forecast', fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,rmse,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature RMSE (deg C)')

# ---- set plot title

plot_title = 't2m_RMSE_NovDec2018_f'+clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now plot the random error

clevs = [0.0,0.5,1.0,1.5,2.0,2.5,3,3.5,4.0,5.0,6.0] 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
fig = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.1,0.96,0.82]
ax1 = fig.add_axes(axloc)
ax1.set_title('ECMWF forecast random error, Nov-Dec 2018, '+clead+'-hour forecast', fontsize=16,color='Black')
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,random_error,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature RMSE (deg C)')

# ---- set plot title

plot_title = 't2m_random_NovDec2018_f'+clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')
