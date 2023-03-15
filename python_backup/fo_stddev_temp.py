#!/usr/bin/env python

"""
marginal_distributions_temp.py 
"""
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import sys
import _pickle as cPickle
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
rcParams['legend.fontsize']='small'

cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'

clead = sys.argv[1]
cseason = sys.argv[2]
       
models = ['raw', 'dav', 'qmap', 'mos', 'mvmos']  

# --- read land-surface mask, lons, lats

infilename = 'lsmask_0p5.cPick'
print (infilename)
inf = open(infilename, 'rb')
lsmask = cPickle.load(inf)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
inf.close()

# ---- populate the histogram of correlations, std devs


for imodel, model in enumerate(models):
          
    # ---- read in the standardized anomaly and correlation from files
    #      prepared for Taylor diagram generation
    
    if model == 'raw':
        ftype = 'raw_forecast_errorstats_2019_'
    elif model ==  'dav':
        ftype = 'decayavg_forecast_errorstats_2019_'  
    elif model == 'qmap':
        ftype = 'quantile_mapping_errorstats_2019_'  
    elif model == 'mos':
        ftype = 'MOS_GEFSv12_2019_'
    elif model == 'mvmos':
        ftype = 'MOS_mvr_GEFSv12_2019_'

    statsfile = cpath_errorstats + ftype + clead+'h_taylor.cPick'
    inf = open(statsfile, 'rb')
    corr_winter = cPickle.load(inf)
    corr_spring = cPickle.load(inf)
    corr_summer = cPickle.load(inf)
    corr_fall = cPickle.load(inf)
    std_normalized_winter = cPickle.load(inf)
    std_normalized_spring = cPickle.load(inf)
    std_normalized_summer = cPickle.load(inf)
    std_normalized_fall = cPickle.load(inf)
    std_obs_winter = cPickle.load(inf)
    std_obs_spring = cPickle.load(inf)
    std_obs_summer = cPickle.load(inf)
    std_obs_fall = cPickle.load(inf)
    
    if model == 'raw':
        nlats, nlons = np.shape(corr_winter)
        corr_winter_multimodel = np.zeros((nlats, nlons,5), dtype=np.float32)
        corr_spring_multimodel = np.zeros((nlats, nlons,5), dtype=np.float32)
        corr_summer_multimodel = np.zeros((nlats, nlons,5), dtype=np.float32)
        corr_fall_multimodel = np.zeros((nlats, nlons,5), dtype=np.float32)
        std_normalized_winter_multimodel = np.zeros((nlats, nlons,5), dtype=np.float32)
        std_normalized_spring_multimodel = np.zeros((nlats, nlons,5), dtype=np.float32)
        std_normalized_summer_multimodel = np.zeros((nlats, nlons,5), dtype=np.float32)
        std_normalized_fall_multimodel = np.zeros((nlats, nlons,5), dtype=np.float32)
    inf.close()
    
    corr_winter_multimodel[:,:,imodel] = corr_winter[:,:]
    corr_spring_multimodel[:,:,imodel] = corr_spring[:,:] 
    corr_summer_multimodel[:,:,imodel] = corr_summer[:,:] 
    corr_fall_multimodel[:,:,imodel] = corr_fall[:,:]
    std_normalized_winter_multimodel[:,:,imodel] = std_normalized_winter[:,:] 
    std_normalized_spring_multimodel[:,:,imodel] = std_normalized_spring[:,:] 
    std_normalized_summer_multimodel[:,:,imodel] = std_normalized_summer[:,:] 
    std_normalized_fall_multimodel[:,:,imodel] = std_normalized_fall[:,:] 
   
   
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],\
    urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],\
    resolution='l', projection='mill')
x, y = m(lons, lats)
clevs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.83,0.86,0.9,0.92,0.94,0.96,0.98]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
    

# ---- spatial plot of the F-A correlations 
fig1 = plt.figure(figsize=(9.,6.))

fig1.suptitle('Forecast vs. analyzed correlation, '+cseason+' '+\
    clead+' h', size='x-large')
for imodel, model in enumerate(models):
        
    if cseason == 'Jan-Feb-Mar':
        corrcoef_multimodel = corr_winter_multimodel
    elif cseason == 'Apr-May-Jun':
        corrcoef_multimodel = corr_spring_multimodel
    elif cseason == 'Jul-Aug-Sep':
        corrcoef_multimodel = corr_summer_multimodel
    elif cseason == 'Oct-Nov-Dec':
        corrcoef_multimodel = corr_fall_multimodel
        
    if model == 'raw':
        title = '(a) Raw model output'
        axlocs = [0.02,0.54,0.46,0.35]
        corrcoef = corrcoef_multimodel[:,:,0]
    elif model ==  'dav': 
        title = '(b) Decaying-average bias correction'
        cletter = '(b) '
        axlocs = [0.52,0.54,0.46,0.35]
        corrcoef = corrcoef_multimodel[:,:,1]
    elif model == 'qmap':
        title = '(c) Quantile mapping'
        axlocs = [0.02,0.13,0.46,0.35]
        corrcoef = corrcoef_multimodel[:,:,2]
    elif model == 'mvmos':
        title = '(d) Multi-variate MOS'
        axlocs = [0.52,0.13,0.46,0.35]   
        corrcoef = corrcoef_multimodel[:,:,4]  
    
    # --- add to correlation plots
    
    ax1 = fig1.add_axes(axlocs)
    ax1.set_title(title, fontsize=15)
    CS2 = m.contourf(x,y,corrcoef,clevs,cmap=None,\
        colors=colorst,extend='both')        
    ax1.set_title(title, fontsize=11,color='Black')
    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')

ax1 = [0.02,0.09,0.96,0.02]
cax = fig1.add_axes(ax1)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Forecast vs. analyzed correlation',fontsize=9)
    
# ---- set plot title

plot_title = 'fo_correlation_multicorrection_lead'+\
    clead+'h_season='+cseason+'.png'
fig1.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



# ---- spatial plot of the normalized standard deviations

clevs = [0.5,0.7,0.8,0.9,0.95, 1.05,1.1,1.2,1.5,2.0]

colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']

fig1 = plt.figure(figsize=(9.,6.))
fig1.suptitle(r'$\sigma_f$ / $\sigma_a$ '+cseason+' '+\
    clead+' h', size='x-large')
for imodel, model in enumerate(models):
        
    if cseason == 'Jan-Feb-Mar':
        std_normalized_multimodel = std_normalized_winter_multimodel
    elif cseason == 'Apr-May-Jun':
        std_normalized_multimodel = std_normalized_spring_multimodel
    elif cseason == 'Jul-Aug-Sep':
        std_normalized_multimodel = std_normalized_summer_multimodel
    elif cseason == 'Oct-Nov-Dec':
        std_normalized_multimodel = std_normalized_fall_multimodel
        
    if model == 'raw':
        title = '(a) Raw model output'
        axlocs = [0.02,0.54,0.46,0.35]
        std_normalized = std_normalized_multimodel[:,:,0]
    elif model ==  'dav': 
        title = '(b) Decaying-average bias correction'
        cletter = '(b) '
        axlocs = [0.52,0.54,0.46,0.35]
        std_normalized = std_normalized_multimodel[:,:,1]
    elif model == 'qmap':
        title = '(c) Quantile mapping'
        axlocs = [0.02,0.13,0.46,0.35]
        std_normalized = std_normalized_multimodel[:,:,2]
    elif model == 'mvmos':
        title = '(d) Multi-variate MOS'
        axlocs = [0.52,0.13,0.46,0.35]   
        std_normalized = std_normalized_multimodel[:,:,4]  
    
    # --- add to correlation plots
    
    ax1 = fig1.add_axes(axlocs)
    ax1.set_title(title, fontsize=15)
    CS2 = m.contourf(x,y,std_normalized,clevs,cmap=None,\
        colors=colorst,extend='both')        
    ax1.set_title(title, fontsize=11,color='Black')
    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')

ax1 = [0.02,0.09,0.96,0.02]
cax = fig1.add_axes(ax1)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label(r'$\sigma_f$ / $\sigma_a$',fontsize=9)
    
# ---- set plot title

plot_title = 'sigmaf_sigmao_ratio_multicorrection_lead'+\
    clead+'h_season='+cseason+'.png'
fig1.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


