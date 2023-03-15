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
rcParams['legend.fontsize']='small'

cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'

clead = sys.argv[1]

# Reference std
stdrefs = {"Jan-Feb-Mar":1.0,   "Apr-May-Jun":1.0,\
    "Jul-Aug-Sep":1.0, "Oct-Nov-Dec":1.0}
       
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

xstd = np.arange(0.5,1.501,0.01)
xcorr = np.arange(0.0,1.001,0.01)
nxstd = len(xstd)
nxcorr = len(xcorr)
corr_winter_hist = ma.zeros((nxcorr,5), dtype=np.float64)
corr_spring_hist = ma.zeros((nxcorr,5), dtype=np.float64)
corr_summer_hist = ma.zeros((nxcorr,5), dtype=np.float64)
corr_fall_hist   = ma.zeros((nxcorr,5), dtype=np.float64)
std_winter_hist = ma.zeros((nxstd,5), dtype=np.float64)
std_spring_hist = ma.zeros((nxstd,5), dtype=np.float64)
std_summer_hist = ma.zeros((nxstd,5), dtype=np.float64)
std_fall_hist   = ma.zeros((nxstd,5), dtype=np.float64)

for imodel, model in enumerate(models):
          
    # ---- read in the standardized anomaly and correlation from files
    #      prepared for Taylor diagram generation
    
    if model == 'raw':
        ftype = 'raw_forecast_errorstats_2019_'
        title = 'Raw\nmodel output\n+'+clead+' h'
        color = 'Red'
    elif model ==  'dav':
        ftype = 'decayavg_forecast_errorstats_2019_'  
        title = 'Decaying-average\nbias correction\n+'+clead+' h'
        color = 'Peru'
    elif model == 'qmap':
        ftype = 'quantile_mapping_errorstats_2019_'  
        title = 'Quantile\nmapping\n+'+clead+' h'
        color = 'Blue'
    elif model == 'mos':
        ftype = 'MOS_GEFSv12_2019_'
        title = 'Univariate\nMOS\n+'+clead+' h'
        color = 'Green'
    elif model == 'mvmos':
        ftype = 'MOS_mvr_GEFSv12_2019_'
        title = 'Multi-variate\nMOS\n+'+clead+' h'
        color = 'Gray'

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
    inf.close()
    nlats, nlons = np.shape(corr_fall)

    for jy in range(nlats):
        for ix in range(nlons):
            if lsmask[jy,ix] == 1:
                corr_idx = np.int(corr_winter[jy,ix]*100.0)
                corr_winter_hist[corr_idx,imodel] = corr_winter_hist[corr_idx,imodel] + 1.0
                corr_idx = np.int(corr_spring[jy,ix]*100.0)
                corr_spring_hist[corr_idx,imodel] = corr_spring_hist[corr_idx,imodel] + 1.0
                corr_idx = np.int(corr_summer[jy,ix]*100.0)
                corr_summer_hist[corr_idx,imodel] = corr_summer_hist[corr_idx,imodel] + 1.0
                corr_idx = np.int(corr_winter[jy,ix]*100.0)
                corr_fall_hist[corr_idx,imodel] = corr_fall_hist[corr_idx,imodel] + 1.0
        
                std_idx = np.amin([nxstd-1, np.int(std_normalized_winter[jy,ix]*100.0 - 50.0)])
                std_winter_hist[std_idx,imodel] = std_winter_hist[std_idx,imodel] + 1.0
                std_idx = np.amin([nxstd-1, np.int(std_normalized_spring[jy,ix]*100.0 - 50.0)])
                std_spring_hist[std_idx,imodel] = std_spring_hist[std_idx,imodel] + 1.0
                std_idx = np.amin([nxstd-1, np.int(std_normalized_summer[jy,ix]*100.0 - 50.0)])
                std_summer_hist[std_idx,imodel] = std_summer_hist[std_idx,imodel] + 1.0
                std_idx = np.amin([nxstd-1, np.int(std_normalized_fall[jy,ix]*100.0 - 50.0)])
                std_fall_hist[std_idx,imodel] = std_fall_hist[std_idx,imodel] + 1.0


fig1 = plt.figure(figsize=(8.5,8))
fig1.suptitle('Histograms of forecast vs. analyzed correlation, +'+clead+' h', size='xx-large')

# ---- plot the histograms in fig1

for iseason, season in enumerate(['Jan-Feb-Mar','Apr-May-Jun','Jul-Aug-Sep','Oct-Nov-Dec']):
    for imodel, model in enumerate(models):
        
        if season == 'Jan-Feb-Mar':
            corrcoef = corr_winter_hist[:,imodel]
            cletter = '(a) '
            axlocs = [0.08,0.54,0.4,0.36]
        elif season == 'Apr-May-Jun':
            corrcoef = corr_spring_hist[:,imodel]
            cletter = '(b) '
            axlocs = [0.57,0.54,0.4,0.36]
        elif season == 'Jul-Aug-Sep':
            corrcoef = corr_summer_hist[:,imodel]
            cletter = '(c) '
            axlocs = [0.08,0.06,0.4,0.36]
        elif season == 'Oct-Nov-Dec':
            corrcoef = corr_fall_hist[:,imodel]
            cletter = '(d) '  
            axlocs = [0.57,0.06,0.4,0.36]
            print ('iseason, imodel, corr = ', iseason, imodel, corrcoef)
         
        if model == 'raw':
            ftype = 'raw_forecast_errorstats_2019_'
            label = 'Raw model output'
            color = 'Red'
            linestyle = '--'
        elif model ==  'dav':
            ftype = 'decayavg_forecast_errorstats_2019_'  
            label = 'Decaying-average'
            color = 'Peru'
            linestyle = '-'
        elif model == 'qmap':
            ftype = 'quantile_mapping_errorstats_2019_'  
            label = 'Quantile mapping'
            color = 'Blue'
            linestyle = '-.'
        elif model == 'mos':
            ftype = 'MOS_GEFSv12_2019_'
            label = 'Univariate MOS'
            color = 'Green'
            linestyle = '-'
        elif model == 'mvmos':
            ftype = 'MOS_mvr_GEFSv12_2019_'
            label = 'Multi-variate MOS'
            color = 'Gray'  
            linestyle = ':'       
    
        # --- add to correlation plots
        
        
        corrcoef = corrcoef / np.sum(corrcoef)
        corrcoef = ma.masked_where(corrcoef == 0.0, corrcoef)
        if imodel == 0: ax1 = fig1.add_axes(axlocs)
        if imodel == 0: ax1.set_title(cletter+season, fontsize=15)
        if season == 'Jul-Aug-Sep':
            ax1.plot(xcorr, corrcoef, color=color,linestyle=linestyle,lw=1.5,label=label)
        else:
            ax1.plot(xcorr, corrcoef, color=color,linestyle=linestyle,lw=1.5)
        if imodel == 0: ax1.set_xlim(0.7,1.0)
        if imodel == 0: ax1.set_xticks([0.7,0.75,0.8,0.85,0.9,0.95,1.0])
        #if imodel == 0: ax1.set_ylim(0.0,0.2)
        if imodel == 0: ax1.set_xlabel('Correlation',fontsize=12)
        if imodel == 0: ax1.set_ylabel('Relative frequency',fontsize=12)
        #if imodel == 0: ax1.set_yticks([0.001,0.005, 0.01,0.05, 0.1, 0.5])
        #if imodel == 0: ax1.set_yticks([0.0,0.04,0.08,0.12,0.16,0.2])
        if imodel == 0: ax1.grid(True, lw=0.25)
        
        if season == 'Jul-Aug-Sep' and model == 'mvmos':
            ax1.legend(loc=0)
    
# ---- set plot titles and write to files

outfile1 = 'corr_histogram_4panel_'+clead+'h.pdf'
print ('saving to ', outfile1)
plt.savefig(outfile1)

fig2 = plt.figure(figsize=(8.5,8))
fig2.suptitle('Histograms of forecast vs. analyzed normalized standard deviation, +'+clead+' h', size='x-large')
  
# ---- plot the histograms in fig1

for iseason, season in enumerate(['Jan-Feb-Mar','Apr-May-Jun','Jul-Aug-Sep','Oct-Nov-Dec']):

    for imodel, model in enumerate(models):
        
        if season == 'Jan-Feb-Mar':
            stddev = std_winter_hist[:,imodel]
            cletter = '(a) '
            axlocs = [0.08,0.54,0.4,0.37]
        elif season == 'Apr-May-Jun':
            stddev = std_spring_hist[:,imodel]
            cletter = '(b) '
            axlocs = [0.58,0.54,0.4,0.37]
        elif season == 'Jul-Aug-Sep':
            stddev = std_summer_hist[:,imodel]
            cletter = '(c) '
            axlocs = [0.08,0.06,0.4,0.37]
        elif season == 'Oct-Nov-Dec':
            stddev = std_fall_hist[:,imodel] 
            cletter = '(d) '  
            axlocs = [0.58,0.06,0.4,0.37]
            
            
        if model == 'raw':
            ftype = 'raw_forecast_errorstats_2019_'
            label = 'Raw model output +'+clead+' h'
            color = 'Red'
            linestyle = '--'
        elif model ==  'dav':
            ftype = 'decayavg_forecast_errorstats_2019_'  
            label = 'Decaying-average'
            color = 'Peru'
            linestyle = '-'
        elif model == 'qmap':
            ftype = 'quantile_mapping_errorstats_2019_'  
            label = 'Quantile mapping '+clead+' h'
            color = 'Blue'
            linestyle = '-.'
        elif model == 'mos':
            ftype = 'MOS_GEFSv12_2019_'
            label = 'Univariate MOS '+clead+' h'
            color = 'Green'
            linestyle = '-'
        elif model == 'mvmos':
            ftype = 'MOS_mvr_GEFSv12_2019_'
            label = 'Multi-variate MOS '+clead+' h'
            color = 'Gray'  
            linestyle = ':'
        
        # --- add to stddev plots
    
        stddev = stddev / np.sum(stddev)
        stddev = ma.masked_where(stddev == 0.0, stddev)
        if imodel == 0: ax2 = fig2.add_axes(axlocs)
        if imodel == 0: ax2.set_title(cletter+season, fontsize=14)
        if season == 'Jan-Feb-Mar':
            ax2.semilogy(xstd, stddev, color=color,linestyle=linestyle,lw=1.5,label=label)
        else:
            ax2.semilogy(xstd, stddev, color=color,linestyle=linestyle,lw=1.5)
        if imodel == 0: ax2.set_xlim(0.5,1.5)
        if imodel == 0: ax2.set_ylim(0.001,0.1)
        if imodel == 0: ax2.set_xlabel(r'$\sigma_f$ / $\sigma_a$',fontsize=12)
        if imodel == 0: ax2.set_ylabel('Relative frequency',fontsize=12)
        if imodel == 0: ax2.set_xticks([0.5,0.75,1.0,1.25,1.5])
        if imodel == 0: ax2.set_yticks([0.001,0.005, 0.01,0.05, 0.1])
        if imodel == 0: ax2.grid(True, lw=0.25)
        
        if season == 'Jan-Feb-Mar' and model == 'mvmos':
            ax2.legend(loc=3)
    
outfile2 = 'stddev_histogram_4panel_'+clead+'h.pdf'
print ('saving to ', outfile2)
plt.savefig(outfile2)
