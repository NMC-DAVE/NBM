#!/usr/bin/env python
# taylor_4panel_multimethod.py 

__version__ = "Time-stamp: <2018-12-06 11:55:22 ycopin>"
__author__ = "Yannick Copin <yannick.copin@laposte.net>"

"""
Example of use of TaylorDiagram. Illustration dataset courtesy of Michael
Rawlins.
Rawlins, M. A., R. S. Bradley, H. F. Diaz, 2012. Assessment of regional climate
model simulation estimates over the Northeast United States, Journal of
Geophysical Research (2012JGRD..11723112R).
"""

from taylorDiagram import TaylorDiagram
import numpy as NP
import matplotlib.pyplot as PLT
import sys
import _pickle as cPickle
from matplotlib import rcParams
import numpy as np
import numpy.ma as ma
rcParams['legend.fontsize']='small'

cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'

clead = sys.argv[1]
season = 'Jan-Feb-Mar' #'Jul-Aug-Sep'
# Reference std
stdrefs = {"Jan-Feb-Mar":1.0,   "Apr-May-Jun":1.0,\
    "Jul-Aug-Sep":1.0, "Oct-Nov-Dec":1.0}

fig = PLT.figure(figsize=(10,8))
fig.suptitle(season+'\n+'+clead+' h', size='xx-large')

# --- read land-surface mask, lons, lats

infilename = 'lsmask_0p5.cPick'
print (infilename)
inf = open(infilename, 'rb')
lsmask = cPickle.load(inf)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
inf.close()
      
for imodel, model in enumerate(['raw','dav', 'qmap','mvmos']):    
               
    if model == 'raw':
        ftype = 'raw_forecast_errorstats_2019_'
        title = 'Raw model output'
        cletter = '(a) '
        rect = 221
    elif model ==  'dav':
        ftype = 'decayavg_forecast_errorstats_2019_'  
        title = 'DAV'
        cletter = '(b) '
        rect = 222        
    elif model == 'qmap':
        ftype = 'quantile_mapping_errorstats_2019_'  
        title = 'QM'
        cletter = '(c) '
        rect = 223
    elif model == 'mos':
        ftype = 'MOS_GEFSv12_2019_'
        title = 'uMOS'
    elif model == 'mvmos':
        ftype = 'MOS_mvr_GEFSv12_2019_'
        title = 'mvMOS'
        cletter = '(d) '  
        rect = 224
        
    statsfile = cpath_errorstats + ftype + clead+'h_taylor.cPick'
    inf = open(statsfile, 'rb')
    corr_winter = cPickle.load(inf)
    corr_spring = cPickle.load(inf)
    corr_summer = cPickle.load(inf)
    corr_fall = cPickle.load(inf)
    print ('max, min corr_fall  = ', NP.max(corr_fall ), NP.min(corr_fall ))
    std_normalized_winter = cPickle.load(inf)
    std_normalized_spring = cPickle.load(inf)
    std_normalized_summer = cPickle.load(inf)
    std_normalized_fall = cPickle.load(inf)
    std_obs_winter = cPickle.load(inf)
    std_obs_spring = cPickle.load(inf)
    std_obs_summer = cPickle.load(inf)
    std_obs_fall = cPickle.load(inf)

    inf.close()
    nlats, nlons = NP.shape(corr_fall)

    if season == 'Jan-Feb-Mar':
        corrcoef = corr_winter
        stddev = std_normalized_winter
        stdobs = std_obs_winter
    elif season == 'Apr-May-Jun':
        corrcoef = corr_spring
        stddev = std_normalized_spring
        stdobs = std_obs_spring
    elif season == 'Jul-Aug-Sep':
        corrcoef = corr_summer
        stddev = std_normalized_summer
        stdobs = std_obs_summer
    elif season == 'Oct-Nov-Dec':
        corrcoef = corr_fall
        stddev = std_normalized_fall 
        stdobs = std_obs_fall  
        
    stdobs_masked = ma.masked_where(lsmask==0, stdobs)
    stdobs_sort = ma.sort(stdobs_masked.flatten())
    print (np.shape(stdobs_sort))
    nstd = np.sum(lsmask).astype(int)
    print (nstd//4, nstd//2, 3*nstd//4, nstd)
    stdobs_quarter = stdobs_sort[nstd//4]
    stdobs_half = stdobs_sort[nstd//2]
    stdobs_threequarter = stdobs_sort[3*nstd//4]
    print ('nstd stdobs_quarter stdobs_half stdobs_threequarter = ',\
        nstd, stdobs_quarter, stdobs_half, stdobs_threequarter )

    meancoef = NP.mean(corrcoef)
    meanstddev = NP.mean(stddev)
    dia = TaylorDiagram(1.0, fig=fig, rect=rect, \
        label=season, srange=(0.0, 2.0))
    
    # Add samples to Taylor diagram
    
    for jlat in range(nlats):
        for ilon in range(nlons):
            if lsmask[jlat,ilon] == 1:
                if stdobs[jlat,ilon] <= stdobs_quarter: # 2.0:
                    dia.add_sample(stddev[jlat,ilon], corrcoef[jlat,ilon],\
                        marker='o',fillstyle='full',markersize=0.5, \
                        markerfacecolor='Red', markeredgecolor='Red',zorder=4) #, label=name)
                #elif stdobs[jlat,ilon] > 2.0 and stdobs[jlat,ilon] <= 4.0:
                elif stdobs[jlat,ilon] > stdobs_quarter and stdobs[jlat,ilon] <= stdobs_half:
                    dia.add_sample(stddev[jlat,ilon], corrcoef[jlat,ilon],\
                        marker='o',fillstyle='full',markersize=0.5, \
                        markerfacecolor='RoyalBlue', markeredgecolor='RoyalBlue',\
                        zorder=5) #, label=name)
                #elif stdobs[jlat,ilon] > 4.0 and stdobs[jlat,ilon] <= 6.0:
                elif stdobs[jlat,ilon] >= stdobs_half and stdobs[jlat,ilon] < stdobs_threequarter:
                    dia.add_sample(stddev[jlat,ilon], corrcoef[jlat,ilon],\
                        marker='o',fillstyle='full',markersize=0.5, \
                        markerfacecolor='Violet', markeredgecolor='Violet',\
                        zorder=6) #, label=name)
                else:
                    dia.add_sample(stddev[jlat,ilon], corrcoef[jlat,ilon],\
                        marker='o',fillstyle='full',markersize=0.5, \
                        markerfacecolor='DarkGoldenRod', markeredgecolor='DarkGoldenRod',\
                        zorder=7) #, label=name)
    
    dia.add_sample(meanstddev, meancoef,marker='o',\
        markersize=3.5, markerfacecolor='Blue',\
        markeredgecolor='Blue',zorder=20) #, label=name)
        
    dia.add_sample([0.0,2.0], [0.2,0.2], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.4,0.4], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.6,0.6], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.7,0.7], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.8,0.8], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.9,0.9], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.95,0.95], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.99,0.99], 'k-',lw=0.25)
    
    if model == 'raw':  
        dia.add_sample(1.0, 1.0, marker='o',markersize=2, markerfacecolor='Red', \
            #markeredgecolor='Red',zorder=1,lw=0, label=r'$\sigma_{a} < 2$') #, label=name) 
            markeredgecolor='Red',zorder=1,lw=0, label=r'$\sigma_{a} < q_{0.25}$') #, label=name) 
        dia.add_sample(1.0, 1.0, marker='o',markersize=2, markerfacecolor='RoyalBlue', \
            #markeredgecolor='RoyalBlue',zorder=1,lw=0,label=r'$2 \leq \sigma_{a} < 4$') #, label=name) 
            markeredgecolor='RoyalBlue',zorder=1,lw=0,label=r'$q_{0.25} \leq \sigma_{a} < q_{0.5}$') #, label=name)  
        dia.add_sample(1.0, 1.0, marker='o',markersize=2, markerfacecolor='Violet', \
            #markeredgecolor='Violet',zorder=1,lw=0,label=r'4 $\leq \sigma_{a} < 6$') #, label=name)
            markeredgecolor='Violet',zorder=1,lw=0,label=r'q_{0.5} $\leq \sigma_{a} < q_{0.75}$') #, label=name)
        dia.add_sample(1.0, 1.0, marker='o',markersize=2, markerfacecolor='DarkGoldenRod', \
            #markeredgecolor='DarkGoldenRod',zorder=1,lw=0,label=r'$\sigma_{a} \geq 6 $') #, label=name)
            markeredgecolor='DarkGoldenRod',zorder=1,lw=0,label=r'$\sigma_{a} \geq q_{0.75}$') #, label=name)
        fig.legend(loc='center')
        
    # Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5') # 5 levels
    dia.ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Tricky: ax is the polar ax (used for plots), _ax is the
    # container (used for layout)
    dia._ax.set_title(cletter+title, fontsize=12)

fig.tight_layout()

# ---- set plot title

outfile = 'taylor_multimodel_'+season+'_'+clead+'h.png'
print ('saving to ', outfile)
PLT.savefig(outfile,dpi=400)
#PLT.show()