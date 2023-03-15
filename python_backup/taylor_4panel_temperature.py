#!/usr/bin/env python
# taylor_4panel_temperature.py 

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
rcParams['legend.fontsize']='small'

cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'

model = sys.argv[1] # raw, dav, qmap, mos, mvmos
clead = sys.argv[2]

# Reference std
stdrefs = {"Jan-Feb-Mar":1.0,   "Apr-May-Jun":1.0,\
    "Jul-Aug-Sep":1.0, "Oct-Nov-Dec":1.0}
               
if model == 'raw':
    ftype = 'raw_forecast_errorstats_2019_'
    title = 'Raw\nmodel output\n+'+clead+' h'
elif model ==  'dav':
    ftype = 'decayavg_forecast_errorstats_2019_'  
    title = 'Decaying-average\nbias correction\n+'+clead+' h'
elif model == 'qmap':
    ftype = 'quantile_mapping_errorstats_2019_'  
    title = 'Quantile\nmapping\n+'+clead+' h'
elif model == 'mos':
    ftype = 'MOS_GEFSv12_2019_'
    title = 'Univariate\nMOS\n+'+clead+' h'
elif model == 'mvmos':
    ftype = 'MOS_mvr_GEFSv12_2019_'
    title = 'Multi-variate\nMOS\n+'+clead+' h'
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

print ('max, min std_normalized_fall  = ', NP.max(std_normalized_fall ), NP.min(std_normalized_fall ))
inf.close()
nlats, nlons = NP.shape(corr_fall)


rects = dict(JFM=221, AMJ=222, JAS=223, OND=224)

fig = PLT.figure(figsize=(10,8))
fig.suptitle(title, size='large')

for season in ['Jan-Feb-Mar','Apr-May-Jun','Jul-Aug-Sep','Oct-Nov-Dec']:

    if season == 'Jan-Feb-Mar':
        corrcoef = corr_winter
        stddev = std_normalized_winter
        stdobs = std_obs_winter
        cletter = '(a) '
        rect = 221
    elif season == 'Apr-May-Jun':
        corrcoef = corr_spring
        stddev = std_normalized_spring
        stdobs = std_obs_spring
        cletter = '(b) '
        rect = 222
    elif season == 'Jul-Aug-Sep':
        corrcoef = corr_summer
        stddev = std_normalized_summer
        stdobs = std_obs_summer
        cletter = '(c) '
        rect = 223
    elif season == 'Oct-Nov-Dec':
        corrcoef = corr_fall
        stddev = std_normalized_fall 
        stdobs = std_obs_fall  
        cletter = '(d) '  
        rect = 224

    meancoef = NP.mean(corrcoef)
    meanstddev = NP.mean(stddev)
    dia = TaylorDiagram(1.0, fig=fig, rect=rect, \
        label=season, srange=(0.0, 2.0))
        
    print ('max, min corrcoef = ', NP.max(corrcoef), NP.min(corrcoef))
    print ('max, min stddev = ', NP.max(stddev), NP.min(stddev))
    
    # Add samples to Taylor diagram
    
    print (NP.shape(stddev))
    print (NP.shape(corrcoef))
    for jlat in range(nlats):
        for ilon in range(nlons):
            if stdobs[jlat,ilon] < 2.0:
                dia.add_sample(stddev[jlat,ilon], corrcoef[jlat,ilon],marker='o',fillstyle='full',\
                    markersize=0.3, markerfacecolor='Red', markeredgecolor='Red',zorder=4) #, label=name)
            elif stdobs[jlat,ilon] > 2.0 and stdobs[jlat,ilon] <= 4.0:
                dia.add_sample(stddev[jlat,ilon], corrcoef[jlat,ilon],marker='o',fillstyle='full',\
                    markersize=0.3, markerfacecolor='RoyalBlue', markeredgecolor='RoyalBlue',zorder=5) #, label=name)
            elif stdobs[jlat,ilon] > 4.0 and stdobs[jlat,ilon] <= 6.0:
                dia.add_sample(stddev[jlat,ilon], corrcoef[jlat,ilon],marker='o',fillstyle='full',\
                    markersize=0.3, markerfacecolor='Violet', markeredgecolor='Violet',zorder=6) #, label=name)
            else:
                dia.add_sample(stddev[jlat,ilon], corrcoef[jlat,ilon],marker='o',fillstyle='full',\
                    markersize=0.3, markerfacecolor='DarkGoldenRod', \
                    zorder=7) #, label=name)
    
    dia.add_sample(meanstddev, meancoef,marker='o',\
        markersize=3.5, markerfacecolor='Blue',markeredgecolor='Blue',zorder=20) #, label=name)
        
    dia.add_sample([0.0,2.0], [0.2,0.2], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.4,0.4], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.6,0.6], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.7,0.7], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.8,0.8], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.9,0.9], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.95,0.95], 'k-',lw=0.25)
    dia.add_sample([0.0,2.0], [0.99,0.99], 'k-',lw=0.25)
    
      
    if season == 'Oct-Nov-Dec':  
        dia.add_sample(1.0, 1.0, marker='o',markersize=1.5, markerfacecolor='Red', \
            markeredgecolor='Red',zorder=1,lw=0, label=r'$\sigma_{obs} < 2$') #, label=name) 
        dia.add_sample(1.0, 1.0, marker='o',markersize=1.5, markerfacecolor='RoyalBlue', \
            markeredgecolor='RoyalBlue',zorder=1,lw=0,label=r'$2 \leq \sigma_{obs} < 4$') #, label=name)  
        dia.add_sample(1.0, 1.0, marker='o',markersize=1.5, markerfacecolor='Violet', \
            markeredgecolor='Violet',zorder=1,lw=0,label=r'4 $\leq \sigma_{obs} < 6$') #, label=name)
        dia.add_sample(1.0, 1.0, marker='o',markersize=1.5, markerfacecolor='DarkGoldenRod', \
            markeredgecolor='DarkGoldenRod',zorder=1,lw=0,label=r'$\sigma_{obs} \geq 6 $') #, label=name)
        fig.legend(loc='center')
        
    # Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5') # 5 levels
    dia.ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Tricky: ax is the polar ax (used for plots), _ax is the
    # container (used for layout)
    dia._ax.set_title(cletter+season, fontsize=12)

    


fig.tight_layout()

# ---- set plot title

outfile = 'taylor_4panel_'+model+'_'+clead+'h.png'
print ('saving to ', outfile)
PLT.savefig(outfile,dpi=300)
#PLT.show()