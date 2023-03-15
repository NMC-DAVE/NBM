from netCDF4 import Dataset
import numpy as np
from dateutils import daterange, datetohrs
import sys
import os
import os.path
from os import path
import numpy.ma as ma
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams 
from mpl_toolkits.basemap import Basemap
import scipy.stats as stats

rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'
rcParams['legend.fontsize']='xx-small'
rcParams['legend.fancybox']=True

# ---- loop through domains, read in the time series of soil moistures

domains = ['Great_Plains','NEq_Africa', 'NIndia']

swmean_domain1 = ma.zeros((0),dtype=np.float32)
swmean_domain2 = ma.zeros((0),dtype=np.float32)
swmean_domain3 = ma.zeros((0),dtype=np.float32)
date_list_domain1 = ma.zeros((0),dtype=np.float32)
date_list_domain2 = ma.zeros((0),dtype=np.float32)
date_list_domain3 = ma.zeros((0),dtype=np.float32)


fig = plt.figure(figsize=(9.,6.))

for idomain, domain in enumerate(domains):
    
    # ==============================================================
    # ---- first the left-hand plots of time series of soil moisture
    # ===============================================================
    
    if idomain == 0:
        llbox = [255.0, 27.0, 265.0, 40.0]
        axloc = [0.06,0.73,0.44,0.21]
        title = '(a) Southern Great Plains time series'
        lims = [0.1,0.3]
    elif idomain == 1:
        llbox = [0.0, 5.0, 30.0, 12.0]
        axloc = [0.06,0.41,0.44,0.21]
        title = '(c) N. Equatorial Africa time series'
        lims = [0.05,0.3]
    else: # idomain == 2
        llbox = [73.0, 18.0, 84.0, 29.0]
        axloc = [0.06,0.08,0.44,0.21]
        lims = [0.1,0.4]
        title = '(e) India time series'
    loncenter = (llbox[0] + llbox[2]) / 2.
    latcenter = (llbox[1] + llbox[3]) / 2.
        
    for istream, cstream in enumerate(['1999','2003','2007','2011','2015']):
        
        infile = '../gefsv12/'+domain+'_'+cstream+'_soilmoisture.dump'
        swmean_masked = np.load(infile, allow_pickle=True)
        infile = '../gefsv12/'+domain+'_'+cstream+'_datelist.dump'
        date_list_vec = np.load(infile, allow_pickle=True)
        
        if idomain == 0:
            swmean_domain1 = ma.concatenate([swmean_domain1, swmean_masked]) 
            date_list_domain1 = ma.concatenate([date_list_domain1, date_list_vec])
            hourssince1CE = np.zeros(len(swmean_domain1), dtype=np.float32)
            for i in range(len(swmean_domain1)):
                hourssince1CE[i] = datetohrs(date_list_domain1[i])
        elif idomain == 1:
            swmean_domain2 = ma.concatenate([swmean_domain2, swmean_masked]) 
            date_list_domain2 = ma.concatenate([date_list_domain2, date_list_vec]) 
        else:
            swmean_domain3 = ma.concatenate([swmean_domain3, swmean_masked]) 
            date_list_domain3 = ma.concatenate([date_list_domain3, date_list_vec])
            

    if idomain == 0:
        swmean = swmean_domain1
    elif idomain == 1:
        swmean = swmean_domain2
    else: # idomain == 2
        swmean = swmean_domain3

    # ---- plot the time series of mean soil moisture

    ndays = 365*4+1
    a1 = fig.add_axes(axloc)
    a1.set_title(title,fontsize=10)
    a1.plot(2000.+np.arange(ndays)/365.,swmean[0:ndays], color='Black',linewidth=1.,label='1999 stream')
    a1.plot(2004.+np.arange(ndays)/365.,swmean[ndays:2*ndays], color='RoyalBlue',linewidth=1.,label='2003 stream')
    a1.plot(2008.+np.arange(ndays)/365.,swmean[2*ndays:3*ndays], color='Red',linewidth=1.,label='2007 stream')
    a1.plot(2012.+np.arange(ndays)/365.,swmean[3*ndays:4*ndays], color='LimeGreen',linewidth=1.,label='2011 stream')
    ndays2 = len(date_list_domain1) - 4*ndays
    a1.plot(2016.+np.arange(ndays2)/365.,swmean[4*ndays:4*ndays+ndays2], color='Violet',linewidth=1.,label='2015 stream')
    a1.set_ylim(lims[0],lims[1])
    a1.set_ylabel('Soil Moisture', fontsize=7)
    a1.set_xlim(2000,2020)
    a1.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
    a1.grid(color='LightGray')
    if idomain==0: a1.legend(loc=0)
    
    # --- next to this, add basemap
    
    if idomain == 0:
        axloc = [0.51,0.73,0.2,0.21]
    elif idomain == 1:
        axloc = [0.51,0.41,0.2,0.21]
    else: # idomain == 2
        axloc = [0.51,0.08,0.2,0.21]
    
    a1 = fig.add_axes(axloc)
    print (loncenter, latcenter)
    m = Basemap(llcrnrlon=loncenter-30.,llcrnrlat=latcenter - 25.,\
        urcrnrlon=loncenter+30.,urcrnrlat=latcenter+25.,\
        projection='mill',resolution='l')
    m.drawcoastlines(linewidth=0.5,color='Gray')
    m.drawcountries(linewidth=0.3,color='Gray')
    m.drawstates(linewidth=0.2,color='Gray')
    x,y = m([llbox[0],llbox[0],llbox[2],llbox[2],llbox[0]],\
        [llbox[1],llbox[3],llbox[3],llbox[1],llbox[1]])
    m.plot(x,y,color='Red',lw=1)
    
    
    llbox = [255.0, 25.0, 265.0, 40.0]
    
    # ==============================================================
    # ---- second the right-hand scatterplots for begin/end of streams
    # ===============================================================   
    
    if idomain == 0:
        axloc = [0.78,0.73,0.15,0.21]
        title = '(b) Southern Great Plains scatter'
        infile = '../gefsv12/Great_Plains_streamboundary_soilq.dump'
        lims = [0.1,0.3]
    elif idomain == 1:
        axloc = [0.78,0.41,0.15,0.21]
        title = '(d) N. Equatorial Africa scatter'
        infile = '../gefsv12/NEq_Africa_streamboundary_soilq.dump'
        lims = [0.05,0.3]
    else: # idomain == 2
        axloc = [0.78,0.08,0.15,0.21]
        title = '(f) India scatter'
        infile = '../gefsv12/NIndia_streamboundary_soilq.dump'
        lims = [0.1,0.4]
    
    inf = open(infile,'rb')
    sw1_save = cPickle.load(inf)
    sw2_save = cPickle.load(inf)
    sw3_save = cPickle.load(inf)
    inf.close()   
    
    a1 = fig.add_axes(axloc)
    a1.plot([0.0,0.5], [0.0,0.5], color='LightGray',lw=1,zorder=10)
    a1.set_title(title,fontsize=8)
    a1.scatter(sw1_save,sw2_save, c='Red',s=0.1,zorder=11)
    slope, intercept, r_value, p_value, std_err = \
        stats.linregress(sw1_save.flatten(),sw2_save.flatten())
    r = stats.pearsonr(sw1_save.flatten(),sw2_save.flatten())
    print (r[0])
    cr = 'r = '+'{0:7.4f}'.format(r[0])
    a1.text(0.03,0.43,cr,fontsize=8)
    a1.set_ylim(0.0,0.5)
    a1.set_xlim(0.0,0.5)
    a1.plot([0.0,0.5],[intercept,intercept+slope*0.5],'r-',lw=0.5,zorder=19)
    a1.set_xticks([0.0,0.1,0.2,0.3,0.4,0.5])
    a1.set_ylabel('End of first stream\nsoil moisture', fontsize=7)
    a1.set_xlabel('Beginning of second stream\nsoil moisture', fontsize=6)
    a1.grid(color='LightGray',lw=0.5)

    
plot_title = 'soilwater.png'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title,dpi=600)
print ('Plot done')


        
        
