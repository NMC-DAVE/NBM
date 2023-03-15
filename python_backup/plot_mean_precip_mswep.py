"""
plot_mean_precip_mswep.py cmonth clead 

"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

import _pickle as cPickle
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

# =====================================================================

def fraczero_possamps(nsamps, precip):
    """
    from the vector input sample precip_ens, define the fraction of
    samples with zero precipitation.   For the positive samples, add
    a small random number to deal with the fact that the data was 
    discretized to 0.1 mm, so that when later creating CDFs we don't 
    have values with lots of tied amounts.   Sort the nonzero amounts 
    and return.
    """
    number_zeros = 0
    precip_nonzero = np.delete(precip, \
        np.where(precip <= 0.0))  # censor at 0.1 mm
    nz = len(precip_nonzero)
    # data discretized, so add random component of this magnitude
    precip_nonzero = precip_nonzero + \
        np.random.uniform(low=-0.005,high=0.005,size=nz) 
    precip_nonzero = np.sort(precip_nonzero)  
    #print (precip_ens_nonzero[0:10]) 
    ntotal = len(precip)
    nzero = ntotal - len(precip_nonzero)
    fraction_zero = float(nzero) / float(ntotal)
    return fraction_zero, precip_nonzero, nz

# =====================================================================
    
# ---- inputs from command line

cmonth = sys.argv[1] # 01, etc
imonth = int(cmonth)-1
clead = sys.argv[2] # 03, 06, 12, etc.
ilead = int(clead)
master_directory = '/Volumes/Backup Plus/mswep/'
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ccmonth = cmonths[imonth]
process_it = False
mswep_directory = '/Volumes/Backup Plus/mswep/'


if process_it == True:
    
    # ---- read in the previously generated netCDF file with precipitation
    #      for this month and lead time.  All members, dates for this 
    #      month have been smushed into one leading index, dimension
    #      nsamps, since the date of the forecast within the month and 
    #      the member number is irrelevant for the distribution fitting.

    if imonth == 1:  # Feb
        nsamps = ndaysomo[imonth]*20
    else:
        nsamps = 6*ndaysomo_leap[imonth] + 14*ndaysomo[imonth]
    print ('nsamps = ', nsamps)
    precip = np.zeros((nsamps,1597,2345), dtype=np.float32)

    ktr = 0
    for iyear in range(2000,2020):
        if iyear%4 == 0:
            ndays = ndaysomo_leap[imonth]
        else:
            ndays = ndaysomo[imonth]
        
        if ilead == 6: # first time slice is precip from 0 to 6 UTC.
            ioff = 0
        elif ilead == 12:
            ioff = 1
        elif ilead == 18:
            ioff = 2
        else:
            ioff = 3
        
        isamplist = range(ioff,ioff+ndays*4,4)
        print (isamplist)
        cyear = str(iyear)
        infile = mswep_directory + cyear + cmonth + '_on_ndfd_grid_6hourly.nc'
        print (iyear, infile)
        nc = Dataset(infile)
        precip_in = nc.variables['apcp_anal'][isamplist,:,:]
        precip[ktr:ktr+ndays,:,:] = precip_in[:,:,:]
        if iyear == 2000:
            lons = nc.variables['lons'][:,:]
            lats = nc.variables['lats'][:,:]
            nlats, nlons = np.shape(lons)
            print ('nlats, nlons = ', nlats, nlons)
        ktr = ktr+ndays
        nc.close()

    # ---- now, a long loop to calculate fraction zero and mean

    print ('now processing to get fraction zero, mean')
    fraction_zero = np.zeros((nlats,nlons), dtype=np.float32)
    precip_mean = np.zeros((nlats,nlons), dtype=np.float32)
    for ix in range(nlons):
        if ix%10 == 0: print ('processing ix = ',ix,' of ',nlons)
        for jy in range(nlats):
            precip_1d = precip[:,jy,ix]
            fz, precip_nonzero, nz = fraczero_possamps(nsamps, precip_1d)
            fraction_zero[jy,ix] = fz
            precip_mean[jy,ix] = np.mean(precip_nonzero)
        
    # ---- save to cPickle

    outfile = 'MSWEP_'+cmonth+'_on_ndfd_grid_'+clead+'UTC.cPick'
    ouf = open(outfile, 'wb')
    cPickle.dump(fraction_zero, ouf)
    cPickle.dump(precip_mean, ouf)
    ouf.close()
    
else:

    infile = mswep_directory + '200001_on_ndfd_grid_6hourly.nc'
    nc = Dataset(infile)
    lons = nc.variables['lons'][:,:]
    lats = nc.variables['lats'][:,:]
    nlats, nlons = np.shape(lons)
    print ('nlats, nlons = ', nlats, nlons)
    nc.close()

    # ---- load from cPickle

    infile = 'MSWEP_'+cmonth+'_on_ndfd_grid_'+clead+'UTC.cPick'
    inf = open(infile, 'rb')
    fraction_zero = cPickle.load(inf)
    precip_mean = cPickle.load(inf)
    inf.close()


# ----- plot!

m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons, lats)

# ---- make plots of fraction positive precip

clevs = [0.0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

fig = plt.figure(figsize=(8,6.5))
axloc = [0.08,0.1,0.87,0.88]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'MSWEP climatological probability of non-zero precipitation, '+\
    ccmonth+', 6 h ending '+clead+' UTC'
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, 1.0-fraction_zero, clevs,\
    cmap=None, colors=colorst, extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
parallels = np.arange(-20.,80,10.)
m.drawparallels(parallels,labels=[1,0,0,0],color='LightGray')
meridians = np.arange(0.,360.,20.)
m.drawmeridians(meridians,labels=[0,0,0,1],color='LightGray')
    
# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
#cax = divider.append_axes("bottom", size="3%", pad=0.28)
cax = fig.add_axes([0.08, 0.09, 0.88, 0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('1 - fraction of samples with zero precipitation',fontsize=9)

# ---- set plot title

plot_title = 'MSWEP_fraction_zero_'+ccmonth+'_'+clead+'UTC.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')




# --- plot the ensemble mean

clevs = [0.0,0.3,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10]
fig = plt.figure(figsize=(8.,6.5))
axloc = [0.08,0.1,0.87,0.88]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'MSWEP mean of non-zero precipitation amounts, '+\
    ccmonth+' 6 h ending '+clead+' UTC'
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, precip_mean, clevs,\
    cmap=None, colors=colorst, extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
parallels = np.arange(-20.,80,10.)
m.drawparallels(parallels,labels=[1,0,0,0],color='LightGray')
meridians = np.arange(0.,360.,20.)
m.drawmeridians(meridians,labels=[0,0,0,1],color='LightGray')
    
# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
#cax = divider.append_axes("bottom", size="3%", pad=0.28)
cax = fig.add_axes([0.08, 0.09, 0.88, 0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('mean of non-zero precipitation amounts (mm)',fontsize=9)

# ---- set plot title

plot_title = 'MSWEP_mean_precip_'+ccmonth+'_'+clead+'UTC.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


