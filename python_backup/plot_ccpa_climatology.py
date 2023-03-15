"""
plot_ccpa_climatology.py cmonth cend_hour

computes and plots mean precip and the 95th and 98th percentiles 
of the precipitation climatology.  Easily modified to plot other
percentiles.

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
rcParams['legend.fontsize']='medium'

# =====================================================================

# ---- inputs from command line

cmonth = sys.argv[1] # '01', '02' etc.
cend_hour = sys.argv[2] # 06, 12, 18, 00 -- end hour of 6-h period
imonth = int(cmonth) - 1
nstride = 1

# ---- set parameters

pflag = False # for print statements
master_directory = '/Volumes/Backup Plus/ccpa/'
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
cmonths = ['Jan','Feb','Mar','Apr','May',\
    'Jun','Jul','Aug','Sep','Oct','Nov','Dec']

cmonths_early = ['12','01','02','03','04','05',\
    '06','07','08','09','10','11']
cmonths_late =  ['02','03','04','05','06',\
    '07','08','09','10','11','12','01']

# ---- determine the overall number of daily precipitation 
#      samples across all years for this month

iearly = int(cmonths_early[imonth])-1
ilate = int(cmonths_late[imonth])-1

if imonth != 1:  # not Feb
    nsamps_mid = ndaysomo[imonth]*18
else:
    nsamps_mid = 4*ndaysomo_leap[imonth] + 14*ndaysomo[imonth]
    
if iearly != 1:  # not Feb    
    nsamps_early = ndaysomo[iearly]*20
else:
    nsamps_early = 4*ndaysomo_leap[iearly] + 14*ndaysomo[iearly]
if ilate != 1:  # not Feb    
    nsamps_late = ndaysomo[ilate]*20
else:
    nsamps_late = 4*ndaysomo_leap[ilate] + 14*ndaysomo[ilate]
nsamps = nsamps_mid + nsamps_early + nsamps_late


# --- now either read in the previously generated mean and percentiles, or 
#     compute them and write the output to a cPickle file.  Reading them
#     back in is much faster than recomputing.

infile = 'quantiles_ccpa_'+cmonths[imonth]+'_'+cend_hour+'UTC.cPick'
fexist = False
fexist = os.path.exists(infile)
if fexist:

    inf = open(infile, 'rb')
    precip_mean = cPickle.load(inf)
    precip_q95 = cPickle.load(inf)
    precip_q98 = cPickle.load(inf)
    lons = cPickle.load(inf)
    lats = cPickle.load(inf)
    inf.close()
   
else:   
   
    # ---- read in the previously generated netCDF file with precipitation
    #      for this month and lead time as well as the surrounding
    #      two months.  All dates for this month have
    #      been smushed into one leading index, dimension nsamps,
    #      since the date of the forecast within the month and 
    #      the member number is irrelevant for the distribution 
    #      fitting.
    
    ktr = 0
    for iyear in range(2002,2020):
        print (iyear)
        for cmo in [cmonth, cmonths_early[imonth], cmonths_late[imonth]]:
            imo = int(cmo)-1
            if iyear%4 == 0:
                ndays = ndaysomo_leap[imo]
            else:
                ndays = ndaysomo[imo]
            cyear = str(iyear)    
            infile = master_directory + cyear + cmo + \
                '_ccpa_on_ndfd_grid_6hourly.nc'
            nc = Dataset(infile)
            yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
            for iday in range(1,ndays+1):
                if iday < 10:
                    cday = '0'+str(iday)
                else:
                    cday = str(iday)
                iyyyymmddhh = int(str(iyear)+cmo+cday+cend_hour)
                idx = np.where(yyyymmddhh_end == iyyyymmddhh)[0]
                precip_in = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
                if iyear == 2002 and iday == 1 and cmo == cmonth:
                    nyin, nxin = np.shape(precip_in)
                    precip_tseries = np.zeros((nsamps,nyin,nxin), \
                        dtype=np.float64)
                    missingv = -99.99*np.ones((nyin, nxin), dtype=np.float64)
                    lons = nc.variables['lons'][:,:]
                    lats = nc.variables['lats'][:,:]
                precip_in = np.where(precip_in < 500., precip_in, missingv)
                precip_tseries[ktr,:,:] = precip_in[:,:]
                ktr = ktr+1
            nc.close()
    
    print ('calculating mean')
    precip_mean = np.mean(precip_tseries, axis=0)
    print ('sorting')
    precip_tseries = np.sort(precip_tseries, axis=0)
    print ('calculating q95, q98')
    
    # -- if you want to compute other percentiles, modify here.
    
    precip_q95 = precip_tseries[int(ktr*0.95),:,:]
    precip_q98 = precip_tseries[int(ktr*0.98),:,:]

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
clevs = [0,0.05,0.1,0.2,0.3,0.5,0.7,1.0,\
    1.3,1.7,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
    urcrnrlon = 300.95782, urcrnrlat = 54.37,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons, lats)

fig = plt.figure(figsize=(8.,6.5))
axloc = [0.02,0.1,0.96,0.81]
ax1 = fig.add_axes(axloc)
title = cmonths[imonth]+\
    ' CCPA/MSWEP mean precipitation amount (mm) for 6 h period ending '+\
    cend_hour+' UTC'
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, precip_mean, clevs,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Mean 6-hourly precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'ccpa_mean_precip_6h_'+cmonths[imonth-1]+\
    '_'+cend_hour+'UTC.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- plot 95th percentile

clevs = [0,0.1,0.3,0.5,1.0,1.3,1.6,2.0,3.0,4.0,6.0,8.0,10.0,15.0,20.0]
m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
    urcrnrlon = 300.95782, urcrnrlat = 54.37,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons, lats)

fig = plt.figure(figsize=(8.,6.5))
axloc = [0.02,0.1,0.96,0.81]
ax1 = fig.add_axes(axloc)
title = cmonths[imonth]+\
    ' CCPA/MSWEP 95th percentile precipitation amount (mm)\nfor the 6 h period ending '+\
    cend_hour+' UTC'
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, precip_q95, clevs,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('95th percentile of 6-hourly precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'ccpa_q95_precip_6h_'+cmonths[imonth-1]+\
    '_'+cend_hour+'UTC.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- plot 98th percentile

clevs = [0,0.1,0.3,0.5,1.0,1.3,1.6,2.0,3.0,4.0,6.0,8.0,10.0,15.0,20.0, 25.0]
m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
    urcrnrlon = 300.95782, urcrnrlat = 54.37,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons, lats)

fig = plt.figure(figsize=(8.,6.5))
axloc = [0.02,0.1,0.96,0.81]
ax1 = fig.add_axes(axloc)
title = cmonths[imonth]+\
    ' CCPA/MSWEP 98th percentile precipitation amount (mm)\nfor the 6 h period ending '+\
    cend_hour+' UTC'
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, precip_q98, clevs,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('98th percentile of 6-hourly precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'ccpa_q98_precip_6h_'+cmonths[imonth-1]+\
    '_'+cend_hour+'UTC.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now write the output to cPickle file.

if fexist == False:
    outfile = 'quantiles_ccpa_'+cmonths[imonth]+'_'+cend_hour+'UTC.cPick'
    ouf = open(outfile, 'wb')
    cPickle.dump(precip_mean, ouf)
    cPickle.dump(precip_q95, ouf)
    cPickle.dump(precip_q98, ouf)
    cPickle.dump(lons, ouf)
    cPickle.dump(lats, ouf)
    ouf.close()

