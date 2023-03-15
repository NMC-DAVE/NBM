"""
plot_gefsv12_climatology.py cmonth clead

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
clead = sys.argv[2] # 006, 012, 018
imonth = int(cmonth) - 1
nstride = 1

# ---- set parameters

pflag = False # for print statements 
master_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
cmonths_late = ['Feb','Mar','Apr','May','Jun',\
        'Jul','Aug','Sep','Oct','Nov','Dec','Jan']    
cmonths_early = ['Dec','Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov',]

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time as well as the surrounding
#      two months.  All dates for this month have
#      been smushed into one leading index, dimension nsamps,
#      since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution 
#      fitting.
   
infile = 'quantiles_gefsv12_'+cmonths[imonth]+'_lead'+clead+'.cPick'
fexist = False
fexist = os.path.exists(infile)
if fexist:
    inf = open(infile, 'rb')
    precip_mean = cPickle.load(inf)
    precip_q95 = cPickle.load(inf)
    precip_q98 = cPickle.load(inf)
    lons_1d = cPickle.load(inf)
    lats_1d = cPickle.load(inf)
    lons = cPickle.load(inf)
    lats = cPickle.load(inf)
    inf.close()
else:

    ktr = 0
    ncfile = master_directory + cmonths[imonth] + '_apcp_sfc_h' + clead + '.nc'
    print (ncfile)
    nc = Dataset(ncfile)
    precip_middle = nc.variables['apcp_fcst'][:,:,:]
    nsamps_middle, ny_gefsv12, nx_gefsv12 = np.shape(precip_middle)
    print ('nsamps_middle, ny_gefsv12, nx_gefsv12 = ', nsamps_middle, ny_gefsv12, nx_gefsv12)
    lons_1d = nc.variables['lons_fcst'][:]
    lats_1d = nc.variables['lats_fcst'][:]
    print ('lons_1d = ', lons_1d)
    print ('lats_1d = ', lats_1d)
    nc.close()

    ncfile = master_directory + cmonths_early[imonth] + '_apcp_sfc_h' + clead + '.nc'
    print (ncfile)
    nc = Dataset(ncfile)
    precip_early = nc.variables['apcp_fcst'][:,:,:]
    nsamps_early, ny_gefsv12, nx_gefsv12 = np.shape(precip_early)
    nc.close()

    ncfile = master_directory + cmonths_late[imonth] + '_apcp_sfc_h' + clead + '.nc'
    print (ncfile)
    nc = Dataset(ncfile)
    precip_late = nc.variables['apcp_fcst'][:,:,:]
    nsamps_late, ny_gefsv12, nx_gefsv12 = np.shape(precip_late)
    nc.close()

    nsamps = nsamps_middle + nsamps_early + nsamps_late
        
    precip_mean = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.float32)
    precip_q95 = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.float32)
    precip_q98 = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.float32)
            
    for jy in range(0,ny_gefsv12):             
        print ('jy = ', jy) 
        for ix in range(0,nx_gefsv12,nstride):

            # ---- there is a grib round-off error that can give negative
            #      values slightly smaller than teeny precip.  to make sure 
            #      that we don't have either negative values or lots of the 
            #      same tiny values, so subtract teeny_precip

            precip_ens_1d_middle = precip_middle[:,jy,ix]
            precip_ens_1d_early = precip_early[:,jy,ix]
            precip_ens_1d_late = precip_late[:,jy,ix]
            precip_ens_1d = np.concatenate((precip_ens_1d_middle, \
                precip_ens_1d_early, precip_ens_1d_late))
            precip_sorted = np.sort(precip_ens_1d)
            precip_mean[jy,ix] = np.mean(precip_sorted)
            precip_q95[jy,ix] = precip_sorted[int(.95*nsamps)]   
            precip_q98[jy,ix] = precip_sorted[int(.98*nsamps)]        

lons, lats = np.meshgrid(lons_1d, lats_1d)
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
    ' GEFSv12 mean precipitation amount (mm) for '+\
    clead+' h forecast'
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

plot_title = 'GEFSv12_mean_precip_'+clead+'h_'+\
    cmonths[imonth-1]+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# ---- now 95th percentile

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
    ' GEFSv12 95th percentile precipitation amount (mm) for '+\
    clead+' h forecast'
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

plot_title = 'GEFSv12_q95_precip_'+clead+'h_'+\
    cmonths[imonth-1]+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



# ---- now 98th percentile

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
    ' GEFSv12 98th percentile precipitation amount (mm) for '+\
    clead+' h forecast'
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

plot_title = 'GEFSv12_q98_precip_'+clead+'h_'+\
    cmonths[imonth-1]+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now write the output to cPickle file.

if fexist == False:
    outfile = 'quantiles_gefsv12_'+cmonths[imonth]+'_lead'+clead+'.cPick'
    ouf = open(outfile, 'wb')
    cPickle.dump(precip_mean, ouf)
    cPickle.dump(precip_q95, ouf)
    cPickle.dump(precip_q98, ouf)
    cPickle.dump(lons_1d, ouf)
    cPickle.dump(lats_1d, ouf)
    cPickle.dump(lons, ouf)
    cPickle.dump(lats, ouf)
    ouf.close()



