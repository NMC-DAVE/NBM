"""
plot_reforecast_precipitation_grib.py

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
from netCDF4 import Dataset

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

infile = '/Volumes/Backup Plus/gefsv12/precip/apcp_sfc_2019011000_p01.grib2'
print (infile)
flatlon = pygrib.open(infile)
fcst = flatlon.select(shortName='tp', dataDate=20190110, validityTime=600)[0]

flatlon.close()

# ---- define bounding coordinates on Gaussian grid for MDL domain
   
nib1 = 518 # ~lon 220E
nie1 = 1440 # up to ~lon 310E
nib2 = 0 # ~lon 220E
nie2 = 45 # up to ~lon 310E
njb = 38 # lat ~ 80.5
nje = 483 # down to lat ~ -30.5
nj = nje - njb 
ni1 = nie1 - nib1
ni2 = nie2 - nib2
ni = ni1+ni2 
print ('nj, ni = ', nj, ni)



infile = '/Volumes/Backup Plus/gefsv12/precip/apcp_sfc_2000013000_p02.grib2'
validityDate1 = 20000130
validityDate2 = 20000130
validityTime1 = 300
validityTime2 = 600

fcstfile = pygrib.open(infile)
grb_late = fcstfile.select(shortName='tp',\
    validityDate=validityDate2, \
    validityTime=validityTime2)[0]
grb_early = fcstfile.select(shortName='tp',\
    validityDate=validityDate1, \
    validityTime=validityTime1)[0]
lats, lons = grb_late.latlons()
nlats, nlons = lons.shape
lats_1d = np.zeros((nlats),dtype=np.float32)
lons_1d = np.zeros((nlons),dtype=np.float32)
lats_1d[:] = lats[:,0]
lons_1d[:] = lons[0,:]
lons_1d = np.where(lons_1d > 90.0, \
    lons_1d - 360., lons_1d)  
fcstfile.close()
lons_1d_out = \
    np.hstack((lons_1d[nib1:nie1], \
    lons_1d[nib2:nie2]))
lats_1d_out = lats_1d[njb:nje]
                
grb_global = grb_late.values - grb_early.values

work_array = np.zeros((nj,ni),dtype=np.float)
work_array1 = np.zeros((nj,ni),dtype=np.float)
work_array2 = np.zeros((nj,ni),dtype=np.float)
    
work_array2[:,:] = np.hstack((grb_late.values\
    [njb:nje,nib1:nie1], grb_late.values[njb:nje,nib2:nie2]))
work_array1[:,:] = np.hstack((grb_early.values\
    [njb:nje,nib1:nie1], grb_early.values[njb:nje,nib2:nie2]))
work_array[:,:] = work_array2[:,:] - work_array1[:,:]
print ('minimum of 6 - 3 = ', np.min(work_array))


print ('work_array2[194,582] = ', work_array2[194,582])
print ('work_array1[194,582] = ', work_array1[194,582])
print ('work_array[194,582]  = ', work_array[194,582])
#sys.exit()

lats_1d_out = np.flipud(lats_1d_out)
work_array = np.flipud(work_array)
work_array1 = np.flipud(work_array1)
work_array2 = np.flipud(work_array2)
lons_2d, lats_2d = np.meshgrid(lons_1d_out,lats_1d_out)

# --- now plot the precipitation

clevs = [0.0,0.1,0.5,0.7,1.0,2.0,3.0,4.0,5.0,\
    6.0,8.0,10.0,12.0,15.0,20.0,25.0,30.0,50.0,100.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','Gray','Black']

fig = plt.figure(figsize=(9.,6.5))
axloc = [0.02,0.11,0.96,0.8]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 precipitation hour +6-h values, apcp_sfc_2000013000_p02.grib2'
ax1.set_title(title, fontsize=14,color='Black')
m = Basemap(llcrnrlon=lons_2d[0,0],llcrnrlat=lats_2d[0,0],\
    urcrnrlon=lons_2d[-1,-1],urcrnrlat=lats_2d[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons_2d, lats_2d)
CS2 = m.contourf(x,y,work_array2, clevs,alpha=0.6, \
    cmap=None, colors=colorst, extend='both') #
m.drawcoastlines(linewidth=0.5,color='Gray',zorder=32)
m.drawcountries(linewidth=0.5,color='Gray',zorder=32)
m.drawstates(linewidth=0.5,color='Gray',zorder=32)

# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.02,0.08,0.96,0.03])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,format='%g') # ticks=clevs,
cb.ax.tick_params(labelsize=7)
cb.ax.set_xticks(clevs)
#cb.ax.set_xticklabels(['0.0','0.1','0.5','0.7','1.0','2.0','3.0','4.0','5.0',\
#    '6.0','8.0','10.0','12.0','15.0','20.0','25.0','30.0','50.0','100.0'])  # horizontal colorbar
cb.set_label('Ensemble precipitation amount (mm)')

# ---- set plot title

plot_title = 'precip_ens_mean_test_hour6.png'
fig.savefig(plot_title, dpi=300,fontsize=9)

print ('saving plot to file = ',plot_title)
print ('Done!')




fig = plt.figure(figsize=(9.,6.5))
axloc = [0.02,0.11,0.96,0.8]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 precipitation hour +3-h values, apcp_sfc_2000013000_p02.grib2'
ax1.set_title(title, fontsize=14,color='Black')
m = Basemap(llcrnrlon=lons_2d[0,0],llcrnrlat=lats_2d[0,0],\
    urcrnrlon=lons_2d[-1,-1],urcrnrlat=lats_2d[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons_2d, lats_2d)
CS2 = m.contourf(x,y,work_array1, clevs,alpha=0.6, \
    cmap=None, colors=colorst, extend='both') #
m.drawcoastlines(linewidth=0.5,color='Gray',zorder=32)
m.drawcountries(linewidth=0.5,color='Gray',zorder=32)
m.drawstates(linewidth=0.5,color='Gray',zorder=32)

# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.02,0.08,0.96,0.03])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,format='%g') # ticks=clevs,
cb.ax.tick_params(labelsize=7)
cb.ax.set_xticks(clevs)
#cb.ax.set_xticklabels(['0.0','0.1','0.5','0.7','1.0','2.0','3.0','4.0','5.0',\
#    '6.0','8.0','10.0','12.0','15.0','20.0','25.0','30.0','50.0','100.0'])  # horizontal colorbar
cb.set_label('Ensemble precipitation amount (mm)')

# ---- set plot title

plot_title = 'precip_ens_mean_test_hour3.png'
fig.savefig(plot_title, dpi=300,fontsize=9)

print ('saving plot to file = ',plot_title)
print ('Done!')



clevs = [-0.1,0.0,0.1,0.5,0.7,1.0,2.0,3.0,4.0,5.0,\
    6.0,8.0,10.0,12.0,15.0,20.0,25.0,30.0,50.0,100.0]
colorst = ['Black','White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','Gray','Black']

fig = plt.figure(figsize=(9.,6.5))
axloc = [0.02,0.11,0.96,0.8]
ax1 = fig.add_axes(axloc)
title = r'GEFSv12 precipitation hour +6h - +3h values, apcp_sfc_2000013000_p02.grib2'
ax1.set_title(title, fontsize=14,color='Black')
m = Basemap(llcrnrlon=lons_2d[0,0],llcrnrlat=lats_2d[0,0],\
    urcrnrlon=lons_2d[-1,-1],urcrnrlat=lats_2d[-1,-1],\
    resolution='l', projection='mill')
x, y = m(lons_2d, lats_2d)
CS2 = m.contourf(x,y,work_array1, clevs,alpha=0.6, \
    cmap=None, colors=colorst, extend='both') #
m.drawcoastlines(linewidth=0.5,color='Gray',zorder=32)
m.drawcountries(linewidth=0.5,color='Gray',zorder=32)
m.drawstates(linewidth=0.5,color='Gray',zorder=32)

# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.02,0.08,0.96,0.03])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,format='%g') # ticks=clevs,
cb.ax.tick_params(labelsize=7)
cb.ax.set_xticks(clevs)
#cb.ax.set_xticklabels(['0.0','0.1','0.5','0.7','1.0','2.0','3.0','4.0','5.0',\
#    '6.0','8.0','10.0','12.0','15.0','20.0','25.0','30.0','50.0','100.0'])  # horizontal colorbar
cb.set_label('Ensemble precipitation amount (mm)')

# ---- set plot title

plot_title = 'precip_ens_mean_test_hour6-3.png'
fig.savefig(plot_title, dpi=300,fontsize=9)

print ('saving plot to file = ',plot_title)
print ('Done!')
