""" python compare_skint_parallel.py """

from netCDF4 import Dataset
import numpy as np
from dateutils import daterange, datetohrs, dayofyear
import sys
import pygrib
import os
import os.path
from os import path
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PIL import Image

rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'

date_list = daterange('2017120100', '2019113000', 24)
date_list_winter = daterange('2017120100', '2018043000', 24) + \
     daterange('2018100100', '2019043000', 24) + \
     daterange('2019100100', '2019113000', 24)
date_list_summer = daterange('2018050100', '2018093000', 24) + \
     daterange('2019050100', '2019093000', 24) 
ndates = len(date_list)
ndates_summer = len(date_list_summer)
ndates_winter = len(date_list_winter)
print ('ndates, ndates_summer, ndates_winter = ',ndates, ndates_summer, ndates_winter )
#sys.exit()


# --- make a land-sea mask

infile = 'hgt_sfc.grib2'
fexist = os.path.exists(infile)
if fexist == True:
    grbfile = pygrib.open(infile)
    grb = grbfile.select()[0] # forecastTime
    height = grb.values
    lats, lons = grb.latlons()
    nlats, nlons = np.shape(lats)
    print (nlats, nlons)
    print ('height[nlats//4,0:nlons:5] = ',height[nlats//4,0:nlons:5])
    print ('lats[nlats//4,0:nlons:5] = ',lats[nlats//4,0:nlons:5])
    print ('lons[nlats//4,0:nlons:5] = ',lons[nlats//4,0:nlons:5])
    grbfile.close()
    

land_sea_mask_0p25 = np.ones((nlats,nlons), dtype=np.int32)
zeros = np.zeros((nlats,nlons), dtype=np.int32)
land_sea_mask_0p25 = np.where(np.logical_and(height  > 0.0, height  < 0.03), zeros, land_sea_mask_0p25)
print ('land_sea_mask_0p25[180,0::12] = ',land_sea_mask_0p25[180,0::12])

im = Image.fromarray(land_sea_mask_0p25)
imu = im.resize(size=(181,360),resample=Image.BOX)
print ('np.shape(imu) = ', np.shape(imu))
land_sea_mask_1p0 = np.transpose(np.asarray(imu))
ny1p0, nx1p0 = np.shape(land_sea_mask_1p0)
print ('land_sea_mask_1p0[ny1p0//4, 0:nx1p0:4] = ', land_sea_mask_1p0[ny1p0//4, 0:nx1p0:4])


im = Image.fromarray(lats)
imu = im.resize(size=(360,181),resample=Image.BOX)
lats_1p0 = np.transpose(np.asarray(imu))
print ('lats_1p0[ny1p0//4, 0:nx1p0:4] = ', lats_1p0[ny1p0//4, 0:nx1p0:4])

im = Image.fromarray(lons)
imu = im.resize(size=(360,181),resample=Image.BOX)
lons_1p0 = np.transpose(np.asarray(imu))
print ('np.shape(land_sea_mask_1p0) = ',np.shape(land_sea_mask_1p0))
print ('lons_1p0[ny1p0//4, 0:nx1p0:4] = ', lons_1p0[ny1p0//4, 0:nx1p0:4])

sys.exit()

# ---- loop over dates.

ktrsum = 0
ktrwin = 0
for idate, date in enumerate(date_list):
        
    yy = date[0:4]
    yyyymmdd = date[0:8]
    mm = date[4:6]
    dd = date[6:8]
    
    infile1 = '/Volumes/Backup Plus/gefsv12/skint/'+date+'_skint.grib2'
    infile2 = '/Volumes/Backup Plus/bfg/2015/bfg_'+yyyymmdd+'06_fhr00_control2.nc4'
    print (infile1, infile2)
    fexist1 = path.exists(infile1)
    fexist2 = path.exists(infile2)
    if fexist1 == True and fexist2 == True:
        
        # ---- read in the skin temperatures from the pre-operational parallel
        
        sstfile = pygrib.open(infile1)
        grb = sstfile.select(shortName = 't')[0]
        sst = grb.values
        if idate == 0:
            lats, lons = grb.latlons()
            nlats, nlons = np.shape(lats)
            nsst_save = ma.zeros((ndates,nlats,nlons), dtype=np.float32)
            reanal_save = ma.zeros((ndates,256,512), dtype=np.float32)
            nsst_save_summer = ma.zeros((ndates_summer,nlats,nlons), dtype=np.float32)
            reanal_save_summer = ma.zeros((ndates_summer,256,512), dtype=np.float32)
            nsst_save_winter = ma.zeros((ndates_winter,nlats,nlons), dtype=np.float32)
            reanal_save_winter = ma.zeros((ndates_winter,256,512), dtype=np.float32)
        sstfile.close()
        nsst_save[idate,:,:] = sst[:,:]
        
        # --- read the associated reanalysis skin temperature
                   
        nc = Dataset(infile2)
        tmpsfc_in = nc.variables['tmpsfc'][0,:,:]
        #orography = nc.variables['orogsfc'][0,:,:]
        #land_sea_mask = np.ones((nlats,nlons), dtype=np.int32)
        #zeros = np.zeros((nlats,nlons), dtype=np.int32)
        #land_sea_mask = np.where(np.logical_and(orography  > 0.0, orography  < 0.03), zeros, land_sea_mask)
        
        #tmpsfc_in = nc.variables['tmp2m'][0,:,:]
        if idate == 0:
            lon_reanal = nc.variables['lon'][:]
            lat_reanal = nc.variables['lat'][:]
            print ('lat_reanal = ',lat_reanal)
            landsfc = nc.variables['landsfc'][0,:,:] # land/water mask
        nc.close()
        reanal_save[idate,:,:] = tmpsfc_in[:,:]
        if date_list_summer.count(date) > 0:
             nsst_save_summer[ktrsum,:,:] = sst[:,:]
             reanal_save_summer[ktrsum,:,:] = tmpsfc_in[:,:]
             ktrsum = ktrsum + 1
        if date_list_winter.count(date) > 0:
             nsst_save_winter[ktrwin,:,:] = sst[:,:]
             reanal_save_winter[ktrwin,:,:] = tmpsfc_in[:,:]
             ktrwin = ktrwin + 1     
        

    else:
        if fexist1 == False:
            print ('unable to read ', infile1)
        if fexist2 == False:
            print ('unable to read ', infile2)
        nsst_save[idate,:,:] = ma.masked
        reanal_save[idate,:,:] = ma.masked
        if date_list_summer.count(date) > 0:
             nsst_save_summer[ktrsum,:,:] = ma.masked
             reanal_save_summer[ktrsum,:,:] = ma.masked
             ktrsum = ktrsum + 1
        if date_list_winter.count(date) > 0:
             nsst_save_winter[ktrwin,:,:] = ma.masked
             reanal_save_winter[ktrwin,:,:] = ma.masked
             ktrwin = ktrwin + 1
    
        
print ('nsst_save_summer[:,60,60]', nsst_save_summer[:,60,60])
print ('nsst_save_winter[:,60,60]', nsst_save_winter[:,60,60]) 
print ('reanal_save_summer[:,60,60]', reanal_save_summer[:,60,60])
print ('reanal_save_winter[:,60,60]', reanal_save_winter[:,60,60])        
        
nsst_mean = ma.mean(nsst_save,axis=0)
nsst_mean_summer = ma.mean(nsst_save_summer,axis=0)
nsst_mean_winter = ma.mean(nsst_save_winter,axis=0)

# --- flip nsst and reanalysis upside down so lats oriented S to N, increasing.
        
nsst_mean = np.flipud(nsst_mean)  # need latitudes ascending order 
nsst_mean_summer = np.flipud(nsst_mean_summer)  
nsst_mean_winter = np.flipud(nsst_mean_winter)  
lons = np.flipud(lons)
lats = np.flipud(lats)

reanal_mean = ma.mean(reanal_save,axis=0)
reanal_mean_summer = ma.mean(reanal_save_summer,axis=0)
reanal_mean_winter = ma.mean(reanal_save_winter,axis=0)

reanal_mean = np.flipud(reanal_mean)
reanal_mean_summer = np.flipud(reanal_mean_summer)
reanal_mean_winter = np.flipud(reanal_mean_winter)
lat_reanal = np.flipud(lat_reanal)

# ---- interpolate the reanalysis data to the NSST grid.

reanal_mean_NSSTgrid = interp(reanal_mean, lon_reanal, lat_reanal, \
    lons, lats, checkbounds=False, masked=False, order=1) # - 273.15    
reanal_mean_NSSTgrid_summer = interp(reanal_mean_summer, lon_reanal, lat_reanal, \
    lons, lats, checkbounds=False, masked=False, order=1) # - 273.15
reanal_mean_NSSTgrid_winter = interp(reanal_mean_winter, lon_reanal, lat_reanal, \
    lons, lats, checkbounds=False, masked=False, order=1) # - 273.15
print ('max, min reanal_mean_NSSTgrid = ', ma.max(reanal_mean_NSSTgrid), \
    ma.min(reanal_mean_NSSTgrid))
    
# ---- plot difference

sst_difference = nsst_mean - reanal_mean_NSSTgrid
sst_difference_summer = nsst_mean_summer - reanal_mean_NSSTgrid_summer
sst_difference_winter = nsst_mean_winter - reanal_mean_NSSTgrid_winter

# ---- code to use for plotting

fig = plt.figure(figsize=(9.,5.6))
axloc = [0.07,0.075,0.9,0.82]
ax = fig.add_axes(axloc)
ax.set_title('00 UTC global skin temperature differences,\npre-production parallel minus reanalysis, 1 Dec 2017 to 30 Nov 2019',fontsize=16)
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
#clevs = [-2.0,-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0]
clevs = [-6,-4,-2,-1,-0.5,0.5,1,2,4,6]
colorstblack='Black'
parallels = np.arange(-80.,90,10.)
meridians = np.arange(0.,360.,20.)
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=-75,\
    urcrnrlon=lons[-1,-1],urcrnrlat=75.,\
    projection='mill',resolution='l')
x, y = m(lons, lats)

print ('np.shape(sst_difference) = ', np.shape(sst_difference))
print ('np.shape(land_sea_mask_1p0) = ',np.shape(land_sea_mask_1p0))
print ('land_sea_mask_1p0[45,0::3] = ',land_sea_mask_1p0[45,0::3])
CS2 = m.contourf(x,y,sst_difference,clevs,cmap=None,colors=colorst,extend='both')

m.drawcoastlines(linewidth=0.5,color='Gray')
m.drawcountries(linewidth=0.3,color='Gray')
m.drawparallels(parallels,labels=[1,0,0],linewidth=0.15,fontsize=8)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.15,fontsize=8)

# ---- use axes_grid toolkit to make colorbar axes.

#ax = fig.add_axes([0.,0.,1.,1.])
divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Difference (deg C)')

# ---- set plot title

plot_title = 'skint_difference.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# ---- code to use for plotting

fig = plt.figure(figsize=(9.,9.))
axloc = [0.07,0.435,0.9,0.52]
ax = fig.add_axes(axloc)
ax.set_title('(a) 00 UTC global skin temperature differences,\npre-production parallel minus reanalysis, 1 Dec 2017 to 30 Nov 2019',fontsize=16)
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
#clevs = [-2.0,-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0]
clevs = [-6,-4,-2,-1,-0.5,0.5,1,2,4,6]
colorstblack='Black'
parallels = np.arange(-80.,90,10.)
meridians = np.arange(0.,360.,20.)
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=-75,\
    urcrnrlon=lons[-1,-1],urcrnrlat=75.,\
    projection='mill',resolution='l')
x, y = m(lons, lats)

CS2 = m.contourf(x,y,sst_difference*(1.-land_sea_mask_1p0),clevs,cmap=None,colors=colorst,extend='both')

m.drawcoastlines(linewidth=0.7,color='Gray')
m.drawcountries(linewidth=0.5,color='Gray')
m.drawparallels(parallels,labels=[1,0,0],linewidth=0.15,fontsize=8)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.15,fontsize=8)



axloc = [0.07,0.095,0.42,0.32]
ax = fig.add_axes(axloc)
ax.set_title('(b) Nov-Apr differences',fontsize=16)
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
#clevs = [-2.0,-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0]
clevs = [-6,-4,-2,-1,-0.5,0.5,1,2,4,6]
colorstblack='Black'
parallels = np.arange(20.,60.,10.)
meridians = np.arange(60.,140.,20.)
m = Basemap(llcrnrlon=60.,llcrnrlat=20,\
    urcrnrlon=140.,urcrnrlat=60.,\
    projection='mill',resolution='l')
x, y = m(lons, lats)
print ('max, min sst_difference_winter ', np.max(sst_difference_winter), np.min(sst_difference_winter))
CS2 = m.contourf(x,y,sst_difference_winter*land_sea_mask_1p0,clevs,cmap=None,colors=colorst,extend='both')

m.drawcoastlines(linewidth=0.7,color='Gray')
m.drawcountries(linewidth=0.7,color='Gray')
m.drawparallels(parallels,labels=[1,0,0],linewidth=0.15,fontsize=8)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.15,fontsize=8)


axloc = [0.55,0.095,0.42,0.32]
ax = fig.add_axes(axloc)
ax.set_title('(c) May-Oct differences',fontsize=16)
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
#clevs = [-2.0,-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0]
clevs = [-6,-4,-2,-1,-0.5,0.5,1,2,4,6]
colorstblack='Black'
parallels = np.arange(20.,60.,10.)
meridians = np.arange(60.,140.,20.)
m = Basemap(llcrnrlon=60.,llcrnrlat=20,\
    urcrnrlon=140.,urcrnrlat=60.,\
    projection='mill',resolution='l')
x, y = m(lons, lats)
print ('max, min sst_difference_summer ', np.max(sst_difference_summer), np.min(sst_difference_summer))
CS2 = m.contourf(x,y,sst_difference_summer*land_sea_mask_1p0,clevs,cmap=None,colors=colorst,extend='both')

m.drawcoastlines(linewidth=0.7,color='Gray')
m.drawcountries(linewidth=0.7,color='Gray')
m.drawparallels(parallels,labels=[1,0,0],linewidth=0.15,fontsize=8)
m.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.15,fontsize=8)


# ---- use axes_grid toolkit to make colorbar axes.

axloc = [0.02, 0.06, 0.96, 0.02]
cax = fig.add_axes(axloc)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Difference (deg C)')

# ---- set plot title

plot_title = 'skint_difference_3panel.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



