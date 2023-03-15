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

rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'

date_list = daterange('2019090100', '2019123100', 24)
ndates = len(date_list)

# ---- loop over dates.

for idate, date in enumerate(date_list):
        
    yy = date[0:4]
    yyyymmdd = date[0:8]
    mm = date[4:6]
    dd = date[6:8]
    infile = 'nsst/'+date[0:8]+'_gdas.t00z.pgrb2.0p25.f000'
    print (infile)
    fexist = path.exists(infile)
    if fexist == True:
        
        # ---- read in the skin temperatures from the operational
        
        sstfile = pygrib.open(infile)
        grb = sstfile.select(shortName = 't')[0]
        sst = grb.values
        if idate == 0:
            lats, lons = grb.latlons()
            nlats, nlons = np.shape(lats)
            nsst_save = ma.zeros((ndates,nlats,nlons), dtype=np.float32)
            reanal_save = ma.zeros((ndates,256,512), dtype=np.float32)
        sstfile.close()
        nsst_save[idate,:,:] = sst[:,:]
        
        # --- read the associated reanalysis skin temperature
            
        infile ='/Users/Tom/python/gefsv12/2015/bfg_'+\
            yyyymmdd+'00_fhr00_control2.nc4'  
        print (infile)          
        nc = Dataset(infile)
        tmpsfc_in = nc.variables['tmpsfc'][0,:,:]
        #tmpsfc_in = nc.variables['tmp2m'][0,:,:]
        if idate == 0:
            lon_reanal = nc.variables['lon'][:]
            lat_reanal = nc.variables['lat'][:]
            print ('lat_reanal = ',lat_reanal)
            landsfc = nc.variables['landsfc'][0,:,:]
        nc.close()
        reanal_save[idate,:,:] = tmpsfc_in[:,:]

    else:
        print ('unable to read ', infile)
        nsst_save[idate,:,:] = ma.masked
        reanal_save[idate,:,:] = ma.masked
    
        
nsst_mean = ma.mean(nsst_save,axis=0)

# --- flip nsst and reanalysis upside down so lats oriented S to N, increasing.
        
nsst_mean = np.flipud(nsst_mean)  # need latitudes ascending order 
lons = np.flipud(lons)
lats = np.flipud(lats)

reanal_mean = ma.mean(reanal_save,axis=0)
reanal_mean = np.flipud(reanal_mean)
lat_reanal = np.flipud(lat_reanal)

print ('reanal_mean max, min = ', ma.max(reanal_mean), ma.min(reanal_mean))
print ('nsst_mean max, min = ', ma.max(nsst_mean), ma.min(nsst_mean))

# ---- interpolate the reanalysis data to the NSST grid.

reanal_mean_NSSTgrid = interp(reanal_mean, lon_reanal, lat_reanal, \
    lons, lats, checkbounds=False, masked=False, order=1) # - 273.15
print ('max, min reanal_mean_NSSTgrid = ', ma.max(reanal_mean_NSSTgrid), \
    ma.min(reanal_mean_NSSTgrid))
    
# ---- plot difference

sst_difference = nsst_mean - reanal_mean_NSSTgrid

# ---- code to use for plotting

fig = plt.figure(figsize=(9.,5.6))
axloc = [0.07,0.11,0.9,0.82]
ax = fig.add_axes(axloc)
ax.set_title('00 UTC global skin temperature differences, operational minus reanalysis, 1 Sep 2019 to 31 Dec 2019')
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
#clevs = [-2.0,-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0]
clevs = [-6,-4,-2,-1,-0.5,0.5,1,2,4,6]
colorstblack='Black'
parallels = np.arange(-80.,90,20.)
meridians = np.arange(0.,360.,20.)
m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=-75,\
    urcrnrlon=lons[-1,-1],urcrnrlat=75.,\
    projection='mill',resolution='l')
x, y = m(lons, lats)

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

