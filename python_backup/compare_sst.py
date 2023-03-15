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

date_list_nsst_00 = daterange('2019090100', '2020042300', 24)
date_list_nsst_06 = daterange('2019090106', '2020042306', 24)
date_list_nsst_12 = daterange('2019090112', '2020042312', 24)
date_list_nsst_18 = daterange('2019090118', '2020042318', 24)
ndates = len(date_list_nsst_00)
infile_2019 = 'nsst/sst.day.mean.2019.nc'
infile_2020 = 'nsst/sst.day.mean.2020.nc'
infile_ice_2019 = 'nsst/icec.day.mean.2019.nc'
infile_ice_2020 = 'nsst/icec.day.mean.2020.nc'
missing_list_yyyymmdd = []
missing_list_hour = []

# ---- read in the grib files of the NSST analyses.

for iset in range(4):
    if iset == 0:
        date_list = date_list_nsst_00
        cycle = '00'
    elif iset == 1:
        date_list = date_list_nsst_06
        cycle = '06'
    elif iset == 2:
        date_list = date_list_nsst_12
        cycle = '12'
    else:
        date_list = date_list_nsst_18
        cycle = '18'
      
    for idate, date in enumerate(date_list):
        
        yy = date[0:4]
        yyyymmdd = date[0:8]
        infile = 'nsst/'+date[0:8]+'_gdas.t'+cycle+'z.pgrb2.0p25.f000'
        #print (infile)
        fexist = path.exists(infile)
        if fexist == True:
            sstfile = pygrib.open(infile)
            grb = sstfile.select(shortName = 't')[0]
            sst = grb.values
            if idate == 0 and iset == 0:
                lats, lons = grb.latlons()
                #print (lats)
                #sys.exit()
                nlats, nlons = np.shape(lats)
                #print ('nlats, nlons = ', nlats, nlons)
                #sys.exit()
                nsst_save = ma.zeros((ndates,nlats,nlons), dtype=np.float32)
                OI_save = ma.zeros((ndates,nlats-1,nlons), dtype=np.float32)
                icec_save = ma.zeros((ndates,nlats-1,nlons), dtype=np.float32)
                nsst_mean_allcycles = np.zeros((4, nlats, nlons), dtype=np.float32)
            sstfile.close()
            nsst_save[idate,:,:] = sst[:,:]
            #print ('after writing to nsst_save')
            
            if yy == '2019':            
                yyyy = int(date[0:4])
                mm = int(date[4:6])
                dd = int(date[6:8])
                julday = dayofyear(yyyy,mm,dd) - 1 # so day 1 = 0 for python index
                nc = Dataset(infile_2019)
                sstin = nc.variables['sst'][julday,:,:]
                if idate == 0:
                    lon_OI = nc.variables['lon'][:]
                    lat_OI = nc.variables['lat'][:]
                nc.close()
                OI_save[idate,:,:] = sstin[:,:]
                
                nc = Dataset(infile_ice_2019)
                icein = nc.variables['icec'][julday,:,:]
                icec_save[idate,:,:] = icein[:,:]
                nc.close()
                
            else:
                yyyy = int(date[0:4])
                mm = int(date[4:6])
                dd = int(date[6:8])
                julday = dayofyear(yyyy,mm,dd) - 1 # so day 1 = 0 for python index
                nc = Dataset(infile_2020)
                sstin = nc.variables['sst'][julday,:,:]
                nc.close()
                OI_save[idate,:,:] = sstin[:,:]  
                
                nc = Dataset(infile_ice_2020)
                icein = nc.variables['icec'][julday,:,:]
                icec_save[idate,:,:] = icein[:,:]  
                nc.close()               
        else:
            missing_list_yyyymmdd.append(str(yyyymmdd))
            missing_list_hour.append(cycle)
            print ('unable to read ', infile)
            nsst_save[idate,:,:] = ma.masked
            OI_save[idate,:,:] = ma.masked
    
        #sys.exit()
    nsst_mean_allcycles[iset,:,:] = ma.mean(nsst_save,axis=0)
        
nsst_mean_overall = np.mean(nsst_mean_allcycles, axis=0)
sst_OI_mean = ma.mean(OI_save,axis=0)
icec_OI_mean = ma.mean(icec_save,axis=0)
sst_OI_mean = ma.masked_where(ma.logical_or(icec_OI_mean>0, sst_OI_mean.mask == True), sst_OI_mean)

nsst_mean_overall = np.flipud(nsst_mean_overall)  # need latitudes ascending order 
lons = np.flipud(lons)
lats = np.flipud(lats)
print ('OI max, min = ', ma.max(sst_OI_mean), ma.min(sst_OI_mean))
print ('nsst_mean_overall max, min = ', ma.max(nsst_mean_overall), ma.min(nsst_mean_overall))
# ---- interpolate the nsst data to the OI grid.

lon_1D = lons[0,:]
lat_1D = lats[:,0]
lon_OI_2D, lat_OI_2D = np.meshgrid(lon_OI, lat_OI)

nsst_mean_OIgrid = interp(nsst_mean_overall, lon_1D, lat_1D, \
    lon_OI_2D, lat_OI_2D, checkbounds=False, masked=False, order=1) - 273.15

# ---- plot difference

sst_difference = nsst_mean_OIgrid - sst_OI_mean

# ---- code to use for plotting

fig = plt.figure(figsize=(9.,5.6))
axloc = [0.07,0.11,0.9,0.82]
ax = fig.add_axes(axloc)
ax.set_title('Global differences, NSST - OI SST, 1 Sep 2019 to 22 April 2020')
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
clevs = [-2.0,-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0]
colorstblack='Black'
parallels = np.arange(-80.,90,20.)
meridians = np.arange(0.,360.,20.)
m = Basemap(llcrnrlon=lon_OI[0],llcrnrlat=-75,\
    urcrnrlon=lon_OI[-1],urcrnrlat=75.,\
    projection='mill',resolution='l')
x, y = m(lon_OI_2D, lat_OI_2D)

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

plot_title = 'NSST_vs_OI_difference.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')

ouf = open('nsst/missing_dates_and_cycles.txt', 'w') 
for listitem1, listitem2 in zip(missing_list_yyyymmdd, missing_list_hour):
    line="'{0}', '{1}'\n".format(listitem1,listitem2)
    ouf.write(line)
ouf.close()
        
print ('wrote to nsst/missing_dates_and_cycles.txt')


