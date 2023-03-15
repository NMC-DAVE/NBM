"""
stencil_example_ith.py

coded by: Tom Hamill, Aug 2021, tom.hamill@noaa.gov 

"""

import os, sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, interp
import matplotlib.pyplot as plt
from matplotlib import rcParams

# =====================================================================

def read_gribdata(gribfilename, endStep):
    
    """ read grib data"""
    
    istat = -1
    fexist_grib = False
    fexist_grib = os.path.exists(gribfilename)
    print ('   reading ',gribfilename, fexist_grib)
    if fexist_grib:
        try:
            fcstfile = pygrib.open(gribfilename)
            grb = fcstfile.select(shortName='tp',endStep=endStep)[0]
            precip_realtime = grb.values
            lats_full, lons_full = grb.latlons()
            istat = 0
            fcstfile.close()
        except IOError:
            print ('   IOError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1
        except ValueError:
            print ('   ValueError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1
        except RuntimeError:
            print ('   RuntimeError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1
    return istat, precip_realtime, lats_full, lons_full

# =====================================================================
    
def get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead):

    """ read in the global forecast.  Subset to CONUS. """

    nib = 886  # precomputed for set of lat-lons that surrounds
    nie = nib+318 # the NDFD CONUS 2.5-km Lambert conformal domain.
    njb = 131 
    nje = njb+153 
    
    # ---- read in forecast grid covering the whole globe.
    
    cycle = cyyyymmddhh[8:10]
    #input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
    input_directory = '/Volumes/NBM/gefsv12/2018_retro/'
    #infile = input_directory + cyyyymmddhh + \
    #    '_ge'+cmem+'.t'+cycle+'z.pgrb2s.0p25.f' + clead  
    infile = input_directory + 'apcp_'+\
        cyyyymmddhh+'_'+cmem+'_.f'+clead+'.grib2'
        
    print ('   reading from ', infile)
    endStep = int(clead)
    istat, precip_realtime, lats_full, lons_full = \
        read_gribdata(infile, endStep)
    
    # ---- subset for CONUS.
    
    precip_realtime = precip_realtime[njb:nje,nib:nie]
    lons_1D_realtime = lons_full[0,nib:nie]-360.
    lats_1D_realtime = lats_full[njb:nje,0]

    return precip_realtime, lons_1D_realtime, lats_1D_realtime    
                
    
# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = '2018021200'    
cmem = 'c00'       
clead = '024'
nstride_qmap = 15 

# ---- various initialization and carving out month information.

master_directory = '/Volumes/NBM/conus_panal/'
infile = master_directory + '200201_ccpa_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile)
lons = nc.variables['lons'][:,:]
lats = nc.variables['lats'][:,:]
ny, nx = np.shape(lons)
nc.close()

londiff = lons + 76.5
latdiff = lats - 42.5
diff = np.sqrt(londiff**2 + latdiff**2)

ind = np.unravel_index(np.argmin(diff, axis=None), diff.shape)
jymin = ind[0]
ixmin = ind[1]
print ('minimum at ',jymin, ixmin, lons[jymin,ixmin], lats[jymin,ixmin])

title = '(a) Stencil at 24 h'
clead = '024'
nstr = 12 + int(clead)//18
m = Basemap(llcrnrlon=-79,llcrnrlat=40,
    urcrnrlon = -74, urcrnrlat = 45,\
    projection='lcc',lat_1=42.5,lat_2=42.5,lon_0=285.,\
    resolution ='i',area_thresh=1000.)
x, y = m(lons, lats)
fig = plt.figure(figsize=(7,4.5))
axloc = [0.03,0.05,0.43,0.88]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=16,color='Black')
for ix in range(-2,3):
    for jy in range(-2,3):
        xdot, ydot = m(lons[jymin+jy*nstr,ixmin+ix*nstr],\
            lats[jymin+jy*nstr,ixmin+ix*nstr])
        if ix == 0 and jy == 0:
            m.plot(xdot,ydot,marker='.',markersize=7,color='Red') 
        else:
            m.plot(xdot,ydot,marker='.',markersize=4,color='Black')
m.drawstates(linewidth=0.6,color='Gray')
m.drawcounties(linewidth=0.2,color='LightGray')
m.drawcoastlines(linewidth=0.9,color='Gray')


title = '(b) Stencil at 240 h'
clead = '240'
nstr = 12 + int(clead)//18
axloc = [0.53,0.05,0.43,0.88]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=16,color='Black')
for ix in range(-2,3):
    for jy in range(-2,3):
        xdot, ydot = m(lons[jymin+jy*nstr,ixmin+ix*nstr],\
            lats[jymin+jy*nstr,ixmin+ix*nstr])
        if ix == 0 and jy == 0:
            m.plot(xdot,ydot,marker='.',markersize=7,color='Red') 
        else:
            m.plot(xdot,ydot,marker='.',markersize=4,color='Black')
m.drawstates(linewidth=0.6,color='Gray')
m.drawcounties(linewidth=0.2,color='LightGray')
m.drawcoastlines(linewidth=0.9,color='Gray')


# ---- set plot title

plot_title = 'stencil_ith.png'
fig.savefig(plot_title, dpi=400)
print ('saving plot to file = ',plot_title)
print ('Done!')

