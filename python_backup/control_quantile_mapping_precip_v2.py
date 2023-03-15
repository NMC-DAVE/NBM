"""
control_quantile_mapping_precip_v2.py cyyyymmddhh clead

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
import pygrib
from quantile_mapping_gamma_mixture_v2_f90 import \
    quantile_mapping_gamma_mixture_v2_f90
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

import _pickle as cPickle
import scipy.stats as stats


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

def set_domain_boundaries(cdomain):

    """ used grib file of 2.5-km blend output grid to determine bounding
        lat and lon, and from that, the domain bounding indices for the
        0.25 GEFSv12 reforecast data that will encompass the domain.
    """
    if cdomain == 'conus':
        jmin = 93
        jmax = 246
        imin = 368
        imax = 686
    elif cdomain == 'pr':
        jmin = 243
        jmax = 256
        imin = 649
        imax = 667
    elif cdomain == 'ak':
        jmin = 19
        jmax = 161
        imin = 201
        imax = 967
    else:
        print ('invalid domain.  Exiting.')
        sys.exit()

    return jmin, jmax, imin, imax
    
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

# --- read grib data on a single level

def read_gribdata(gribfilename, endStep):
    istat = -1
    fexist_grib = False
    fexist_grib = os.path.exists(gribfilename)
    print (gribfilename, endStep)
    if fexist_grib:
        try:
            fcstfile = pygrib.open(gribfilename)
            #print ('opened ', gribfilename)
            #grb = fcstfile.select(shortName='tp',\
            #    validityDate=validityDate, \
            #    validityTime=validityTime)[0]
            grb = fcstfile.select(shortName='tp',endStep=endStep)[0]
            #print ('selected grb')
            precip_realtime = grb.values
            #print ('read values')
            lats_full, lons_full = grb.latlons()
            #print ('got lat/lon ',lats_full[0,0], lons_full[0,0])
            istat = 0
            fcstfile.close()
            #print ('closed file ')
        except IOError:
            print ('   IOError in read_gribdata reading ', \
                gribfilename, validityDate, validityTime)
            istat = -1
        except ValueError:
            print ('   ValueError in read_gribdata reading ', \
                gribfilename, validityDate, validityTime)
            istat = -1
        except RuntimeError:
            print ('   RuntimeError in read_gribdata reading ', \
                gribfilename, validityDate, validityTime)
            istat = -1
    return istat, precip_realtime, lats_full, lons_full


    
# =====================================================================

# ---- inputs from command line

nstride = 1
cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
mswep_directory = '/Volumes/Backup Plus/mswep/'

# ---- read in the previously generated netCDF file with precipitation
#      to get lat/lon of MSWEP grid

infile = mswep_directory + '200001_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile)
lons_mswep = nc.variables['lons'][:,:]
lons_mswep = lons_mswep - 360.0
lats_mswep = nc.variables['lats'][:,:]
#print ('lats_mswep[0], [-1] = ', lats_mswep[0], lats_mswep[-1])
#print ('mswep lons dtype = ',lons_mswep.dtype)
nc.close()
   
    
# ---- read the MSWEP fitted Gamma parameters from cPickle file

data_directory = '/Volumes/Backup Plus/mswep/'
infile = data_directory + cmonth+'_conus'+\
    '_MSWEP_apcp_gamma_parameters_h'+clead+'.cPick'       
    
print ('reading from ', infile)
inf = open(infile, 'rb')
weights_mswep = cPickle.load(inf)
alpha_mswep = cPickle.load(inf)
beta_mswep = cPickle.load(inf)
fzero_mswep = cPickle.load(inf)

Dn = cPickle.load(inf)
nmixture_mswep = cPickle.load(inf)
ny_mswep, nx_mswep = np.shape(nmixture_mswep)
ncomponents = 3
print ('ny_nswep, nx_mswep = ',ny_mswep, nx_mswep)
inf.close()

# ---- read the GEFSv12 fitted Gamma parameters from cPickle file

cdomain = 'conus'
gefs_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'
infile = gefs_directory + ccmonth+'_'+cdomain+\
    '_apcp_gamma_parameters_h' + clead + '.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
weights_gefsv12 = cPickle.load(inf)
alpha_gefsv12 = cPickle.load(inf)
beta_gefsv12 = cPickle.load(inf)
fzero_gefsv12 = cPickle.load(inf)

Dnstat1a_gefsv12 = cPickle.load(inf)
Dnstat2a_gefsv12 = cPickle.load(inf)
Dnstat3a_gefsv12 = cPickle.load(inf)
nmixture_gefsv12 = cPickle.load(inf)
inf.close()

# ---- get the GEFSv12 domain subgrid lat/lon

jmin, jmax, imin, imax = set_domain_boundaries(cdomain)

ncfile = gefs_directory + ccmonth + '_apcp_h' + clead + '.nc'
nc = Dataset(ncfile)
lons_1d_gefsv12 = nc.variables['lons_fcst'][imin:imax]
lats_1d_gefsv12 = nc.variables['lats_fcst'][jmin:jmax]
lons_1d_gefsv12_full = nc.variables['lons_fcst'][:]
lats_1d_gefsv12_full = nc.variables['lats_fcst'][:]
print ('lons_1d_gefsv12 dtype = ', lons_1d_gefsv12.dtype)
nx_gefsv12 = len(lons_1d_gefsv12)
ny_gefsv12 = len(lats_1d_gefsv12)
print ('ny_gefsv12, nx_gefsv12 = ', ny_gefsv12, nx_gefsv12)
nc.close()

lons_fcst_2d, lats_fcst_2d = np.meshgrid(lons_1d_gefsv12,lats_1d_gefsv12)

# ---- get the desired 2021 GEFSv12 forecast as grib file downloaded
#      from NOMADS server

input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
infile = input_directory + cyyyymmddhh + \
    '_gec00.t00z.pgrb2s.0p25.f0' + clead 
    
endStep = int(clead)
istat, precip_realtime, lats_full, lons_full = \
    read_gribdata(infile, endStep)
    
if lats_full[0,0] > lats_full[-1,0]: 
    lats_full = np.flipud(lats_full)
    precip_realtime = np.flipud(precip_realtime)   
lons_full = lons_full - 360.
    
print (lons_full[0,:])  
print (lats_full[:,0])  
#sys.exit() 
    
precip_gefsv12_on_mswep = interp(precip_realtime, \
    lons_full[0,:], lats_full[:,0], \
    lons_mswep, lats_mswep, checkbounds=False, \
    masked=False, order=1)    

# ---- now call the fortran routine to perform the quantile mapping 
#      more quickly


print ('np.shape(weights_mswep) = ', np.shape(weights_mswep) )
print ('np.shape(alpha_mswep) = ', np.shape(alpha_mswep) )
print ('np.shape(beta_mswep) = ', np.shape(beta_mswep) )
print ('np.shape(fzero_mswep) = ', np.shape(fzero_mswep) )

print ('np.shape(weights_gefsv12) = ', np.shape(weights_gefsv12) )
print ('np.shape(alpha_gefsv12) = ', np.shape(alpha_gefsv12) )
print ('np.shape(beta_gefsv12) = ', np.shape(beta_gefsv12) )
print ('np.shape(fzero_gefsv12) = ', np.shape(fzero_gefsv12) )

print (quantile_mapping_gamma_mixture_v2_f90.__doc__)
qmapped_precip = np.zeros((ny_mswep, nx_mswep), dtype=np.float64)
ncomponents = 3
qmapped_precip = quantile_mapping_gamma_mixture_v2_f90( \
    weights_mswep, alpha_mswep, beta_mswep, fzero_mswep, \
    lons_mswep, lats_mswep, weights_gefsv12, alpha_gefsv12, \
    beta_gefsv12, fzero_gefsv12, precip_gefsv12_on_mswep, \
    lons_1d_gefsv12, lats_1d_gefsv12, ncomponents, \
    ny_mswep, nx_mswep, ny_gefsv12, nx_gefsv12 )
        
# ---- plot the raw forecast.

m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_mswep, lats_mswep)

clevs = [0.0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Quantile mapped forecast, IC = '+cyyyymmddhh+' lead = '+clead+' h'
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, qmapped_precip, clevs,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'qmapped_precip_'+cyyyymmddhh+'_lead'+clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')






# ---- first plot GEFS control forecast amount

m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_mswep, lats_mswep)

# ---- make plots of fraction positive precip

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

clevs = [0.0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20]
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Control forecast precipitation amount (mm) for '+clead+\
    '-h IC = '+cyyyymmddhh
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, precip_gefsv12_on_mswep, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Mean precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'forecast_precip_'+clead+'_h_IC'+cyyyymmddhh+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


