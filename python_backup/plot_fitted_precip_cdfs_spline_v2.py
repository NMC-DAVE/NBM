"""
plot_fitted_precip_cdfs_allmixtures.py cmonth clead jy ix

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
from scipy.interpolate import LSQUnivariateSpline, splrep, splev

import _pickle as cPickle
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

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

def find_nearest(vec, value):
    
    """ given a vector vec and a particular value, find the index in vec
    that is nearest to value"""
    
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

def get_surrounding_months(cmonth):

    if cmonth == 'Jan':
        cmonth_early = 'Dec'
        cmonth_late = 'Feb'
    elif cmonth == 'Feb':
        cmonth_early = 'Jan'
        cmonth_late = 'Mar'
    elif cmonth == 'Mar':
        cmonth_early = 'Feb'
        cmonth_late = 'Apr'
    elif cmonth == 'Apr':
        cmonth_early = 'Mar'
        cmonth_late = 'May'
    elif cmonth == 'May':
        cmonth_early = 'Apr'
        cmonth_late = 'Jun'
    elif cmonth == 'Jun':
        cmonth_early = 'May'
        cmonth_late = 'Jul'
    elif cmonth == 'Jul':
        cmonth_early = 'Jun'
        cmonth_late = 'Aug'
    elif cmonth == 'Aug':
        cmonth_early = 'Jul'
        cmonth_late = 'Sep'
    elif cmonth == 'Sep':
        cmonth_early = 'Aug'
        cmonth_late = 'Oct'
    elif cmonth == 'Oct':
        cmonth_early = 'Sep'
        cmonth_late = 'Nov'
    elif cmonth == 'Nov':
        cmonth_early = 'Oct'
        cmonth_late = 'Dec'
    elif cmonth == 'Dec':
        cmonth_early = 'Nov'
        cmonth_late = 'Jan'
    else:
        print ('invalid month')
        sys.exit()
    
    return cmonth_early, cmonth_late

# =====================================================================

def fraczero_possamps(nsamps, precip_ens):
    
    """
    
    from the vector input sample precip_ens, define the fraction of
    samples with zero precipitation.   For the positive samples, add
    a small random number to deal with the fact that the data was 
    discretized to ~0.1 mm, so that when later creating CDFs we don't 
    have empirical values with lots of tied amounts.  Also, sort the 
    nonzero amounts and return.
    
    """
    number_zeros = 0
    
    # data discretized, so add random component of this magnitude
    #precip_ens = np.where(precip_ens < 4.0, precip_ens* \
    #    np.random.uniform(low=-0.5,high=1.5,size=len(precip_ens)), precip_ens)
    
    precip_ens_nonzero = np.delete(precip_ens, \
        np.where(precip_ens <= 0.0))  # censor at 0.0 mm
    precip_ens_nonzero = precip_ens_nonzero + \
        np.random.uniform(low=-0.1,high=0.1,size=len(precip_ens_nonzero))
    precip_ens_nonzero = np.delete(precip_ens_nonzero, \
        np.where(precip_ens_nonzero <= 0.0))  # censor at 0.0 mm
    nz = len(precip_ens_nonzero)    
        
    #precip_ens_nonzero = precip_ens_nonzero + \
    #    np.random.uniform(low=-0.005,high=0.005,size=nz) 
    precip_ens_nonzero = np.sort(precip_ens_nonzero)  
    ntotal = len(precip_ens)
    nzero = ntotal - len(precip_ens_nonzero)
    fraction_zero = float(nzero) / float(ntotal)
    return fraction_zero, precip_ens_nonzero, nz

# =====================================================================
    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clead = sys.argv[2] # 03, 06, 12, etc.
clat = sys.argv[3] 
clon = sys.argv[4] 
cdomain = 'conus'
cmonth_early, cmonth_late = get_surrounding_months(cmonth)
jmin, jmax, imin, imax = set_domain_boundaries(cdomain)
pflag = True

# --- load from cPickle file the spline parameters and Dn statistics
#     for the GEFSv12 forecasts

master_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'

infile = master_directory + cmonth+'_'+cdomain+\
    '_GEFSv12_spline_info_h' + clead + '.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
spline_info = cPickle.load(inf)
indices_to_query = cPickle.load(inf)
print ('np.shape(spline_info) = ', np.shape(spline_info))
print ('np.shape(indices_to_query) = ', np.shape(indices_to_query))

inf.close()


infile = master_directory + cmonth+'_'+cdomain+\
    '_GEFSv12_Dnstat_h' + clead + '.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
Dnstat = cPickle.load(inf)
lons_1d = cPickle.load(inf)
lats_1d = cPickle.load(inf)
print ('max, min lons_1d = ', np.max(lons_1d), np.min(lons_1d))
print ('max Dnstat = ', np.max(Dnstat))
inf.close()


# ---- read in the precipitation forecasts used to build CDF
#      using the current month and surrounding two months

ncfile = master_directory + cmonth + '_apcp_h' + clead + '.nc'
if pflag == True: print (ncfile)
nc = Dataset(ncfile)
precip_middle = nc.variables['apcp_fcst'][:,jmin:jmax,imin:imax]
lons_fcst = nc.variables['lons_fcst'][imin:imax]
print ('max, min lons_fcst = ', np.max(lons_fcst), np.min(lons_fcst))
lats_fcst = nc.variables['lats_fcst'][jmin:jmax]
nsamps_middle = len(precip_middle)
nc.close()

ncfile = master_directory + cmonth_early + '_apcp_h' + clead + '.nc'
if pflag == True: print (ncfile)
nc = Dataset(ncfile)
precip_early = nc.variables['apcp_fcst'][:,jmin:jmax,imin:imax]
nsamps_early = len(precip_early)
nc.close()

ncfile = master_directory + cmonth_late + '_apcp_h' + clead + '.nc'
if pflag == True: print (ncfile)
nc = Dataset(ncfile)
precip_late = nc.variables['apcp_fcst'][:,jmin:jmax,imin:imax]
nsamps_late = len(precip_late)
nc.close()
nsamps = nsamps_middle + nsamps_early + nsamps_late


# ---- determine the nearest grid point given the input lat, lon

rlat = float(clat)
rlon = float(clon)
print ('rlon, rlat = ', rlon, rlat)
ilon = find_nearest(lons_1d, rlon)
jlat = find_nearest(lats_1d, rlat)
print ('lats_1d = ', lats_1d)
print ('lons_1d = ', lons_1d)
print ('lons_1d[ilon] = ', lons_1d[ilon])
print ('lats_1d[jlat] = ', lats_1d[jlat])
print ('nearest grid point  = ', jlat, ilon)

precip_middle_1d = precip_middle[:,jlat,ilon]
precip_early_1d = precip_early[:,jlat,ilon]
precip_late_1d = precip_late[:,jlat,ilon]

precip_ens_1d = np.concatenate((precip_middle_1d, \
    precip_early_1d, precip_late_1d))

fraction_zero, precip_ens_nonzero, nz = \
    fraczero_possamps(nsamps, precip_ens_1d) # return sorted

# ---- build empirical CDF for this grid point's data

x = np.arange(0.0,50.01,0.1)
nx = len(x)
len_nonzero = len(precip_ens_nonzero)
fnzero = float(len_nonzero)
nz = int(fnzero)
empirical_cdf = 1.0 / (2.0*nz) + np.arange(nz)/nz  
#cdf_empirical = np.zeros((nx),dtype=np.float64)
#for i, xi in enumerate(x):
#    nbelow = (precip_ens_nonzero < xi).sum() 
#    cdf_empirical[i] = fraction_zero + \
#        (1.-fraction_zero)*float(nbelow)/fnzero

splines_tuple = (spline_info[jlat,ilon,0,:], spline_info[jlat,ilon,1,:], 3)
spline_hazard = splev(precip_ens_nonzero, splines_tuple)
spline_cdf = 1.0 - np.exp(-spline_hazard)

# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('(a) Spline fit, '+cmonth+', '+clead+'-h forecast, '+\
    clon+' W '+clat+' N',fontsize=13)
ax.plot(precip_ens_nonzero, spline_cdf,color='Blue',lw=2,label='Cubic spline')
ax.plot(precip_ens_nonzero, empirical_cdf,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.6,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_spline_GEFSv12.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 
              
                










