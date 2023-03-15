"""
python plot_splines_forecast_analysis.py cmonth cend_hour clat clon 

where cend_hour is 00, 06, 12, 18

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
from qmapping_spline_1pt import qmapping_spline_1pt # from f2py of fortran90 code.
from control_splev import control_splev

import _pickle as cPickle
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

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
    
    precip_ens_nonzero = np.delete(precip_ens, \
        np.where(precip_ens <= 0.0))  # censor at 0.0 mm
    precip_ens_nonzero = precip_ens_nonzero + \
        np.random.uniform(low=-0.1,high=0.1,size=len(precip_ens_nonzero))
    precip_ens_nonzero = np.delete(precip_ens_nonzero, \
        np.where(precip_ens_nonzero <= 0.0))  # censor at 0.0 mm
    nz = len(precip_ens_nonzero)    
        
    precip_ens_nonzero = np.sort(precip_ens_nonzero)  
    ntotal = len(precip_ens)
    nzero = ntotal - len(precip_ens_nonzero)
    fraction_zero = float(nzero) / float(ntotal)
    return fraction_zero, precip_ens_nonzero, nz

# =====================================================================

def read_analysis_spline_info (clead_use, master_directory_panal_spline,\
    cmonth):

    """ read netCDF file here for inverse spline parameters for
    the combined CCPA/MSWEP precip analysis CDFs on the 
    Lambert-conformal NDFD 2.5 km grid surrounding the CONUS.  
  
    **** Note that if we are applying to cycles other than 
    forecasts with 00UTC init, we'll need adjust the dates 
    of the spline files that we use to sync up. 
    """ 
   
    ndays = int(clead_use) // 24
    ilead = int(clead_use)-ndays*24
    if ilead == 0:
        cleada = '00'
    elif ilead == 6:
        cleada = '06'
    elif ilead == 12:
        cleada = '12'
    elif ilead == 18:
        cleada = '18'
        
    infile = master_directory_panal_spline + cmonth + \
        '_conus_CCPA_spline_info_h' + cleada + 'UTC.nc'
    #infile = master_directory_panal_spline + cmonth+\
    #    '_conus_MSWEP_spline_info_h' + cleada + 'UTC.nc' 
    print ('reading from ', infile)
    nc = Dataset(infile)
    spline_info_inv = nc.variables['spline_info_inv'][:,:,:,:]
    fraction_zero_ndfd = nc.variables['fzero'][:,:]
    number_knots = nc.variables['number_knots'][:,:]
    usegamma_ndfd = nc.variables['usegamma'][:,:]
    lons_ndfd = nc.variables['lons'][:,:]
    lons_ndfd = lons_ndfd
    lats_ndfd = nc.variables['lats'][:,:]
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
    nc.close()     
    
    return spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
        number_knots, lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd
        
# ===================================================================== 

def read_forecast_spline_info(ccmonth, \
    master_directory_forecast_spline):

    """ load spline information for forecast netCDF file.
        For initial times other than 00 UTC, we need to find the
        appropriate spline file to read in based on 00 UTC
        forecast initialization times.   Assuming that
        the CDF characteristics are primarily a function of the
        diurnal cycle, make this adjustment. """

    ccycle = '00'
    icycle = int(ccycle)
    ilead = int(clead)
    ilead_use = icycle + ilead # diurnal cycle adjust
    if ilead_use > 240: ilead_use = ilead_use - 24
    if ilead_use < 10: # make 3-character variable of ilead_use
        clead_use = '00' + str(ilead_use)
    elif ilead_use >= 10 and ilead_use < 100:
        clead_use = '0' + str(ilead_use)
    else:
        clead_use = str(ilead_use)

    infile = master_directory_forecast_spline + \
        ccmonth + '_conus_GEFSv12_spline_info_h' + clead_use + '.nc'
    print ('reading forecast spline information from ', infile)
    nc = Dataset(infile)
    lons_spline_gefsv12_1d = nc.variables['lons'][:]
    lats_spline_gefsv12_1d = nc.variables['lats'][:]
    spline_info_gefsv12 = nc.variables['spline_info'][:,:,:,:]
    fraction_zero_gefsv12 = nc.variables['fzero'][:,:]
    usegamma_gefsv12 = nc.variables['usegamma'][:,:]
    quantile_99_gefsv12 = nc.variables['quantile_99'][:,:]
    number_knots_gefsv12 = nc.variables['number_knots'][:,:]
    ny_gefsv12, nx_gefsv12 = np.shape(usegamma_gefsv12)
    nc.close()
    lons_fcst_2d, lats_fcst_2d = \
        np.meshgrid(lons_spline_gefsv12_1d,lats_spline_gefsv12_1d)

    return lons_spline_gefsv12_1d, lats_spline_gefsv12_1d, \
        spline_info_gefsv12, fraction_zero_gefsv12, \
        usegamma_gefsv12, quantile_99_gefsv12, \
        number_knots_gefsv12, ny_gefsv12, nx_gefsv12, \
        lons_fcst_2d, lats_fcst_2d, clead_use

# =====================================================================

def find_nearest(vec, value):
    
    """ given a vector vec and a particular value, find the index in vec
    that is nearest to value"""
    
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

def get_surrounding_months(cmonth):

    if cmonth == 'Jan':
        cmonth_middle = 'Jan'
        cmonth_early = 'Dec'
        cmonth_late = 'Feb'
    elif cmonth == 'Feb':
        cmonth_middle = 'Feb'
        cmonth_early = 'Jan'
        cmonth_late = 'Mar'
    elif cmonth == 'Mar':
        cmonth_middle = 'Mar'
        cmonth_early = 'Feb'
        cmonth_late = 'Apr'
    elif cmonth == 'Apr':
        cmonth_middle = 'Apr'
        cmonth_early = 'Mar'
        cmonth_late = 'May'
    elif cmonth == 'May':
        cmonth_middle = 'May'
        cmonth_early = 'Apr'
        cmonth_late = 'Jun'
    elif cmonth == 'Jun':
        cmonth_middle = 'Jun'
        cmonth_early = 'May'
        cmonth_late = 'Jul'
    elif cmonth == 'Jul':
        cmonth_middle = 'Jul'
        cmonth_early = 'Jun'
        cmonth_late = 'Aug'
    elif cmonth == 'Aug':
        cmonth_middle = 'Aug'
        cmonth_early = 'Jul'
        cmonth_late = 'Sep'
    elif cmonth == 'Sep':
        cmonth_middle = 'Sep'
        cmonth_early = 'Aug'
        cmonth_late = 'Oct'
    elif cmonth == 'Oct':
        cmonth_middle = 'Oct'
        cmonth_early = 'Sep'
        cmonth_late = 'Nov'
    elif cmonth == 'Nov':
        cmonth_middle = 'Nov'
        cmonth_early = 'Oct'
        cmonth_late = 'Dec'
    elif cmonth == 'Dec':
        cmonth_middle = 'Dec'
        cmonth_early = 'Nov'
        cmonth_late = 'Jan'
    else:
        print ('invalid month')
        sys.exit()
    
    return cmonth_middle, cmonth_early, cmonth_late       
        
# =====================================================================


# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
cend_hour = sys.argv[2] # 06, 12, 18, 00
clat = sys.argv[3] 
clon = sys.argv[4] 
clead = '0'+cend_hour

# ---- define constants.

rlon_input = float(clon)
rlat_input = float(clat)
cdomain = 'conus'
pflag = True

# ================================================================
# --- load analysis spline info, compute for selected point
# ================================================================

master_directory_panal_spline = \
    '/Volumes/NBM/conus_panal/CDF_spline/'
spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, number_knots, \
    lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd = \
    read_analysis_spline_info (cend_hour, \
    master_directory_panal_spline, cmonth)
if np.max(lons_ndfd) > 180.: lons_ndfd = lons_ndfd - 180.
    
# ---- find index of nearest ndfd point

distx = np.abs(lons_ndfd - rlon_input)
disty = np.abs(lats_ndfd - rlat_input)
dist = distx**2 + disty**2
i = np.argmin(dist)
print (i)
jndfd, indfd = np.unravel_index(i, shape=(ny_ndfd, nx_ndfd), order='C')   

# ---- get the spline fit.

cdf = 0.0005 + np.arange(1000.)/1000.
hazard_fn = -np.log(1.0 - cdf)
numknots = number_knots[jndfd, indfd]  # now this number includes boundary
knots_in = spline_info_inv[jndfd, indfd,0,0:numknots]
bspline_coef_in = spline_info_inv[jndfd, indfd,1,0:numknots]
nthree = 3
ier = 0
quantile_mapped_precip = np.zeros((1000), dtype=np.float64)
ier, quantile_mapped_precip = control_splev(knots_in, \
    bspline_coef_in, hazard_fn, numknots, 1000)
cdf_with_zeros = fraction_zero_ndfd[jndfd, indfd] + \
    (1.-fraction_zero_ndfd[jndfd, indfd])*cdf
    
# ================================================================
# --- load forecast spline info, compute for selected point.
# ================================================================

master_directory_forecast_spline = '/Volumes/NBM/conus_gefsv12/CDF_spline/'
print ('reading forecast spline parameters')
lons_1d_realtime, lats_1d_realtime, \
    spline_info_gefsv12, fraction_zero_gefsv12, \
    usegamma_gefsv12, quantile_99_gefsv12, \
    number_knots_gefsv12, ny_gefsv12, nx_gefsv12, \
    lons_fcst_2d, lats_fcst_2d, clead_use = \
    read_forecast_spline_info( \
    cmonth, master_directory_forecast_spline)
if np.max(lons_1d_realtime) > 180.: \
    lons_1d_realtime = lons_1d_realtime - 180.
lons_2d, lats_2d = np.meshgrid(lons_1d_realtime, lats_1d_realtime)
    
# ---- find index of nearest grid point for spline data

print ('rlon_input, rlat_input = ', rlon_input, rlat_input)
distx = np.abs(lons_2d - rlon_input)
disty = np.abs(lats_2d - rlat_input)
dist = distx**2 + disty**2
i = np.argmin(dist)
jspline, ispline = np.unravel_index(i, shape=(ny_gefsv12, nx_gefsv12), order='C')    
    
# ---- get the spline fit.

numknots = number_knots_gefsv12[jspline, ispline]  # now this number includes boundary
knots_in = spline_info_gefsv12[jspline, ispline,0,0:numknots]
bspline_coef_in = spline_info_gefsv12[jspline, ispline,1,0:numknots]
precip_regular = 0.0001+np.arange(0.,251.,1.)/10.
nz = len(precip_regular)
hazard_fn = -np.log(1.0 - precip_regular)
npr = len(precip_regular)
nthree = 3
ier = 0
cdf_spline = np.zeros((nz), dtype=np.float64)
splines_tuple = (spline_info_gefsv12[jspline, ispline,0,0:numknots], \
    spline_info_gefsv12[jspline, ispline,1,0:numknots], 3)
spline_hazard = splev(precip_regular, splines_tuple)
qpositive = 1.0 - np.exp(-spline_hazard)    
cdf_spline_with_zeros = fraction_zero_gefsv12[jspline, ispline] + \
    (1.-fraction_zero_gefsv12[jspline, ispline])*qpositive
    
# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('Spline fit to analysis & forecast,\n'+\
    cmonth+', '+cend_hour+' UTC, '+\
    clon+' W '+clat+' N',fontsize=13)

ax.plot(quantile_mapped_precip, cdf_with_zeros,\
    color='Blue',lw=2,label='Analysis spline')
ax.plot(precip_regular, cdf_spline_with_zeros,\
    color='Red',lw=2,label='Forecast spline')

plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.8,1.0)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,30)
ax.set_yticks(np.arange(0.8,1.01,0.01))
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_forecast_analysis_example_spline.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 
              
              

