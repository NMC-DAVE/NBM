"""
python plot_fcst_analyzed_cdfs.py cmonth clead clat clon
    where cmonth = 01 to 12
    clead = 006 to 240 by 006
    clat = latitude
    clon = longitude (west is negative)

    This simultaneously plots the forecast and analyzed CDFs for 
    chosen point

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
from qmapping_spline_1pt import \
    qmapping_spline_1pt # from f2py of fortran90 code.
from control_splev import control_splev

import _pickle as cPickle
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

# =====================================================================

def find_nearest(vec, value):
    
    """ given a vector vec and a particular value, find the index in vec
    that is nearest to value"""
    
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

def get_surrounding_months(cmonth):

    if cmonth == 'Jan':
        cmonth_middle = '01'
        cmonth_early = '12'
        cmonth_late = '02'
    elif cmonth == 'Feb':
        cmonth_middle = '02'
        cmonth_early = '01'
        cmonth_late = '03'
    elif cmonth == 'Mar':
        cmonth_middle = '03'
        cmonth_early = '02'
        cmonth_late = '04'
    elif cmonth == 'Apr':
        cmonth_middle = '04'
        cmonth_early = '03'
        cmonth_late = '05'
    elif cmonth == 'May':
        cmonth_middle = '05'
        cmonth_early = '04'
        cmonth_late = '06'
    elif cmonth == 'Jun':
        cmonth_middle = '06'
        cmonth_early = '05'
        cmonth_late = '07'
    elif cmonth == 'Jul':
        cmonth_middle = '07'
        cmonth_early = '06'
        cmonth_late = '08'
    elif cmonth == 'Aug':
        cmonth_middle = '08'
        cmonth_early = '07'
        cmonth_late = '09'
    elif cmonth == 'Sep':
        cmonth_middle = '09'
        cmonth_early = '08'
        cmonth_late = '10'
    elif cmonth == 'Oct':
        cmonth_middle = '10'
        cmonth_early = '09'
        cmonth_late = '11'
    elif cmonth == 'Nov':
        cmonth_middle = '11'
        cmonth_early = '10'
        cmonth_late = '12'
    elif cmonth == 'Dec':
        cmonth_middle = '12'
        cmonth_early = '11'
        cmonth_late = '01'
    else:
        print ('invalid month')
        sys.exit()
    
    return cmonth_middle, cmonth_early, cmonth_late

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
    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clead = sys.argv[2] # 006, 012, etc
clat = sys.argv[3] 
clon = sys.argv[4] 
ileadb = int(clead)-6
if ileadb < 10:
    cleadb = '00'+str(ileadb)
else:
    cleadb = '0'+str(ileadb)

ilead_anal = int(clead)%24
if ilead_anal < 10:
    clead_anal = '0'+str(ilead_anal)
else:
    clead_anal = str(ilead_anal)

# ---- define constants.

rlon_input = float(clon)
rlat_input = float(clat)
cdomain = 'conus'
pflag = True
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
cmonths_early = ['12','01','02','03','04','05','06','07','08','09','10','11']
cmonths_late =  ['02','03','04','05','06','07','08','09','10','11','12','01']
cmonths_middle =  ['01','02','03','04','05','06','07','08','09','10','11','12']
cmonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
yearstart = 2002 # CCPA only available starting 2002
yearend = 2020 # companion reforecasts end at end of 2019
imonth = cmonths.index(cmonth)

# ----------------------------------------------------------------
# --- load analysis spline info.
# ----------------------------------------------------------------

master_directory_panal_spline = \
    '/Volumes/NBM/conus_panal/CDF_spline/'
spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, number_knots, \
    lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd = \
    read_analysis_spline_info (clead_anal, \
    master_directory_panal_spline, cmonth)
if np.max(lons_ndfd) > 180.: lons_ndfd = lons_ndfd - 180.
    
# ---- find index of nearest ndfd point

distx = np.abs(lons_ndfd - rlon_input)
disty = np.abs(lats_ndfd - rlat_input)
dist = distx**2 + disty**2
i = np.argmin(dist)
jndfd, indfd = np.unravel_index(i, shape=(ny_ndfd, nx_ndfd), order='C')   

# ---- for diagnostic purposes, also generate a spline fit to CDF
#      with many more percentiles.

nz5000 = 5000
empirical_cdf_hires = 1.0 / (2.0*nz5000) + np.arange(nz5000)/nz5000
hazard_fn_hires = -np.log(1.0 - empirical_cdf_hires) 

# ---- get the spline fit.

numknots = number_knots[jndfd, indfd]  # now this number includes boundary
knots_in = spline_info_inv[jndfd, indfd,0,0:numknots]
bspline_coef_in = spline_info_inv[jndfd, indfd,1,0:numknots]

nthree = 3
ier = 0
quantile_mapped_precip_hires = np.zeros((nz5000), dtype=np.float64)
ier, quantile_mapped_precip_hires = control_splev(knots_in, \
    bspline_coef_in, hazard_fn_hires, numknots, nz5000)

# ----------------------------------------------------------------
# --- load forecast spline info.
# ----------------------------------------------------------------

master_directory_forecast_spline = '/Volumes/NBM/conus_gefsv12/CDF_spline/'
    # where forecast spline info is stored

# ---- read forecast spline parameters from netCDF file

print ('reading forecast spline parameters')
lons_1d_realtime, lats_1d_realtime, \
    spline_info_gefsv12, fraction_zero_gefsv12, \
    usegamma_gefsv12, quantile_99_gefsv12, \
    number_knots_gefsv12, ny_gefsv12, nx_gefsv12, \
    lons_fcst_2d, lats_fcst_2d, clead_anal = \
    read_forecast_spline_info( \
    cmonth, master_directory_forecast_spline)
if np.max(lons_1d_realtime) > 180.: lons_1d_realtime = lons_1d_realtime - 180.
lons_2d, lats_2d = np.meshgrid(lons_1d_realtime, lats_1d_realtime)

# ---- find index of nearest grid point for spline data

distx = np.abs(lons_2d - rlon_input)
disty = np.abs(lats_2d - rlat_input)
dist = distx**2 + disty**2
i = np.argmin(dist)
jspline, ispline = np.unravel_index(i, shape=(ny_gefsv12, nx_gefsv12), order='C')

# ---- get the spline fit.

numknots = number_knots_gefsv12[jspline, ispline]  # now this number includes boundary
knots_in = spline_info_gefsv12[jspline, ispline,0,0:numknots]
bspline_coef_in = spline_info_gefsv12[jspline, ispline,1,0:numknots]
precip_regular = 0.0001+np.arange(0.,501.,1.)/10.
hazard_fn = -np.log(1.0 - precip_regular)
npr = len(precip_regular)
nthree = 3
ier = 0
splines_tuple = (spline_info_gefsv12[jspline, ispline,0,0:numknots], \
    spline_info_gefsv12[jspline, ispline,1,0:numknots], 3)
spline_hazard = splev(precip_regular, splines_tuple)
qpositive = 1.0 - np.exp(-spline_hazard)
cdf_spline_with_zeros = fraction_zero_gefsv12[jspline, ispline] + \
    (1.-fraction_zero_gefsv12[jspline, ispline])*qpositive

# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('Spline-fitted GEFSv12 reforecast and analyzed CDFs,\n'+\
    cmonth+', '+cleadb+'-'+clead+' h forecast, '+\
    clon+' W '+clat+' N',fontsize=13)

ax.plot(precip_regular, cdf_spline_with_zeros,\
    color='Red',lw=2,label='Forecast')
ax.plot(quantile_mapped_precip_hires, fraction_zero_ndfd[jndfd, indfd] +\
    (1.-fraction_zero_ndfd[jndfd, indfd])*empirical_cdf_hires,\
    color='Blue',lw=2,label='Analyzed')

plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.4,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,75)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_'+cmonth+'_lead'+clead+'_lon'+clon+'_lat'+clat+'.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 


# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(4,4.2))#
ax = f.add_axes([.14,.14,.82,.75])
ax.set_title('Jul-Aug-Sep fitted CDFs',fontsize=13)

ax.plot(precip_regular, cdf_spline_with_zeros,\
    color='Red',lw=2,label='Forecast')
ax.plot(quantile_mapped_precip_hires, fraction_zero_ndfd[jndfd, indfd] +\
    (1.-fraction_zero_ndfd[jndfd, indfd])*empirical_cdf_hires,\
    color='Blue',lw=2,label='Analyzed')

plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.4,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_'+cmonth+'_lead'+clead+'_lon'+clon+'_lat'+clat+'_ITH.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 
            