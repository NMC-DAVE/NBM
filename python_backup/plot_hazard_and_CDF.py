""" python plot_hazard_and_CDF.py cyyyymmddhh clead cjy cix"""
import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
from netCDF4 import Dataset
import _pickle as cPickle
import scipy.stats as stats
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

# =====================================================================

def read_forecast_spline_info(cyyyymmddhh, ccmonth, master_directory_forecast_spline):

    """ load spline information for forecast netCDF file.  
        For initial times other than 00 UTC, we need to find the
        appropriate 00 UTC spline file to read in.   Assuming that
        the CDF characteristics are primarily a function of the 
        diurnal cycle, make this adjustment."""

    ccycle = cyyyymmddhh[8:10]
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
    ny_gefsv12, nx_gefsv12 = np.shape(usegamma_gefsv12)
    nc.close()
    lons_fcst_2d, lats_fcst_2d = \
        np.meshgrid(lons_spline_gefsv12_1d,lats_spline_gefsv12_1d)

    return lons_spline_gefsv12_1d, lats_spline_gefsv12_1d, \
        spline_info_gefsv12, fraction_zero_gefsv12, \
        usegamma_gefsv12, ny_gefsv12, nx_gefsv12, \
        lons_fcst_2d, lats_fcst_2d, clead_use

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
        
    infile = master_directory_panal_spline + cmonth+\
        '_conus_CCPA_spline_info_h' + cleada + 'UTC.nc' 
    print ('reading from ', infile)
    nc = Dataset(infile)
    spline_info_inv = nc.variables['spline_info_inv'][:,:,:,:]
    spline_info = nc.variables['spline_info_inv'][:,:,:,:]
    fraction_zero_ndfd = nc.variables['fzero'][:,:]
    usegamma_ndfd = nc.variables['usegamma'][:,:]
    lons_ndfd = nc.variables['lons'][:,:]
    lons_ndfd = lons_ndfd
    lats_ndfd = nc.variables['lats'][:,:]
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
    nc.close() 
    
    return spline_info_inv, spline_info, fraction_zero_ndfd, \
        usegamma_ndfd, lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd

# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # the initial year/month/day/hour
clead = sys.argv[2]  # lead time as 3-digit number, e.g., 
                     # forecast ending at +18h is 018.
cjy  = sys.argv[3]
cix = sys.argv[4]
jy = int(cjy)
ix = int(cix)

# ---- various initialization and carving out month information.

cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cmembers = ['c00','p01', 'p02','p03','p04',\
    'p05','p06','p07','p08','p09','p10',\
    'p11', 'p12','p13','p14','p15','p16','p17','p18','p19','p20',\
    'p21', 'p22','p23','p24','p25','p26','p27','p28','p29','p30']
nmembers = len(cmembers)
master_directory_panal_spline = '/Volumes/NBM/conus_panal/CDF_spline/'
master_directory_forecast_spline = '/Volumes/NBM/conus_gefsv12/CDF_spline/'
master_directory_probability_output = '/Volumes/NBM/conus_gefsv12/probabilities/'
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
nthresholds = len(thresholds)

# ---- read forecast spline parameters from netCDF file

now = datetime.now()
time = now.strftime("%H:%M:%S")
print ('starting quantile mapping ', time)
lons_1d_realtime, lats_1d_realtime, \
    spline_info_gefsv12, fraction_zero_gefsv12, \
    usegamma_gefsv12, ny_gefsv12, nx_gefsv12, \
    lons_fcst_2d, lats_fcst_2d, clead_use = \
    read_forecast_spline_info(cyyyymmddhh, \
    ccmonth, master_directory_forecast_spline)
    
# ---- read precipitation analysis inverse CDF 
#      spline parameters from netCDF file    

spline_info_inv, spline_info, fraction_zero_ndfd, usegamma_ndfd, \
    lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd = \
    read_analysis_spline_info (clead_use, \
    master_directory_panal_spline, ccmonth)

# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
title = 'Spline fit of hazard function, jy, ix = '+str(jy)+','+str(ix)
ax.set_title(title,fontsize=13)
ax.plot(precip_ens_nonzero, spline_hazard,color='Blue',\
    lw=2,label='Cubic spline')
ax.plot(precip_ens_nonzero, hazard_function_empirical,\
    color='Red',lw=2,label='Empirical')
plt.ylabel('Hazard function',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(-0.2,6.0)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'hazard_function_precip_example_spline'+str(jy)+'.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()


f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
title = 'Spline fit of CDF via hazard function, jy, ix = '+\
    str(jy)+','+str(ix)
ax.set_title(title,fontsize=13)
ax.plot(precip_ens_nonzero, spline_cdf,color='Blue',\
    lw=2,label='Cubic spline')
ax.plot(precip_ens_nonzero, empirical_cdf,color='Red',\
    lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_spline'+str(jy)+'.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()