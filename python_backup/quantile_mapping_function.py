"""
quantile_mapping_function.py cyyyymmddhh clead clon clat

coded by: Tom Hamill Feb 2021, tom.hamill@noaa.gov 

"""

import os, sys
from datetime import datetime # Jeff Whitaker's datetime library
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from dateutils import daterange, dateshift
import pygrib # grib-reading routine
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
#from qmapping_spline_onept import qmapping_spline_onept # from f2py of fortran90 code.
from qmapping_spline_1d import qmapping_spline_1d
import scipy.stats as stats
from mpl_toolkits.basemap import Basemap, interp
import _pickle as cPickle
from PIL import Image
    
# =====================================================================

def find_nearest(vec, value):

    """ given a vector vec and a particular value, find the index in vec
        that is nearest to value"""

    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

def read_forecast_spline_info(cyyyymmddhh, ccmonth, \
    master_directory_forecast_spline):

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
    quantile_98 = nc.variables['quantile_98'][:,:]
    ny_gefsv12, nx_gefsv12 = np.shape(usegamma_gefsv12)
    nc.close()
    lons_fcst_2d, lats_fcst_2d = \
        np.meshgrid(lons_spline_gefsv12_1d,lats_spline_gefsv12_1d)
        
    
    return lons_spline_gefsv12_1d, lats_spline_gefsv12_1d, \
        spline_info_gefsv12, fraction_zero_gefsv12, \
        usegamma_gefsv12, quantile_98, ny_gefsv12, nx_gefsv12, \
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
        
    infile = master_directory_panal_spline + cmonth + \
        '_conus_CCPA_spline_info_h' + cleada + 'UTC.nc'
    #infile = master_directory_panal_spline + cmonth+\
    #    '_conus_MSWEP_spline_info_h' + cleada + 'UTC.nc' 
    print ('reading from ', infile)
    nc = Dataset(infile)
    spline_info_inv = nc.variables['spline_info_inv'][:,:,:,:]
    fraction_zero_ndfd = nc.variables['fzero'][:,:]
    usegamma_ndfd = nc.variables['usegamma'][:,:]
    lons_ndfd = nc.variables['lons'][:,:]
    lons_ndfd = lons_ndfd
    lats_ndfd = nc.variables['lats'][:,:]
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
    nc.close() 
    
    return spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
        lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd

    
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

def get_quantile_gefsv12(pamt, jy, ix, \
    spline_info_gefsv12, fraction_zero_gefsv12,\
    usegamma_gefsv12):

    """ this gets the quantile associated with a given precipitation 
    amount for GEFSv12 data, this month and 6-hour period. """
    
    offset_out = 0.0
    if pamt == 0.0:
        
        # ---- arbitrarily assign the CDF to zero if precip is zero.
        
	    quantile = 0.0
    else:	
	
        if usegamma_gefsv12[jy,ix] == 0:
            
	        # ---- flagged as a wet-enough point to estimate the CDF with 
            #      the spline fit to a hazard function. 
            
            splines_tuple = (spline_info_gefsv12[jy,ix,0,:], \
                spline_info_gefsv12[jy,ix,1,:], 3)
            spline_hazard = splev(pamt, splines_tuple)
            qpositive = 1.0 - np.exp(-spline_hazard)
            quantile = fraction_zero_gefsv12[jy,ix] + \
                (1.0 - fraction_zero_gefsv12[jy,ix])*qpositive
        else:
            
            if usegamma_gefsv12[jy,ix] == -1:
                # --- flagged as basically no training data.
                quantile = 0.0
            else:  # --- flagged as minimal training data - use Gamma
                alpha_hat = spline_info_fcst[jy,ix,0,0] 
                beta_hat = spline_info_fcst[jy,ix,1,0] 
                y0 = pamt / beta_hat
                qpositive = stats.gamma.cdf(y0, alpha_hat)
                quantile = fraction_zero_gefsv12[jy,ix] + \
                    (1.0 - fraction_zero_gefsv12[jy,ix])*qpositive 
                    
            
    if quantile >= 1.0: quantile = 0.9999
    if quantile <= 0.0: quantile = 0.0

    return quantile
    
# =====================================================================
    
def get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead):

    """ read in the global forecast.  Subset to CONUS. """

    nib = 886  # precomputed for set of lat-lons that surrounds
    nie = nib+318 # the NDFD CONUS 2.5-km Lambert conformal domain.
    njb = 131 
    nje = njb+153 
    
    # ---- read in forecast grid covering the whole globe.
    
    cycle = cyyyymmddhh[8:10]
    input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
    infile = input_directory + cyyyymmddhh + \
        '_ge'+cmem+'.t'+cycle+'z.pgrb2s.0p25.f' + clead 
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

cyyyymmddhh_begin = sys.argv[1] # the initial year/month/day/hour
clead = sys.argv[2]  # lead time as 3-digit number, e.g., 
                     # forecast ending at +18h is 018.
clon = sys.argv[3]
clat = sys.argv[4]
rlon_input = float(clon)
rlat_input = float(clat)

use98 = True # apply fixed offset beyond 98th percentile

# ---- various initialization and carving out month information.

iyear = int(cyyyymmddhh_begin[0:4])
cmonth = cyyyymmddhh_begin[4:6]
cyyyy = cyyyymmddhh_begin[0:4]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
master_directory_panal_spline = '/Volumes/NBM/conus_panal/CDF_spline/'
master_directory_forecast_spline = '/Volumes/NBM/conus_gefsv12/CDF_spline/'
master_directory_forecast_qmapped = '/Volumes/NBM/conus_gefsv12/qmapped/'

# ---- read forecast spline parameters from netCDF file

print ('reading forecast spline parameters')
lons_1d_realtime, lats_1d_realtime, \
    spline_info_gefsv12, fraction_zero_gefsv12, \
    usegamma_gefsv12, quantile_98, ny_gefsv12, nx_gefsv12, \
    lons_fcst_2d, lats_fcst_2d, clead_use = \
    read_forecast_spline_info(cyyyymmddhh_begin, \
    ccmonth, master_directory_forecast_spline)
print ('min, max lons_1d_realtime = ', \
    np.min(lons_1d_realtime), np.max(lons_1d_realtime))
if np.max(lons_1d_realtime) > 180.:
      lons_1d_realtime = lons_1d_realtime - 360.
    
print ('max, min quantile_98 = ', \
    np.max(quantile_98), np.min(quantile_98))
ixgefs = find_nearest(lons_1d_realtime, rlon_input)
jygefs = find_nearest(lats_1d_realtime, rlat_input)
print ('ixgefs, jygefs = ', ixgefs, jygefs)
print ('lons_1d_realtime = ', lons_1d_realtime[ixgefs])
print ('lats_1d_realtime = ', lats_1d_realtime[jygefs])
print ('fraction_zero_gefs = ', fraction_zero_gefsv12[jygefs,ixgefs])
print ('quantile_98[jygefs,ixgefs] = ', quantile_98[jygefs,ixgefs])
#sys.exit()


# ---- read precipitation analysis inverse CDF 
#      spline parameters from netCDF file    

print ('reading analyzed spline parameters')
spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
    lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd = \
    read_analysis_spline_info (clead_use, \
    master_directory_panal_spline, ccmonth)
if np.max(lons_ndfd) > 180.: lons_ndfd = lons_ndfd - 180.

    
# ---- find index of nearest ndfd point

distx = np.abs(lons_ndfd - rlon_input)
disty = np.abs(lats_ndfd - rlat_input)
dist = distx**2 + disty**2
i = np.argmin(dist)
print (i)
jndfd, indfd = np.unravel_index(i, shape=(ny_ndfd, nx_ndfd), order='C')    
print ('jndfd, indfd, ny_ndfd, nx_ndfd, lons_ndfd, lats_ndfd = ',\
    jndfd, indfd, ny_ndfd, nx_ndfd, lons_ndfd[jndfd,indfd], lats_ndfd[jndfd,indfd])


print ('spline_info_inv[jndfd, indfd, 0,:] = ', spline_info_inv[jndfd, indfd, 0,:])
print ('spline_info_inv[jndfd, indfd, 1,:] = ', spline_info_inv[jndfd, indfd, 1,:])

# ---- process range of precip amts and get quantile and then 
#      quantile-mapped value    
    
precip_amts_forecast = np.arange(401, dtype=np.float64)/4.
percentiles = np.zeros((401), dtype=np.float64)
quantile_mapped_amounts = np.zeros((401), dtype=np.float64)

# ---- get the percentile associated with the 98th percentile
#      of positive precipitation values

percentile_at_q98plus = get_quantile_gefsv12(\
    quantile_98[jygefs,ixgefs], \
    jygefs, ixgefs, spline_info_gefsv12, \
    fraction_zero_gefsv12, usegamma_gefsv12)

# ---- get the quantile mapped value associated with 
#      this percentile    

sinv = spline_info_inv[jndfd,indfd,:,:]
usegamma = usegamma_ndfd[jndfd,indfd]
fzero = fraction_zero_ndfd[jndfd,indfd]
qmapped_value_at_quantile_98 = \
    qmapping_spline_1d(percentile_at_q98plus, \
    quantile_98[jygefs,ixgefs], sinv, fzero, usegamma)
    
for iamt, pamt in enumerate(precip_amts_forecast):
    
    # ---- get the percentile in the overall (including zeros)
    #      distribution associated with this forecast
    #      precipitation amount
    
    percentiles[iamt] = \
        get_quantile_gefsv12(pamt, jygefs, ixgefs, \
        spline_info_gefsv12, fraction_zero_gefsv12, \
        usegamma_gefsv12)   

    # ---- using spline function that maps back from a 
    #      percentile to quantile-mapped analyzed 
    #      precipitation amount.  Note that pamt 
    #      passed in to handle case of input as zero
    #      precip.
    
    quantile_mapped_precip = qmapping_spline_1d \
        (percentiles[iamt], pamt, sinv, fzero, usegamma)
    print ('%ile, pamt, qmp = ', percentiles[iamt], \
        pamt, quantile_mapped_precip)
    
    #quantile_mapped_amounts[iamt] = \
    #    quantile_mapped_precip
        
    #if pamt > quantile_98[jygefs,ixgefs]:
    #    #offset = pamt - quantile_98[jygefs,ixgefs]
    #    offset = 0.0 #pamt - quantile_98[jygefs,ixgefs]
    #    quantile_mapped_amounts[iamt] = \
    #        qmapped_value_at_quantile_98 + offset
    #else:
    #    quantile_mapped_amounts[iamt] = \
    #        quantile_mapped_precip
    quantile_mapped_amounts[iamt] = quantile_mapped_precip
    
    

f = plt.figure(figsize=(6.5,4.5))#
ax = f.add_axes([.15,.14,.8,.77])
ax.set_title('CDFs',fontsize=11)
ax.plot(precip_amts_forecast,percentiles,color='Blue',lw=2,label='Forecast')
ax.plot(quantile_mapped_amounts,percentiles,color='Red',lw=2,label='Analyzed')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.7,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_use98.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()


f = plt.figure(figsize=(6.5,6.7))#
ax = f.add_axes([.15,.14,.8,.77])
ax.set_title('Quantile mapping function',fontsize=11)
ax.plot(precip_amts_forecast,quantile_mapped_amounts,color='Blue',lw=2)
ax.plot([0,40],[0,40],lw=0.5,color='Gray')
plt.ylabel('6-hourly analyzed total precip. (mm)',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0,40)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
ax.set_xlabel('6-hourly forecast precip. (mm)',fontsize=11)
figname = 'CDF_mapping_function.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()


print ('DONE!') 