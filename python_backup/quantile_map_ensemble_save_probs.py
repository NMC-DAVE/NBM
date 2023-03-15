"""
quantile_map_ensemble_save_probs.py cyyyymmddhh clead 

This python routine reads in pre-generated information on the 
fitted CDF for the forecast, the inverse of the CDF for the
analyzed, the forecast in question.   It then applies
quantile mapping, generates probabilities, and writes these
to a netCDF-formatted output file

coded by: Tom Hamill Feb 2021, tom.hamill@noaa.gov 

"""

import os, sys
from datetime import datetime # Jeff Whitaker's datetime library
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import pygrib # grib-reading routine
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
from qmapping_spline import qmapping_spline # from f2py of fortran90 code.
import scipy.stats as stats
from mpl_toolkits.basemap import Basemap, interp
    
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
        
    infile = master_directory_panal_spline + cmonth+\
        '_conus_CCPA_spline_info_h' + cleada + 'UTC.nc' 
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
    #print (gribfilename, endStep)
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

def get_quantile_gefsv12(precip_amount, jy, ix, \
    spline_info_gefsv12, fraction_zero_gefsv12,\
    usegamma_gefsv12):

    """ this gets the quantile associated with a given precipitation 
    amount for GEFSv12 data, this month and 6-hour period. """
    
    if precip_amount[jy,ix] == 0.0:
        
        # ---- arbitrarily assign the CDF to zero if precip is zero.
        
	    quantile = 0.0
    else:	
	
        if usegamma_gefsv12[jy,ix] == 0:
            
	        # ---- flagged as a wet-enough point to estimate the CDF with 
            #      the spline fit to a hazard function. 
            
            splines_tuple = (spline_info_gefsv12[jy,ix,0,:], \
                spline_info_gefsv12[jy,ix,1,:], 3)
            spline_hazard = splev(precip_amount[jy,ix], splines_tuple)
            spline_cdf = 1.0 - np.exp(-spline_hazard)
            quantile = fraction_zero_gefsv12[jy,ix] + \
                (1.0 - fraction_zero_gefsv12[jy,ix])*spline_cdf
        else:
            
            if usegamma_gefsv12[jy,ix] == -1:
                # --- flagged as basically no training data.
                quantile = 0.0
            else:  # --- flagged as minimal training data - use Gamma
                alpha_hat = spline_info_fcst[jy,ix,0,0] 
                beta_hat = spline_info_fcst[jy,ix,1,0] 
                y0 = precip_amount[jy,ix] / beta_hat
                quantile = stats.gamma.cdf(y0, alpha_hat)            

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

def generate_probabilities(nmembers, precip_ens, \
    zeros, ones, thresh, prob, ny_ndfd, nx_ndfd):
    
    """ determine the ensemble probability of exceeding the threshold """
    
    prob[:,:] = 0.0
    for imem in range(nmembers):
        onezero = np.where(precip_ens[imem,:,:] >= thresh, ones, zeros)
        prob = prob + onezero
    prob = prob / float(nmembers)
    return prob

# =====================================================================    
    
def write_probabilities_to_netcdf(master_directory_probability_output,\
    cyyyymmddhh, clead, ny_ndfd, nx_ndfd,  nthresholds, \
    thresholds, prob_raw_out, prob_qmapped_out):    
    
    """ write the raw and quantile-mapped ensemble probabilities to 
    a netCDF file """
    
    outfile = master_directory_probability_output+\
        cyyyymmddhh+'_'+clead+'h_conus_probs.nc'
    print ('writing to ', outfile)
    ncout = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

    xf = ncout.createDimension('xf',nx_ndfd)
    xvf = ncout.createVariable('xf','f4',('xf',))
    xvf.long_name = "eastward grid point number on NDFD grid"
    xvf.units = "n/a"

    yf = ncout.createDimension('yf',ny_ndfd)
    yvf = ncout.createVariable('yf','f4',('yf',))
    yvf.long_name = "northward grid point number on 1/4-degree lat-lon grid"
    yvf.units = "n/a"

    xthresh = ncout.createDimension('xthresh',nthresholds)
    xvthresh = ncout.createVariable('xthresh','f4',('xthresh',))
    xvthresh.long_name = "threshold values (mm)"
    xvthresh.units = "mm"

    lonsa = ncout.createVariable('lons','f4',('yf','xf',))
    lonsa.long_name = "longitude"
    lonsa.units = "degrees (neg = west)"

    latsa = ncout.createVariable('lats','f4',('yf','xf',))
    latsa.long_name = "latitude"
    latsa.units = "degrees_north"

    prob_raw = ncout.createVariable('prob_raw','f4',('xthresh','yf','xf',),
        zlib=True,least_significant_digit=3)
    prob_raw.units = "n/a"
    prob_raw.long_name = \
        "Raw ensemble probability of 6-h "+\
        "accumulated precipitation exceeding threshold"
    prob_raw.valid_range = [0.,1.]
    prob_raw.missing_value = np.array(-99.99,dtype=np.float32)

    prob_qmapped = ncout.createVariable('prob_qmapped','f4',('xthresh','yf','xf',),
        zlib=True,least_significant_digit=3)
    prob_qmapped.units = "n/a"
    prob_qmapped.long_name = \
        "Quantile-mapped ensemble probability of 6-h "+\
        "accumulated precipitation exceeding threshold"
    prob_qmapped.valid_range = [0.,1.]
    prob_qmapped.missing_value = np.array(-99.99,dtype=np.float32)

    # ---- metadata

    ncout.title = "NDFD CONUS grid probabilities of precipitation exceeding '+\
        'various thresholds, raw GEFSv12 and quantile mapped"
    ncout.history = "GEFSv12 implemented at NCEP/EMC Sep 2020-"
    ncout.institution =  "NCEP/EMC and PSL"
    ncout.platform = ""
    ncout.references = ""

    # ---- initialize dimensions, thresholds, lat/lon

    xvf[:] = np.arange(nx_ndfd)
    yvf[:] = np.arange(ny_ndfd)
    xvthresh[:] = thresholds[:]
    lonsa[:] = lons_ndfd[:,:]
    latsa[:] = lats_ndfd[:,:]

    # ---- write probabilities and close

    prob_raw[:] = prob_raw_out[:,:,:]
    prob_qmapped[:] = prob_qmapped_out[:,:,:]

    ncout.close()  
    istat = 0
    return istat
    
# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # the initial year/month/day/hour
clead = sys.argv[2]  # lead time as 3-digit number, e.g., 
                     # forecast ending at +18h is 018.

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
use98 = False # apply fixed offset beyond 98th percentile
   
# ---- read forecast spline parameters from netCDF file

now = datetime.now()
time = now.strftime("%H:%M:%S")
print ('starting quantile mapping ', time)
lons_1d_realtime, lats_1d_realtime, \
    spline_info_gefsv12, fraction_zero_gefsv12, \
    usegamma_gefsv12, quantile_98, ny_gefsv12, nx_gefsv12, \
    lons_fcst_2d, lats_fcst_2d, clead_use = \
    read_forecast_spline_info(cyyyymmddhh, \
    ccmonth, master_directory_forecast_spline)
    
print ('spline_info_gefsv12[-115,240,:,;] = ', spline_info_gefsv12[-115,240,:,:])
print ('usegamma_gefsv12[-115,240] = ', usegamma_gefsv12[-115,240])
print ('fraction_zero_gefsv12[-115,240] = ', fraction_zero_gefsv12[-115,240])
    
# ---- read precipitation analysis inverse CDF 
#      spline parameters from netCDF file    

spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
    lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd = \
    read_analysis_spline_info (clead_use, \
    master_directory_panal_spline, ccmonth)

print ('spline_info_inv[-115,240,:,;] = ', spline_info_inv[-115,240,:,:])
print ('usegamma_ndfd[-115,240] = ', usegamma_ndfd[-115,240])
print ('fraction_zero_ndfd[-115,240] = ', fraction_zero_ndfd[-115,240])    

# ---- flip and then interpolate the GEFSv12 climatological fraction zeros 
#      and latitude array to NDFD grid. Flipping is necessary as
#      Basemap.interp requires input lats and data be oriented S to N.

print ('interpolating fraction_zero of GEFSv12 to NDFD grid. ')  
fraction_zero_gefsv12_flipped = np.flipud(fraction_zero_gefsv12)
lats_1d_realtime_flipud = np.flipud(lats_1d_realtime)
fraction_zero_gefsv12_on_ndfd = interp(fraction_zero_gefsv12_flipped, \
    lons_1d_realtime, lats_1d_realtime_flipud, \
    lons_ndfd, lats_ndfd, checkbounds=False, \
    masked=False, order=1)

# ---- set up output grids and working arrays

precip_ens_raw_ndfd = np.zeros((nmembers, ny_ndfd, nx_ndfd), \
    dtype=np.float32)
precip_ens_qmapped_ndfd = np.zeros((nmembers, ny_ndfd, nx_ndfd), \
    dtype=np.float32)
prob_raw_out = np.zeros((nthresholds, ny_ndfd, nx_ndfd), \
    dtype=np.float32)
prob_qmapped_out = np.zeros((nthresholds, ny_ndfd, nx_ndfd), \
    dtype=np.float32)
prob = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
zeros = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
ones = np.ones((ny_ndfd, nx_ndfd), dtype=np.float64)

# ---- loop over members, Read in this member of GEFS, 
#      interpolate to NDFD grid and quantile map

for imem, cmem in enumerate(cmembers): 
    
    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print ('------ processing member = ',cmem,time)
    
    # ---- read in & extract the real-time precipitation for this domain.
    #      Set GEFSv12 grid dimensions and 1/4-deg lon and lat.

    precip_realtime, lons_1d_realtime, lats_1d_realtime = \
        get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead) 

    # ---- flip upside down, as subsequent interpolation requires
    #      S to N with increasing j index.  Then
    #      interpolate precipitation to NDFD grid. 

    print ('   interpolating and flipping real-time forecast')
    precip_realtime_flipud = np.flipud(precip_realtime) 
    precip_gefsv12_on_ndfd = interp(precip_realtime_flipud, \
        lons_1d_realtime, lats_1d_realtime_flipud,  \
        lons_ndfd, lats_ndfd, checkbounds=False, \
        masked=False, order=1)  
    precip_ens_raw_ndfd[imem,:,:] = precip_gefsv12_on_ndfd[:,:]
        
    # ---- now loop over grid points and obtain the forecast quantile
    #      associated with this 0.25 degree forecast grid point

    print ('   getting quantiles of forecast')
    gefsv12_quantiles = np.zeros((ny_gefsv12, nx_gefsv12), \
        dtype=np.float64)
    for jy in range(ny_gefsv12):
        for ix in range(nx_gefsv12):
            gefsv12_quantiles[jy,ix] = get_quantile_gefsv12( \
                precip_realtime, jy, ix, spline_info_gefsv12, \
                fraction_zero_gefsv12, usegamma_gefsv12)

    # ---- interpolate the forecast quantile to the NDFD grid.

    print ('   interpolating quantiles to NDFD grid')
    gefsv12_quantiles_flipud = np.flipud(gefsv12_quantiles) 
    gefsv12_quantiles_on_ndfd = interp(gefsv12_quantiles_flipud, \
        lons_1d_realtime, lats_1d_realtime_flipud,  \
        lons_ndfd, lats_ndfd, checkbounds=False, \
        masked=False, order=1)  
        
    # ---- apply quantile mapping procedure to this member

    print ('   applying quantile mapping')
    qmapped_precip = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
    qmapped_precip = qmapping_spline(gefsv12_quantiles_on_ndfd, \
    	precip_gefsv12_on_ndfd, spline_info_inv, \
    	fraction_zero_ndfd, usegamma_ndfd, ny_ndfd, nx_ndfd)
       
    precip_ens_qmapped_ndfd[imem,:,:] = qmapped_precip[:,:]

    
# ---- generate probabilities, raw and quantile mapped for each
#      specified precipitation threshold

for ithresh, thresh in enumerate(thresholds):
    print ('   generating probabilities for threshold = ', thresh)
    prob = generate_probabilities(nmembers, precip_ens_raw_ndfd, \
        zeros, ones, thresh, prob, ny_ndfd, nx_ndfd)
    prob_raw_out[ithresh,:,:] = prob[:,:]
    prob = generate_probabilities(nmembers, precip_ens_qmapped_ndfd, \
        zeros, ones, thresh, prob, ny_ndfd, nx_ndfd)
    prob_qmapped_out[ithresh,:,:] = prob[:,:]

# ---- write the probabilities to netCDF file.

now = datetime.now()
time = now.strftime("%H:%M:%S")
print ('writing probabilities to file. ',time)
istat = write_probabilities_to_netcdf( \
    master_directory_probability_output, \
    cyyyymmddhh, clead, ny_ndfd, nx_ndfd, \
    nthresholds, thresholds, prob_raw_out, \
    prob_qmapped_out)

print ('writing completed.  Done.')



