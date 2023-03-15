"""
quantile_map_ensemble_bymonth_andlead_v3.py cyyyymmddhh clead 

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
from dateutils import daterange, dateshift
import pygrib # grib-reading routine
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
from qmapping_spline import qmapping_spline # from f2py of fortran90 code.
import scipy.stats as stats
from mpl_toolkits.basemap import Basemap, interp
import _pickle as cPickle
from PIL import Image
    
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
    
    
    print ('read_analysis_spline_info: ')
    print ('spline_info_inv[809,923,:,:] = ', spline_info_inv[809,923,:,:])
    print ('fraction_zero_ndfd = ', fraction_zero_ndfd[809,923])
    print ('usegamma_ndfd = ', usegamma_ndfd[809,923])
    
    return spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
        lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd

# =====================================================================

def initialize_netCDF_probfiles(outfile, ncname, ny, nx, \
    lats_in, lons_in, nthresholds, thresholds):
        
    """ initialize the netCDF files for writing probability output"""

    print ('initializing netCDF file ',outfile)
    if ncname == 'nc_fullfield':
        nc_fullfield = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
        ncout = nc_fullfield
    elif ncname == 'nc_upscaled':
        nc_upscaled = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
        ncout = nc_upscaled
    elif ncname == 'nc_thinned':
        nc_thinned = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
        ncout = nc_thinned
    else:
        print ('invalid ncname in initalize_netCDF_probfiles = ', ncname)
        sys.exit()
        
    xf = ncout.createDimension('xf',nx)
    xvf = ncout.createVariable('xf','f4',('xf',))
    xvf.long_name = "eastward grid point number on 1/4-degree lat-lon grid"
    xvf.units = "n/a"

    yf = ncout.createDimension('yf',ny)
    yvf = ncout.createVariable('yf','f4',('yf',))
    yvf.long_name = "northward grid point number on 1/4-degree lat-lon grid"
    yvf.units = "n/a"

    sample = ncout.createDimension('sample',None)
    samplev = ncout.createVariable('samplev','i4',('sample',))
    samplev.units = "n/a"

    lons_out = ncout.createVariable('lons','f4',('yf','xf',))
    lons_out.long_name = "longitude"
    lons_out.units = "degrees_east"

    lats_out = ncout.createVariable('lats','f4',('yf','xf',))
    lats_out.long_name = "latitude"
    lats_out.units = "degrees_north"
    
    threshold = ncout.createDimension('threshold',nthresholds)
    thresholdv = ncout.createVariable('thresholdv','i4',('threshold',))
    thresholdv.units = "mm"

    yyyymmddhh_init = ncout.createVariable('yyyymmddhh_init','i4',('sample',))
    yyyymmddhh_init.longname = "Initial condition date/time in yyyymmddhh format"

    yyyymmddhh_fcst = ncout.createVariable('yyyymmddhh_fcst','i4',('sample',))
    yyyymmddhh_fcst.longname = "Forecast valid date/time in yyyymmddhh format"

    # --- declare the single-level variable information on lat-lon grid

    probability_raw = ncout.createVariable('probability_raw','f4',\
        ('sample','threshold','yf','xf',),zlib=True,least_significant_digit=3)
    probability_raw.units = "n/a"
    probability_raw.long_name = 'Probability of precipitation '+\
        'exceeding threshold amount in raw ens.'
    probability_raw.valid_range = [0.,1.]
    probability_raw.missing_value = np.array(-99.99,dtype=np.float32)
    
    probability_qmapped = ncout.createVariable('probability_qmapped',\
        'f4',('sample','threshold','yf','xf',),
        zlib=True,least_significant_digit=3)
    probability_qmapped.units = "n/a"
    probability_qmapped.long_name = 'Probability of precipitation '+\
        'exceeding threshold amount in qmapped ensemble'
    probability_qmapped.valid_range = [0.,1.]
    probability_qmapped.missing_value = np.array(-99.99,dtype=np.float32)
    
    
    ensmean_raw = ncout.createVariable('ensmean_raw','f4',\
        ('sample','yf','xf',),zlib=True,least_significant_digit=3)
    ensmean_raw.units = "mm"
    ensmean_raw.long_name = 'Raw ensemble-mean precicipitation amount (mm)'
    ensmean_raw.valid_range = [0.,200.]
    ensmean_raw.missing_value = np.array(-99.99,dtype=np.float32)
    
    ensmean_qmapped = ncout.createVariable('ensmean_qmapped',\
        'f4',('sample','yf','xf',),
        zlib=True,least_significant_digit=3)
    ensmean_qmapped.units = "mm"
    ensmean_qmapped.long_name = 'Quantile-mapped ensemble-mean '+\
        'precicipitation amount (mm)'
    ensmean_qmapped.valid_range = [0.,200.]
    ensmean_qmapped.missing_value = np.array(-99.99,dtype=np.float32)
    
    # ---- metadata

    ncout.title = 'NDFD CONUS gridded probabilities of precipitation exceeding '+\
        'various thresholds, raw GEFSv12 and quantile mapped'
    ncout.history = "GEFSv12 implemented at NCEP/EMC Sep 2020"
    ncout.institution =  "NCEP/EMC and PSL"
    ncout.platform = ""
    ncout.references = ""

    # ---- initialize

    xvf[:] = np.arange(nx)
    yvf[:] = np.arange(ny)
    lons_out[:] = lons_in[:,:]
    lats_out[:] = lats_in[:,:]
    thresholdv[:] = thresholds[:]
    
    if ncname == 'nc_fullfield':
        nc_fullfield.close()
    elif ncname == 'nc_upscaled':
        nc_upscaled.close() 
    elif ncname == 'nc_thinned':
        nc_thinned.close() 
    
    istat = 0
    print ('returning')
    return istat
    
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
    
def write_probabilities_to_netcdf (outfile, ncname, idate, cyyyymmddhh, \
    cyyymmddhh_fcst, prob_raw_out, prob_qmapped_out, ensmean_raw_out, \
    ensmean_qmapped_out):    
    
    """ write the raw and quantile-mapped ensemble probabilities to 
    a netCDF file """
    
    # ---- write probabilities and close

    print ('writing probs to netCDF file ',outfile)
    if ncname == 'nc_fullfield':
        nc_fullfield = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        nc_fullfield['yyyymmddhh_init'][idate] = int(cyyyymmddhh)
        nc_fullfield['yyyymmddhh_fcst'][idate] = int(cyyyymmddhh_fcst)
        nc_fullfield['probability_raw'][idate] = prob_raw_out[:,:,:]
        nc_fullfield['probability_qmapped'][idate] = prob_qmapped_out[:,:,:]
        nc_fullfield['ensmean_raw'][idate] = ensmean_raw_out[:,:]
        nc_fullfield['ensmean_qmapped'][idate] = ensmean_qmapped_out[:,:]
        nc_fullfield.close()
    elif ncname == 'nc_upscaled':
        nc_upscaled = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        nc_upscaled['yyyymmddhh_init'][idate] = int(cyyyymmddhh)
        nc_upscaled['yyyymmddhh_fcst'][idate] = int(cyyyymmddhh_fcst)
        nc_upscaled['probability_raw'][idate] = prob_raw_out[:,:,:]
        nc_upscaled['probability_qmapped'][idate] = prob_qmapped_out[:,:,:]
        nc_upscaled['ensmean_raw'][idate] = ensmean_raw_out[:,:]
        nc_upscaled['ensmean_qmapped'][idate] = ensmean_qmapped_out[:,:]
        nc_upscaled.close()
    elif ncname == 'nc_thinned':
        nc_thinned = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        nc_thinned['yyyymmddhh_init'][idate] = int(cyyyymmddhh)
        nc_thinned['yyyymmddhh_fcst'][idate] = int(cyyyymmddhh_fcst)
        nc_thinned['probability_raw'][idate] = prob_raw_out[:,:,:]
        nc_thinned['probability_qmapped'][idate] = prob_qmapped_out[:,:,:]
        nc_thinned['ensmean_raw'][idate] = ensmean_raw_out[:,:]
        nc_thinned['ensmean_qmapped'][idate] = ensmean_qmapped_out[:,:]
        nc_thinned.close()
    else:
        print ('invalid ncname in initalize_netCDF_probfiles = ', ncname)
        sys.exit()
        

    #if ncname == 'nc_fullfield':
    #    nc_fullfield.close()
    #elif ncname == 'nc_upscaled':
    #    nc_upscaled.close() 
    #elif ncname == 'nc_thinned':
    #    nc_thinned.close()
    istat = 0
    return istat
    
# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh_begin = sys.argv[1] # the initial year/month/day/hour
clead = sys.argv[2]  # lead time as 3-digit number, e.g., 
                     # forecast ending at +18h is 018.

# ---- various initialization and carving out month information.

iyear = int(cyyyymmddhh_begin[0:4])
cmonth = cyyyymmddhh_begin[4:6]
cyyyy = cyyyymmddhh_begin[0:4]
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
master_directory_forecast_qmapped = '/Volumes/NBM/conus_gefsv12/qmapped/'
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'
master_directory_upscaled_output = '/Volumes/NBM/conus_gefsv12/upscaled/'
master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/thinned/'
master_directory_forecast_spline = '/Volumes/NBM/conus_gefsv12/CDF_spline/'
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
nthresholds = len(thresholds)
firstdate = True
output_qmapped = True
netcdf_fullfield = True
netcdf_thinned = True
netcdf_upscaled = True
nstride = 10

if iyear%4 == 0:
    ndays = ndaysomo_leap[imonth]
else:
    ndays = ndaysomo[imonth]

#cyyyymmddhh_end = dateshift(cyyyymmddhh_begin,24*ndays)
cyyyymmddhh_end = dateshift(cyyyymmddhh_begin,0)

date_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 24)
   
# ---- loop over dates to process

for idate, cyyyymmddhh in enumerate(date_list):
    
    cyyyymmddhh_fcst = dateshift(cyyyymmddhh, int(clead))
    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print ('starting quantile mapping for date = ', \
        cyyyymmddhh,' at time = ',time)
    if firstdate == True:
    
        # ---- read forecast spline parameters from netCDF file
    
        print ('reading forecast spline parameters')
        lons_1d_realtime, lats_1d_realtime, \
            spline_info_gefsv12, fraction_zero_gefsv12, \
            usegamma_gefsv12, ny_gefsv12, nx_gefsv12, \
            lons_fcst_2d, lats_fcst_2d, clead_use = \
            read_forecast_spline_info(cyyyymmddhh, \
            ccmonth, master_directory_forecast_spline)
    
        # ---- read precipitation analysis inverse CDF 
        #      spline parameters from netCDF file    

        print ('reading analyzed spline parameters')
        spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
            lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd = \
            read_analysis_spline_info (clead_use, \
            master_directory_panal_spline, ccmonth)
            
        # ---- flip and then interpolate the GEFSv12  
        #      climatological fraction zerosand latitude 
        #      array to NDFD grid. Flipping is necessary as
        #      Basemap.interp requires input lats and data 
        #      be oriented S to N.

        print ('interpolating fraction_zero of GEFSv12 to NDFD grid. ')  
        fraction_zero_gefsv12_flipped = np.flipud(fraction_zero_gefsv12)
        lats_1d_realtime_flipud = np.flipud(lats_1d_realtime)
        fraction_zero_gefsv12_on_ndfd = interp(fraction_zero_gefsv12_flipped, \
            lons_1d_realtime, lats_1d_realtime_flipud, \
            lons_ndfd, lats_ndfd, checkbounds=False, \
            masked=False, order=1)
            
        # ---- initialize netCDF file for full probabilistic fields
        
        if netcdf_fullfield == True:
            outfile = master_directory_fullfield_output + ccmonth + cyyyy + \
                '_lead'+clead+'_probabilities_fullfield.nc'
            ncname = 'nc_fullfield'
            istat = initialize_netCDF_probfiles(outfile, ncname, \
                ny_ndfd, nx_ndfd, lats_ndfd, lons_ndfd, nthresholds, \
                thresholds)
            
        # ---- set up full output grids and working arrays

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
            
        # ---- Upscale lat/lon array 10-fold.  Determine the grid 
        #      size of an upscaled field and pass these in to initialize
        #      netCDF file for the upscaled
        
        if netcdf_upscaled == True:
            print ('performing upscaling of lat/lon arrays, initializing work arrays.')
            im = Image.fromarray(lons_ndfd)
            imu = im.resize(size=(ny_ndfd//10, nx_ndfd//10),resample=Image.BOX)
            lons_upscaled = np.transpose(np.asarray(imu))
            im = Image.fromarray(lats_ndfd)
            imu = im.resize(size=(ny_ndfd//10, nx_ndfd//10),resample=Image.BOX)
            lats_upscaled = np.transpose(np.asarray(imu))
            ny_upscaled, nx_upscaled = np.shape(lons_upscaled)
            
            precip_ens_raw_upscaled = np.zeros((nmembers, ny_upscaled, nx_upscaled), \
                dtype=np.float32)
            precip_ens_qmapped_upscaled = np.zeros((nmembers, ny_upscaled, nx_upscaled), \
                dtype=np.float32)
            prob_raw_upscaled = np.zeros((nthresholds, ny_upscaled, nx_upscaled), \
                dtype=np.float32)
            prob_qmapped_upscaled = np.zeros((nthresholds, ny_upscaled, nx_upscaled), \
                dtype=np.float32)
            prob_upscaled = np.zeros((ny_upscaled, nx_upscaled), dtype=np.float64)
            zeros_upscaled = np.zeros((ny_upscaled, nx_upscaled), dtype=np.float64)
            ones_upscaled = np.ones((ny_upscaled, nx_upscaled), dtype=np.float64)
            
            outfile = master_directory_upscaled_output + ccmonth + cyyyy + \
                '_lead'+clead+'_probabilities_upscaled.nc'
            ncname = 'nc_upscaled'
            istat = initialize_netCDF_probfiles(outfile, ncname, \
                ny_upscaled, nx_upscaled, lats_upscaled, lons_upscaled, \
                nthresholds, thresholds)
        
        # ---- Take every 10th pixel in lat/lon to subsample starting at
        #      the sixth pixel.  Determine the grid size of an upscaled
        #      field and pass these in to initialize netCDF file 
        #      for the thinned
        
        if netcdf_thinned == True:
            lons_thinned = lons_ndfd[5::nstride, 5::nstride]
            lats_thinned = lats_ndfd[5::nstride, 5::nstride]        
            ny_thinned, nx_thinned = np.shape(lons_thinned)
            precip_ens_raw_thinned = \
                np.zeros((nmembers, ny_thinned, nx_thinned), dtype=np.float32)
            precip_ens_qmapped_thinned = \
                np.zeros((nmembers, ny_thinned, nx_thinned), dtype=np.float32)
            outfile = master_directory_thinned_output + ccmonth + cyyyy + \
                '_lead'+clead+'_probabilities_thinned.nc'
            ncname = 'nc_thinned'
            istat = initialize_netCDF_probfiles(outfile, ncname, \
                ny_thinned, nx_thinned, lats_thinned, lons_thinned, \
                nthresholds, thresholds)
                
            prob_raw_thinned = np.zeros((nthresholds, ny_thinned, nx_thinned), \
                dtype=np.float32)
            prob_qmapped_thinned = np.zeros((nthresholds, ny_thinned, nx_thinned), \
                dtype=np.float32)
            prob_thinned = np.zeros((ny_thinned, nx_thinned), dtype=np.float64)
            zeros_thinned = np.zeros((ny_thinned, nx_thinned), dtype=np.float64)
            ones_thinned = np.ones((ny_thinned, nx_thinned), dtype=np.float64)
        
        # ---- having done the upfront work, set a flag so this isn't
        #      repeated as we loop thru other dates
        
        firstdate = False
        
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
        print ('   gefsv12_quantiles_on_ndfd[809,923] = ', \
            gefsv12_quantiles_on_ndfd[809,923])
        
        # ---- apply quantile mapping procedure to this member

        print ('   applying quantile mapping')
        qmapped_precip = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
        a = np.where(fraction_zero_ndfd < 0)
        aunraveled = np.unravel_index(a,shape=(ny_ndfd,nx_ndfd))
        qmapped_precip = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
        qmapped_precip = qmapping_spline(gefsv12_quantiles_on_ndfd, \
    	    precip_gefsv12_on_ndfd, spline_info_inv, \
    	    fraction_zero_ndfd, usegamma_ndfd, ny_ndfd, nx_ndfd)
        precip_ens_qmapped_ndfd[imem,:,:] = qmapped_precip[:,:]   
        print ('   precip_ens_qmapped_ndfd[imem, 809,923] = ', \
             precip_ens_qmapped_ndfd[imem,809,923])     

        # ---- if upscaling to a 10x coarser grid is desired, 
        #      perform that now and save.
        
        if netcdf_upscaled == True:
            print ('   performing upscaling')
            im = Image.fromarray(precip_gefsv12_on_ndfd)
            imu = im.resize(size=(ny_ndfd//10, nx_ndfd//10),resample=Image.BOX)
            precip_ens_raw_upscaled[imem,:,:] = np.transpose(np.asarray(imu))
            im = Image.fromarray(qmapped_precip)
            imu = im.resize(size=(ny_ndfd//10, nx_ndfd//10),resample=Image.BOX)
            precip_ens_qmapped_upscaled[imem,:,:] = np.transpose(np.asarray(imu))
                
        # ---- if subsampling every 10th pixel, 
        #      perform that now and save.
        
        if netcdf_thinned == True:
            print ('   performing thinning')
            precip_ens_raw_thinned[imem,:,:] = \
                 precip_gefsv12_on_ndfd[5::nstride, 5::nstride]
            precip_ens_qmapped_thinned[imem,:,:] = \
                 qmapped_precip[5::nstride, 5::nstride]
               
    # ---- generate probabilities, raw and quantile mapped for each
    #      specified precipitation threshold

    for ithresh, thresh in enumerate(thresholds):
        print ('   generating probabilities for threshold = ', thresh)
    
        if netcdf_fullfield == True:
            print ('   generating probabilities for full field')
            prob = generate_probabilities(nmembers, precip_ens_raw_ndfd, \
                zeros, ones, thresh, prob, ny_ndfd, nx_ndfd)
            prob_raw_out[ithresh,:,:] = prob[:,:]
            
            prob = generate_probabilities(nmembers, precip_ens_qmapped_ndfd, \
                zeros, ones, thresh, prob, ny_ndfd, nx_ndfd)
            prob_qmapped_out[ithresh,:,:] = prob[:,:]
            
            ensmean_raw_ndfd = np.squeeze(np.mean(precip_ens_raw_ndfd, axis=0))
            ensmean_qmapped_ndfd = np.squeeze(np.mean(precip_ens_qmapped_ndfd, axis=0))
    
        if netcdf_upscaled == True:
            print ('   generating probabilities for upscaled')
            prob_upscaled = generate_probabilities(nmembers, precip_ens_raw_upscaled, \
                zeros_upscaled, ones_upscaled, thresh, prob_upscaled, \
                ny_upscaled, nx_upscaled)
            prob_raw_upscaled[ithresh,:,:] = prob_upscaled[:,:]
            
            prob_upscaled = generate_probabilities(nmembers, precip_ens_qmapped_upscaled, \
                zeros_upscaled, ones_upscaled, thresh, prob_upscaled, \
                ny_upscaled, nx_upscaled)
            prob_qmapped_upscaled[ithresh,:,:] = prob_upscaled[:,:]
            
            ensmean_raw_upscaled = np.squeeze(np.mean(precip_ens_raw_upscaled, \
                axis=0))
            ensmean_qmapped_upscaled = np.squeeze(np.mean(\
                precip_ens_qmapped_upscaled, axis=0))
        
        if netcdf_thinned == True:
            print ('   generating probabilities for thinned')
            prob_thinned = generate_probabilities(nmembers, precip_ens_raw_thinned, \
                zeros_thinned, ones_thinned, thresh, prob_thinned, \
                ny_thinned, nx_thinned)
            prob_raw_thinned[ithresh,:,:] = prob_thinned[:,:]        
        
            prob_thinned = generate_probabilities(nmembers, precip_ens_qmapped_thinned, \
                zeros_thinned, ones_thinned, thresh, prob_thinned, \
                ny_upscaled, nx_upscaled)
            prob_qmapped_thinned[ithresh,:,:] = prob_thinned[:,:]
            
            ensmean_raw_thinned = np.squeeze(np.mean(precip_ens_raw_thinned, \
                axis=0))
            ensmean_qmapped_thinned = np.squeeze(np.mean(\
                precip_ens_qmapped_thinned, axis=0))
    
    print('b4 prob precip_ens_qmapped_ndfd[:,809,923] = ', \
        precip_ens_qmapped_ndfd[:,809,923])
        
    # ---- write the probabilities over full set of thresholds to netCDF file.

    if netcdf_fullfield == True:
        now = datetime.now()
        time = now.strftime("%H:%M:%S")
        ncname = 'nc_fullfield'
        outfile = master_directory_fullfield_output + ccmonth + cyyyy + \
            '_lead'+clead+'_probabilities_fullfield.nc'
        print ('   writing full field probabilities to file. ',time)
        istat = write_probabilities_to_netcdf(outfile, ncname, idate, \
            cyyyymmddhh, cyyyymmddhh_fcst, prob_raw_out, \
            prob_qmapped_out, ensmean_raw_ndfd, ensmean_qmapped_ndfd)
        print ('   writing full field probabilities and ens means completed.')

    if netcdf_thinned == True:
        now = datetime.now()
        time = now.strftime("%H:%M:%S")
        ncname = 'nc_thinned'
        outfile = master_directory_thinned_output + ccmonth + cyyyy + \
            '_lead'+clead+'_probabilities_thinned.nc'
        print ('   writing thinned field probabilities to file. ',time)
        istat = write_probabilities_to_netcdf(outfile, ncname, idate, \
            cyyyymmddhh, cyyyymmddhh_fcst, prob_raw_thinned, \
            prob_qmapped_thinned, ensmean_raw_thinned, \
            ensmean_qmapped_thinned)
        print ('   writing thinned field probabilities and ens means completed.')

    if netcdf_upscaled == True:
        now = datetime.now()
        time = now.strftime("%H:%M:%S")
        ncname = 'nc_upscaled'
        outfile = master_directory_upscaled_output + ccmonth + cyyyy + \
            '_lead'+clead+'_probabilities_upscaled.nc'
        print ('   writing upscaled field probabilities to file. ',time)
        istat = write_probabilities_to_netcdf(outfile, ncname, idate, \
            cyyyymmddhh, cyyyymmddhh_fcst, prob_raw_upscaled, \
            prob_qmapped_upscaled, ensmean_raw_upscaled, \
            ensmean_qmapped_upscaled)
        print ('   writing upscaled field probabilities and ens means completed.') 
           
    if output_qmapped == True:
        now = datetime.now()
        time = now.strftime("%H:%M:%S")
        outfile = master_directory_forecast_qmapped + \
            cyyyymmddhh+'_lead='+clead+'.cPick'
        print ('writing raw and qmapped ens to ', outfile) 
        ouf = open(outfile,'wb')
        cPickle.dump(precip_ens_raw_ndfd, ouf)
        cPickle.dump(precip_ens_qmapped_ndfd, ouf)
        print('precip_ens_qmapped_ndfd[:,809,923] = ', precip_ens_qmapped_ndfd[:,809,923])
        cPickle.dump(lons_ndfd, ouf)
        cPickle.dump(lats_ndfd, ouf)
        ouf.close()
        
print ('DONE!') 