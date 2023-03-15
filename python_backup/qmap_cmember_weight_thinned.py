"""
qmap_cmember_weight_thinned.py cyyyymmddhh clead 

This python routine reads in pre-generated spline coefficient information 
on the cumulative hazard function of the fitted CDF for the forecast, 
the inverse of the hazard function CDF for the
analyzed, and the forecast in question.   It then applies
quantile mapping, generates probabilities, and applies closest-member
histogram weighting to the output on the thinned grid for this 
date and all other dates in the same month. It writes these
to a netCDF-formatted output file.   Because probabilities on the
high-resolution grid are too large, it permits saving every 10th point
or upscaling by a factor of 10.


coded by: Tom Hamill, Aug 2021, tom.hamill@noaa.gov 

"""

import os, sys
from datetime import datetime # Jeff Whitaker's datetime library
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from dateutils import daterange, dateshift
import pygrib # grib-reading routine
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
from qmapping_spline_fivexfive_thinned import \
    qmapping_spline_fivexfive_thinned # from f2py  
from weighted_probability_f90 import weighted_probability_f90 # from f2py 
from weighted_probability_dressed_f90 import weighted_probability_dressed_f90 # from f2py 
import scipy.stats as stats
from mpl_toolkits.basemap import Basemap, interp
import _pickle as cPickle
from PIL import Image
    
# =====================================================================

def read_forecast_spline_info(cyyyymmddhh, ccmonth, \
    master_directory_forecast_spline):

    """ load spline information for forecast netCDF file.  
        For initial times other than 00 UTC, we need to find the
        appropriate spline file to read in based on 00 UTC 
        forecast initialization times.   Assuming that
        the CDF characteristics are primarily a function of the 
        diurnal cycle, make this adjustment. """

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
    quantile_99 = nc.variables['quantile_99'][:,:]
    number_knots_gefsv12 = nc.variables['number_knots'][:,:]  
    ny_gefsv12, nx_gefsv12 = np.shape(usegamma_gefsv12)
    print ('forecast usegamma, total number of points = ', \
        np.sum(usegamma_gefsv12), ny_gefsv12*nx_gefsv12 )
    nc.close()
    lons_fcst_2d, lats_fcst_2d = \
        np.meshgrid(lons_spline_gefsv12_1d,lats_spline_gefsv12_1d)
        
    
    return lons_spline_gefsv12_1d, lats_spline_gefsv12_1d, \
        spline_info_gefsv12, fraction_zero_gefsv12, \
        usegamma_gefsv12, quantile_99, number_knots_gefsv12, \
        ny_gefsv12, nx_gefsv12, lons_fcst_2d, lats_fcst_2d, clead_use

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
    usegamma_ndfd = nc.variables['usegamma'][:,:]
    number_knots_ndfd = nc.variables['number_knots'][:,:] 
    lons_ndfd = nc.variables['lons'][:,:]
    lons_ndfd = lons_ndfd
    lats_ndfd = nc.variables['lats'][:,:]
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
    nc.close() 
    
    return spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
        number_knots_ndfd, lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd

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
    
    probability_qmapped = ncout.createVariable('probability_qmapped',\
        'f4',('sample','threshold','yf','xf',),
        zlib=True,least_significant_digit=3)
    probability_qmapped.units = "n/a"
    probability_qmapped.long_name = 'Probability of precipitation '+\
        'exceeding threshold amount in qmapped ensemble'
    probability_qmapped.valid_range = [0.,1.]
    probability_qmapped.missing_value = np.array(-99.99,dtype=np.float32)
    
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
        'various thresholds, weighted and quantile mapped'
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

def get_quantile_gefsv12(precip_amount, precip_gefsv12_on_ndfd, \
    spline_info_gefsv12, ny_gefsv12, nx_gefsv12, \
    fraction_zero_gefsv12, fraction_zero_gefsv12_on_ndfd, \
    usegamma_gefsv12, quantile_99_ndfd, use99, \
    number_knots_gefsv12, lons_1d_realtime, \
    lats_1d_realtime_flipud, lons_ndfd, lats_ndfd, jnd, ind):    

    """ this gets the quantile associated with a given precipitation 
    amount for GEFSv12 data, this month and 6-hour period.  This is
    the quantile in the overall distribution, including zeros. """
    
    # ---- first loop thru the GEFSv12 grid points and determine the cumulative
    #      percentile of the forecast in the distribution of positive amounts.
    
    #print ('   precip_gefsv12_on_ndfd[jnd,ind] = ', precip_gefsv12_on_ndfd[jnd,ind])
    gefsv12_quantiles = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.float64)
    print ('   finding quantiles in + distribution')
    for jy in range(ny_gefsv12):
        for ix in range(nx_gefsv12):
            offset_out = 0.0
            if precip_amount[jy,ix] == 0.0:
        
                # ---- arbitrarily assign the CDF to zero if precip is zero.
        
                quantile = 0.0
                qpositive = 0.0
            else:	
        
                if usegamma_gefsv12[jy,ix] == 0: 
            
        	        # ---- flagged as a wet-enough point to estimate the CDF with 
                    #      the spline fit to a hazard function. 
            
                    nk = number_knots_gefsv12[jy,ix]
                    splines_tuple = (spline_info_gefsv12[jy,ix,0,0:nk], \
                        spline_info_gefsv12[jy,ix,1,0:nk], 3)
                    spline_hazard = splev(precip_amount[jy,ix], splines_tuple)
                    qpositive = 1.0 - np.exp(-spline_hazard)
                    
                else:
            
                    if usegamma_gefsv12[jy,ix] == -1:
                        # --- flagged as basically no training data.
                        qpositive = 0.0
                    else:  # --- flagged as minimal training data 
                        #        so use Gamma distribution fit
                        alpha_hat = spline_info_gefsv12[jy,ix,0,0] 
                        beta_hat = spline_info_gefsv12[jy,ix,1,0] 
                        y0 = precip_amount[jy,ix] / beta_hat
                        qpositive = stats.gamma.cdf(y0, alpha_hat)
                
            gefsv12_quantiles[jy,ix] = qpositive
            
    # --- interpolate to the NDFD grid
    
    print ('   interpolating quantiles to NDFD grid')
    gefsv12_quantiles_flipud = np.flipud(gefsv12_quantiles) 
    gefsv12_quantiles_on_ndfd = interp(gefsv12_quantiles_flipud, \
        lons_1d_realtime, lats_1d_realtime_flipud,  \
        lons_ndfd, lats_ndfd, checkbounds=False, \
        masked=False, order=1)  
        
    # ---- if the quantile is above the 99th percentile and use99 == True, 
    #      then set the quantile to 0.99, and determine the offset
    
    print ('   truncating + quantiles above 0.99 to 0.99')
    ninetynines = 0.99*np.ones((ny_ndfd, nx_ndfd), dtype=np.float64)
    gefsv12_quantiles_on_ndfd = np.where( gefsv12_quantiles_on_ndfd < 0.99, 
        gefsv12_quantiles_on_ndfd, ninetynines)

    # ---- for points with quantiles >= 0.99, set the offset

    print ('   determining the offset for points with + quantiles >= 0.99')
    zeros = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
    offset_on_ndfd = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
    offset_on_ndfd = np.where( gefsv12_quantiles_on_ndfd >= 0.99, \
        precip_gefsv12_on_ndfd - quantile_99_ndfd, zeros)

    # ---- change the output quantile from the percentile of positive 
    #      values to the percentile in the distribution including zeros
    
    print ('   changing output quantiles to include zeros.')
    gefsv12_quantiles_on_ndfd = fraction_zero_gefsv12_on_ndfd + \
        (1.0 - fraction_zero_gefsv12_on_ndfd)*gefsv12_quantiles_on_ndfd

    return gefsv12_quantiles_on_ndfd, offset_on_ndfd    

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
    #input_directory = '/Volumes/NBM/gefsv12/2018_retro/'
    input_directory = '/Volumes/NBM/retro/'+clead+'/'
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
            
def generate_probabilities(nmembers, thresholds, precip_mean, precip_ens, \
    prob, ny, nx, hist_thresholds, closest_hist, dressit) : 
    
    """ determine the ensemble probability of exceeding the threshold,
        weighted by closest-member histograms.   Simplify calculations
        and eliminate sums by using CDF of closest-member histograms. 
        Calls a fortran routine for computational efficiency."""
    
    # ---- set weights based on the ensemble-mean amount
    
    nchcats = len(hist_thresholds)+1
    nthresholds = len(thresholds)
    prob = np.zeros((ny, nx), dtype=np.float64)
    if dressit == True:
        prob = weighted_probability_dressed_f90(precip_mean, \
            precip_ens, thresholds, hist_thresholds, closest_hist, \
            ny, nx, nmembers, nthresholds, nchcats)
    else:
        prob = weighted_probability_f90(precip_mean, \
            precip_ens, thresholds, hist_thresholds, closest_hist, \
            ny, nx, nmembers, nthresholds, nchcats)
            
    a = np.where(prob > 1.0)
    if a: prob[a] = 1.0
    return prob
    
# =====================================================================    
    
def write_probabilities_to_netcdf (outfile, ncname, idate, cyyyymmddhh, \
    cyyymmddhh_fcst,  prob_qmapped_out,  ensmean_qmapped_out):    
    
    """ write the raw and quantile-mapped ensemble probabilities to 
    a netCDF file """
    
    # ---- write probabilities and close

    print ('   writing probs to netCDF file ',outfile)
    if ncname == 'nc_fullfield':
        nc_fullfield = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        nc_fullfield['yyyymmddhh_init'][idate] = int(cyyyymmddhh)
        nc_fullfield['yyyymmddhh_fcst'][idate] = int(cyyyymmddhh_fcst)
        nc_fullfield['probability_qmapped'][idate] = prob_qmapped_out[:,:,:]
        nc_fullfield['ensmean_qmapped'][idate] = ensmean_qmapped_out[:,:]
        nc_fullfield.close()
    elif ncname == 'nc_upscaled':
        nc_upscaled = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        nc_upscaled['yyyymmddhh_init'][idate] = int(cyyyymmddhh)
        nc_upscaled['yyyymmddhh_fcst'][idate] = int(cyyyymmddhh_fcst)
        nc_upscaled['probability_qmapped'][idate] = prob_qmapped_out[:,:,:]
        nc_upscaled['ensmean_qmapped'][idate] = ensmean_qmapped_out[:,:]
        nc_upscaled.close()
    elif ncname == 'nc_thinned':
        nc_thinned = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        nc_thinned['yyyymmddhh_init'][idate] = int(cyyyymmddhh)
        nc_thinned['yyyymmddhh_fcst'][idate] = int(cyyyymmddhh_fcst)
        nc_thinned['probability_qmapped'][idate] = prob_qmapped_out[:,:,:]
        nc_thinned['ensmean_qmapped'][idate] = ensmean_qmapped_out[:,:]
        nc_thinned.close()
    else:
        print ('invalid ncname in initalize_netCDF_probfiles = ', ncname)
        sys.exit()
        
    istat = 0
    return istat
    
# =====================================================================

def read_closest_member_histogram(master_directory_histogram, ctype, \
    cmonth, clead):
    
    """ read the closest-member histogram data from netCDF file"""    
    
    infile = master_directory_histogram +'closest_member_histogram_'+\
        ctype+'_month='+cmonth+'_lead='+clead+'h.nc'
    print ('reading from ',infile)
    nc = Dataset(infile,'r')
    histogram_thresholds = nc.variables['thresholds'][:]
    nthresh = len(histogram_thresholds)
    closest_member_histogram = nc.variables['closest_hist'][:,:]
    ns,nt = np.shape(closest_member_histogram)
    closest_member_histogram[:,0] = closest_member_histogram[:,0] / \
        np.sum(closest_member_histogram[:,0])
    closest_member_histogram[:,1] = closest_member_histogram[:,1] / \
        np.sum(closest_member_histogram[:,1])
    closest_member_histogram[:,2] = closest_member_histogram[:,2] / \
        np.sum(closest_member_histogram[:,2])
    closest_member_histogram[:,3] = closest_member_histogram[:,3] / \
        np.sum(closest_member_histogram[:,3])
    closest_member_histogram[:,4] = closest_member_histogram[:,4] / \
        np.sum(closest_member_histogram[:,4])
    closest_member_histogram[:,5] = closest_member_histogram[:,5] / \
        np.sum(closest_member_histogram[:,5])
        
    #print ('closest_member_histograms')
    #for i in range(ns):
    #    ch1 = '%.4f' %(closest_member_histogram[i,0])
    #    ch2 = '%.4f' %(closest_member_histogram[i,1])
    #    ch3 = '%.4f' %(closest_member_histogram[i,2])
    #    ch4 = '%.4f' %(closest_member_histogram[i,3])
    #    print (ch1,' ',ch2,' ', ch3, ' ',ch4)
        
    nmemx25, ncats = np.shape(closest_member_histogram)
    nc.close()
    return histogram_thresholds, closest_member_histogram
        #, closest_member_histogram_cdf
    
# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh_begin = sys.argv[1] # the initial year/month/day/hour
clead = sys.argv[2]  # lead time as 3-digit number, e.g., 
                     # forecast ending at +18h is 018.
                     
use99 = True # apply fixed offset beyond 99th percentile
jnd = 1350
ind = 1110

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
    # where analysis spline info is stored
master_directory_forecast_spline = '/Volumes/NBM/conus_gefsv12/CDF_spline/'
    # where forecast spline info is stored
master_directory_forecast_qmapped = '/Volumes/NBM/conus_gefsv12/qmapped/'
    # where output quantile-mapped members are stored.
master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/thinned/'
    # where the thinned field probabilities are stored, if netcdf_thinned == True
master_directory_panal = '/Volumes/NBM/conus_panal/'    #Lesley, this is new.
    # where the gridded precipitation analyses are archived.
master_directory_histogram_thinned = '/Volumes/NBM/conus_gefsv12/thinned/chist/'

ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31] # days of month, normal year
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31] # leap year
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
nthresholds = len(thresholds)
firstdate = True
output_qmapped = True
netcdf_thinned = True
nstride = 10 # for thinning the output
nstride_qmap = 12 + int(clead)//18 
#nstride_qmap = 15 # for using grid point stencil. 2.5 km pixels, 15 goes out to +/- 75 km
dressit = False

# ---- read the closest-member histograms for thinned

if netcdf_thinned == True:
    histogram_thresholds, closest_member_histogram_thinned = \
        read_closest_member_histogram(master_directory_histogram_thinned, \
        'thinned', cmonth, clead)
        
if iyear%4 == 0:
    ndays = ndaysomo_leap[imonth]
else:
    ndays = ndaysomo[imonth]
        
cyyyymmddhh_end = dateshift(cyyyymmddhh_begin,24*(ndays-1)) 
    # set this to the first date of the month if you want to run a 
    # full month of forecasts
#cyyyymmddhh_end = dateshift(cyyyymmddhh_begin,0)

date_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 24)
   
# ========================================================================
# ---- loop over dates to process
# ========================================================================

for idate, cyyyymmddhh in enumerate(date_list):
    
    cyyyymmddhh_fcst = dateshift(cyyyymmddhh, int(clead))
    iyyyymmddhh_veriftime = int(cyyyymmddhh_fcst)
    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print ('starting quantile mapping for date = ', \
        cyyyymmddhh,' at time = ',time)
    if firstdate == True:
    
        # ---- read forecast spline parameters from netCDF file
    
        print ('reading forecast spline parameters')
        lons_1d_realtime, lats_1d_realtime, \
            spline_info_gefsv12, fraction_zero_gefsv12, \
            usegamma_gefsv12, quantile_99, \
            number_knots_gefsv12, ny_gefsv12, nx_gefsv12, \
            lons_fcst_2d, lats_fcst_2d, clead_use = \
            read_forecast_spline_info(cyyyymmddhh, \
            ccmonth, master_directory_forecast_spline)
    
        # ---- read precipitation analysis inverse CDF 
        #      spline parameters from netCDF file    

        print ('reading analyzed spline parameters')
        spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
            number_knots_ndfd, lons_ndfd, lats_ndfd, ny_ndfd, \
            nx_ndfd = read_analysis_spline_info (clead_use, \
            master_directory_panal_spline, ccmonth)
            
        # ---- flip and then interpolate the GEFSv12  
        #      climatological fraction zero and latitude 
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
            
        quantile_99_flipped = np.flipud(quantile_99)
        quantile_99_ndfd = interp(quantile_99_flipped, \
            lons_1d_realtime, lats_1d_realtime_flipud, \
            lons_ndfd, lats_ndfd, checkbounds=False, \
            masked=False, order=1)
            
        # ---- set up full output grids and working arrays

        precip_ens_raw_ndfd = np.zeros((nmembers, ny_ndfd, nx_ndfd), \
            dtype=np.float32)
        precip_ens_qmapped_ndfd = np.zeros((nmembers*25, ny_ndfd, nx_ndfd), \
            dtype=np.float32)
        offset_ens_ndfd = np.zeros((nmembers, ny_ndfd, nx_ndfd), \
            dtype=np.float32)
        
        # ---- Take every 10th pixel in lat/lon to subsample starting at
        #      the sixth pixel.  Determine the grid size of an upscaled
        #      field and pass these in to initialize netCDF file 
        #      for the thinned
        
        lons_thinned = lons_ndfd[5::nstride, 5::nstride]
        lats_thinned = lats_ndfd[5::nstride, 5::nstride]        
        ny_thinned, nx_thinned = np.shape(lons_thinned)
        precip_ens_raw_thinned = \
            np.zeros((nmembers*25, ny_thinned, nx_thinned), dtype=np.float32)
        precip_ens_qmapped_thinned = \
            np.zeros((nmembers*25, ny_thinned, nx_thinned), dtype=np.float32)
        if use99 == False:
            if dressit == True:
                outfile = master_directory_thinned_output + ccmonth + cyyyy + \
                    '_lead'+clead+'_probabilities_weighted_dressed_thinned.nc'
            else:
                outfile = master_directory_thinned_output + ccmonth + cyyyy + \
                        '_lead'+clead+'_probabilities_weighted_thinned.nc'
        else:
            if dressit == True:
                outfile = master_directory_thinned_output + ccmonth + cyyyy + \
                    '_use99_lead'+clead+'_probabilities_weighted_dressed_thinned.nc'
            else:
                outfile = master_directory_thinned_output + ccmonth + cyyyy + \
                    '_use99_lead'+clead+'_probabilities_weighted_thinned.nc'
                        
        ncname = 'nc_thinned'
        istat = initialize_netCDF_probfiles(outfile, ncname, \
            ny_thinned, nx_thinned, lats_thinned, lons_thinned, \
            nthresholds, thresholds)
                
        prob_raw_thinned = np.zeros((nthresholds, ny_thinned, nx_thinned), \
            dtype=np.float32)
        prob_qmapped_thinned = np.zeros((nthresholds, ny_thinned, nx_thinned), \
            dtype=np.float32)
        prob_qmapped_weighted_thinned = np.zeros((nthresholds, \
            ny_thinned, nx_thinned), dtype=np.float32)
        prob_thinned = np.zeros((nthresholds, ny_thinned, nx_thinned), dtype=np.float64)
        zeros_thinned = np.zeros((ny_thinned, nx_thinned), dtype=np.float64)
        ones_thinned = np.ones((ny_thinned, nx_thinned), dtype=np.float64)
        missingv_thinned = -99.99*np.ones((ny_thinned, nx_thinned), dtype=np.float64)
        closest_member_stats_thinned = np.zeros((25*nmembers,4), dtype=np.int64)
        
        # ---- having done the upfront work, set a flag so this isn't
        #      repeated as we loop thru other dates
        
        firstdate = False
        
    # ========================================================================
    # ---- loop over members, read in this member of GEFS, 
    #      interpolate to NDFD grid and quantile map
    # ========================================================================

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
        print ('   member max precip = ', np.max(precip_gefsv12_on_ndfd))
        precip_ens_raw_ndfd[imem,:,:] = precip_gefsv12_on_ndfd[:,:]
        
        # ---- now loop over grid points and obtain the forecast quantile
        #      associated with this 0.25 degree forecast grid point.  This is the
        #      quantile in the overall distribution, including zeros.

        print ('   getting quantiles of forecast')
        gefsv12_quantiles_on_ndfd, offset_on_ndfd = get_quantile_gefsv12(\
            precip_realtime, precip_gefsv12_on_ndfd, spline_info_gefsv12, \
            ny_gefsv12, nx_gefsv12, fraction_zero_gefsv12, \
            fraction_zero_gefsv12_on_ndfd, \
            usegamma_gefsv12, quantile_99_ndfd, use99, \
            number_knots_gefsv12, lons_1d_realtime, \
            lats_1d_realtime_flipud, lons_ndfd, lats_ndfd, jnd, ind)
        
        # ---- apply quantile mapping procedure to this member.  Now with
        #      5x5 stencil, so 25x more qmapped members.
        
        print ('   applying quantile mapping')
        qmapped_precip_thinned = np.zeros((25,ny_thinned, nx_thinned), dtype=np.float64)    
        qmapped_precip_thinned = \
            qmapping_spline_fivexfive_thinned(gefsv12_quantiles_on_ndfd, \
    	    precip_gefsv12_on_ndfd, spline_info_inv, \
            fraction_zero_gefsv12_on_ndfd, fraction_zero_ndfd, \
            usegamma_ndfd, use99, offset_on_ndfd, \
            number_knots_ndfd, nstride_qmap, nstride, \
            ny_thinned, nx_thinned, ny_ndfd, nx_ndfd )         
            
        print ('   min, max qmapped_precip = ', \
            np.min(qmapped_precip_thinned), np.max(qmapped_precip_thinned))
        
        # ---- store these quantile-mapped members
        
        # ---- store these quantile-mapped members
        
        istart = imem*25                
        precip_ens_raw_thinned[imem,:,:] = \
            precip_gefsv12_on_ndfd[0:-1:nstride, 0:-9:nstride]                 
        istart = imem*25
        for istencil in range(25):
            precip_ens_qmapped_thinned[istart+istencil,:,:] = \
                qmapped_precip_thinned[istencil,:,:]
     
    # ========================================================================        
    # ---- generate probabilities, raw and quantile mapped for each
    #      specified precipitation threshold, and for full field,
    #      thinned, and upscaled output
    # ========================================================================

    ensmean_qmapped_thinned = np.squeeze(np.mean(\
        precip_ens_qmapped_thinned, axis=0))
    print ('   generating probabilities for thinned') 
    prob_qmapped_weighted_thinned = \
        generate_probabilities(nmembers*25, thresholds, \
        ensmean_qmapped_thinned, precip_ens_qmapped_thinned, \
        prob_thinned, ny_thinned, nx_thinned, \
        histogram_thresholds, closest_member_histogram_thinned, dressit)
    print ('   min, max, mean prob thinned = ',\
        np.min(prob_qmapped_weighted_thinned), \
        np.max(prob_qmapped_weighted_thinned), \
        np.mean(prob_qmapped_weighted_thinned))

    # ========================================================================
    # ---- write the probabilities over full set of thresholds to netCDF file.
    #      for the thinned forecasts if desired.
    #      Also read in thinned verification precipitation analysis and
    #      get closest-member histogram statistics, write out to file
    # ========================================================================

    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    ncname = 'nc_thinned'
    if use99 == False:
        if dressit == True:
            outfile = master_directory_thinned_output + ccmonth + cyyyy + \
                '_lead'+clead+'_probabilities_weighted_dressed_thinned.nc'
        else:
            outfile = master_directory_thinned_output + ccmonth + cyyyy + \
                '_lead'+clead+'_probabilities_weighted_thinned.nc'
    else:
        if dressit == True:
            outfile = master_directory_thinned_output + ccmonth + cyyyy + \
                '_use99_lead'+clead+'_probabilities_weighted_dressed_thinned.nc'
        else:
            outfile = master_directory_thinned_output + ccmonth + cyyyy + \
                '_use99_lead'+clead+'_probabilities_weighted_thinned.nc'
    print ('   writing thinned field probabilities to file. ',time)
    print (outfile)
    istat = write_probabilities_to_netcdf(outfile, ncname, idate, \
        cyyyymmddhh, cyyyymmddhh_fcst, prob_qmapped_weighted_thinned,  \
        ensmean_qmapped_thinned)
    print ('   writing thinned field probabilities and ens means completed.')
    
print ('DONE!') 