"""
qmap_reforecast_for_histogram_stats_thinned_x2.py cyyyymmddhh clead

This python routine quantile maps GECFSv12 forecasts and then
compares the quantile-mapped forecast information to analyzed data
in order to generate the statistics that are required to generate
estimates of the closest-member histograms and dressing statistics.
It does this by reading pre-generated spline coefficient information
on the cumulative hazard function of the fitted CDF for the forecast,
the inverse of the hazard function CDF for the
analyzed, and the forecast in question.   It then applies
quantile mapping, and compares to observations, saving statistics.

coded by: Tom Hamill, Jun 2021, tom.hamill@noaa.gov

"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from dateutils import daterange, dateshift # Jeff Whitaker's utility
import pygrib
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
from qmapping_spline_fivexfive_thinned import \
    qmapping_spline_fivexfive_thinned # from f2py
import scipy.stats as stats
from mpl_toolkits.basemap import Basemap, interp
import _pickle as cPickle
from compute_closest_member_dress_f90 \
    import compute_closest_member_dress_f90
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
    print ('    reading forecast spline information from ', infile)
    nc = Dataset(infile)
    lons_spline_gefsv12_1d = nc.variables['lons'][:]
    lats_spline_gefsv12_1d = nc.variables['lats'][:]
    spline_info_gefsv12 = nc.variables['spline_info'][:,:,:,:]
    fraction_zero_gefsv12 = nc.variables['fzero'][:,:]
    usegamma_gefsv12 = nc.variables['usegamma'][:,:]
    quantile_99 = nc.variables['quantile_99'][:,:]
    number_knots_gefsv12 = nc.variables['number_knots'][:,:]
    ny_gefsv12, nx_gefsv12 = np.shape(usegamma_gefsv12)
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
    print ('    reading from ', infile)
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
    print ('    reading ',gribfilename, fexist_grib)
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
    usegamma_gefsv12, quantile_99, number_knots_gefsv12, use99):

    """ this gets the quantile associated with a given precipitation
    amount for GEFSv12 data, this month and 6-hour period.  This is
    the quantile in the overall distribution, including zeros """

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
            quantile = fraction_zero_gefsv12[jy,ix] + \
                (1.0 - fraction_zero_gefsv12[jy,ix])*qpositive
        else:

            if usegamma_gefsv12[jy,ix] == -1:
                # --- flagged as basically no training data.
                quantile = 0.0
            else:  # --- flagged as minimal training data - use Gamma
                alpha_hat = spline_info_gefsv12[jy,ix,0,0]
                beta_hat = spline_info_gefsv12[jy,ix,1,0]
                y0 = precip_amount[jy,ix] / beta_hat
                qpositive = stats.gamma.cdf(y0, alpha_hat)
                quantile = fraction_zero_gefsv12[jy,ix] + \
                    (1.0 - fraction_zero_gefsv12[jy,ix])*qpositive

        if use99 == True and qpositive > 0.99:
            offset_out = precip_amount[jy,ix] - quantile_99[jy,ix]
            if qpositive > 0.99:
                quantile = fraction_zero_gefsv12[jy,ix] + \
                    (1.0 - fraction_zero_gefsv12[jy,ix])*0.99

    if quantile >= 1.0: quantile = 0.9999
    if quantile <= 0.0: quantile = 0.0

    return quantile, offset_out

# =====================================================================

def get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead):

    """ read in the global forecast.  Subset to CONUS. """

    nib = 886  # precomputed for set of lat-lons that surrounds
    nie = nib+318 # the NDFD CONUS 2.5-km Lambert conformal domain.
    njb = 131
    nje = njb+153

    # ---- read in forecast grid covering the whole globe.

    cycle = cyyyymmddhh[8:10]

    input_directory = '/Volumes/NBM/retro/'+clead+'/'
    infile = input_directory + 'apcp_'+cyyyymmddhh + \
        '_'+cmem+'_.f'+clead+'.grib2'
    print ('    reading from ', infile)
    endStep = int(clead)
    istat, precip_realtime, lats_full, lons_full = \
        read_gribdata(infile, endStep)

    # ---- subset for CONUS.

    precip_realtime = precip_realtime[njb:nje,nib:nie]
    lons_1D_realtime = lons_full[0,nib:nie]-360.
    lats_1D_realtime = lats_full[njb:nje,0]

    return precip_realtime, lons_1D_realtime, lats_1D_realtime

# =====================================================================

def read_precipitation_analyses(master_directory, cyear, \
    cmonth, ctype, iyyyymmddhh_veriftime, missingv, nstride):

    """ read the merged CCPA/MSWEP precipitation analysis
        appropriate to this type ('','_thinned','_upscaled'
        and appropriate for the input verification time)"""

    ktr = 0
    infile = master_directory + cyear + cmonth + \
        '_ccpa_on_ndfd_grid_6hourly'+ctype+'.nc'
    print ('    trying to read',infile)
    fexist = os.path.exists(infile)
    print ('    fexist = ',fexist)
    if fexist == True:
        nc = Dataset(infile)
        yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
        idx = np.where(yyyymmddhh_end == iyyyymmddhh_veriftime)[0]
        precip_anal = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
        precip_anal = np.where(precip_anal < 500., precip_anal, missingv)
        nc.close()
        istat = 1
    else:
        istat = -1
        precip_anal = np.copy(missingv)

    return precip_anal,istat

# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh_begin = sys.argv[1] # the initial year/month/day/hour
clead = sys.argv[2]  # lead time as 3-digit number, e.g.,
                     # forecast ending at +18h is 018.

# --- other parameters

use99 = True # apply fixed offset beyond 99th percentile
nstride = 2 # for thinning the output and speeding processing
nstride_qmap = 12 + int(clead)//18 # determines spacing between stencil pts.
closest_member_thresh_low = 0.01  # thresholds for closest-member histograms
closest_member_thresh_lowmod = 0.1
closest_member_thresh_mod = 0.5
closest_member_thresh_modhigh = 2.0
closest_member_thresh_high = 6.0
closest_member_thresh_superhigh = 15.0
save_chist = True # save the closest-member histogram data

# ---- various initialization and carving out month information.

firstdate = True # used to flag the first date through, read in stuff
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
master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/thinned2/'
    # where the thinned field probabilities are stored, if netcdf_thinned == True
master_directory_panal = '/Volumes/NBM/conus_panal/'
    # where the gridded precipitation analyses are archived.

ndaysomo = [31, 28, 31,  30, 31, 30,  \
    31, 31, 30,  31, 30, 31] # days of month, normal year
ndaysomo_leap = [31, 28, 31,  30, 31, \
    30,  31, 31, 30,  31, 30, 31] # leap year
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
nthresholds = len(thresholds)


if iyear%4 == 0:
    ndays = ndaysomo_leap[imonth]
else:
    ndays = ndaysomo[imonth]

cyyyymmddhh_end = dateshift(cyyyymmddhh_begin,24*(ndays-1))

date_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 24)

# ========================================================================
# ---- loop over dates to process
# ========================================================================

for idate, cyyyymmddhh in enumerate(date_list):

    cyyyymmddhh_fcst = dateshift(cyyyymmddhh, int(clead))
    iyyyymmddhh_veriftime = int(cyyyymmddhh_fcst)
    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print ('--- starting quantile mapping for date = ', \
        cyyyymmddhh,' at time = ',time)
    if firstdate == True:

        # ---- read forecast spline parameters from netCDF file

        print ('--- reading forecast spline parameters')
        lons_1d_realtime, lats_1d_realtime, \
            spline_info_gefsv12, fraction_zero_gefsv12, \
            usegamma_gefsv12, quantile_99, \
            number_knots_gefsv12, ny_gefsv12, nx_gefsv12, \
            lons_fcst_2d, lats_fcst_2d, clead_use = \
            read_forecast_spline_info(cyyyymmddhh, \
            ccmonth, master_directory_forecast_spline)

        # ---- read precipitation analysis inverse CDF
        #      spline parameters from netCDF file

        print ('--- reading analyzed spline parameters')
        spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
            number_knots_ndfd, lons_ndfd, lats_ndfd, ny_ndfd, \
            nx_ndfd = read_analysis_spline_info (clead_use, \
            master_directory_panal_spline, ccmonth)

        # ---- flip and then interpolate the GEFSv12
        #      climatological fraction zero and latitude
        #      array to NDFD grid. Flipping is necessary as
        #      Basemap.interp requires input lats and data
        #      be oriented S to N.

        print ('--- interpolating fraction_zero of GEFSv12 to NDFD grid. ')
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
        offset_ens_ndfd = np.zeros((nmembers, ny_ndfd, nx_ndfd), \
            dtype=np.float32)
        precip_ens_qmapped_ndfd = np.zeros((nmembers*25, ny_ndfd, nx_ndfd), \
            dtype=np.float32)

        # ---- Take every 2nd pixel in lat/lon to subsample.
        #      Determine the grid size of an upscaled
        #      field and pass these in to initialize netCDF file
        #      for the thinned

        lons_thinned = lons_ndfd[nstride//2::nstride, nstride//2::nstride]
        lats_thinned = lats_ndfd[nstride//2::nstride, nstride//2::nstride]
        ny_thinned, nx_thinned = np.shape(lons_thinned)
        precip_ens_raw_thinned = \
            np.zeros((nmembers, ny_thinned, nx_thinned), dtype=np.float32)
        precip_ens_qmapped_thinned = \
            np.zeros((nmembers*25, ny_thinned, nx_thinned), dtype=np.float32)
        zeros_thinned = np.zeros((ny_thinned, nx_thinned), dtype=np.float64)
        ones_thinned = np.ones((ny_thinned, nx_thinned), dtype=np.float64)
        missingv_thinned = -99.99*np.ones((ny_thinned, nx_thinned), dtype=np.float64)

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
        print ('--- processing member = ',cmem,time)

        # ---- read in & extract the real-time precipitation for this domain.
        #      Set GEFSv12 grid dimensions and 1/4-deg lon and lat.

        precip_realtime, lons_1d_realtime, lats_1d_realtime = \
            get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead)

        # ---- flip upside down, as subsequent interpolation requires
        #      S to N with increasing j index.  Then
        #      interpolate precipitation to NDFD grid.

        print ('    interpolating and flipping real-time forecast')
        precip_realtime_flipud = np.flipud(precip_realtime)
        precip_gefsv12_on_ndfd = interp(precip_realtime_flipud, \
            lons_1d_realtime, lats_1d_realtime_flipud,  \
            lons_ndfd, lats_ndfd, checkbounds=False, \
            masked=False, order=1)
        precip_ens_raw_ndfd[imem,:,:] = precip_gefsv12_on_ndfd[:,:]

        # ---- now loop over grid points and obtain the forecast quantile
        #      associated with this 0.25 degree forecast grid point.  This is the
        #      quantile in the overall distribution, including zeros.

        print ('    getting quantiles of forecast')
        gefsv12_quantiles = np.zeros((ny_gefsv12, nx_gefsv12), \
            dtype=np.float64)
        offset = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.float64)
        for jy in range(ny_gefsv12):
            for ix in range(nx_gefsv12):
                gefsv12_quantiles[jy,ix], offset[jy,ix] = get_quantile_gefsv12( \
                    precip_realtime, jy, ix, spline_info_gefsv12, \
                    fraction_zero_gefsv12, usegamma_gefsv12, \
                    quantile_99, number_knots_gefsv12, use99)

        # ---- interpolate the forecast quantile to the NDFD grid.

        print ('    interpolating quantiles to NDFD grid')
        gefsv12_quantiles_flipud = np.flipud(gefsv12_quantiles)
        gefsv12_quantiles_on_ndfd = interp(gefsv12_quantiles_flipud, \
            lons_1d_realtime, lats_1d_realtime_flipud,  \
            lons_ndfd, lats_ndfd, checkbounds=False, \
            masked=False, order=1)
        offset_flipud = np.flipud(offset)
        offset_on_ndfd = interp(offset_flipud, \
            lons_1d_realtime, lats_1d_realtime_flipud,  \
            lons_ndfd, lats_ndfd, checkbounds=False, \
            masked=False, order=1)

        # ---- apply quantile mapping procedure to this member.  Now with
        #      5x5 stencil, so 25x more qmapped members.

        print ('    applying quantile mapping')
        qmapped_precip_thinned = np.zeros((25,ny_thinned, nx_thinned), \
            dtype=np.float64)
        qmapped_precip_thinned = \
            qmapping_spline_fivexfive_thinned(gefsv12_quantiles_on_ndfd, \
    	    precip_gefsv12_on_ndfd, spline_info_inv, \
            fraction_zero_gefsv12_on_ndfd, fraction_zero_ndfd, \
            usegamma_ndfd, use99, offset_on_ndfd, \
            number_knots_ndfd, nstride_qmap, nstride, \
            ny_thinned, nx_thinned, ny_ndfd, nx_ndfd )
        print ('    min, max qmapped_precip = ', \
            np.min(qmapped_precip_thinned), \
            np.max(qmapped_precip_thinned))

        # ---- store these quantile-mapped members

        istart = imem*25
        precip_ens_raw_thinned[imem,:,:] = \
            precip_gefsv12_on_ndfd[0:-1:nstride, 0:-1:nstride]
        istart = imem*25
        for istencil in range(25):
            precip_ens_qmapped_thinned[istart+istencil,:,:] = \
                qmapped_precip_thinned[istencil,:,:]

    ensmean_raw_thinned = np.squeeze(np.mean(\
        precip_ens_raw_thinned, axis=0))
    ensmean_qmapped_thinned = np.squeeze(np.mean(\
        precip_ens_qmapped_thinned, axis=0))

    # ---- compute the closest-member histogram and the ensemble
    #      mean of the quantile-mapped precipitation amount.

    if save_chist == True:

        # ---- read in the merged CCPA/MSWEP precipitation analyses
        #      needed for histogram statistics

        now = datetime.now()
        time = now.strftime("%H:%M:%S")
        print ('--- computing closest-member histogram statistics. ', time)
        ctype = '_thinned2'
        cmonth_verif = cyyyymmddhh_fcst[4:6]
        cyear_verif = cyyyymmddhh_fcst[0:4]
        precip_anal_thinned, istat = read_precipitation_analyses( \
            master_directory_panal, cyear_verif, cmonth_verif, ctype, \
            iyyyymmddhh_veriftime, missingv_thinned, nstride)
        if istat > 0:

            closest_histogram = np.zeros((25*nmembers, 7), dtype=int)
            nseven = 7
            nmembers_x25 = nmembers*25
            n251 = 251
            sum_fracrank = np.zeros((n251), dtype=np.float64)
            sum_fracrank_squared = np.zeros((n251), dtype=np.float64)
            nsamps_fracrank = np.zeros((n251), dtype=np.int32)
            sumxi_low = np.zeros((n251), dtype=np.float64)
            sumxi2_low = np.zeros((n251), dtype=np.float64)
            nsamps_low = np.zeros((n251), dtype=np.int32)
            sumxi_mid = np.zeros((n251), dtype=np.float64)
            sumxi2_mid = np.zeros((n251), dtype=np.float64)
            nsamps_mid = np.zeros((n251), dtype=np.int32)
            sumxi_high = np.zeros((n251), dtype=np.float64)
            sumxi2_high = np.zeros((n251), dtype=np.float64)
            nsamps_high = np.zeros((n251), dtype=np.int32)

            closest_histogram, sumxi_low, sumxi2_low, nsamps_low, \
                sumxi_mid, sumxi2_mid, nsamps_mid,\
                sumxi_high, sumxi2_high, nsamps_high,\
                istat = compute_closest_member_dress_f90( \
                nmembers_x25, ny_thinned, nx_thinned, nseven, n251, \
                closest_member_thresh_low, closest_member_thresh_lowmod, \
                closest_member_thresh_mod, closest_member_thresh_modhigh, \
                closest_member_thresh_high, closest_member_thresh_superhigh, \
                ensmean_qmapped_thinned, \
                precip_ens_qmapped_thinned, precip_anal_thinned)

            # --- write c-m histogram for day of interest.   These are
            #     later synthesized into statistics and plotted with
            #     plot_closest_member_histogram.py

            outfile = master_directory_thinned_output +\
                'closest_histogram'+cyyyymmddhh_begin+'_'+clead+'.cPick'
            ouf = open(outfile, 'wb')
            cPickle.dump(closest_histogram, ouf)
            cPickle.dump(sumxi_low, ouf)
            cPickle.dump(sumxi2_low, ouf)
            cPickle.dump(nsamps_low, ouf)
            cPickle.dump(sumxi_mid, ouf)
            cPickle.dump(sumxi2_mid, ouf)
            cPickle.dump(nsamps_mid, ouf)
            cPickle.dump(sumxi_high, ouf)
            cPickle.dump(sumxi2_high, ouf)
            cPickle.dump(nsamps_high, ouf)
            ouf.close()

        else:
            print ('there was a problem with verification data, not using this date.')

    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print ('--- Finished processing this date. ', time)

print ('DONE!')