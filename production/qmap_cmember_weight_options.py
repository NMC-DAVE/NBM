"""
qmap_cmember_weight_options.py cyyyymmddhh clead 

This python routine reads in pre-generated spline coefficient information 
on the cumulative hazard function of the fitted CDF for the forecast, 
the inverse of the hazard function CDF for the
analyzed, and the forecast in question.   It then applies
quantile mapping with and without closest-member histogram weighting
and dressing, generates probabilities. It writes these
to a netCDF-formatted output file.   Because probabilities on the
high-resolution grid are too large, it permits saving every 10th point
or upscaling by a factor of 10.  (netcdf_thinned = True)

coded by: Tom Hamill, Dec 2021, tom.hamill@noaa.gov 

"""

import os, sys
from datetime import datetime 
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from dateutils import daterange, dateshift # Whitaker library routinez
import pygrib # grib-reading routine
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
from qmapping_spline_fivexfive_flexiknot import \
    qmapping_spline_fivexfive_flexiknot # from f2py  

import scipy.stats as stats
from mpl_toolkits.basemap import Basemap, interp
import _pickle as cPickle
from PIL import Image


# ---- locally developed procedures

from read_forecast_spline_info import read_forecast_spline_info
from read_analysis_spline_info import read_analysis_spline_info
from fraczero_possamps import fraczero_possamps
from initialize_netCDF_probfiles import initialize_netCDF_probfiles
from get_domain_subset_of_gefs import get_domain_subset_of_gefs
from generate_probabilities import generate_probabilities
from get_quantile_gefsv12 import get_quantile_gefsv12
from weighted_probability_f90 import weighted_probability_f90 # from f2py 
from weighted_probability_dressed_f90 import weighted_probability_dressed_f90 # from f2py 
from read_closest_member_histogram import read_closest_member_histogram
from write_probabilities_to_netcdf import write_probabilities_to_netcdf

# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh_begin = sys.argv[1] # the initial year/month/day/hour

iyear = int(cyyyymmddhh_begin[0:4])
cmonth = cyyyymmddhh_begin[4:6]
cyyyy = cyyyymmddhh_begin[0:4]
cyyyymm = cyyyymmddhh_begin[0:6]
chh = cyyyymmddhh_begin[8:10]
imonth = int(cmonth)-1

clead = sys.argv[2]  # lead time as 3-digit number, e.g., 
                     # forecast ending at +18h is 018.
                     

# ---- various initialization and carving out month information.


cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cmembers = ['c00','p01', 'p02','p03','p04',\
    'p05','p06','p07','p08','p09','p10',\
    'p11', 'p12','p13','p14','p15','p16','p17','p18','p19','p20',\
    'p21', 'p22','p23','p24','p25','p26','p27','p28','p29','p30']
nmembers = len(cmembers)
#gefs_data_directory = '/Volumes/NBM/retro/'+clead+'/'
gefs_data_directory = '/Volumes/NBM/gefsv12/2021/'
master_directory_panal_spline = '/Volumes/NBM/conus_panal/CDF_spline/'
    # where analysis spline info is stored
master_directory_forecast_spline = '/Volumes/NBM/conus_gefsv12/CDF_spline/'
    # where forecast spline info is stored
master_directory_forecast_qmapped = '/Volumes/NBM/conus_gefsv12/qmapped/'
    # where output quantile-mapped members are stored.
master_directory_fullfield_output = '/Volumes/NBM/conus_gefsv12/fullfield/'
    # where the full field probabilities are stored, if netcdf_fullfield == True
master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/thinned/'
    # where the thinned field probabilities are stored, if netcdf_thinned == True
master_directory_panal = '/Volumes/NBM/conus_panal/'    #Lesley, this is new.
    # where the gridded precipitation analyses are archived.
#master_directory_histogram_thinned = '/Volumes/NBM/conus_gefsv12/thinned2/chist/'
master_directory_histogram_thinned = '/Volumes/NBM/chist/'

ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31] # days of month, normal year
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31] # leap year
thresholds = np.array([0.254, 1.0, 5.0, 10.0, 25.0]) # for probabilistic forecasts
nthresholds = len(thresholds)
firstdate = True
output_qmapped = True
netcdf_thinned = True
netcdf_fullfield = True
nstride = 10 # for thinning the output
nstride_qmap = 12 + int(clead)//18 # spacing between stencil points
use99 = True # apply fixed offset beyond 99th percentile
jnd = 1350 # grid point for diagnostics output
ind = 1110

# ---- read the closest-member histograms for thinned.  Use these
#      statistics for fullfield as well.   If not 00 UTC, we'll adjust
#      the lead time of the closest-member histogram at the same
#      time of the diurnal cycle as for 00 UTC forecast, but 
#      a later forecast lead.

print ('chh = ', chh)
if chh == '00':
    ilead_hist = int(clead)
elif chh == '06':
    ilead_hist = ilead_hist + 6
elif chh == '12':
    ilead_hist = ilead_hist + 12
elif chh == '18':
    ilead_hist = ilead_hist + 18
else:
    print ('invalid initial hour.  Stopping')
    sys.exit()
if ilead_hist < 10:
    clead_hist = '00'+str(ilead_hist)
elif ilead_hist >= 10 and ilead_hist < 100:
    clead_hist = '0'+str(ilead_hist)
else:
    clead_hist = str(ilead_hist)

histogram_thresholds, closest_member_histogram, \
    nmemx25, ncats,  nthresh, b0_mean_lowrank, \
    b1_mean_lowrank, b0_std_lowrank, b1_std_lowrank, \
    b0_mean_midrank, b1_mean_midrank, b0_std_midrank, \
    b1_std_midrank, b0_mean_highrank, b1_mean_highrank, \
    b0_std_highrank, b1_std_highrank = \
    read_closest_member_histogram( \
    master_directory_histogram_thinned, \
    'thinned', cmonth, clead_hist)
        
if iyear%4 == 0:
    ndays = ndaysomo_leap[imonth]
else:
    ndays = ndaysomo[imonth]
        
#cyyyymmddhh_end = dateshift(cyyyymmddhh_begin,24*(ndays-1)) 
    # set this to the first date of the month if you want to run a 
    # full month of forecasts
cyyyymmddhh_end = dateshift(cyyyymmddhh_begin,0)

date_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 24)
   
# ========================================================================
# ---- loop over dates to process
# ========================================================================

for idate, cyyyymmddhh in enumerate(date_list):
    
    cyyyymmddhh_fcst = dateshift(cyyyymmddhh, int(clead))
    print ('cyyyymmddhh_fcst = ', cyyyymmddhh_fcst)
    iyyyymmddhh_veriftime = int(cyyyymmddhh_fcst)
    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print ('starting quantile mapping for date = ', \
        cyyyymmddhh,' at time = ',time)
    if firstdate == True:
    
        # ---- read forecast spline parameters from netCDF file.  Also
        #      returns "clead_use" which indicates which cycle of 
        #      analysis data to read in the following step.
    
        print ('reading forecast spline parameters')
        lons_1d_realtime, lats_1d_realtime, \
            spline_info_gefsv12, fraction_zero_gefsv12, \
            usegamma_gefsv12, quantile_99, \
            number_knots_gefsv12, ny_gefsv12, nx_gefsv12, \
            lons_fcst_2d, lats_fcst_2d, clead_use = \
            read_forecast_spline_info(cyyyymmddhh, \
            ccmonth, clead, master_directory_forecast_spline)
    
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
            
        # ---- initialize netCDF file for full probabilistic fields
        
        outfile_fullfield = master_directory_fullfield_output + ccmonth + cyyyy + \
            '_lead'+clead+'_probabilities_all_fullfield.nc'
        ncname = 'nc_fullfield'
        istat = initialize_netCDF_probfiles(outfile_fullfield, ncname, \
            ny_ndfd, nx_ndfd, lats_ndfd, lons_ndfd, nthresholds, \
            thresholds)
            
        # ---- set up full output grids and working arrays

        precip_ens_raw_ndfd = np.zeros((nmembers, ny_ndfd, nx_ndfd), \
            dtype=np.float32)
        precip_ens_qmapped_ndfd = np.zeros((nmembers*25, ny_ndfd, nx_ndfd), \
            dtype=np.float32)
        prob = np.zeros((nthresholds, ny_ndfd, nx_ndfd), dtype=np.float64)
        prob_raw_out = np.zeros((nthresholds, ny_ndfd, nx_ndfd), \
            dtype=np.float64)
        prob_qmapped_out = np.zeros((nthresholds, ny_ndfd, nx_ndfd), \
            dtype=np.float64)
        prob_qmapped_weighted_out = np.zeros((nthresholds, ny_ndfd, nx_ndfd), \
            dtype=np.float64)
        prob_qmapped_weighted_dressed_out = np.zeros((nthresholds, ny_ndfd, nx_ndfd), \
            dtype=np.float64)
        
        # ---- Take every nth pixel in lat/lon to subsample starting at
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
                np.zeros((nmembers*25, ny_thinned, nx_thinned), dtype=np.float32)
            outfile_thinned = master_directory_thinned_output + ccmonth + cyyyy + \
                '_lead'+clead+'_probabilities_all_thinned.nc'   
            ncname = 'nc_thinned'
            istat = initialize_netCDF_probfiles(outfile_thinned, ncname, \
                ny_thinned, nx_thinned, lats_thinned, lons_thinned, \
                nthresholds, thresholds)    
            prob_thinned = np.zeros((nthresholds, \
                ny_thinned, nx_thinned), dtype=np.float32)            
            prob_raw_thinned_out = np.zeros((nthresholds, \
                ny_thinned, nx_thinned), dtype=np.float32)    
            prob_qmapped_thinned_out = np.zeros((nthresholds, \
                ny_thinned, nx_thinned), dtype=np.float32)
            prob_qmapped_weighted_thinned_out = np.zeros((nthresholds, \
                ny_thinned, nx_thinned), dtype=np.float32)
            prob_qmapped_weighted_dressed_thinned_out = np.zeros((nthresholds, \
                ny_thinned, nx_thinned), dtype=np.float32)
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
            get_domain_subset_of_gefs(gefs_data_directory, cyyyymmddhh, cmem, clead) 

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
        a = np.where(fraction_zero_ndfd < 0)
        aunraveled = np.unravel_index(a,shape=(ny_ndfd,nx_ndfd))
        qmapped_precip = np.zeros((25,ny_ndfd, nx_ndfd), dtype=np.float64)    
        qmapped_precip = \
            qmapping_spline_fivexfive_flexiknot(gefsv12_quantiles_on_ndfd, \
    	    precip_gefsv12_on_ndfd, spline_info_inv, fraction_zero_gefsv12_on_ndfd, \
            fraction_zero_ndfd, usegamma_ndfd, use99, offset_on_ndfd, \
            number_knots_ndfd, nstride_qmap, ny_ndfd, nx_ndfd)
            
        print ('   max qmapped_precip = ',  np.max(qmapped_precip))
        
        # ---- store these quantile-mapped members
        
        istart = imem*25
        for istencil in range(25):
            precip_ens_qmapped_ndfd[istart+istencil,:,:] = qmapped_precip[istencil,:,:]    

        # ---- if subsampling every 10th pixel, perform that now and save.
        
        if netcdf_thinned == True:
            precip_ens_raw_thinned[imem,:,:] = precip_ens_raw_ndfd[imem,5::nstride, 5::nstride]
            print ('   performing thinning')
            istart = imem*25
            for istencil in range(25):
                precip_ens_qmapped_thinned[istart+istencil,:,:] = \
                    qmapped_precip[istencil,5::nstride, 5::nstride]
     
    # ========================================================================        
    # ---- generate probabilities, raw and quantile mapped for each
    #      specified precipitation threshold, and for full field,
    #      thinned, and upscaled output
    # ========================================================================
    
    print ('before netcdf_fullfield, cyyyymmddhh_fcst = ', cyyyymmddhh_fcst)
    if netcdf_fullfield == True:
        ensmean_raw_ndfd = np.squeeze(np.mean(precip_ens_raw_ndfd, axis=0))
        ensmean_qmapped_ndfd = np.squeeze(np.mean(precip_ens_qmapped_ndfd, axis=0))
        print ('---- generating probabilities for full field')
        
        # ---- generate probabilities
        
        isitunweighted = True
        dressit = False
        print ('   thresholds = ', thresholds)
    
        # --- raw
        print ('   generating raw probability')
        prob_raw_out = \
            generate_probabilities(nmembers, thresholds, \
            ensmean_raw_ndfd, precip_ens_raw_ndfd, \
            prob, ny_ndfd, nx_ndfd, histogram_thresholds, \
            closest_member_histogram, dressit, isitunweighted, \
            b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, \
            b1_std_lowrank, b0_mean_midrank, b1_mean_midrank, \
            b0_std_midrank, b1_std_midrank, b0_mean_highrank, \
            b1_mean_highrank, b0_std_highrank, b1_std_highrank)
            
        # --- quantile mapped, unweighted, undressed
        print ('   generating quantile-mapped, unweighted probability')
        prob_qmapped_out = \
            generate_probabilities(nmembers*25, thresholds, \
            ensmean_qmapped_ndfd, precip_ens_qmapped_ndfd, \
            prob, ny_ndfd, nx_ndfd, histogram_thresholds, \
            closest_member_histogram, dressit, isitunweighted, \
            b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, \
            b1_std_lowrank, b0_mean_midrank, b1_mean_midrank, \
            b0_std_midrank, b1_std_midrank, b0_mean_highrank, \
            b1_mean_highrank, b0_std_highrank, b1_std_highrank)  
            
        # --- quantile mapped, weighted, undressed
        isitunweighted = False  
        print ('   generating quantile-mapped, weighted probability')
        prob_qmapped_weighted_out = \
            generate_probabilities(nmembers*25, thresholds, \
            ensmean_qmapped_ndfd, precip_ens_qmapped_ndfd, \
            prob, ny_ndfd, nx_ndfd, histogram_thresholds, \
            closest_member_histogram, dressit, isitunweighted, \
            b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, \
            b1_std_lowrank, b0_mean_midrank, b1_mean_midrank, \
            b0_std_midrank, b1_std_midrank, b0_mean_highrank, \
            b1_mean_highrank, b0_std_highrank, b1_std_highrank) 
            
        # --- quantile mapped, weighted, dressed
        dressit = True
        print ('   generating quantile-mapped, weighted, dressed probability')
        prob_qmapped_weighted_dressed_out = \
            generate_probabilities(nmembers*25, thresholds, \
            ensmean_qmapped_ndfd, precip_ens_qmapped_ndfd, \
            prob, ny_ndfd, nx_ndfd, histogram_thresholds, \
            closest_member_histogram, dressit, isitunweighted, \
            b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, \
            b1_std_lowrank, b0_mean_midrank, b1_mean_midrank, \
            b0_std_midrank, b1_std_midrank, b0_mean_highrank, \
            b1_mean_highrank, b0_std_highrank, b1_std_highrank)
    
        print ('   min, max, mean prob raw fullfield = ',\
            np.min(prob_raw_out), \
            np.max(prob_raw_out), \
            np.mean(prob_raw_out))
        print ('   min, max, mean prob qmapped fullfield unweighted, undressed = ',\
            np.min(prob_qmapped_out), \
            np.max(prob_qmapped_out), \
            np.mean(prob_qmapped_out))
        print ('   min, max, mean prob qmapped fullfield weighted, undressed = ',\
            np.min(prob_qmapped_weighted_out), \
            np.max(prob_qmapped_weighted_out), \
            np.mean(prob_qmapped_weighted_out))
        print ('   min, max, mean prob qmapped fullfield weighted, dressed = ',\
            np.min(prob_qmapped_weighted_dressed_out), \
            np.max(prob_qmapped_weighted_dressed_out), \
            np.mean(prob_qmapped_weighted_dressed_out))
    
    print ('before netcdf_thinned, cyyyymmddhh_fcst = ', cyyyymmddhh_fcst)
    if netcdf_thinned == True:
        ensmean_raw_thinned = np.squeeze(np.mean(\
            precip_ens_raw_thinned, axis=0))
        ensmean_qmapped_thinned = np.squeeze(np.mean(\
            precip_ens_qmapped_thinned, axis=0))
        print ('---- generating probabilities for thinned') 
        isitunweighted = True
        dressit = False
        
        # --- raw
        print ('   generating raw probability')
        prob_raw_thinned_out = \
            generate_probabilities(nmembers, thresholds, \
            ensmean_raw_thinned, precip_ens_raw_thinned, \
            prob_thinned, ny_thinned, nx_thinned, \
            histogram_thresholds, closest_member_histogram, \
            dressit, isitunweighted, \
            b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, \
            b1_std_lowrank, b0_mean_midrank, b1_mean_midrank, \
            b0_std_midrank, b1_std_midrank, b0_mean_highrank, \
            b1_mean_highrank, b0_std_highrank, b1_std_highrank)
            
        # --- quantile mapped, unweighted, undressed
        print ('   generating quantile-mapped, unweighted probability')
        prob_qmapped_thinned_out = \
            generate_probabilities(nmembers*25, thresholds, \
            ensmean_qmapped_thinned, precip_ens_qmapped_thinned, \
            prob_thinned, ny_thinned, nx_thinned, \
            histogram_thresholds, closest_member_histogram, \
            dressit, isitunweighted, \
            b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, \
            b1_std_lowrank, b0_mean_midrank, b1_mean_midrank, \
            b0_std_midrank, b1_std_midrank, b0_mean_highrank, \
            b1_mean_highrank, b0_std_highrank, b1_std_highrank)  

        # --- quantile mapped, weighted, undressed
        
        isitunweighted = False
        print ('   generating quantile-mapped, weighted probability')
        prob_qmapped_weighted_thinned_out = \
            generate_probabilities(nmembers*25, thresholds, \
            ensmean_qmapped_thinned, precip_ens_qmapped_thinned, \
            prob_thinned, ny_thinned, nx_thinned, \
            histogram_thresholds, closest_member_histogram, \
            dressit, isitunweighted, \
            b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, \
            b1_std_lowrank, b0_mean_midrank, b1_mean_midrank, \
            b0_std_midrank, b1_std_midrank, b0_mean_highrank, \
            b1_mean_highrank, b0_std_highrank, b1_std_highrank)         
            
        # --- quantile mapped, weighted, dressed
        
        dressit = True
        print ('   generating quantile-mapped, weighted, dressed probability')
        prob_qmapped_weighted_dressed_thinned_out = \
            generate_probabilities(nmembers*25, thresholds, \
            ensmean_qmapped_thinned, precip_ens_qmapped_thinned, \
            prob_thinned, ny_thinned, nx_thinned, \
            histogram_thresholds, closest_member_histogram, \
            dressit, isitunweighted, \
            b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, \
            b1_std_lowrank, b0_mean_midrank, b1_mean_midrank, \
            b0_std_midrank, b1_std_midrank, b0_mean_highrank, \
            b1_mean_highrank, b0_std_highrank, b1_std_highrank)
            
        print ('   min, max, mean prob raw thinned = ',\
            np.min(prob_raw_thinned_out), \
            np.max(prob_raw_thinned_out), \
            np.mean(prob_raw_thinned_out))
        print ('   min, max, mean prob qmapped thinned unweighted, undressed = ',\
            np.min(prob_qmapped_thinned_out), \
            np.max(prob_qmapped_thinned_out), \
            np.mean(prob_qmapped_thinned_out))
        print ('   min, max, mean prob qmapped thinned weighted, undressed = ',\
            np.min(prob_qmapped_weighted_thinned_out), \
            np.max(prob_qmapped_weighted_thinned_out), \
            np.mean(prob_qmapped_weighted_thinned_out))
        print ('   min, max, mean prob qmapped thinned weighted, dressed = ',\
            np.min(prob_qmapped_weighted_dressed_thinned_out), \
            np.max(prob_qmapped_weighted_dressed_thinned_out), \
            np.mean(prob_qmapped_weighted_dressed_thinned_out))
            
    # ========================================================================
    # ---- write the probabilities over full set of thresholds to netCDF file.
    #      for the full-field forecasts (every point) if desired
    # ========================================================================

    print ('before netcdf_fullfield write, cyyyymmddhh_fcst = ', cyyyymmddhh_fcst)
    if netcdf_fullfield == True:
        now = datetime.now()
        time = now.strftime("%H:%M:%S")
        ncname = 'nc_fullfield'
        print ('---- writing full field probabilities to file. ',time)
        print ('   cyyyymmddhh_fcst = ', cyyyymmddhh_fcst)
        istat = write_probabilities_to_netcdf(outfile_fullfield, ncname, idate, \
            cyyyymmddhh, cyyyymmddhh_fcst, prob_raw_out, prob_qmapped_out, \
            prob_qmapped_weighted_out, prob_qmapped_weighted_dressed_out, \
            ensmean_raw_ndfd, ensmean_qmapped_ndfd)         
        print ('   writing full field probabilities and ens means completed.')

    # ========================================================================
    # ---- write the probabilities over full set of thresholds to netCDF file.
    #      for the thinned forecasts if desired.
    #      Also read in thinned verification precipitation analysis and
    #      get closest-member histogram statistics, write out to file
    # ========================================================================

    if netcdf_thinned == True:
        now = datetime.now()
        time = now.strftime("%H:%M:%S")
        ncname = 'nc_thinned'
        print ('---- writing thinned field probabilities to file. ',time)
        print ('   cyyyymmddhh_fcst = ', cyyyymmddhh_fcst)
    
        istat = write_probabilities_to_netcdf(outfile_thinned, ncname, idate, \
            cyyyymmddhh, cyyyymmddhh_fcst, prob_raw_thinned_out, \
            prob_qmapped_thinned_out, prob_qmapped_weighted_thinned_out, \
            prob_qmapped_weighted_dressed_thinned_out, \
            ensmean_raw_thinned, ensmean_qmapped_thinned)
        print ('   writing thinned field probabilities and ens means completed.')
    
        
print ('******* DONE! ********') 