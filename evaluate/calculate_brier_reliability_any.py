""" calculate_brier_reliability_any.py cyyyymm clead ctype
"""

import os, sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import _pickle as cPickle

# ============================================================

def compute_brier_score_and_contingency_tables( idate, ithresh, \
    brier_score_overall, brier_score_overall_west, \
    brier_score_overall_east, brier_score_daily, \
    brier_score_gridded, prob, contab, contab_daily, \
    precip_anal, thresh, ones, zeros, weight):

    # ---- compute addition to the Brier Score for this case to that
    #      tallied for other case.   Populate the contingency tables
    #      that are used to calculate the reliability and frequency
    #      of usage.   

    ny, nx = np.shape(precip_anal)
    binary_anal = np.where(precip_anal > thresh, ones, zeros)
    binary_anal[np.where(precip_anal < 0)] = -1
    brier_score_overall[ithresh] = brier_score_overall[ithresh] + \
        np.sum(weight*(binary_anal-prob)**2)
    brier_score_overall_west[ithresh] = brier_score_overall_west[ithresh] + \
        np.sum(weight[:,0:3*nx//7]*(binary_anal[:,0:3*nx//7]- \
        prob[:,0:3*nx//7])**2)
    brier_score_overall_east[ithresh] = brier_score_overall_east[ithresh] + \
        np.sum(weight[:,3*nx//7:-1]*(binary_anal[:,3*nx//7:-1]- \
        prob[:,3*nx//7:-1])**2)
    brier_score_gridded[:,:] = brier_score_gridded[:,:] + \
        weight[:,:]*(binary_anal[:,:]-prob[:,:])**2            
    brier_score_daily[ithresh,idate] = \
        brier_score_daily[ithresh,idate] + \
        np.sum(weight*(binary_anal-prob)**2)
    
    # ---- compute increment to contingency table array
    
    for icat in range(32):
        
        # ---- saved for 31 categories given 31 ens mbrs.
        
        pmin = np.max([0.0,float(icat)/31. - 1./62.])
        pmax = np.min([1.0,float(icat)/31. + 1./62.])
        
        a = np.where(np.logical_and(np.logical_and( \
            prob >= pmin, prob < pmax), binary_anal == 1)  )
        if len(a) > 0:  # a[0] != -1:
            contab[ithresh,icat,1] = \
                contab[ithresh,icat,1] + \
                np.sum(weight[a])
            contab_daily[idate,ithresh,icat,1] = \
                contab_daily[idate,ithresh,icat,1] + \
                np.sum(weight[a])
        a = np.where(np.logical_and(np.logical_and( \
            prob >= pmin, prob < pmax), binary_anal == 0)  )
        if len(a) > 0:   # a[0] != -1:
            contab[ithresh,icat,0] = \
                contab[ithresh,icat,0] + \
                np.sum(weight[a])
            contab_daily[idate,ithresh,icat,0] = \
                contab_daily[idate,ithresh,icat,0] + \
                np.sum(weight[a])
    
    return brier_score_overall, brier_score_overall_west, \
        brier_score_overall_east,  brier_score_daily, \
        brier_score_gridded, contab, contab_daily

# ============================================================

# ---- get the month and end time from the commmand line.  The first 00
#      hour analysis of the month will need to access the data from
#      the previous month.

cyyyymm = sys.argv[1] # 202001 etc
clead = sys.argv[2] # 018 etc
ctype = sys.argv[3] 
    # probability_raw, probability_qmapped, probability_qmapped_weighted, 
    # probability_qmapped_weighted_dressed
cmo = cyyyymm[4:6]
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
imonth = int(cyyyymm[4:6])-1
cmonth = cmonths[imonth]
cyyyy = cyyyymm[0:4]
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
nthresh = len(thresholds)

# ---- read the probability forecasts for this lead time and month
#      saved by the quantile-mapping routine software.

forecast_directory = '/Volumes/NBM/conus_gefsv12/all_thinned/'
#forecast_directory = '/Volumes/NBM/conus_gefsv12/thinned/'
infile = forecast_directory + cmonth + cyyyy + \
     '_lead'+clead+'_probabilities_all_thinned.nc'
print ('reading ', infile)
nc = Dataset(infile,'r')
yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
yyyymmddhh_fcst = nc.variables['yyyymmddhh_fcst'][:]
probability = nc.variables[ctype][:,:,:,:]
ndates, nthresh, ny_fcst, nx_fcst = np.shape(probability)
nc.close()

# --- set up various working arrays needed

ones = np.ones((ny_fcst, nx_fcst), dtype = np.int32)
zeros = np.zeros((ny_fcst, nx_fcst), dtype = np.int32)

# ---- need to develop and read in the climatological probability 
#      by month, and for thinned data with the appropriate 
#      time of day (00, 06, 12, 18)

ilead = int(clead)
iend_hour = ilead%24 # remainder of division
if iend_hour < 10:
    cend_hour = '0'+str(iend_hour)
else:
    cend_hour = str(iend_hour)
climo_directory = '/Volumes/NBM/conus_panal/'
infile = climo_directory + cmo + '_'+cend_hour+\
    'UTC_ccpa_on_ndfd_grid_6hourly_climo_probability_thinned.cPick'

print ('reading from ', infile) 
inf = open(infile, 'rb')
probability_climo = cPickle.load(inf)
lons_climo = cPickle.load(inf)
lats_climo = cPickle.load(inf)
inf.close()

# --- set up work arrays.

contab_forecast = np.zeros((nthresh,32,2), dtype=np.float64)
contab_climo = np.zeros((nthresh,32,2), dtype=np.float64)
contab_forecast_daily = np.zeros((ndates,nthresh,32,2), dtype=np.float64)
contab_climo_daily = np.zeros((ndates, nthresh,32,2), dtype=np.float64)
brier_score_overall = np.zeros((nthresh), dtype=np.float64)
brier_score_climo_overall = np.zeros((nthresh), dtype=np.float64)
brier_score_west = np.zeros((nthresh), dtype=np.float64)
brier_score_climo_west = np.zeros((nthresh), dtype=np.float64)
brier_score_east = np.zeros((nthresh), dtype=np.float64)
brier_score_climo_east = np.zeros((nthresh), dtype=np.float64)
brier_score_daily = np.zeros((nthresh, ndates), dtype=np.float64)
brier_score_climo_daily = np.zeros((nthresh, ndates), dtype=np.float64)
brier_score_gridded = \
    np.zeros((nthresh, ny_fcst, nx_fcst), dtype=np.float64)
brier_score_climo_gridded = \
    np.zeros((nthresh, ny_fcst, nx_fcst), dtype=np.float64)

for idate, fcst_date in enumerate(yyyymmddhh_fcst):
    
    # ---- read the precipitation analyses for the chosen date

    cyyyymm_anal = str(yyyymmddhh_fcst[idate])[0:6]
    master_directory = '/Volumes/NBM/conus_panal/'
    infile = master_directory + cyyyymm_anal + \
        '_ccpa_on_ndfd_grid_6hourly_thinned.nc'
    nc = Dataset(infile)
    yyyymmddhh_end_in = nc.variables['yyyymmddhh_end'][:]
    if idate == 0:
        conusmask_in = nc.variables['conusmask'][:,:]
        lons_in = nc.variables['lons'][:,:]
        lats_in = nc.variables['lats'][:,:]        
        weight = conusmask_in # np.ones((ny_fcst, nx_fcst), dtype=np.float64)
        #np.where(conusmask_in == 1.0, \
        #    ones*np.cos(lats_in*3.1415926/180.), zeros)

    idx = np.where(yyyymmddhh_end_in == fcst_date)[0]
    precip_anal = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
    ny_anal, nx_anal = np.shape(precip_anal)
    nc.close()    

    # ---- loop thru thresholds, compute scores for raw, qmapped forecast, 
    #      and climatology
    
    for ithresh, thresh in enumerate(thresholds):
        
        prob = probability[idate,ithresh,:,:]
            
        # ---- forecast
        
        brier_score_overall, brier_score_west, brier_score_east, brier_score_daily, \
            brier_score_gridded, contab_forecast, contab_forecast_daily = \
            compute_brier_score_and_contingency_tables( idate, ithresh, \
            brier_score_overall, brier_score_west, brier_score_east, \
            brier_score_daily, brier_score_gridded, prob, contab_forecast, \
            contab_forecast_daily, precip_anal, thresh, ones, zeros, weight)
            
        # --- climatology

        prob = probability_climo[ithresh,:,:]
        brier_score_climo_overall, brier_score_climo_west, \
            brier_score_climo_east, brier_score_climo_daily,\
            brier_score_climo_gridded, contab_climo, contab_climo_daily = \
            compute_brier_score_and_contingency_tables( idate, ithresh, \
            brier_score_climo_overall, brier_score_climo_west, \
            brier_score_climo_east, brier_score_climo_daily, \
            brier_score_climo_gridded, prob, \
            contab_climo, contab_climo_daily, precip_anal, \
            thresh, ones, zeros, weight)
            
# --- save data via cPickle file for later computation of skill scores,  
#     reliability diagrams

print ('************* Scores for ',cyyyymm, clead, ctype,' **************')
print ('BSS_overall = ', 1. - brier_score_overall/brier_score_climo_overall)

outfile = forecast_directory + cyyyymm + \
    '_lead'+clead+'_'+ctype+'_Brier_contingency_table.cPick'
print ('writing to ', outfile)
ouf = open(outfile, 'wb')

cPickle.dump(brier_score_overall, ouf)
cPickle.dump(brier_score_west, ouf)
cPickle.dump(brier_score_east, ouf)
cPickle.dump(brier_score_daily, ouf)
cPickle.dump(brier_score_gridded, ouf)
cPickle.dump(contab_forecast, ouf)
cPickle.dump(contab_forecast_daily, ouf)

cPickle.dump(brier_score_climo_overall, ouf)
cPickle.dump(brier_score_climo_west, ouf)
cPickle.dump(brier_score_climo_east, ouf)
cPickle.dump(brier_score_climo_daily, ouf)
cPickle.dump(brier_score_climo_gridded, ouf)
cPickle.dump(contab_climo, ouf)
cPickle.dump(contab_climo_daily, ouf)

cPickle.dump(yyyymmddhh_init, ouf)
ouf.close()

