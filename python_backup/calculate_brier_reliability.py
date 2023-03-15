""" calculate_brier_reliability.py cyyyymm clead ctype
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
    brier_score_gridded, probability, \
    contab, contab_daily, precip_anal, thresh, ones, zeros, weight):

    # ---- compute addition to the Brier Score for this case to that
    #      tallied for other case.   Populate the contingency tables
    #      that are used to calculate the reliability and frequency
    #      of usage.   

    ny, nx = np.shape(precip_anal)
    binary_anal = np.where(precip_anal > thresh, ones, zeros)
    binary_anal[np.where(precip_anal < 0)] = -1
    brier_score_overall[ithresh] = brier_score_overall[ithresh] + \
        np.sum(weight*(binary_anal-probability)**2)
    brier_score_overall_west[ithresh] = brier_score_overall_west[ithresh] + \
        np.sum(weight[:,0:3*nx//7]*(binary_anal[:,0:3*nx//7]- \
        probability[:,0:3*nx//7])**2)
    brier_score_overall_east[ithresh] = brier_score_overall_east[ithresh] + \
        np.sum(weight[:,3*nx//7:-1]*(binary_anal[:,3*nx//7:-1]- \
        probability[:,3*nx//7:-1])**2)
    brier_score_gridded[:,:] = brier_score_gridded[:,:] + \
        weight[:,:]*(binary_anal[:,:]-probability[:,:])**2            
        
    brier_score_daily[ithresh,idate] = \
        brier_score_daily[ithresh,idate] + \
        np.sum(weight*(binary_anal-probability)**2)
    
    # ---- compute increment to contingency table array
    
    for icat in range(32):
        
        # ---- saved for 31 categories given 31 ens mbrs.
        
        pmin = np.max([0.0,float(icat)/31. - 1./62.])
        pmax = np.min([1.0,float(icat)/31. + 1./62.])
        
        a = np.where(np.logical_and(np.logical_and( \
            probability >= pmin, probability < pmax), binary_anal == 1)  )
        if len(a) > 0:  # a[0] != -1:
            contab[ithresh,icat,1] = \
                contab[ithresh,icat,1] + \
                np.sum(weight[a])
            contab_daily[idate,ithresh,icat,1] = \
                contab_daily[idate,ithresh,icat,1] + \
                np.sum(weight[a])
        a = np.where(np.logical_and(np.logical_and( \
            probability >= pmin, probability < pmax), binary_anal == 0)  )
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
ctype = sys.argv[3] # thinned or upscaled
cmo = cyyyymm[4:6]
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#cmonths = ['Mar']
imonth = int(cyyyymm[4:6])-1
cmonth = cmonths[imonth]
cyyyy = cyyyymm[0:4]
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
nthresh = len(thresholds)

# ---- read the probability forecasts for this lead time and month
#      saved by the quantile-mapping routine software.

forecast_directory = '/Volumes/NBM/conus_gefsv12/'+ctype+'/'
infile = forecast_directory + cmonth + cyyyy + \
     '_use99_lead'+clead+'_probabilities_'+ctype+'.nc'
print ('reading ', infile)
nc = Dataset(infile)
#yyyymmddhh_init_in = nc.variables['yyyymmddhh_init'][0:-1]
yyyymmddhh_init_in = nc.variables['yyyymmddhh_init'][:]
print (yyyymmddhh_init_in)
#yyyymmddhh_fcst_in = nc.variables['yyyymmddhh_fcst'][0:-1]
yyyymmddhh_fcst_in = nc.variables['yyyymmddhh_fcst'][:]
lats_fcst = nc.variables['lats'][:,:]
lons_fcst = nc.variables['lons'][:,:]
#probability_raw = nc.variables['probability_raw'][0:-1,:,:,:]
#yyyymmddhh_init_in = nc.variables['yyyymmddhh_init'][0:-1]
#probability_qmapped = nc.variables['probability_qmapped'][0:-1,:,:,:]
probability_raw = nc.variables['probability_raw'][:,:,:,:]
yyyymmddhh_init_in = nc.variables['yyyymmddhh_init'][:]
probability_qmapped = nc.variables['probability_qmapped'][:,:,:,:]
print ('np.shape(probability_raw) = ', np.shape(probability_raw))
ndates, nthresh, ny_fcst, nx_fcst = np.shape(probability_qmapped )
lons_in = nc.variables['lons'][:,:]
lats_in = nc.variables['lats'][:,:]
nc.close()

# --- set up various working arrays needed

ones = np.ones((ny_fcst, nx_fcst), dtype = np.int32)
zeros = np.zeros((ny_fcst, nx_fcst), dtype = np.int32)

# ---- need to develop and read in the climatological probability 
#      by month

climo_directory = '/Volumes/NBM/conus_panal/'
infile = climo_directory + cmo + \
    '_ccpa_on_ndfd_grid_6hourly_climo_probability_'+ctype+'.cPick'
print ('reading from ', infile) 
inf = open(infile, 'rb')
probability_climo = cPickle.load(inf)
lons_climo = cPickle.load(inf)
lats_climo = cPickle.load(inf)
inf.close()

# --- set up work arrays.

contab_forecast_raw = np.zeros((nthresh,32,2), dtype=np.float64)
contab_forecast_qmapped = np.zeros((nthresh,32,2), dtype=np.float64)
contab_climo = np.zeros((nthresh,32,2), dtype=np.float64)

contab_forecast_raw_daily = np.zeros((ndates,nthresh,32,2), dtype=np.float64)
contab_forecast_qmapped_daily = np.zeros((ndates,nthresh, 32,2), dtype=np.float64)
contab_climo_daily = np.zeros((ndates, nthresh,32,2), dtype=np.float64)

brier_score_raw_overall = np.zeros((nthresh), dtype=np.float64)
brier_score_qmapped_overall = np.zeros((nthresh), dtype=np.float64)
brier_score_climo_overall = np.zeros((nthresh), dtype=np.float64)

brier_score_raw_west = np.zeros((nthresh), dtype=np.float64)
brier_score_qmapped_west = np.zeros((nthresh), dtype=np.float64)
brier_score_climo_west = np.zeros((nthresh), dtype=np.float64)

brier_score_raw_east = np.zeros((nthresh), dtype=np.float64)
brier_score_qmapped_east = np.zeros((nthresh), dtype=np.float64)
brier_score_climo_east = np.zeros((nthresh), dtype=np.float64)

brier_score_raw_daily = np.zeros((nthresh, ndates), dtype=np.float64)
brier_score_qmapped_daily = np.zeros((nthresh, ndates), dtype=np.float64)
brier_score_climo_daily = np.zeros((nthresh, ndates), dtype=np.float64)

brier_score_raw_gridded = \
    np.zeros((nthresh, ny_fcst, nx_fcst), dtype=np.float64)
brier_score_qmapped_gridded = \
    np.zeros((nthresh, ny_fcst, nx_fcst), dtype=np.float64)
brier_score_climo_gridded = \
    np.zeros((nthresh, ny_fcst, nx_fcst), dtype=np.float64)

for idate, fcst_date in enumerate(yyyymmddhh_fcst_in):
    
    # ---- read the precipitation analyses for the chosen date

    cyyyymm_anal = str(yyyymmddhh_fcst_in[idate])[0:6]
    master_directory = '/Volumes/NBM/conus_panal/'
    infile = master_directory + cyyyymm_anal + \
        '_ccpa_on_ndfd_grid_6hourly_'+ctype+'.nc'
    nc = Dataset(infile)
    yyyymmddhh_end_in = nc.variables['yyyymmddhh_end'][:]
    if idate == 0:
        conusmask_in = nc.variables['conusmask'][:,:]
        lons_in = nc.variables['lons'][:,:]
        lats_in = nc.variables['lats'][:,:]
        weight = np.where(conusmask_in == 1.0, \
            ones*np.cos(lats_in*3.1415926/180.), zeros)
    idx = np.where(yyyymmddhh_end_in == fcst_date)[0]
    precip_anal = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
    ny_anal, nx_anal = np.shape(precip_anal)
    nc.close()    

    # ---- loop thru thresholds, compute scores for raw, qmapped forecast, 
    #      and climatology
    
    for ithresh, thresh in enumerate(thresholds):
        print ('thresh = ', thresh)
        
        # --- raw
        
        probability = probability_raw[idate,ithresh,:,:]
        print ('  max, min probability raw = ', np.max(probability))
        brier_score_raw_overall, brier_score_raw_west, \
            brier_score_raw_east, brier_score_raw_daily, \
            brier_score_raw_gridded, contab_forecast_raw, \
            contab_forecast_raw_daily = \
            compute_brier_score_and_contingency_tables( idate, ithresh, \
            brier_score_raw_overall, brier_score_raw_west, \
            brier_score_raw_east, brier_score_raw_daily,\
            brier_score_raw_gridded, probability, \
            contab_forecast_raw, contab_forecast_raw_daily, \
            precip_anal, thresh, ones, zeros, weight)
        
        # --- quantile mapped
        
        probability = probability_qmapped[idate,ithresh,:,:]
        print ('  max, min probability qmapped = ', np.max(probability))
        brier_score_qmapped_overall, brier_score_qmapped_west, \
            brier_score_qmapped_east, brier_score_qmapped_daily, \
            brier_score_qmapped_gridded, contab_forecast_qmapped, \
            contab_forecast_qmapped_daily = \
            compute_brier_score_and_contingency_tables( idate, ithresh, \
            brier_score_qmapped_overall, brier_score_qmapped_west, \
            brier_score_qmapped_east, brier_score_qmapped_daily, \
            brier_score_qmapped_gridded, probability, \
            contab_forecast_qmapped, contab_forecast_qmapped_daily, \
            precip_anal, thresh, ones, zeros, weight)
            
        # --- climatology
        
        probability = probability_climo[ithresh,:,:]
        print ('  max, min probability weighted qmapped = ', np.max(probability))
        brier_score_climo_overall, brier_score_climo_west, \
            brier_score_climo_east, brier_score_climo_daily,\
            brier_score_climo_gridded, contab_climo, contab_climo_daily = \
            compute_brier_score_and_contingency_tables( idate, ithresh, \
            brier_score_climo_overall, brier_score_climo_west, \
            brier_score_climo_east, brier_score_climo_daily, \
            brier_score_climo_gridded, probability, \
            contab_climo, contab_climo_daily, precip_anal, \
            thresh, ones, zeros, weight)            
            

# --- save data via cPickle file for later computation of skill scores,  
#     reliability diagrams

print ('--------------- ', cyyyymm,' ',clead,' h, ',ctype,' -----------------')
print ('BSS_raw_overall = ', 1. - brier_score_raw_overall/brier_score_climo_overall)
print ('BSS_qmapped_overall = ', 1.-brier_score_qmapped_overall/brier_score_climo_overall)
print ('BSS_raw_west = ', 1. - brier_score_raw_west/brier_score_climo_west)
print ('BSS_qmapped_west = ', 1. - brier_score_qmapped_west/brier_score_climo_west)
print ('BSS_raw_east = ', 1. - brier_score_raw_east/brier_score_climo_east)
print ('BSS_qmapped_east = ', 1. - brier_score_qmapped_east/brier_score_climo_east)


outfile = forecast_directory + cyyyymm + \
    '_lead'+clead+'_'+ctype+'_Brier_contingency_table.cPick'
    
ouf = open(outfile, 'wb')

cPickle.dump(brier_score_raw_overall, ouf)
cPickle.dump(brier_score_raw_west, ouf)
cPickle.dump(brier_score_raw_east, ouf)
cPickle.dump(brier_score_raw_daily, ouf)
cPickle.dump(brier_score_raw_gridded, ouf)
cPickle.dump(contab_forecast_raw, ouf)
cPickle.dump(contab_forecast_raw_daily, ouf)

cPickle.dump(brier_score_qmapped_overall, ouf)
cPickle.dump(brier_score_qmapped_west, ouf)
cPickle.dump(brier_score_qmapped_east, ouf)
cPickle.dump(brier_score_qmapped_daily, ouf)
cPickle.dump(brier_score_qmapped_gridded, ouf)
cPickle.dump(contab_forecast_qmapped, ouf)
cPickle.dump(contab_forecast_qmapped_daily, ouf)

cPickle.dump(brier_score_climo_overall, ouf)
cPickle.dump(brier_score_climo_west, ouf)
cPickle.dump(brier_score_climo_east, ouf)
cPickle.dump(brier_score_climo_daily, ouf)
cPickle.dump(brier_score_climo_gridded, ouf)
cPickle.dump(contab_climo, ouf)
cPickle.dump(contab_climo_daily, ouf)

ouf.close()

