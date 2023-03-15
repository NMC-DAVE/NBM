"""
plot_bss_reliability_precip_all.py cseason 
    where cseason == 'cool' or 'warm'

"""

import os, sys
import _pickle as cPickle
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ====================================================================

def calc_BSS_relia(brier_score_overall, brier_score_climo_overall, \
    brier_score_west, brier_score_east, \
    brier_score_climo_west, brier_score_climo_east, \
    brier_score_daily, brier_score_climo_daily, brier_score_gridded, \
    brier_score_climo_gridded, contab_forecast,\
    contab_forecast_daily, ndays_total, nthresh, \
    ncats, ny, nx, ilead, BSS, BSS_05, \
    BSS_95, BSS_west, BSS_east, BSS_gridded,\
    relia, relia_05, relia_95, frequse, \
    yyyymmddhh_init_daily):

    # ---- iterate over precip thresholds.
    
    for ithresh in range(nthresh): 
        
        # ---- first calculate the Brier skill score, overall, west/east, by gridpt
           
        BSS[ilead,ithresh] = 1. - brier_score_overall[ithresh] / \
            brier_score_climo_overall[ithresh]
        BSS_west[ilead,ithresh] = 1. - brier_score_west[ithresh] / \
            brier_score_climo_west[ithresh]
        BSS_east[ilead,ithresh] = 1. - brier_score_east[ithresh] / \
            brier_score_climo_east[ithresh]
            
        for jy in range(ny):
            for ix in range(nx):
                if brier_score_climo_gridded[ithresh,jy,ix] > 0.0:
                    BSS_gridded[ilead,ithresh,jy,ix] = \
                        1.0 - brier_score_gridded[ithresh,jy,ix] / \
                        brier_score_climo_gridded[ithresh,jy,ix]
                else:
                    BSS_gridded[ilead,ithresh,jy,ix] = -99.99
  
        # ---- calculate the frequency of use and reliability if enough samples
        
        nsamps_total = np.sum(contab_forecast[ithresh,:,:])
        for icat in range(ncats):
            frequse[ilead,ithresh,icat] = \
                np.sum(contab_forecast[ithresh,icat,:]) / float(nsamps_total)
            if np.sum(contab_forecast[ithresh,icat,:]) > 50:
                relia[ilead,ithresh,icat] = contab_forecast[ithresh,icat,1] / \
                    np.sum(contab_forecast[ithresh,icat,:])

            
            else:
                relia[ilead,ithresh,icat] = -99.99
        
        # ---- resampled confidence intervals for reliability. 
        
        relia_05[ilead,ithresh,:], relia_95[ilead,ithresh,:] = \
            bootstrap_reliability(contab_forecast_daily, \
            ndays_total, ithresh, ncats)
 
    print ('contab_forecast[nthresh-1,:,0] = ', contab_forecast[nthresh-1,:,0])
    print ('contab_forecast[nthresh-1,:,1] = ', contab_forecast[nthresh-1,:,1])
    print ('relia[nthresh-1,:,1] = ', relia[ilead,nthresh-1,:])

    BSSes = np.zeros((nthresh), dtype=np.float32)
    for idate in range(ndays_total):
        BSSes[:] = 1. -  brier_score_daily[idate,:] / \
            brier_score_climo_daily[idate,:]
        #print ("{0:4},{1:10},{2:8.2}".format(idate, \
        #    yyyymmddhh_init_daily[idate], BSSes[3]))
    
    return BSS, BSS_05, BSS_95, \
        BSS_west, BSS_east, BSS_gridded, \
        relia, relia_05, relia_95, frequse \

# ====================================================================

def bootstrap_reliability(contab_forecast_daily, ndays_total, \
    ithresh, ncats):
    
    relia_05 = np.zeros((ncats), dtype=np.float64)
    relia_95 = np.zeros((ncats), dtype=np.float64)
    contab = contab_forecast_daily[:,ithresh,:,:]
    nresa = 1000
    relia_resamped = np.zeros((nresa, ncats), dtype=np.float64)
    for iresa in range(nresa):
        
        # ---- generate random indices of which days to sample
        
        ix = np.random.uniform(low=0.0, high=float(ncats), \
            size=ndays_total).astype(int)
        
        # ---- using these indices, calculate a resampled reliability
        
        for icat in range(ncats):
            numer = np.sum(contab_forecast_daily[ix,ithresh,icat,1])
            denom = np.sum(contab_forecast_daily[ix,ithresh,icat,:])
            if numer > 100.0:
                relia_resamped[iresa,icat] = \
                    np.sum(contab_forecast_daily[ix,ithresh,icat,1]) / \
                    np.sum(contab_forecast_daily[ix,ithresh,icat,:])
            else:
                relia_resamped[iresa,icat] = -99.99
                
    # ---- get 5th and 95th percentiles and return these
    
    for icat in range(ncats):
        rsamps = relia_resamped[:,icat]
        rsort = np.sort(rsamps)
        rsort_positive = rsort[np.where(rsort > 0.0)]
        nresa_pos = rsort_positive.size
        if nresa_pos > 10:
            relia_05[icat] = rsort_positive[int(.05*float(nresa_pos))]
            relia_95[icat] = rsort_positive[int(.95*float(nresa_pos))]
        else:
            relia_05[icat] = -99.99
            relia_95[icat] = -99.99
        
    return relia_05, relia_95

# ==================================================================== 

def paired_bootstrap(x, y, zclim):
    
    """ given BS vectors x, y, and zclim return a paired 
    bootstrap estimate of the 5th and 95th percentiles
    of the resampled BSS distribution  Follows Hamill, T. M., 
    1999: Hypothesis tests for evaluating numerical 
    precipitation forecasts. Wea. Forecasting, 14, 155-167. 
    """
    import numpy as np

    nresa = 100   # 10000 resamplings
    nelts = np.size(x)  # what is the size of the vector x
    sum1 = np.zeros((nresa),dtype=np.float64)
    sum2 = np.zeros((nresa),dtype=np.float64)
    sumclim = np.zeros((nresa),dtype=np.float64)
    dist = np.zeros((nresa),dtype=np.float64)
    ones = np.ones((nelts),dtype=np.float64)
    zeros = np.zeros((nelts),dtype=np.float64)
    
    for i in range(nresa):
        x0 = np.random.rand(nelts)
        iusex = np.where(x0 < 0.5, ones, zeros)
        iusey = 1.0 - iusex
        sum1[i] = np.sum(x*iusey + y*iusex)
        sum2[i] = np.sum(x*iusex + y*iusey)
        sumclim[i] = np.sum(zclim)
        BSS1 = 1. - sum1[i] / sumclim[i]
        BSS2 = 1. - sum2[i] / sumclim[i]
        dist[i] = BSS1 - BSS2

    dsort = np.sort(dist)
    d05 = dsort[int(.05*float(nresa))]
    d95 = dsort[int(.95*float(nresa))]

    return d05, d95
    
# ====================================================================    

cseason_in = sys.argv[1] # warm, cool

cleads = ['012','024','036','048','060','072',\
    '084','108','120','132','144','156',\
    '168','180','192','204','216','228','240']    
   
    
nleads = len(cleads)
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
cthresholds = ['POP', '> 1 mm', '> 5 mm', '> 10 mm', '> 25 mm']
ccthresholds = ['POP', '1mm', '5mm', '10mm', '25mm']
nthresh = len(thresholds)
ncats = 32
probs = np.arange(ncats)/31.
printit = True

if cseason_in == 'warm':
    cseason = 'Warm season, '
    ccseason = 'warm'
    cyyyymm_list = ['201804','201805','201806','201807','201808','201809',\
        '201904','201905','201906','201907','201908','201909']
    ndaysomo = [30,31,30,31,31,30, 30,31,30,31,31,30]
else:
    cseason = 'Cool season, '
    ccseason = 'cool'
    cyyyymm_list = ['201712','201801','201802','201803','201810',\
        '201811','201812','201901','201902','201903','201910','201911']
    ndaysomo = [31,31,28,31,31,30, 31,31,28,31,31,30]
ndays_total = sum(ndaysomo)

forecast_directory = '/Volumes/NBM/conus_gefsv12/all_thinned/'

# ---- read data needed to form contingency tables, brier scores

for ilead, clead in enumerate(cleads):
    
    print ('processing lead = ',ilead, clead)
    
    # --- set up arrays needed to compute Brier skill scores, reliability
    #     for various test options.

    brier_score_raw_daily = np.zeros((ndays_total, nthresh), dtype=np.float64)
    brier_score_qmapped_daily = np.zeros((ndays_total, nthresh), dtype=np.float64)
    brier_score_qmapped_weighted_daily = np.zeros((ndays_total, nthresh), dtype=np.float64)
    brier_score_qmapped_weighted_dressed_daily = np.zeros((ndays_total, nthresh), dtype=np.float64)
    brier_score_climo_daily = np.zeros((ndays_total, nthresh), dtype=np.float64)
    
    contab_forecast_raw_daily = np.zeros((ndays_total, nthresh, ncats, 2), dtype=np.float64)
    contab_forecast_qmapped_daily = np.zeros((ndays_total, nthresh, ncats, 2), dtype=np.float64)
    contab_forecast_qmapped_weighted_daily = \
        np.zeros((ndays_total, nthresh, ncats, 2), dtype=np.float64)
    contab_forecast_qmapped_weighted_dressed_daily = \
        np.zeros((ndays_total, nthresh, ncats, 2), dtype=np.float64)
    contab_climo_daily = np.zeros((ndays_total, nthresh, ncats, 2), dtype=np.float64)
    
    yyyymmddhh_init_daily = np.zeros((ndays_total), dtype=int)
    
    # ---- loop over months to read in
    
    ktr = 0
    for item, cyyyymm in enumerate(cyyyymm_list):
        
        ktrend = ktr + ndaysomo[item]      
        
        # ---- read raw
        
        infile = forecast_directory + cyyyymm + \
            '_lead' + clead+'_probability_raw_Brier_contingency_table.cPick'
        #print ('reading ', infile)
        inf = open(infile, 'rb')
        brier_score_raw_overall_in = cPickle.load(inf)
        brier_score_raw_west_in = cPickle.load(inf)
        brier_score_raw_east_in = cPickle.load(inf)
        brier_score_raw_daily_in = cPickle.load(inf)
        brier_score_raw_gridded_in = cPickle.load(inf)
        contab_forecast_raw_in = cPickle.load(inf)
        contab_forecast_raw_daily_in = cPickle.load(inf)
        
        brier_score_climo_overall_in = cPickle.load(inf)
        brier_score_climo_west_in = cPickle.load(inf)
        brier_score_climo_east_in = cPickle.load(inf)
        brier_score_climo_daily_in = cPickle.load(inf)
        brier_score_climo_gridded_in = cPickle.load(inf)
        contab_climo_in = cPickle.load(inf)
        contab_climo_daily_in = cPickle.load(inf)
        
        yyyymmddhh_init = cPickle.load(inf)

        # ---- read quantile mapped, unweighted, undressed

        infile = forecast_directory + cyyyymm + \
            '_lead' + clead+'_probability_qmapped_Brier_contingency_table.cPick'
        print ('reading ', infile)
        inf = open(infile, 'rb')
        brier_score_qmapped_overall_in = cPickle.load(inf)
        brier_score_qmapped_west_in = cPickle.load(inf)
        brier_score_qmapped_east_in = cPickle.load(inf)
        brier_score_qmapped_daily_in = cPickle.load(inf)
        brier_score_qmapped_gridded_in = cPickle.load(inf)
        contab_forecast_qmapped_in = cPickle.load(inf)
        contab_forecast_qmapped_daily_in = cPickle.load(inf)
        inf.close()

        # ---- read quantile mapped, weighted, undressed
                
        infile = forecast_directory + cyyyymm + \
            '_lead' + clead+'_probability_qmapped_weighted_Brier_contingency_table.cPick'
        print ('reading ', infile)
        inf = open(infile, 'rb')        
        brier_score_qmapped_weighted_overall_in = cPickle.load(inf)
        brier_score_qmapped_weighted_west_in = cPickle.load(inf)
        brier_score_qmapped_weighted_east_in = cPickle.load(inf)
        brier_score_qmapped_weighted_daily_in = cPickle.load(inf)
        brier_score_qmapped_weighted_gridded_in = cPickle.load(inf)
        contab_forecast_qmapped_weighted_in = cPickle.load(inf)
        contab_forecast_qmapped_weighted_daily_in = cPickle.load(inf)
        inf.close()
        
        # ---- read quantile mapped, weighted, dressed
                
        infile = forecast_directory + cyyyymm + \
            '_lead' + clead+'_probability_qmapped_weighted_dressed_Brier_contingency_table.cPick'
        print ('reading ', infile)
        inf = open(infile, 'rb')        
        brier_score_qmapped_weighted_dressed_overall_in = cPickle.load(inf)
        brier_score_qmapped_weighted_dressed_west_in = cPickle.load(inf)
        brier_score_qmapped_weighted_dressed_east_in = cPickle.load(inf)
        brier_score_qmapped_weighted_dressed_daily_in = cPickle.load(inf)
        brier_score_qmapped_weighted_dressed_gridded_in = cPickle.load(inf)
        contab_forecast_qmapped_weighted_dressed_in = cPickle.load(inf)
        contab_forecast_qmapped_weighted_dressed_daily_in = cPickle.load(inf)
        inf.close()        
        
        if item == 0:
            
            brier_score_raw_overall = brier_score_raw_overall_in
            brier_score_raw_west = brier_score_raw_west_in
            brier_score_raw_east = brier_score_raw_east_in
            brier_score_raw_daily[ktr:ktrend,:] = \
                np.transpose(brier_score_raw_daily_in[:,:])
            brier_score_raw_gridded = brier_score_raw_gridded_in
            contab_forecast_raw = contab_forecast_raw_in
            contab_forecast_raw_daily[ktr:ktrend,:,:,:] = \
                contab_forecast_raw_daily_in[:,:,:,:]

            brier_score_qmapped_overall = brier_score_qmapped_overall_in
            brier_score_qmapped_west = brier_score_qmapped_west_in
            brier_score_qmapped_east = brier_score_qmapped_east_in
            brier_score_qmapped_daily[ktr:ktrend,:] = \
                np.transpose(brier_score_qmapped_daily_in[:,:])
            brier_score_qmapped_gridded = brier_score_qmapped_gridded_in
            contab_forecast_qmapped = contab_forecast_qmapped_in
            contab_forecast_qmapped_daily[ktr:ktrend,:,:,:] = \
                contab_forecast_qmapped_daily_in[:,:,:,:]
                
            brier_score_qmapped_weighted_overall = \
                brier_score_qmapped_weighted_overall_in
            brier_score_qmapped_weighted_west = brier_score_qmapped_weighted_west_in
            brier_score_qmapped_weighted_east = brier_score_qmapped_weighted_east_in
            brier_score_qmapped_weighted_daily[ktr:ktrend,:] = \
                np.transpose(brier_score_qmapped_weighted_daily_in[:,:])
            brier_score_qmapped_weighted_gridded = brier_score_qmapped_weighted_gridded_in
            contab_forecast_qmapped_weighted = contab_forecast_qmapped_weighted_in
            contab_forecast_qmapped_weighted_daily[ktr:ktrend,:,:,:] = \
                contab_forecast_qmapped_weighted_daily_in[:,:,:,:]
                
            brier_score_qmapped_weighted_dressed_overall = \
                brier_score_qmapped_weighted_dressed_overall_in
            brier_score_qmapped_weighted_dressed_west = \
                brier_score_qmapped_weighted_dressed_west_in
            brier_score_qmapped_weighted_dressed_east = \
                brier_score_qmapped_weighted_dressed_east_in
            brier_score_qmapped_weighted_dressed_daily[ktr:ktrend,:] = \
                np.transpose(brier_score_qmapped_weighted_dressed_daily_in[:,:])
            brier_score_qmapped_weighted_dressed_gridded = \
                brier_score_qmapped_weighted_dressed_gridded_in
            contab_forecast_qmapped_weighted_dressed = \
                contab_forecast_qmapped_weighted_dressed_in
            contab_forecast_qmapped_weighted_dressed_daily[ktr:ktrend,:,:,:] = \
                contab_forecast_qmapped_weighted_dressed_daily_in[:,:,:,:]

            brier_score_climo_overall = brier_score_climo_overall_in
            brier_score_climo_west = brier_score_climo_west_in
            brier_score_climo_east = brier_score_climo_east_in
            brier_score_climo_daily[ktr:ktrend,:] = \
                np.transpose(brier_score_climo_daily_in[:,:])
            brier_score_climo_gridded = brier_score_climo_gridded_in
            contab_climo = contab_climo_in
            contab_climo_daily[ktr:ktrend,:,:,:] = contab_climo_daily_in[:,:,:,:]
            
            yyyymmddhh_init_daily[ktr:ktrend] = yyyymmddhh_init[:]
            
            nt, ny, nx = np.shape(brier_score_raw_gridded)
            
            # ---- set up variables to store final output variables.
            
            BSS_raw = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_raw_05 = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_raw_95 = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_raw_west = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_raw_east = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_raw_gridded = np.zeros((nleads, nthresh, ny, nx), dtype=np.float64)

            relia_raw = np.zeros((nleads,nthresh, ncats), dtype=np.float64)
            relia_raw_05 = np.zeros((nleads, nthresh,ncats), dtype=np.float64)
            relia_raw_95 = np.zeros((nleads, nthresh,ncats), dtype=np.float64)
            frequse_raw = np.zeros((nleads, nthresh,ncats), dtype=np.float64)
            
            BSS_qmapped = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_05 = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_95 = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_west = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_east = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_gridded = np.zeros((nleads, nthresh, ny, nx), dtype=np.float64)
            
            relia_qmapped = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            relia_qmapped_05 = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            relia_qmapped_95 = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            frequse_qmapped = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            
            BSS_qmapped_weighted = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_05 = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_95 = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_west = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_east = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_gridded = \
                np.zeros((nleads, nthresh, ny, nx), dtype=np.float64)

            relia_qmapped_weighted = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            relia_qmapped_weighted_05 = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            relia_qmapped_weighted_95 = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            frequse_qmapped_weighted = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            
            BSS_qmapped_weighted_dressed = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_dressed_05 = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_dressed_95 = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_dressed_west = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_dressed_east = np.zeros((nleads, nthresh), dtype=np.float64)
            BSS_qmapped_weighted_dressed_gridded = \
                np.zeros((nleads, nthresh, ny, nx), dtype=np.float64)

            relia_qmapped_weighted_dressed = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            relia_qmapped_weighted_dressed_05 = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            relia_qmapped_weighted_dressed_95 = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            frequse_qmapped_weighted_dressed = np.zeros((nleads, nthresh, ncats), dtype=np.float64)
            
        else:
            
            brier_score_raw_overall = brier_score_raw_overall + brier_score_raw_overall_in
            brier_score_raw_west = brier_score_raw_west + brier_score_raw_west_in
            brier_score_raw_east = brier_score_raw_east + brier_score_raw_east_in
            brier_score_raw_daily[ktr:ktrend,:] = np.transpose(brier_score_raw_daily_in[:,:])
            brier_score_raw_gridded = brier_score_raw_gridded + brier_score_raw_gridded_in
            contab_forecast_raw = contab_forecast_raw + contab_forecast_raw_in
            contab_forecast_raw_daily[ktr:ktrend,:,:,:] = \
                contab_forecast_raw_daily_in[:,:,:,:]

            brier_score_qmapped_overall = brier_score_qmapped_overall + \
                brier_score_qmapped_overall_in
            brier_score_qmapped_west = brier_score_qmapped_west + \
                brier_score_qmapped_west_in
            brier_score_qmapped_east = brier_score_qmapped_east + \
                brier_score_qmapped_east_in
            brier_score_qmapped_daily[ktr:ktrend,:] = \
                np.transpose(brier_score_qmapped_daily_in[:,:])
            brier_score_qmapped_gridded = brier_score_qmapped_gridded + \
                brier_score_qmapped_gridded_in
            contab_forecast_qmapped = contab_forecast_qmapped + \
                contab_forecast_qmapped_in
            contab_forecast_qmapped_daily[ktr:ktrend,:,:,:] = \
                contab_forecast_qmapped_daily_in[:,:,:,:]
                
            brier_score_qmapped_weighted_overall = brier_score_qmapped_weighted_overall + \
                brier_score_qmapped_weighted_overall_in
            brier_score_qmapped_weighted_west = brier_score_qmapped_weighted_west + \
                brier_score_qmapped_weighted_west_in
            brier_score_qmapped_weighted_east = brier_score_qmapped_weighted_east + \
                brier_score_qmapped_weighted_east_in
            brier_score_qmapped_weighted_daily[ktr:ktrend,:] = \
                np.transpose(brier_score_qmapped_weighted_daily_in[:,:])
            brier_score_qmapped_weighted_gridded = brier_score_qmapped_weighted_gridded + \
                brier_score_qmapped_weighted_gridded_in
            contab_forecast_qmapped_weighted = contab_forecast_qmapped_weighted + \
                contab_forecast_qmapped_weighted_in
            contab_forecast_qmapped_weighted_daily[ktr:ktrend,:,:,:] = \
                contab_forecast_qmapped_weighted_daily_in[:,:,:,:]
                
            brier_score_qmapped_weighted_dressed_overall = \
                brier_score_qmapped_weighted_dressed_overall + \
                brier_score_qmapped_weighted_dressed_overall_in
            brier_score_qmapped_weighted_dressed_west = \
                brier_score_qmapped_weighted_dressed_west + \
                brier_score_qmapped_weighted_dressed_west_in
            brier_score_qmapped_weighted_dressed_east = \
                brier_score_qmapped_weighted_dressed_east + \
                brier_score_qmapped_weighted_dressed_east_in
            brier_score_qmapped_weighted_dressed_daily[ktr:ktrend,:] = \
                np.transpose(brier_score_qmapped_weighted_dressed_daily_in[:,:])
            brier_score_qmapped_weighted_dressed_gridded = \
                brier_score_qmapped_weighted_dressed_gridded + \
                brier_score_qmapped_weighted_dressed_gridded_in
            contab_forecast_qmapped_weighted_dressed = \
                contab_forecast_qmapped_weighted_dressed + \
                contab_forecast_qmapped_weighted_dressed_in
            contab_forecast_qmapped_weighted_dressed_daily[ktr:ktrend,:,:,:] = \
                contab_forecast_qmapped_weighted_dressed_daily_in[:,:,:,:]

            brier_score_climo_overall = brier_score_climo_overall + brier_score_climo_overall_in
            brier_score_climo_west = brier_score_climo_west + brier_score_climo_west_in
            brier_score_climo_east = brier_score_climo_east + brier_score_climo_east_in
            brier_score_climo_daily[ktr:ktrend,:] = np.transpose(brier_score_climo_daily_in[:,:])
            brier_score_climo_gridded = brier_score_climo_gridded + brier_score_climo_gridded_in
            contab_climo = contab_climo + contab_climo_in
            contab_climo_daily[ktr:ktrend,:,:,:] = contab_climo_daily_in[:,:,:,:]
            
        ktr = ktrend
        
    # --- now process to get raw Brier Skill Scores, reliabilities and 5/95 confidence intervals.
    
    print ('calling calc_BSS_relia for raw')
    BSS_raw, BSS_raw_05, BSS_raw_95, BSS_raw_west, BSS_raw_east, BSS_raw_gridded, \
        relia_raw, relia_05_raw, relia_95_raw, frequse_raw = \
        calc_BSS_relia(brier_score_raw_overall, brier_score_climo_overall, \
        brier_score_raw_west, brier_score_raw_east, brier_score_climo_west, \
        brier_score_climo_east, brier_score_raw_daily, \
        brier_score_climo_daily, brier_score_raw_gridded, \
        brier_score_climo_gridded, contab_forecast_raw, contab_forecast_raw_daily, \
        ndays_total, nthresh, ncats, ny, nx, ilead, BSS_raw, BSS_raw_05, \
        BSS_raw_95, BSS_raw_west, BSS_raw_east, BSS_raw_gridded,\
        relia_raw, relia_raw_05, relia_raw_95, frequse_raw, yyyymmddhh_init_daily)
        
    # --- now process to get quantile mapped Brier Skill Scores, 
    #     reliabilities and 5/95 confidence intervals.
    
    print ('calling calc_BSS_relia for qmapped')
    BSS_qmapped, BSS_qmapped_05, BSS_qmapped_95, \
        BSS_qmapped_west, BSS_qmapped_east, BSS_qmapped_gridded, \
        relia_qmapped, relia_05_qmapped, relia_95_qmapped, frequse_qmapped = \
        calc_BSS_relia(brier_score_qmapped_overall, brier_score_climo_overall, \
        brier_score_qmapped_west, brier_score_qmapped_east, brier_score_climo_west, \
        brier_score_climo_east, brier_score_qmapped_daily, \
        brier_score_climo_daily, brier_score_qmapped_gridded, \
        brier_score_climo_gridded, contab_forecast_qmapped, contab_forecast_qmapped_daily, \
        ndays_total, nthresh, ncats, ny, nx, ilead, BSS_qmapped, BSS_qmapped_05, \
        BSS_qmapped_95, BSS_qmapped_west, BSS_qmapped_east, BSS_qmapped_gridded,\
        relia_qmapped, relia_qmapped_05, relia_qmapped_95, \
        frequse_qmapped, yyyymmddhh_init_daily)
        
    # --- now process to get raw Brier Skill Scores, reliabilities and 5/95 confidence intervals.
    
    print ('calling calc_BSS_relia for qmapped, weighted')
    BSS_qmapped_weighted, BSS_qmapped_weighted_05, BSS_qmapped_weighted_95, \
        BSS_qmapped_weighted_west, BSS_qmapped_weighted_east, BSS_qmapped_weighted_gridded, \
        relia_qmapped_weighted, relia_05_qmapped_weighted, \
        relia_95_qmapped_weighted, frequse_qmapped_weighted = \
        calc_BSS_relia(brier_score_qmapped_weighted_overall, \
        brier_score_climo_overall, brier_score_qmapped_weighted_west, \
        brier_score_qmapped_weighted_east, brier_score_climo_west, \
        brier_score_climo_east, brier_score_qmapped_weighted_daily, \
        brier_score_climo_daily, brier_score_qmapped_weighted_gridded, \
        brier_score_climo_gridded, contab_forecast_qmapped_weighted, \
        contab_forecast_qmapped_weighted_daily, \
        ndays_total, nthresh, ncats, ny, nx, ilead, \
        BSS_qmapped_weighted, BSS_qmapped_weighted_05, \
        BSS_qmapped_weighted_95, BSS_qmapped_weighted_west, \
        BSS_qmapped_weighted_east, BSS_qmapped_weighted_gridded,\
        relia_qmapped_weighted, relia_qmapped_weighted_05, \
        relia_qmapped_weighted_95, frequse_qmapped_weighted, \
        yyyymmddhh_init_daily)
        
    # --- now process to get raw Brier Skill Scores, reliabilities and 5/95 confidence intervals.
    
    print ('calling calc_BSS_relia for qmapped, weighted, dressed')
    BSS_qmapped_weighted_dressed, BSS_qmapped_weighted_dressed_05, \
        BSS_qmapped_weighted_dressed_95, BSS_qmapped_weighted_dressed_west, \
        BSS_qmapped_weighted_dressed_east, BSS_qmapped_weighted_dressed_gridded, \
        relia_qmapped_weighted_dressed, relia_05_qmapped_weighted_dressed, \
        relia_95_qmapped_weighted_dressed, frequse_qmapped_weighted_dressed = \
        calc_BSS_relia(brier_score_qmapped_weighted_dressed_overall, \
        brier_score_climo_overall, brier_score_qmapped_weighted_dressed_west, \
        brier_score_qmapped_weighted_dressed_east, brier_score_climo_west, \
        brier_score_climo_east, brier_score_qmapped_weighted_dressed_daily, \
        brier_score_climo_daily, brier_score_qmapped_weighted_dressed_gridded, \
        brier_score_climo_gridded, contab_forecast_qmapped_weighted_dressed, \
        contab_forecast_qmapped_weighted_dressed_daily, \
        ndays_total, nthresh, ncats, ny, nx, ilead, \
        BSS_qmapped_weighted_dressed, BSS_qmapped_weighted_dressed_05, \
        BSS_qmapped_weighted_dressed_95, BSS_qmapped_weighted_dressed_west, \
        BSS_qmapped_weighted_dressed_east, BSS_qmapped_weighted_dressed_gridded,\
        relia_qmapped_weighted_dressed, relia_qmapped_weighted_dressed_05, \
        relia_qmapped_weighted_dressed_95, frequse_qmapped_weighted_dressed, \
        yyyymmddhh_init_daily)
        
    print (contab_forecast_qmapped_weighted_dressed)

    if printit == True:
        print ('Brier Skill Score (raw) = ', BSS_raw[ilead,:])
        print ('Brier Skill Score (qmapped) = ', BSS_qmapped[ilead,:])
        print ('Brier Skill Score (qmapped, weighted) = ', \
            BSS_qmapped_weighted[ilead,:])
        print ('Brier Skill Score (qmapped, weighted, dressed) = ', \
            BSS_qmapped_weighted_dressed[ilead,:])
        print ('Brier Skill Score 05 (raw) = ', \
            BSS_raw_05[ilead,:])
        print ('Brier Skill Score 95 (raw) = ', \
            BSS_raw_95[ilead,:])
        print ('Brier Skill Score 05 (qmapped) = ', \
            BSS_qmapped_05[ilead,:])
        print ('Brier Skill Score 95 (qmapped) = ', \
            BSS_qmapped_95[ilead,:])
        print ('Brier Skill Score 05 (qmapped, weighted) = ', \
            BSS_qmapped_weighted_05[ilead,:])
        print ('Brier Skill Score 95 (qmapped, weighted) = ', \
            BSS_qmapped_weighted_95[ilead,:])
        print ('Brier Skill Score 05 (qmapped, weighted, dressed) = ', \
            BSS_qmapped_weighted_dressed_05[ilead,:])
        print ('Brier Skill Score 95 (qmapped, weighted, dressed) = ', \
            BSS_qmapped_weighted_dressed_95[ilead,:])
        print ('Brier Skill Score (raw west) = ', \
            BSS_raw_west[ilead,:])
        print ('Brier Skill Score (qmapped west) = ', \
            BSS_qmapped_west[ilead,:])
        print ('Brier Skill Score (qmapped, weighted west) = ', \
            BSS_qmapped_weighted_west[ilead,:])
        print ('Brier Skill Score (qmapped, weighted, dressed west) = ', \
            BSS_qmapped_weighted_dressed_west[ilead,:])
        print ('Brier Skill Score (raw east) = ', \
            BSS_raw_east[ilead,:])
        print ('Brier Skill Score (qmapped east) = ', \
            BSS_qmapped_east[ilead,:])
        print ('Brier Skill Score (qmapped, weighted east) = ', \
            BSS_qmapped_weighted_east[ilead,:])
        print ('Brier Skill Score (qmapped, weighted, dressed east) = ', \
            BSS_qmapped_weighted_dressed_east[ilead,:])
        print ('Brier Skill Score (raw gridded 1 mm) = ', \
            BSS_raw_gridded[ilead,1,ny//2,0:nx:10])
        print ('Brier Skill Score (qmapped gridded 1 mm) = ', \
            BSS_qmapped_gridded[ilead,1,ny//2,0:nx:10])
        print ('Brier Skill Score (qmapped, weighted gridded 1 mm) = ', \
            BSS_qmapped_weighted_gridded[ilead,1,ny//2,0:nx:10])
        print ('Brier Skill Score (qmapped, weighted, dressed gridded 1 mm) = ', \
            BSS_qmapped_weighted_dressed_gridded[ilead,1,ny//2,0:nx:10])
    
    for ithresh in range(nthresh):
        
        # ---- paired block bootstrap to get Brier Skill Score confidence intervals
        #      for overall.
    
        x = brier_score_raw_daily[:,ithresh]
        y = brier_score_qmapped_daily[:,ithresh]
        zclim = brier_score_climo_daily[:,ithresh]
        q05, q95 = paired_bootstrap(x, y, zclim)
        BSS_raw_05[ilead,ithresh] = BSS_raw[ilead,ithresh] + q05
        BSS_raw_95[ilead,ithresh] = BSS_raw[ilead,ithresh] + q95
        BSS_qmapped_05[ilead,ithresh] = BSS_qmapped[ilead,ithresh] + q05
        BSS_qmapped_95[ilead,ithresh] = BSS_qmapped[ilead,ithresh] + q95

        x = brier_score_qmapped_daily[:,ithresh]
        y = brier_score_qmapped_weighted_daily[:,ithresh]
        zclim = brier_score_climo_daily[:,ithresh]
        q05, q95 = paired_bootstrap(x, y, zclim)
        BSS_qmapped_weighted_05[ilead,ithresh] = \
            BSS_qmapped_weighted[ilead,ithresh] + q05
        BSS_qmapped_weighted_95[ilead,ithresh] = \
            BSS_qmapped_weighted[ilead,ithresh] + q95
        
        x = brier_score_qmapped_weighted_daily[:,ithresh]
        y = brier_score_qmapped_weighted_dressed_daily[:,ithresh]
        zclim = brier_score_climo_daily[:,ithresh]
        q05, q95 = paired_bootstrap(x, y, zclim)
        BSS_qmapped_weighted_dressed_05[ilead,ithresh] = \
            BSS_qmapped_weighted_dressed[ilead,ithresh] + q05
        BSS_qmapped_weighted_dressed_95[ilead,ithresh] = \
            BSS_qmapped_weighted_dressed[ilead,ithresh] + q95
            
        # -- make reliability diagram and frequency of usage
        
        ctitle = cseason+' '+cthresholds[ithresh]+', '+cleads[ilead]+' h lead'
        fig = plt.figure(figsize=(6.,6.))
        a1 = fig.add_axes([.12,.1,.83,.82])
        a1.set_title(ctitle,fontsize=15)
        
        rcParams['legend.fontsize']='xx-small'
        rcParams['legend.fancybox']=True
        rcParams['xtick.labelsize']='medium'
        rcParams['ytick.labelsize']='medium'

        # --- add basic reliability diagram

        relia_m_raw = ma.array(relia_raw[ilead,ithresh,:])
        relia_m_raw = ma.masked_where(relia_m_raw < 0.0, relia_m_raw)
        relia_m_raw_error = \
            ma.array((relia_raw_95[ilead,ithresh,:]-relia_raw_05[ilead,ithresh,:])/2.)
        relia_m_raw_error = ma.masked_where(relia_m_raw_error <= 0.0, relia_m_raw_error)
        
        relia_m_qmapped = ma.array(relia_qmapped[ilead,ithresh,:])
        relia_m_qmapped = ma.masked_where(relia_m_qmapped < 0.0, relia_m_qmapped)
        relia_m_qmapped_error = \
            ma.array((relia_qmapped_95[ilead,ithresh,:]- \
            relia_qmapped_05[ilead,ithresh,:])/2.)
        relia_m_qmapped_error = \
            ma.masked_where(relia_m_qmapped_error <= 0.0, relia_m_qmapped_error)
        
        relia_m_qmapped_weighted = \
            ma.array(relia_qmapped_weighted[ilead,ithresh,:])
        relia_m_qmapped_weighted = \
            ma.masked_where(relia_m_qmapped_weighted < 0.0, \
            relia_m_qmapped_weighted)
        relia_m_qmapped_weighted_error = \
            ma.array((relia_qmapped_weighted_95[ilead,ithresh,:]- \
            relia_qmapped_weighted_05[ilead,ithresh,:])/2.)
        relia_m_qmapped_weighted_error = \
            ma.masked_where(relia_m_qmapped_weighted_error <= 0.0, \
            relia_m_qmapped_weighted_error)
            
        relia_m_qmapped_weighted_dressed = \
            ma.array(relia_qmapped_weighted_dressed[ilead,ithresh,:])
        relia_m_qmapped_weighted_dressed = \
            ma.masked_where(relia_m_qmapped_weighted_dressed < 0.0, \
            relia_m_qmapped_weighted_dressed)
        relia_m_qmapped_weighted_dressed_error = \
            ma.array((relia_qmapped_weighted_dressed_95[ilead,ithresh,:]- \
            relia_qmapped_weighted_dressed_05[ilead,ithresh,:])/2.)
        relia_m_qmapped_weighted_dressed_error = \
            ma.masked_where(relia_m_qmapped_weighted_dressed_error <= 0.0, \
            relia_m_qmapped_weighted_dressed_error)

        probs = np.arange(ncats) * 100./np.real(ncats-1)
        
        strbss_raw = 'BSS = %0.2f' %(BSS_raw[ilead,ithresh])
        strbss_qmapped = 'BSS = %0.2f' %(BSS_qmapped[ilead,ithresh])
        strbss_qmapped_weighted = 'BSS = %0.2f' %(BSS_qmapped_weighted[ilead,ithresh])
        strbss_qmapped_weighted_dressed = \
            'BSS = %0.2f' %(BSS_qmapped_weighted_dressed[ilead,ithresh])
        #a1.text(51,9,strbss_raw,fontsize=16,color='Red')
        #a1.text(51,3,strbss_qmapped,fontsize=16,color='Blue')
        
        a1.plot(probs,100.*relia_m_raw,'o-',\
            color='Red',linewidth=2,label='Raw, '+strbss_raw)
        a1.plot(probs,100.*relia_m_qmapped,'o-',\
            color='RoyalBlue',linewidth=1.5,label='Quantile mapped, '+strbss_qmapped)
        a1.plot(probs,100.*relia_m_qmapped_weighted,'o-',color='LimeGreen',\
            linewidth=1.5,label='Weighted, quantile mapped, '+strbss_qmapped_weighted)
        a1.plot(probs,100.*relia_m_qmapped_weighted_dressed,\
            'o-',color='Orange',linewidth=1.5,\
            label='Weighted, dressed, quantile mapped, '+strbss_qmapped_weighted)
        a1.errorbar(probs,100.*relia_m_raw,\
            100.*relia_m_raw_error,elinewidth=1,color='Red')
        a1.errorbar(probs,100.*relia_m_qmapped,\
            100.*relia_m_qmapped_error,elinewidth=1,color='RoyalBlue')
        a1.errorbar(probs,100.*relia_m_qmapped_weighted,\
            100.*relia_m_qmapped_weighted_error,elinewidth=1,color='LimeGreen')
        a1.errorbar(probs,100.*relia_m_qmapped_weighted_dressed,\
            100.*relia_m_qmapped_weighted_dressed_error,elinewidth=1,color='Orange')
        a1.legend(loc=4, fontsize=10)
        
        a1.plot([0,100],[0,100],'--',color='k')
        a1.set_ylabel('Observed Relative Frequency (%)',fontsize=14)
        a1.set_xlabel('Forecast Probability (%)',fontsize=14)
        a1.set_ylim(-1,101)
        a1.set_xlim(-1,101)

        # --- Frequency of usage inset diagram
        
        a2 = fig.add_axes([.24,.66,.4,.2])
        a2.bar(probs-0.9,frequse_raw[ilead,ithresh,:],width=0.4,bottom=0.0001,\
            log=True,color='red',edgecolor='None',align='center')
        a2.bar(probs-0.3,frequse_qmapped[ilead,ithresh,:],width=0.4,bottom=0.0001,\
            log=True,color='RoyalBlue',edgecolor='None',align='center')
        a2.bar(probs+0.3,frequse_qmapped_weighted[ilead,ithresh,:],width=0.4,bottom=0.0001,\
            log=True,color='LimeGreen',edgecolor='None',align='center')
        a2.bar(probs+0.9,frequse_qmapped_weighted[ilead,ithresh,:],width=0.4,bottom=0.0001,\
            log=True,color='Orange',edgecolor='None',align='center')

        a2.set_xlim(-3,103)
        a2.set_ylim(0.0001,1.)
        a2.set_title('Frequency of usage',fontsize=11)
        a2.set_xlabel('Forecast probability',fontsize=10)
        a2.set_ylabel('Forecast frequency',fontsize=10)
        a2.hlines([.01,.1],0,100,linestyles='dashed',colors='gray',lw=0.5)

        plot_title = 'relia_all_'+ccseason+'_'+ccthresholds[ithresh]+'_'+\
            cleads[ilead]+'h.png'
        print ('saving plot to file = ',plot_title)
        plt.savefig(plot_title, dpi=400)
        print ('Plot done')
