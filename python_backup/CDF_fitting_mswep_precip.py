"""
CDF_fitting_mswep_precip.py cmonth cend_hour

this python script is designed to fit a mixture of Gamma CDFs to closely
match an empirical distribution of precipitation.   Up to 3 Gamma distributions
can be mixed, though 1 or 2 is preferable if they provide a sufficiently
close fit.   The script is tailored to MSWEP precipitation analyses over 
one of the National Digital Forecast Database domains for the National 
Blend of Models

Designed by Tom Hamill, NOAA, with help from Michael Scheuerer, Dec 2020

"""

import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
from netCDF4 import Dataset
import _pickle as cPickle
import scipy.stats as stats
from gammamix import gammamix_em

# =====================================================================

def set_domain_boundaries(cdomain):
    
    """ used grib file of 2.5-km blend output grid to determine bounding 
        lat and lon, and from that, the domain bounding indices for the 
        0.25 GEFSv12 reforecast data that will encompass the domain.    
    """
    if cdomain == 'conus': 
        jmin = 93
        jmax = 246
        imin = 368
        imax = 686
    elif cdomain == 'pr':
        jmin = 243
        jmax = 256
        imin = 649
        imax = 667   
    elif cdomain == 'ak':
        jmin = 19
        jmax = 161
        imin = 201
        imax = 967
    else:
        print ('invalid domain.  Exiting.')     
        sys.exit()    
 
    return jmin, jmax, imin, imax
 
# =====================================================================

def find_nearest(vec, value):
    
    """ given a vector vec and a particular value, find the index in vec
    that is nearest to value"""
    
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

def fraczero_possamps(nsamps, precip_ens):
    
    """
    
    from the vector input sample precip_ens, define the fraction of
    samples with zero precipitation.   For the positive samples, add
    a small random number to deal with the fact that the data was 
    discretized to ~0.1 mm, so that when later creating CDFs we don't 
    have empirical values with lots of tied amounts.  Also, sort the 
    nonzero amounts and return.
    
    """
    
    number_zeros = 0
    precip_ens_nonzero = np.delete(precip_ens, \
        np.where(precip_ens <= 0.006))  # censor at 0.006 mm
    nz = len(precip_ens_nonzero)
    # data discretized, so add random component of this magnitude
    precip_ens_nonzero = precip_ens_nonzero + \
        np.random.uniform(low=-0.005,high=0.005,size=nz) 
    precip_ens_nonzero = np.sort(precip_ens_nonzero)  
    ntotal = len(precip_ens)
    nzero = ntotal - len(precip_ens_nonzero)
    fraction_zero = float(nzero) / float(ntotal)
    return fraction_zero, precip_ens_nonzero, nz

# =====================================================================

def one_parameter_gamma(nz, precip_ens_nonzero, pflag, fraction_zero):

    """ 
    Along with some other ancillary stuff, return the fitted Gamma distribution
    alpha and beta values, along with an evaluation (Dnstat) of how 
    closely the fitted CDF matches the empirical CDF.   Dnstat follows
    the Dn statistics discussed in Wilks (2011) Statistical Methods in
    the Atmospheric Sciences text for distribution fitting. 
    
    inputs:  
    
    nz: number of nonzero samples
    precip_ens_nonzero[nz]:   sorted vector of nonzero precip amounts (mm)
    pflag: print flag, true if diagnostic printing desired
    fraction_zero: fraction of samples that are zero
    
    """

    # ---- define the indices in the previously sorted precip_ens_nonzero vector
    #      that are at the 0.05, 0.10, .... , 0.95 quantiles of the sorted vector, + a few 
    #      others at high quantiles.   Less if there are few samples
    
    query_these_indices = [ nz//20, nz//10, (3*nz)//20, nz//5, nz//4, (3*nz)//10, \
        (7*nz)//20, (2*nz)//5, (9*nz)//20, nz//2, (11*nz)//20, (3*nz)//5, (13*nz)//20, \
        (7*nz)//10, (3*nz)//4, (4*nz)//5, (17*nz)//20, (35*nz)//40,(9*nz)//10, \
        (37*nz)//40, (19*nz)//20, (39*nz)//40, (79*nz)//80]
    empirical_precipvals = precip_ens_nonzero[query_these_indices] 
    if pflag == True: print ('   nz, query_these_indices = ', nz, query_these_indices)
    
    # ---- convert the query_these_indices into the cumulative probability
    
    empirical_CDF = np.array(query_these_indices, dtype=np.float) / float(nz) 
    
    if pflag == True: print ('   empirical_CDF = ',empirical_CDF)
    
    # ---- extract the quantiles at these cumulative probabilities
    
    if pflag == True: print ('   empirical_precipvals = ', empirical_precipvals)
    
    # ---- See Wilks Statistical Meteorology Text, section on Gamma Distribution.   Use the
    #      Thom (1958) method of maximum-likelihood estimator, and from this estimate
    #      the Gamma distribution Parameters.   Ref: Statistical Methods in the Atmospheric
    #      Sciences (3rd Ed), 2011, Daniel S. Wilks (Academic Press)
    
    pmean = np.mean(precip_ens_nonzero)
    lnxbar = np.log(pmean)
    meanlnxi = np.mean(np.log(precip_ens_nonzero))
    D = lnxbar - meanlnxi
    alpha_hat = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
    beta_hat = pmean / alpha_hat
    
    # ---- Now we evaluate how good the fitted distribution matches the empirical one.
    #      Because errors in fit are largely irrelevant for the quantile mapping at very small
    #      values, we will evaluate the goodness of fit only at the minimum of: (a) 
    #      where there is precip > 1.0 mm, or (b) the 90th percentile of the empirical 
    #      distribution, whichever is smaller.   
    
    if pflag == True: print ('   D = ',D,' alpha_hat = ', alpha_hat,' beta_hat = ', beta_hat)
    y0 = empirical_precipvals / beta_hat
    fitted_CDF = stats.gamma.cdf(y0, alpha_hat)
    if pflag == True: print ('   fitted_CDF = ', fitted_CDF)
    if pflag == True: print ('   CDF differences: ', np.abs(fitted_CDF - empirical_CDF))
    
    # --- only bother evaluating Dn statistic either for quantiles associated with precip
    #     > 1.0 mm, or above the 35/40th quantile, whichever is smaller.
    
    i1 = np.min([17,np.argmin(np.abs(empirical_precipvals-1.0))])
    Dnstat = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
    print ('   Dnstat for 1 gamma mixture = ', Dnstat)
    weight_save_1m = 1.0
    alpha_save_1m = alpha_hat
    beta_save_1m = beta_hat
    
    return query_these_indices, empirical_CDF, empirical_precipvals, \
        pmean, alpha_hat, beta_hat, i1, Dnstat, weight_save_1m, \
        alpha_save_1m, beta_save_1m
        
# =====================================================================

def two_parameter_gamma(jy, ix, precip_ens_nonzero, empirical_precipvals, \
    empirical_CDF, weights, alpha, beta, i1, pflag, weight_save_2m, \
    alpha_save_2m, beta_save_2m, nmixture, fraction_zero, nstride, \
    Dnstat2, excessive_threshold):
    
    """
    
    Call an R gamma mixture routine to estimate the weights and parameters of
    a mixture of two Gamma distributions.  Evaluate how well it fits the empirical
    data and return weights, parameters, and a fitting statistic.
    
    inputs:
    jy, ix: grid indices
    precip_ens_nonzero: vector of nonzero precip values
    empirical_precipvals : 0.05, 0.1, ... 0.95 empirical quantiles
    empirical_CDF: roughly 0.05, 0.1, ... 0.95 but subject to number of samples
    weights: input and output array of Gamma mixture weights
    alpha: input and output array of Gamma fitted alphas
    beta: input and output array of Gamma fitted betas
    i1: index of minimum quantile to evaluate when determining goodness of fit.
    pflag: true if printing desired
    weight_save_2m: 2d-vector of saved weights from last call
    alpha_save_2m: 2d-vector of saved alpha estimates from last call
    beta_save_2m: 2d-vector of saved beta estimates from last call
    nmixture: indicates how many parameters are used in the mix
    fraction_zero: fraction of points with zero precip
    nstride: how many grid points to step in each direction
    Dnstat2: previous grid point's goodness of fit statistic
    excessive_threshold: threshold value for Dnstat that triggers action
    
    """

    two_parameter_fail = False # in case routine bombs off
    try:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        
        # --- call the R-inspired routine that estimates the weights and a mixture 
        #     of two Gamma distributions. 
        
        if pflag == True: print ('calling gammamix_em 2 parameter')
        if pflag == True: print ('precip_ens_nonzero[0], [-1] = ', \
            precip_ens_nonzero[0], precip_ens_nonzero[-1] )
        if pflag == True: print ('weight_save_2m = ', weight_save_2m)
        if pflag == True: print ('alpha_save_2m = ', alpha_save_2m)
        if pflag == True: print ('beta_save_2m = ', beta_save_2m)
            
        if ix == 0: 
            best_result = gammamix_em(precip_ens_nonzero, \
                mix_prop=None, alpha=None, invbeta=None,\
                k=2, epsilon=0.003, maxit=60, maxrestarts=20, verb=False)
        elif nmixture[jy,ix-nstride] == 2 and Dnstat2 < excessive_threshold*1.5:
            best_result = gammamix_em(precip_ens_nonzero, \
                mix_prop=weight_save_2m, alpha=alpha_save_2m, invbeta=beta_save_2m,\
                k=2, epsilon=0.003, maxit=60, maxrestarts=20, verb=False)
        else: # better to restart with new values
            best_result = gammamix_em(precip_ens_nonzero, \
                mix_prop=None, alpha=None, invbeta=None,\
                k=2, epsilon=0.003, maxit=60, maxrestarts=20, verb=False)
        
        params = best_result.params
        alpha[0:2,jy,ix] = np.array(params.alpha[:])
        beta[0:2,jy,ix] = np.array(params.invbeta[:])
        weights[0:2,jy,ix] = np.array(params.mix_prop[:])
        weights[2,jy,ix] = 0.0
        alpha[2,jy,ix] = 1.0
        beta[2,jy,ix] = 1.0
        nmixture[jy,ix] = 2
        weight_save_2m = np.array(params.mix_prop[:])
        alpha_save_2m = np.array(params.alpha[:])
        beta_save_2m = np.array(params.invbeta[:])
        
        if pflag == True: print ('   weights = ', weights[:,jy,ix])
        if pflag == True: print ('   alpha = ', alpha[:,jy,ix])
        if pflag == True: print ('   beta = ', beta[:,jy,ix])
                
        #  ---- estimate a CDF from a weighted mixture of the two distributions.  
        
        y0 = empirical_precipvals / beta[0,jy,ix]
        y1 = empirical_precipvals / beta[1,jy,ix]
        fitted_CDF0 = stats.gamma.cdf(y0, alpha[0,jy,ix])
        fitted_CDF1 = stats.gamma.cdf(y1, alpha[1,jy,ix])
        fitted_CDF = weights[0,jy,ix]*fitted_CDF0 + weights[1,jy,ix]*fitted_CDF1
        if pflag == True: print ('   fitted_CDF = ',fitted_CDF)
        
        # ---- calculate a statistic for how far off the fitted CDF is from the 
        #      empirical CDF.   Do this only for either the higher quantiles of the
        #      distribution or where the quantile exceeds 1 mm (the i1 index)
        
        Dnstat2 = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        if pflag == True: print ('   np.abs(fitted_CDF - empirical_CDF) = ', \
            np.abs(fitted_CDF - empirical_CDF))
        print ('   Dnstat for 2 gamma mixture = ', Dnstat2)
    except: # something went wrong ...
        print ('   two_parameter_fail')
        two_parameter_fail = True
        Dnstat2 = 0.10

    return two_parameter_fail, weights, alpha, beta, Dnstat2, weight_save_2m, \
        alpha_save_2m, beta_save_2m, nmixture

    
# =====================================================================

def three_parameter_gamma(jy, ix, precip_ens_nonzero, empirical_precipvals, \
    empirical_CDF, weights, alpha, beta, i1, pflag, weight_save_3m, \
    alpha_save_3m, beta_save_3m, nmixture, fraction_zero, nstride, mop_up,\
    Dnstat3, excessive_threshold):    
    
    """
    
    Call an R-inspired gamma mixture routine to estimate the weights and parameters of
    a mixture of three Gamma distributions.  Evaluate how well it fits the empirical
    data and return weights, parameters, and a fitting statistic.
    
    inputs:
    jy, ix: grid indices
    precip_ens_nonzero: vector of nonzero precip values
    empirical_precipvals : 0.05, 0.1, ... 0.95 empirical quantiles
    empirical_CDF: roughly 0.05, 0.1, ... 0.95 but subject to number of samples
    weights: input and output array of Gamma mixture weights
    alpha: input and output array of Gamma fitted alphas [3, ny, nx]
    beta: input and output array of Gamma fitted betas
    i1: index of minimum quantile to evaluate when determining goodness of fit.
    pflag: true if printing desired
    weight_save_3m: first guess for weights from previous call
    alpha_save_3m: first guess for alpha
    beta_save_3m: first guess for beta
    nmixture: indicates how many parameters are used in the mix
    fraction_zero: fraction of points with zero precip
    nstride: how many grid points to step in each direction
    mop_up: flag for whether we're having a 2nd go at things to address points
        where first time through, no good fit.
    Dnstat3: previous grid point's goodness of fit statistic
    excessive_threshold: threshold value for Dnstat that triggers action
    
    """

    three_parameter_fail = False # in case routine bombs off
    try:    
    
        if pflag == True: print ('   Need to try 3 gamma mixture.')
        
        # --- call the R routine that estimates the weights and a mixture 
        #     of 3 Gamma distributions. Must convert back from R data format.
        
        if pflag == True: print ('   calling gammamix_em 3 parameter')
        if pflag == True: print ('   precip_ens_nonzero[0], [-1] = ', \
            precip_ens_nonzero[0], precip_ens_nonzero[-1] )
        if pflag == True: print ('   weight_save_3m = ', weight_save_3m)
        if pflag == True: print ('   alpha_save_3m = ', alpha_save_3m)
        if pflag == True: print ('   beta_save_3m = ', beta_save_3m)
            
        if ix == 0 or mop_up == True:
            print ('   resetting 3-parameter Gamma mixture')
            best_result = gammamix_em(precip_ens_nonzero, \
                mix_prop=None, alpha=None, invbeta=None,\
                k=3, epsilon=0.001, maxit=100, maxrestarts=20, verb=False)
        elif nmixture[jy,ix-nstride] == 3 and Dnstat3 < excessive_threshold*1.5:
            best_result = gammamix_em(precip_ens_nonzero, \
                mix_prop=weight_save_3m, alpha=alpha_save_3m, invbeta=beta_save_3m,\
                k=3, epsilon=0.003, maxit=60, maxrestarts=20, verb=False)
        else: 
            best_result = gammamix_em(precip_ens_nonzero, \
                mix_prop=None, alpha=None, invbeta=None,\
                k=3, epsilon=0.003, maxit=60, maxrestarts=20, verb=False)   
         
        params = best_result.params
        alpha[:,jy,ix] = np.array(params.alpha[:])
        beta[:,jy,ix] = np.array(params.invbeta[:])
        weights[:,jy,ix] = np.array(params.mix_prop[:])
                
        weight_save_3m = np.array(params.mix_prop[:])
        alpha_save_3m = np.array(params.alpha[:])
        beta_save_3m = np.array(params.invbeta[:])

        if pflag == True: print ('   weights = ', weights[:,jy,ix])
        if pflag == True: print ('   alpha = ', alpha[:,jy,ix])
        if pflag == True: print ('   beta = ', beta[:,jy,ix])
        
        #  ---- estimate a CDF from a weighted mixture of the 3 distributions.  
        
        y0 = empirical_precipvals / beta[0,jy,ix]
        y1 = empirical_precipvals / beta[1,jy,ix]
        y2 = empirical_precipvals / beta[2,jy,ix]
        fitted_CDF0 = stats.gamma.cdf(y0, alpha[0,jy,ix])
        fitted_CDF1 = stats.gamma.cdf(y1, alpha[1,jy,ix])
        fitted_CDF2 = stats.gamma.cdf(y2, alpha[2,jy,ix])
        fitted_CDF = weights[0,jy,ix]*fitted_CDF0 + weights[1,jy,ix]*fitted_CDF1 + \
            weights[2,jy,ix]*fitted_CDF2
        if pflag == True: print ('   fitted_CDF = ',fitted_CDF)
        
        # ---- calculate a statistic for how far off the fitted CDF is from the 
        #      empirical CDF.   Do this only for either the higher quantiles of the
        #      distribution or where the quantile exceeds 1 mm (the i1 index)
        
        Dnstat3 = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        if pflag == True: print ('   np.abs(fitted_CDF - empirical_CDF) = ', \
            np.abs(fitted_CDF - empirical_CDF))
        print ('   Dnstat for 3 gamma mixture = ', Dnstat3)
        nmixture[jy,ix] = 3
    except:
        print ('   three_parameter_fail')
        three_parameter_fail = True
        Dnstat3 = 0.10    
    
    return three_parameter_fail, weights, alpha, beta, Dnstat3, \
        weight_save_3m, alpha_save_3m, beta_save_3m, nmixture
    
# =====================================================================

def decide_which_mixture(jy, ix, weights, alpha, beta, nmixture, pflag, \
    Dnstat1, Dnstat2, Dnstat3, alpha_save_1m, beta_save_1m, \
    weight_save_1m, alpha_save_2m, beta_save_2m, \
    weight_save_2m, alpha_save_3m, beta_save_3m, \
    weight_save_3m):

    """ based on the Dn statistics for 1, 2, 3 Gamma mixtures,
        decide which mixture to use (the one with lowest Dn)
    """

    Dnstats = np.array([Dnstat1, Dnstat2, Dnstat3])
    #if pflag == True: print ('Dnstats 1,2,3,4 = ', Dnstats)
    print ('   Dnstats 1,2,3  = ', Dnstats)
    imin = np.argmin(Dnstats)
    if imin == 0:
        if pflag == True: print ('selected 1-Gamma mixture')
        weights[0,jy,ix] = 1.0
        alpha[0,jy,ix] = alpha_save_1m
        beta[0,jy,ix] = beta_save_1m
        weights[1:,jy,ix] = 0.0
        alpha[1:,jy,ix] = 1.0
        beta[1:,jy,ix] = 1.0
        nmixture[jy,ix] = 1
    elif imin == 1:
        if pflag == True: print ('selected 2-Gamma mixture')
        weights[0:2,jy,ix] = weight_save_2m[:]
        alpha[0:2,jy,ix] = alpha_save_2m[:]
        beta[0:2,jy,ix] = beta_save_2m[:]
        weights[2,jy,ix] = 0.0
        alpha[2,jy,ix] = 1.0
        beta[2,jy,ix] = 1.0
        nmixture[jy,ix] = 2
    elif imin == 2:
        if pflag == True: print ('selected 3-Gamma mixture')
        weights[:,jy,ix] = weight_save_3m[:]
        alpha[:,jy,ix] = alpha_save_3m[:]
        beta[:,jy,ix] = beta_save_3m[:]
        nmixture[jy,ix] = 3
    return weights, alpha, beta, nmixture
    
# =====================================================================
    
# ---- inputs from command line

cmonth = sys.argv[1] # '01', '02' etc.
cend_hour = sys.argv[2] # 06, 12, 18, 00 -- end hour of 6-h period
imonth = int(cmonth) - 1
nstride = 10
cdomain = 'conus'
excessive_threshold = 0.035

# ---- set parameters

pflag = False # for print statements
master_directory = '/Volumes/Backup Plus/mswep/'
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]

# ---- determine the overall number of daily precipitation 
#      samples across all years for this month

if imonth == 1:  # Feb
    nsamps = ndaysomo[imonth]*20
else:
    nsamps = 6*ndaysomo_leap[imonth] + 14*ndaysomo[imonth]
print ('nsamps = ', nsamps)

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time.  All members, dates for this 
#      month have been smushed into one leading index, dimension
#      nsamps, since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution fitting.
   
ktr = 0
for iyear in range(2000,2020):
    if iyear%4 == 0:
        ndays = ndaysomo_leap[imonth]
    else:
        ndays = ndaysomo[imonth]
    cyear = str(iyear)
    infile = master_directory + cyear + cmonth + '_on_ndfd_grid_6hourly.nc'
    print (iyear, infile)
    nc = Dataset(infile)
    yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
    for iday in range(1,ndays+1):
        if iday < 10:
            cday = '0'+str(iday)
        else:
            cday = str(iday)
        iyyyymmddhh = int(str(iyear)+cmonth+cday+cend_hour)
        idx = np.where(yyyymmddhh_end == iyyyymmddhh)[0]
        #print (iyyyymmddhh, idx)
        precip_in = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
        #print ('max min precip_in = ',np.max(precip_in), np.min(precip_in))
        if iyear == 2000 and iday == 1:
            #print ('np.shape(precip_in) = ',np.shape(precip_in))
            #sys.exit()
            nyin, nxin = np.shape(precip_in)
            precip_tseries = np.zeros((nsamps,nyin,nxin), dtype=np.float32)
            lons = nc.variables['lons'][:,:]
            lats = nc.variables['lats'][:,:]
        precip_tseries[ktr,:,:] = precip_in[:,:]
        ktr = ktr+1
    nc.close()

# ---- more initialization of output storage arrays now that 
#      we know the array dimensions

weights = np.zeros((3,nyin, nxin), dtype=np.float)
alpha = np.zeros((3,nyin,nxin), dtype=np.float)
beta = np.zeros((3,nyin,nxin), dtype=np.float)
fzero = np.zeros((nyin, nxin), dtype=np.float)
nmixture = np.zeros((nyin, nxin), dtype=np.int)
Dnstat1a = 0.10*np.ones((nyin, nxin), dtype=np.float)
Dnstat2a = 0.10*np.ones((nyin, nxin), dtype=np.float)
Dnstat3a = 0.10*np.ones((nyin, nxin), dtype=np.float)

weight_save_1m = 1.0
alpha_save_1m = 1.0
beta_save_1m = 1.0
   
# ---- loop over the grid points and estimate the Gamma distributions
#      for each parameter.  First see if a single Gamma distribution
#      is appropriate; if not, try a mixture of two.   If that still
#      doesn't fit well, try a mixture of three.   
        
tktr = 0
weight_save_2m = 0.5*np.ones((2), dtype=np.float)
alpha_save_2m = np.random.uniform(low=0.5, high=1.5, size=2)
beta_save_2m =  np.random.uniform(low=0.5, high=1.5, size=2)
weight_save_3m = 0.33333*np.ones((3), dtype=np.float) # bullshirt values
alpha_save_3m = np.random.uniform(low=0.5, high=1.5, size=3) 
beta_save_3m = np.random.uniform(low=0.5, high=1.5, size=3) 
weight_save_3mop = 0.33333*np.ones((3), dtype=np.float) # bullshirt values
alpha_save_3mop = np.random.uniform(low=0.5, high=1.5, size=3) 
beta_save_3mop = np.random.uniform(low=0.5, high=1.5, size=3) 
Dnstat1 = 0.10
Dnstat2 = 0.10
Dnstat3 = 0.10

now = datetime.now()
begin_time = now.strftime("%H:%M:%S")
for jy in range(0,nyin,nstride):
#for jy in range(60, 61):
        
    for ix in range(0,nxin,nstride):
    #for ix in range (200,251):
        
        mop_up = False
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        tktr = tktr+1  # number of grid points processed
        print ('***** begin time, current time, jy, ix, nyin, nxin, lon, lat = ',\
            begin_time, current_time, jy, ix, nyin, nxin, lons[jy, ix], lats[jy,ix])
            
        # ---- there is a grib round-off error that can give negative
        #      values slightly smaller than teeny precip.  to make sure 
        #      that we don't have either negative values or lots of the 
        #      same tiny values, subtractoff teeny_precip
        
        precip_ens_1d = precip_tseries[:,jy,ix]
        tp = np.min( [np.abs(np.min(precip_ens_1d)), 0.0] )
        teeny_precip = tp*np.ones(nsamps)
        if pflag == True: print ('   tp = ', tp)
        precip_ens_1d = precip_ens_1d - teeny_precip[:]

        fraction_zero, precip_ens_nonzero, nz = \
            fraczero_possamps(nsamps, precip_ens_1d) # return sorted
        if pflag == True: print ('   precip_ens_nonzero[0:-1] = ',\
            precip_ens_nonzero[0:-1])
        fzero[jy,ix] = fraction_zero 
        if nz > 0:
            pmean = np.mean(precip_ens_nonzero)
        else:
            pmean = 0.0
        if pflag == True: print ('   After fraczero_possamps, jy,ix = ', \
            jy,ix, nxin*nyin, current_time)
        if pflag == True: print ('   number of samples with positive precip = ', nz)
        if pflag == True: print ('   precip_ens_nonzero[-10:] = ', \
            precip_ens_nonzero[-20:])
        if nz > 20 and precip_ens_nonzero[-1] > 2.0:
        
            # --- first fit a single Gamma distribution per Wilks; determine 
            #     the quantiles associated with every 1/20th percentile 
            #     (0.05 to 0.95) with a few extras, if there are enough samples.
            #     With smaller number of samples, do the sampling at a fewer
            #     number of quantiles
        
            query_these_indices, empirical_CDF, empirical_precipvals, \
                pmean, alpha_hat, beta_hat, i1, Dnstat1, \
                weight_save_1m, alpha_save_1m, beta_save_1m = \
                one_parameter_gamma(nz, precip_ens_nonzero, pflag, fraction_zero)
            Dnstat1a[jy,ix] = Dnstat1
            
            # ---- decide if Dn excessive.   If so, then try a 2-component 
            #      Gamma mixture model.  excessive_threshold inspired 
            #      by Wilks textbook, Table 5.2, assume for simplicity alpha=1
            #      since there is only a weak dependence on alpha.
            
            excessive_threshold = 0.025
            if pflag == True: print ('excessive_threshold = ', excessive_threshold)  
            if Dnstat1 > excessive_threshold:
                
                two_parameter_fail, weights, alpha, beta, Dnstat2, \
                    weight_save_2m, alpha_save_2m, beta_save_2m, nmixture = \
                    two_parameter_gamma(jy, ix, precip_ens_nonzero, \
                    empirical_precipvals, empirical_CDF, weights, \
                    alpha, beta, i1, pflag, weight_save_2m, alpha_save_2m, \
                    beta_save_2m, nmixture, fraction_zero, nstride, \
                    Dnstat2, excessive_threshold)
                Dnstat2a[jy,ix] = Dnstat2
                if pflag == True:
                    print ('   two_parameter_fail, Dnstat2 = ', two_parameter_fail, Dnstat2)
                    
                if two_parameter_fail == True or Dnstat2 > excessive_threshold:
                
                    # ---- try 3-parameter Gamma
                
                    three_parameter_fail, weights, alpha, beta, Dnstat3,\
                        weight_save_3m, alpha_save_3m, beta_save_3m, nmixture = \
                        three_parameter_gamma(jy, ix, precip_ens_nonzero, \
                        empirical_precipvals, empirical_CDF, weights, \
                        alpha, beta, i1, pflag, weight_save_3m, alpha_save_3m, \
                        beta_save_3m, nmixture, fraction_zero, nstride, mop_up, \
                        Dnstat3, excessive_threshold)
                    Dnstat3a[jy,ix] = Dnstat3
                    
                    if three_parameter_fail == True or Dnstat3 > excessive_threshold: 
                        
                        weights, alpha, beta, nmixture = decide_which_mixture(\
                            jy, ix, weights, alpha, beta, nmixture, pflag, \
                            Dnstat1, Dnstat2, Dnstat3, alpha_save_1m, beta_save_1m, \
                            weight_save_1m, alpha_save_2m, beta_save_2m, \
                            weight_save_2m, alpha_save_3m, beta_save_3m, \
                            weight_save_3m)

            else: # ok to use simple single Gamma
            
                weights[0,jy,ix] = 1.0
                alpha[0,jy,ix] = alpha_hat
                beta[0,jy,ix] = beta_hat
                weights[1:,jy,ix] = 0.0
                alpha[1:,jy,ix] = 1.0
                beta[1:,jy,ix] = 1.0
                nmixture[jy,ix] = 1
            
        else: 
            
            # --- very few positive samples, or light precip; fit an simple maximum
            #      likelihood distribution to single mode Gamma distribution
            
            if pflag == True: print ('   very light precipitation at this grid point. ', pmean)
            if nz < 2:  # not long enough sample to find some nonzero precip amounts
                alpha_hat = 1.0
                beta_hat = 100.0
            else:  # use maximum likelihood method of Wilks text to fit non-mixture distribution.
                pmean = np.mean(precip_ens_nonzero)
                lnxbar = np.log(pmean)
                meanlnxi = np.mean(np.sum(np.log(precip_ens_nonzero)))
                D = lnxbar - meanlnxi
                alpha_hat = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
                beta_hat = pmean / alpha_hat
                if pflag == True: print('   single distribution alpha, beta = ', alpha_hat, beta_hat)
            weights[0,jy,ix] = 1.0
            alpha[0,jy,ix] = alpha_hat
            beta[0,jy,ix] = beta_hat
            weights[1:,jy,ix] = 0.0
            alpha[1:,jy,ix] = 1.0
            beta[1:,jy,ix] = 1.0
            nmixture[jy,ix] = 1
            Dnstat1[jy,ix] = 0.01  # a fictitious value, but precip so low, it doesn't matter
            

        # ----- sometimes, by using the previous grid point's estimate as a prior mixture, it can
        #       go off into non-sensical parameter space.   Check Dn for this, and if so, initiate 
        #       a final 3-parameter mixture *without* using prior grid point's saved values.

        if Dnstat1a[jy,ix] > excessive_threshold and \
        Dnstat2a[jy,ix] > excessive_threshold and \
        Dnstat3a[jy,ix] > excessive_threshold:
        
            mop_up = True
            print ('   mopping up! ')
            three_parameter_fail, weights, alpha, beta, Dnstat3mop,\
                weight_save_3mop, alpha_save_3mop, beta_save_3mop, nmixture = \
                three_parameter_gamma(jy, ix, precip_ens_nonzero, \
                empirical_precipvals, empirical_CDF, weights, \
                alpha, beta, i1, pflag, weight_save_3mop, alpha_save_3mop, \
                beta_save_3mop, nmixture, fraction_zero, nstride, mop_up, \
                Dnstat3, excessive_threshold)

            if Dnstat3mop < Dnstat3a[jy,ix]:
                Dnstat3 = Dnstat3mop
                weight_save_3m = weight_save_3mop
                alpha_save_3m = alpha_save_3mop
                beta_save_3m = beta_save_3mop
                                
            weights, alpha, beta, nmixture = decide_which_mixture(\
                jy, ix, weights, alpha, beta, nmixture, pflag, \
                Dnstat1, Dnstat2, Dnstat3, alpha_save_1m, beta_save_1m, \
                weight_save_1m, alpha_save_2m, beta_save_2m, \
                weight_save_2m, alpha_save_3m, beta_save_3m, \
                weight_save_3m)
            
# --- save to cPickle file

outfile = master_directory + cmonth+'_'+cdomain+'_MSWEP_apcp_gamma_parameters_h' + cend_hour + '.cPick'
#outfile = master_directory + cmonth+'_'+cdomain+'_apcp_gamma_parameters_test_h' + cend_hour + '.cPick'

print ('writing to ', outfile)
ouf = open(outfile, 'wb')
cPickle.dump(weights, ouf)
cPickle.dump(alpha, ouf)
cPickle.dump(beta, ouf)
cPickle.dump(fzero, ouf)
cPickle.dump(Dnstat1a, ouf)
cPickle.dump(Dnstat2a, ouf)
cPickle.dump(Dnstat3a, ouf)
cPickle.dump(nmixture, ouf)
ouf.close()



