"""
CDF_fitting_forecast_precip_v3.py cmonth clead 

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

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

def fraczero_possamps(nsamps, precip_ens):
    """
    from the vector input sample precip_ens, define the fraction of
    samples with zero precipitation.   For the positive samples, add
    a small random number to deal with the fact that the data was 
    discretized to 0.1 mm, so that when later creating CDFs we don't 
    have values with lots of tied amounts.   Sort the nonzero amounts 
    and return.
    """
    number_zeros = 0
    precip_ens_nonzero = np.delete(precip_ens, \
        np.where(precip_ens <= 0.006))  # censor at 0.006 mm
    nz = len(precip_ens_nonzero)
    # data discretized, so add random component of this magnitude
    precip_ens_nonzero = precip_ens_nonzero + \
        np.random.uniform(low=-0.005,high=0.005,size=nz) 
    precip_ens_nonzero = np.sort(precip_ens_nonzero)  
    #print (precip_ens_nonzero[0:10]) 
    ntotal = len(precip_ens)
    nzero = ntotal - len(precip_ens_nonzero)
    fraction_zero = float(nzero) / float(ntotal)
    return fraction_zero, precip_ens_nonzero, nz

# =====================================================================

def one_parameter_gamma(nz, precip_ens_nonzero, pflag, fraction_zero):

    """ 
    Along with some other ancillary stuff, return the fitted Gamma distribution
    alpha and beta values, along with an evaluation (Dnstat) of how 
    closely the fitted CDF matches the empirical CDF 
    
    inputs:  
    
    nz: number of nonzero samples
    precip_ens_nonzero[nz]:   sorted vector of nonzero precip amounts (mm)
    pflag: print flag, true if diagnostic printing desired
    
    """

    # ---- define the indices in the previously sorted precip_ens_nonzero vector
    #      that are at the 0.05, 0.10, .... , 0.95 quantiles of the sorted vector
    
    query_these_indices = [ nz//20, nz//10, (3*nz)//20, nz//5, nz//4, (3*nz)//10, \
        (7*nz)//20, (2*nz)//5, (9*nz)//20, nz//2, (11*nz)//20, (3*nz)//5, (13*nz)//20, \
        (7*nz)//10, (3*nz)//4, (4*nz)//5, (17*nz)//20, (35*nz)//40,(9*nz)//10, \
        (37*nz)//40, (19*nz)//20, (39*nz)//40, (79*nz)//80]
    if pflag == True: print ('   nz, query_these_indices = ', nz, query_these_indices)
    
    # ---- convert the query_these_indices into the cumulative probability
    
    empirical_CDF = np.array(query_these_indices, dtype=np.float) / float(nz) 
    
    if pflag == True: print ('   empirical_CDF = ',empirical_CDF)
    
    # ---- extract the quantiles at these cumulative probabilities
    
    empirical_precipvals = precip_ens_nonzero[query_these_indices] 
    
    if pflag == True: print ('   empirical_precipvals = ', empirical_precipvals)
    #if pflag == True: print ('   precip_ens_nonzero[0:-1:10] = ',precip_ens_nonzero[0:-1:10] )
    
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
    i1 = np.min([17,np.argmin(np.abs(empirical_precipvals-1.0))]) 
    #Dnstat = (1.0 - fraction_zero)*np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
    Dnstat = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
    if pflag == True: print ('   Dnstat for 1 gamma mixture = ', Dnstat)
    weight_save_1m = 1.0
    alpha_save_1m = alpha_hat
    beta_save_1m = beta_hat
    
    
    return query_these_indices, empirical_CDF, empirical_precipvals, \
        pmean, alpha_hat, beta_hat, i1, Dnstat, weight_save_1m, \
        alpha_save_1m, beta_save_1m
        
# =====================================================================

def two_parameter_gamma(jy, ix, precip_ens_nonzero, empirical_precipvals, \
    empirical_CDF, weights, alpha, beta, i1, pflag, weight_save_2m, \
    alpha_save_2m, beta_save_2m, nmixture, fraction_zero):
    
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
    """

    two_parameter_fail = True # in case routine bombs off
    try:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        
        # --- call the R routine that estimates the weights and a mixture 
        #     of two Gamma distributions. Must convert back from R data format..
        
        if pflag == True: print ('calling gammamix_em 2 parameter')
        if pflag == True: print ('precip_ens_nonzero[0], [-1] = ', \
            precip_ens_nonzero[0], precip_ens_nonzero[-1] )
        if pflag == True: print ('weight_save_2m = ', weight_save_2m)
        if pflag == True: print ('alpha_save_2m = ', alpha_save_2m)
        if pflag == True: print ('beta_save_2m = ', beta_save_2m)
        #best_result = gammamix_em(precip_ens_nonzero, \
        #    mix_prop=weight_save_2m, alpha=alpha_save_2m, invbeta=beta_save_2m,
        #    k=2, epsilon=0.001, maxit=50, maxrestarts=20, verb=False)
        best_result = gammamix_em(precip_ens_nonzero, \
            mix_prop=None, alpha=None, invbeta=None,\
            k=2, epsilon=0.002, maxit=100, maxrestarts=20, verb=False)
        
        params = best_result.params
        alpha[0:2,jy,ix] = np.array(params.alpha[:])
        beta[0:2,jy,ix] = np.array(params.invbeta[:])
        weights[0:2,jy,ix] = np.array(params.mix_prop[:])
        weights[2:3,jy,ix] = 0.0
        alpha[2:3,jy,ix] = 1.0
        beta[2:3,jy,ix] = 1.0
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
        
        #Dnstat2 = (1.0-fraction_zeros)*np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        Dnstat2 = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        if pflag == True: print ('   np.abs(fitted_CDF - empirical_CDF) = ', \
            np.abs(fitted_CDF - empirical_CDF))
        if pflag == True: print ('   Dnstat for 2 gamma mixture = ', Dnstat2)
    except:
        print ('two_parameter_fail')
        two_parameter_fail = True
        Dnstat2 = 0.10

    return two_parameter_fail, weights, alpha, beta, Dnstat2, weight_save_2m, \
        alpha_save_2m, beta_save_2m, nmixture

    
# =====================================================================

def three_parameter_gamma(jy, ix, precip_ens_nonzero, empirical_precipvals, \
    empirical_CDF, weights, alpha, beta, i1, pflag, weight_save_3m, alpha_save_3m, \
    beta_save_3m, nmixture, fraction_zero):    
    
    """
    
    Call an R gamma mixture routine to estimate the weights and parameters of
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
    
    """

    three_parameter_fail = False # in case routine bombs off
    try:    
    
        if pflag == True: print ('   Need to try 3 gamma mixture.')
        
        # --- call the R routine that estimates the weights and a mixture 
        #     of 3 Gamma distributions. Must convert back from R data format.
        
        if pflag == True: print ('calling gammamix_em 3 parameter')
        if pflag == True: print ('precip_ens_nonzero[0], [-1] = ', \
            precip_ens_nonzero[0], precip_ens_nonzero[-1] )
        if pflag == True: print ('weight_save_3m = ', weight_save_3m)
        if pflag == True: print ('alpha_save_3m = ', alpha_save_3m)
        if pflag == True: print ('beta_save_3m = ', beta_save_3m)
        
        #best_result = gammamix_em(precip_ens_nonzero, \
        #    mix_prop=weight_save_3m, alpha=alpha_save_3m, invbeta=beta_save_3m,
        #    k=3, epsilon=0.001, maxit=50, maxrestarts=20, verb=False)
            
        best_result = gammamix_em(precip_ens_nonzero, \
            mix_prop=None, alpha=None, invbeta=None,\
            k=3, epsilon=0.002, maxit=100, maxrestarts=20, verb=False)
         
        params = best_result.params
        alpha[0:3,jy,ix] = np.array(params.alpha[:])
        beta[0:3,jy,ix] = np.array(params.invbeta[:])
        weights[0:3,jy,ix] = np.array(params.mix_prop[:])
        alpha[3,jy,ix] = 1.0
        beta[3,jy,ix] = 1.0
        weights[3,jy,ix] = 0.0
                
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
        
        #Dnstat3 = (1.0-fraction_zero)*np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        Dnstat3 = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        if pflag == True: print ('   np.abs(fitted_CDF - empirical_CDF) = ', \
            np.abs(fitted_CDF - empirical_CDF))
        if pflag == True: print ('   Dnstat for 3 gamma mixture = ', Dnstat3)
        nmixture[jy,ix] = 3
    except:
        print ('three_parameter_fail')
        three_parameter_fail = True
        Dnstat3 = 0.10    
    
    return three_parameter_fail, weights, alpha, beta, Dnstat3, \
        weight_save_3m, alpha_save_3m, beta_save_3m, nmixture
        
        
# =====================================================================

def four_parameter_gamma(jy, ix, precip_ens_nonzero, empirical_precipvals, \
    empirical_CDF, weights, alpha, beta, i1, pflag, weight_save_4m, alpha_save_4m, \
    beta_save_4m, nmixture, fraction_zero):    
    
    """
    
    Call an R gamma mixture routine to estimate the weights and parameters of
    a mixture of four Gamma distributions.  Evaluate how well it fits the empirical
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
    weight_save_4m: first guess for weights from previous call
    alpha_save_4m: first guess for alpha
    beta_save_4m: first guess for beta
    nmixture: indicates how many parameters are used in the mix
    fraction_zero: fraction of points with zero precip
    
    """

    four_parameter_fail = False # in case routine bombs off
    try:    
    
        if pflag == True: print ('   Need to try 4 gamma mixture.')
        
        # --- call the R routine that estimates the weights and a mixture 
        #     of 3 Gamma distributions. Must convert back from R data format.
        
        if pflag == True: print ('calling gammamix_em 4 parameter')
        if pflag == True: print ('precip_ens_nonzero[0], [-1] = ', \
            precip_ens_nonzero[0], precip_ens_nonzero[-1] )
        if pflag == True: print ('weight_save = ', weight_save_4m)
        if pflag == True: print ('alpha_save = ', alpha_save_4m)
        if pflag == True: print ('beta_save = ', beta_save_4m)
        
        #best_result = gammamix_em(precip_ens_nonzero, \
        #    mix_prop=weight_save_4m, alpha=alpha_save_4m, invbeta=beta_save_4m,
        #    k=4, epsilon=0.001, maxit=50, maxrestarts=20, verb=True)
        #best_result = gammamix_em(precip_ens_nonzero, \
        #    mix_prop=None, alpha=None, invbeta=None,
        #    k=4, epsilon=0.001, maxit=50, maxrestarts=20, verb=False)
        best_result = gammamix_em(precip_ens_nonzero, \
            mix_prop=None, alpha=None, invbeta=None,\
            k=4, epsilon=0.001, maxit=100, maxrestarts=20, verb=False)
        
        weights[:,jy,ix] = 0.0
        alpha[:,jy,ix] = 0.0
        beta[:,jy,ix] = 0.0    
        
        params = best_result.params
        alpha[:,jy,ix] = np.array(params.alpha[:])
        beta[:,jy,ix] = np.array(params.invbeta[:])
        weights[:,jy,ix] = np.array(params.mix_prop[:])
        
        weight_save_4m = np.array(params.mix_prop[:])
        alpha_save_4m = np.array(params.alpha[:])
        beta_save_4m = np.array(params.invbeta[:])
        

        if pflag == True: print ('   weights = ', weights[:,jy,ix])
        if pflag == True: print ('   alpha = ', alpha[:,jy,ix])
        if pflag == True: print ('   beta = ', beta[:,jy,ix])
        
        #  ---- estimate a CDF from a weighted mixture of the 3 distributions.  
        
        y0 = empirical_precipvals / beta[0,jy,ix]
        y1 = empirical_precipvals / beta[1,jy,ix]
        y2 = empirical_precipvals / beta[2,jy,ix]
        y3 = empirical_precipvals / beta[3,jy,ix]
        fitted_CDF0 = stats.gamma.cdf(y0, alpha[0,jy,ix])
        fitted_CDF1 = stats.gamma.cdf(y1, alpha[1,jy,ix])
        fitted_CDF2 = stats.gamma.cdf(y2, alpha[2,jy,ix])
        fitted_CDF3 = stats.gamma.cdf(y2, alpha[3,jy,ix])
        fitted_CDF = weights[0,jy,ix]*fitted_CDF0 + weights[1,jy,ix]*fitted_CDF1 + \
            weights[2,jy,ix]*fitted_CDF2 + weights[3,jy,ix]*fitted_CDF3
        if pflag == True: print ('   fitted_CDF = ',fitted_CDF)
        
        # ---- calculate a statistic for how far off the fitted CDF is from the 
        #      empirical CDF.   Do this only for either the higher quantiles of the
        #      distribution or where the quantile exceeds 1 mm (the i1 index)
        
        #Dnstat4 = (1.0-fraction_zero)*np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        Dnstat4 = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        if pflag == True: print ('   np.abs(fitted_CDF - empirical_CDF) = ', \
            np.abs(fitted_CDF - empirical_CDF))
        if pflag == True: print ('   Dnstat for 4 gamma mixture = ', Dnstat4)
        nmixture[jy,ix] = 4
    except:
        print ('four_parameter_fail')
        four_parameter_fail = True
        Dnstat4 = 0.10    
    
    return four_parameter_fail, weights, alpha, beta, Dnstat4, \
        weight_save_4m, alpha_save_4m, beta_save_4m, nmixture    
        
# =====================================================================        
        
def decide_which_mixture(jy,ix, weights, alpha, beta, pflag, \
    Dnstat1, Dnstat2, Dnstat3, Dnstat4, alpha_save_1m, beta_save_1m,\
    weight_save_1m, alpha_save_2m, beta_save_2m,\
    weight_save_2m, alpha_save_3m, beta_save_3m,\
    weight_save_3m, alpha_save_4m, beta_save_4m,\
    weight_save_4m, nmixture):
    
    """ based on the Dn statistics for 1, 2, 3, 4 Gamma mixtures, 
        decide which mixture to use (the one with lowest Dn)
    """
    
    Dnstats = np.array([Dnstat1, Dnstat2, Dnstat3, Dnstat4])
    #if pflag == True: print ('Dnstats 1,2,3,4 = ', Dnstats)
    print ('Dnstats 1,2,3,4 = ', Dnstats)
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
        weights[0:2,jy,ix] = weight_save_2m
        alpha[0:2,jy,ix] = alpha_save_2m
        beta[0:2,jy,ix] = beta_save_2m
        weights[2:,jy,ix] = 0.0
        alpha[2:,jy,ix] = 1.0
        beta[2:,jy,ix] = 1.0 
        nmixture[jy,ix] = 2
    elif imin == 2:
        if pflag == True: print ('selected 3-Gamma mixture')
        weights[0:3,jy,ix] = weight_save_3m
        alpha[0:3,jy,ix] = alpha_save_3m
        beta[0:3,jy,ix] = beta_save_3m
        weights[3,jy,ix] = 0.0
        alpha[3,jy,ix] = 1.0
        beta[3,jy,ix] = 1.0 
        nmixture[jy,ix] = 3
    elif imin == 3:
        if pflag == True: print ('selected 4-Gamma mixture')
        weights[:,jy,ix] = weight_save_4m
        alpha[:,jy,ix] = alpha_save_4m
        beta[:,jy,ix] = beta_save_4m
        nmixture[jy,ix] = 1
                                         
    return weights, alpha, beta, nmixture  
    
# =====================================================================

def relax_toward_equal_weight(weight_save_2m, \
    weight_save_3m, weight_save_4m):
    
    w2 = np.array([0.5, 0.5])
    w3 = np.array([0.33333, 0.33333, 0.33333])
    w4 = np.array([0.25, 0.25, 0.25, 0.25])
    weight_save_2m = (weight_save_2m + w2) / 2.
    weight_save_3m = (weight_save_3m + w3) / 2.
    weight_save_4m = (weight_save_4m + w4) / 2.
    
    return weight_save_2m, weight_save_3m, weight_save_4m

# =====================================================================
# =====================================================================

    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clead = sys.argv[2] # 03, 06, 12, etc.

# ---- set parameters

pflag = False # for print statements
master_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'
nmembers = 5

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time.  All members, dates for this 
#      month have been smushed into one leading index, dimension
#      nsamps, since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution fitting.

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
if pflag == True: print ('before read ', current_time)    
ncfile = master_directory + cmonth + '_apcp_h' + clead + '.nc'
if pflag == True: print (ncfile)
nc = Dataset(ncfile)
precip_ens = nc.variables['apcp_fcst'][:,:,:]
nsamps, nyin, nxin = np.shape(precip_ens)
lons_1d = nc.variables['lons_fcst'][:]
lats_1d = nc.variables['lats_fcst'][:]
nc.close()
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
if pflag == True: print ('after read ', current_time) 

# ---- more initialization of output storage arrays now that 
#      we know the array dimensions

weights = np.zeros((4,nyin, nxin), dtype=np.float)
alpha = np.zeros((4,nyin,nxin), dtype=np.float)
beta = np.zeros((4,nyin,nxin), dtype=np.float)
fzero = np.zeros((nyin, nxin), dtype=np.float)
nmixture = np.zeros((nyin, nxin), dtype=np.int)
Dnstat1a = -99.99*np.ones((nyin, nxin), dtype=np.float)
Dnstat2a = -99.99*np.ones((nyin, nxin), dtype=np.float)
Dnstat3a = -99.99*np.ones((nyin, nxin), dtype=np.float)
Dnstat4a = -99.99*np.ones((nyin, nxin), dtype=np.float)
teeny_precip = 0.06*np.ones(nsamps) 

weight_save_1m = 1.0
alpha_save_1m = 1.0
beta_save_1m = 1.0
   
# ---- loop over the grid points and estimate the Gamma distributions
#      for each parameter.  First see if a single Gamma distribution
#      is appropriate; if not, try a mixture of two.   If that still
#      doesn't fit well, try a mixture of three.   
        
tktr = 0
#for jy in range(0,nyin,5):
for jy in range(290,291,5):
    
    
    
    
    
    
    if jy == 290:
        
        
        
        
        
        weight_save_4m = 0.25*np.ones((4), dtype=np.float) # bullshirt values
        alpha_save_4m = np.random.uniform(low=0.5, high=1.5, size=4) 
        beta_save_4m = np.random.uniform(low=0.5, high=1.5, size=4) 
        weight_save_3m = 0.33333*np.ones((3), dtype=np.float) # bullshirt values
        alpha_save_3m = np.random.uniform(low=0.5, high=1.5, size=3) 
        beta_save_3m = np.random.uniform(low=0.5, high=1.5, size=3) 
        weight_save_2m = 0.5*np.ones((2), dtype=np.float)
        alpha_save_2m = np.random.uniform(low=0.5, high=1.5, size=2)
        beta_save_2m =  np.random.uniform(low=0.5, high=1.5, size=2)
    
    #for ix in range(0,nxin,5):
    for ix in range(170,231,5):
        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        tktr = tktr+1  # number of grid points processed
        print ('*** time, jy, ix, lon, lat = ',\
            current_time, jy, ix, tktr, nyin*nxin, lons_1d[ix], lats_1d[jy])
        if pflag == True: print ('   Before fraczero_possamps, jy,ix, \
            nsamps, time = ', jy,ix, nxin*nyin, current_time)
            
        # ---- there is a grib compaction error that can give negative
        #      values slightly smaller than teeny precip.  to make sure 
        #      that we don't have either negative values or lots of the 
        #      same tiny values, subtractoff teeny_precip
        
        precip_ens_1d = precip_ens[:,jy,ix] - teeny_precip[:]
        fraction_zero, precip_ens_nonzero, nz = \
            fraczero_possamps(nsamps, precip_ens_1d) # return sorted
        fzero[jy,ix] = fraction_zero 
        if pflag == True: print ('   After fraczero_possamps, jy,ix = ', \
            jy,ix, nxin*nyin, current_time)
        if pflag == True: print ('   number of samples with positive precip = ', nz)
        #print ('   precip_ens_nonzero[-10:] = ', precip_ens_nonzero[-20:])
        if nz > 30 and precip_ens_nonzero[-1] > 4.0:
        
            # --- first fit a single Gamma distribution per Wilks; determine 
            #     the quantiles associated with every 1/20th percentile 
            #     (0.05 to 0.95)
        
            query_these_indices, empirical_CDF, empirical_precipvals, \
                pmean, alpha_hat, beta_hat, i1, Dnstat1, \
                weight_save_1m, alpha_save_1m, beta_save_1m = \
                one_parameter_gamma(nz, precip_ens_nonzero, pflag, fraction_zero)
            Dnstat1a[jy,ix] = Dnstat1
            
            # ---- decide if Dn excessive.   If so, then try a 2-component 
            #      Gamma mixture model.  excessive_threshold inspired 
            #      by Wilks textbook, Table 5.2, assume for simplicity alpha=1
            #      since there is only a weak dependence on alpha.
            
            excessive_threshold = 0.02
            if pflag == True: print ('excessive_threshold = ', excessive_threshold)  
            if Dnstat1 > excessive_threshold:
                
                two_parameter_fail, weights, alpha, beta, Dnstat2, \
                    weight_save_2m, alpha_save_2m, beta_save_2m, nmixture = \
                    two_parameter_gamma(jy, ix, precip_ens_nonzero, \
                    empirical_precipvals, empirical_CDF, weights, \
                    alpha, beta, i1, pflag, weight_save_2m, alpha_save_2m, \
                    beta_save_2m, nmixture, fraction_zero)
                    
                Dnstat2a[jy,ix] = Dnstat2
                    
                if two_parameter_fail == True or Dnstat2 > excessive_threshold:
                
                    # ---- try 3-parameter Gamma
                
                    three_parameter_fail, weights, alpha, beta, Dnstat3,\
                        weight_save_3m, alpha_save_3m, beta_save_3m, nmixture = \
                        three_parameter_gamma(jy, ix, precip_ens_nonzero, \
                        empirical_precipvals, empirical_CDF, weights, \
                        alpha, beta, i1, pflag, weight_save_3m, alpha_save_3m, \
                        beta_save_3m, nmixture, fraction_zero)
                    Dnstat3a[jy,ix] = Dnstat3
                    
                    if three_parameter_fail == True or Dnstat3 > excessive_threshold: 
                        
                        four_parameter_fail, weights, alpha, beta, Dnstat4,\
                            weight_save_4m, alpha_save_4m, beta_save_4m, nmixture = \
                            four_parameter_gamma(jy, ix, precip_ens_nonzero, \
                            empirical_precipvals, empirical_CDF, weights, \
                            alpha, beta, i1, pflag, weight_save_4m, alpha_save_4m, \
                            beta_save_4m, nmixture, fraction_zero)
                        Dnstat4a[jy,ix] = Dnstat4
                            
                        if four_parameter_fail == True or Dnstat4 > excessive_threshold: 
                            weights, alpha, beta, nmixture = decide_which_mixture(jy,ix,
                                weights, alpha, beta, pflag, Dnstat1, \
                                Dnstat2, Dnstat3, Dnstat4, alpha_save_1m, beta_save_1m,\
                                weight_save_1m, alpha_save_2m, beta_save_2m,\
                                weight_save_2m, alpha_save_3m, beta_save_3m,\
                                weight_save_3m, alpha_save_4m, beta_save_4m,\
                                weight_save_4m, nmixture)
                
                # ---- these weights are saved for the next grid point, but sometimes
                #      the weights are far from equal and not appropriate if the 
                #      next grid point has a much smaller sample.   So nudge toward
                #      equal weighting.
                
                weight_save_2m, weight_save_3m, weight_save_4m = \
                    relax_toward_equal_weight(weight_save_2m, \
                    weight_save_3m, weight_save_4m)
                    
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
                #if pflag == True: print('   lnxbar, meanlnxi = ',lnxbar, meanlnxi)
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

# --- save to cPickle file



#sys.exit()

#outfile = master_directory + cmonth+ '_apcp_gamma_parameters_h' + clead + '.cPick'
#outfile = master_directory + cmonth+ '_apcp_gamma_parameters_v2_h' + clead + '.cPick'
outfile = master_directory + cmonth+ '_apcp_gamma_parameters_v3_h' + clead + '.cPick'
print ('writing to ', outfile)
ouf = open(outfile, 'wb')
cPickle.dump(weights, ouf)
cPickle.dump(alpha, ouf)
cPickle.dump(beta, ouf)
cPickle.dump(fzero, ouf)
cPickle.dump(Dnstat1a, ouf)
cPickle.dump(Dnstat2a, ouf)
cPickle.dump(Dnstat3a, ouf)
cPickle.dump(Dnstat4a, ouf)
cPickle.dump(nmixture, ouf)
ouf.close()



