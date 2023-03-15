"""
CDF_fitting_forecast_precip_v2.py cmonth clead 

"""

import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
from netCDF4 import Dataset
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import _pickle as cPickle
import scipy.stats as stats

base = importr('base')
mixtools = importr('mixtools')

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
        np.where(precip_ens <= 0.0))  # censor at 0.1 mm
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

def one_parameter_gamma(nz, precip_ens_nonzero, pflag):

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
    if pflag == True: print ('   precip_ens_nonzero[0:-1:10] = ',precip_ens_nonzero[0:-1:10] )
    
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
    Dnstat = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
    if pflag == True: print ('   Dnstat for 1 gamma mixture = ', Dnstat)
    
    return query_these_indices, empirical_CDF, empirical_precipvals, \
        pmean, alpha_hat, beta_hat, i1, Dnstat
        
# =====================================================================

def two_parameter_gamma(jy, ix, precip_ens_nonzero, empirical_precipvals, \
    empirical_CDF, weights, alpha, beta, i1, pflag):
    
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
    
    """

    two_parameter_fail = False # in case routine bombs off
    try:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        
        # --- call the R routine that estimates the weights and a mixture 
        #     of two Gamma distributions. Must convert back from R data format..
        
        precip_nonzero_R = robjects.FloatVector(precip_ens_nonzero)
        result_R = mixtools.gammamixEM(precip_nonzero_R, k=2, maxit=1000 ) #, \
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        if pflag == True: print ('   2 Gammas. After mixtools.gammamixEM, jy,ix = ', \
            jy,ix, nxin*nyin, current_time)
        result_np = np.asarray(result_R, dtype=object)
        result_weights_np = np.asarray(result_R[1])
        result_alpha_beta_np = np.asarray(result_R[2])
        weights[0:2,jy,ix] = result_weights_np[:]
        alpha[0:2,jy,ix] = result_alpha_beta_np[0,:]
        beta[0:2,jy,ix] = result_alpha_beta_np[1,:]
        weights[2,jy,ix] = 0.0
        alpha[2,jy,ix] = 1.0
        beta[2,jy,ix] = 1.0
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
        
        Dnstat = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        if pflag == True: print ('   np.abs(fitted_CDF - empirical_CDF) = ', \
            np.abs(fitted_CDF - empirical_CDF))
        if pflag == True: print ('   Dnstat for 2 gamma mixture = ', Dnstat)
        if pflag == True: print ('   weights = ', result_weights_np[:])
        if pflag == True: print ('   result_alpha_beta_np[0,:] = ', \
            result_alpha_beta_np[0,:])
        if pflag == True: print ('   result_alpha_beta_np[1,:] = ', \
            result_alpha_beta_np[1,:])
    except:
        two_parameter_fail = True
        Dnstat = 0.10

    return two_parameter_fail, weights, alpha, beta, Dnstat

    
# =====================================================================

def three_parameter_gamma(jy, ix, precip_ens_nonzero, empirical_precipvals, \
    empirical_CDF, weights, alpha, beta, i1, pflag, weight_save, alpha_save, \
    beta_save):    
    
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
    weight_save: first guess for weights
    alpha_save: first guess for alpha
    beta_save: first guess for beta
    
    """

    three_parameter_fail = False # in case routine bombs off
    try:    
    
        if pflag == True: print ('   Need to try 3 gamma mixture.')
        
        # --- call the R routine that estimates the weights and a mixture 
        #     of 3 Gamma distributions. Must convert back from R data format.
        
        precip_nonzeroR = robjects.FloatVector(precip_ens_nonzero)
        weightsR = robjects.FloatVector(weight_save)
        alphaR = robjects.FloatVector(alpha_save)
        betaR = robjects.FloatVector(beta_save)
        #result_R = mixtools.gammamixEM(precip_nonzero_R, k=3, maxit=1000) #, 
        result_R = mixtools.gammamixEM(precip_nonzeroR,lambda==weightsR,alpha=alphaR,beta=betaR,k=3, maxit=1000) 
        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        result_np = np.asarray(result_R, dtype=object)
        result_weights_np = np.asarray(result_R[1])
        result_alpha_beta_np = np.asarray(result_R[2])
        weights[:,jy,ix] = result_weights_np[:]
        alpha[:,jy,ix] = result_alpha_beta_np[0,:]
        beta[:,jy,ix] = result_alpha_beta_np[1,:]
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
        
        Dnstat = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
        if pflag == True: print ('   np.abs(fitted_CDF - empirical_CDF) = ', \
            np.abs(fitted_CDF - empirical_CDF))
        if pflag == True: print ('   Dnstat for 3 gamma mixture = ', Dnstat)
    except:
        three_parameter_fail = True
        Dnstat = 0.10    
    
    return three_parameter_fail, weights, alpha, beta, Dnstat
    
# =====================================================================
# =====================================================================
# =====================================================================

    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clead = sys.argv[2] # 03, 06, 12, etc.

# ---- set parameters

pflag = True # for print statements
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
ncfile = master_directory + cmonth + '_apcp' '_h' + clead + '.nc'
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

weights = np.zeros((3,nyin, nxin), dtype=np.float)
alpha = np.zeros((3,nyin,nxin), dtype=np.float)
beta = np.zeros((3,nyin,nxin), dtype=np.float)
fzero = np.zeros((nyin, nxin), dtype=np.float)
teeny_precip = 0.06*np.ones(nsamps) 
weight_save = np.zeros((3), dtype=np.float)
alpha_save = np.zeros((3), dtype=np.float)
beta_save = np.zeros((3), dtype=np.float)
   
   
# ---- loop over the grid points and estimate the Gamma distributions
#      for each parameter.  First see if a single Gamma distribution
#      is appropriate; if not, try a mixture of two.   If that still
#      doesn't fit well, try a mixture of three.   
        
tktr = 0
#for jy in range(0,nyin,5):
#    for ix in range(0,nxin,5):
        
for jy in range(200,201):
    for ix in range(600,601):

        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('jy, ix = ', jy, ix, nyin, nxin, current_time)
        tktr = tktr+1  # number of grid points processed
        if pflag == True: print ('*** time, jy, ix, lon, lat = ',\
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
        if pflag == True: print ('   precip_ens_nonzero[-20:] = ', \
            precip_ens_nonzero[-20:])
        if nz > 21 and precip_ens_nonzero[-1] > 3.0:
        
            # --- first fit a single Gamma distribution per Wilks; determine 
            #     the quantiles associated with every 1/20th percentile 
            #     (0.05 to 0.95)
        
            query_these_indices, empirical_CDF, empirical_precipvals, \
                pmean, alpha_hat, beta_hat, i1, Dnstat = \
                one_parameter_gamma(nz, precip_ens_nonzero, pflag)
            
            # ---- decide if Dn excessive.   If so, then try a 2-component 
            #      Gamma mixture model.  excessive_threshold inspired 
            #      by Wilks textbook, Table 5.2, assume for simplicity alpha=1
            #      since there is only a weak dependence on alpha.
            
            excessive_threshold = 0.3 / np.sqrt(np.float(nz))   
            if Dnstat > excessive_threshold:
                
                two_parameter_fail, weights, alpha, beta, Dnstat = \
                    two_parameter_gamma(jy, ix, precip_ens_nonzero, \
                    empirical_precipvals, empirical_CDF, weights, \
                    alpha, beta, i1, pflag)
                    
                if two_parameter_fail == True or Dnstat > excessive_threshold:
                
                    # ---- try 3-parameter Gamma
                
                    three_parameter_fail, weights, alpha, beta, Dnstat = \
                        three_parameter_gamma(jy, ix, precip_ens_nonzero, \
                        empirical_precipvals, empirical_CDF, weights, \
                        alpha, beta, i1, pflag, weight_save, alpha_save, \
                        beta_save)
                    weight_save[:] = weights[:,jy,ix]
                    alpha_save[:] = alpha[:,jy,ix]
                    beta_save[:] = beta[:,jy,ix]
                    if three_parameter_fail == True: # bombed off w. 3-param.  
                        weights[0,jy,ix] = 1.0       # Revert to one parameter
                        alpha[0,jy,ix] = alpha_hat
                        beta[0,jy,ix] = beta_hat
                        weights[1:,jy,ix] = 0.0
                        alpha[1:,jy,ix] = 1.0
                        beta[1:,jy,ix] = 1.
                    
            else: # ok to use simple single Gamma
            
                weights[0,jy,ix] = 1.0
                alpha[0,jy,ix] = alpha_hat
                beta[0,jy,ix] = beta_hat
                weights[1:,jy,ix] = 0.0
                alpha[1:,jy,ix] = 1.0
                beta[1:,jy,ix] = 1.0
            
        else: 
            
            # --- very few positive samples; fit an simple maximum likelihood distribution
            #     to single mode Gamma distribution
            
            if pflag == True: print ('   very light precipitation at this grid point. ', pmean)
            if nz < 2:  # not long enough sample to find some nonzero precip amounts
                alpha_hat = 1.0
                beta_hat = 100.0
            else:  # use maximum likelihood method of Wilks text to fit non-mixture distribution.
                pmean = np.mean(precip_ens_nonzero)
                lnxbar = np.log(pmean)
                meanlnxi = np.mean(np.sum(np.log(precip_ens_nonzero)))
                if pflag == True: print('   lnxbar, meanlnxi = ',lnxbar, meanlnxi)
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

# --- save to cPickle file

#outfile = master_directory + cmonth+ '_apcp_gamma_parameters_h' + clead + '.cPick'
outfile = master_directory + cmonth+ '_apcp_gamma_parameters_v2_h' + clead + '.cPick'
print ('writing to ', outfile)
ouf = open(outfile, 'wb')
cPickle.dump(weights, ouf)
cPickle.dump(alpha, ouf)
cPickle.dump(beta, ouf)
cPickle.dump(fzero, ouf)
ouf.close()



