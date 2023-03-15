"""
CDF_fitting_forecast_precip.py cmonth clead 

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import _pickle as cPickle
import scipy.stats as stats
#from numba import jit
base = importr('base')
mixtools = importr('mixtools')


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

def fraczero_possamps(nsamps, nxin, nyin, precip_ens):
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
    
    along with some other ancillary stuff, return the fitted Gamma distribution
    alpha and beta values, along with an evaluation (Dnstat) of how close the fitted
    CDF matches the empirical CDF 
    
    inputs:  
    
    nz: number of nonzero samples
    precip_ens_nonzero[nz]:   sorted vector of nonzero precip amounts (mm)
    
    """

    # ---- define the indices in the previously sorted precip_ens_nonzero vector
    #      that are at the 0.05, 0.10, .... , 0.95 quantiles of the sorted vector
    
    query_these_indices = [ nz//20, nz//10, (3*nz)//20, nz//5, nz//4, (3*nz)//10, \
        (7*nz)//20, (2*nz)//5, (9*nz)//20, nz//2, (11*nz)//20, (3*nz)//5, (13*nz)//20, \
        (7*nz)//10, (3*nz)//4, (4*nz)//5, (17*nz)//20, (9*nz)//10, (19*nz)//20]
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
    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', etdc
clead = sys.argv[2] # 06, etc.

pflag = True # for print statements
master_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'
nmembers = 5

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('before read ', current_time)    
ncfile = master_directory + cmonth + '_apcp' '_h' + clead + '.nc'
print (ncfile)
nc = Dataset(ncfile)
precip_ens = nc.variables['apcp_fcst'][:,:,:]
nsamps, nyin, nxin = np.shape(precip_ens)
lons_1d = nc.variables['lons_fcst'][:]
lats_1d = nc.variables['lats_fcst'][:]
nc.close()
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('after read ', current_time) 
weights = np.zeros((3,nyin, nxin), dtype=np.float)
alpha = np.zeros((3,nyin,nxin), dtype=np.float)
beta = np.zeros((3,nyin,nxin), dtype=np.float)
fzero = np.zeros((nyin, nxin), dtype=np.float)
teeny_precip = 0.06*np.ones(nsamps)
        
tktr = 0
for jy in range(0,nyin,5):
    for ix in range(0,nxin,5):
        
#for jy in range(5,6):
#    for ix in range(844,845):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        tktr = tktr+1
        print ('*** time, jy, ix, lon, lat = ',current_time, jy, ix, tktr, nyin*nxin, lons_1d[ix], lats_1d[jy])
        print ('   Before fraczero_possamps, jy,ix, nsamps, time = ', jy,ix, nxin*nyin, current_time)
        # there is a grib compaction error that can give negative values slightly smaller than teeny precip.  
        # to make sure that we don't have either negative values or lots of the same tiny values, subtract
        # off teeny_precip
        precip_ens_1d = precip_ens[:,jy,ix] - teeny_precip[:]
        fraction_zero, precip_ens_nonzero, nz = \
            fraczero_possamps(nsamps, nxin, nyin, precip_ens_1d) # return sorted
        fzero[jy,ix] = fraction_zero 
        print ('   After fraczero_possamps, jy,ix = ', jy,ix, nxin*nyin, current_time)
        print ('   number of samples with positive precip = ', nz)
        print ('   precip_ens_nonzero[-20:] = ', precip_ens_nonzero[-20:])
        if nz > 21 and precip_ens_nonzero[-1] > 3.0:
            

        
            # --- first fit a Gamma distribution per Wilks; determine the quantiles every 1/20th
        

            one_parameter_gamma(nz, precip_ens_nonzero, pflag):
            
            # ---- decide if Dn excessive.   If so, then try a Gamma mixture model.
            
            excessive_threshold = 1.0 / np.sqrt(np.float(nz))   # inspired by Wilks table 5.2, assume alpha=1
            if Dnstat > excessive_threshold:
                
                try:
                
                    # ---- try 2-parameter Gamma
                
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    print ('   Dnstat excessive!  try Gamma mixture. ', Dnstat, current_time)
                    precip_nonzero_R = robjects.FloatVector(precip_ens_nonzero)
                    result_R = mixtools.gammamixEM(precip_nonzero_R, k=2, maxit=1000 ) #, \
                        #alpha <- c(alpha_hat, alpha_hat), beta <- c(0.7*beta_hat, 1.3*beta_hat) )
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    #print ('   After mixtools.gammamixEM, jy,ix = ', jy,ix, nxin*nyin, current_time)
                    result_np = np.asarray(result_R, dtype=object)
                    result_weights_np = np.asarray(result_R[1])
                    result_alpha_beta_np = np.asarray(result_R[2])
                    weights[0:2,jy,ix] = result_weights_np[:]
                    alpha[0:2,jy,ix] = result_alpha_beta_np[0,:]
                    beta[0:2,jy,ix] = result_alpha_beta_np[1,:]
                    weights[2,jy,ix] = 0.0
                    alpha[2,jy,ix] = 1.0
                    beta[2,jy,ix] = 1.0
                    print ('   weights = ', weights[:,jy,ix])
                    print ('   alpha = ', alpha[:,jy,ix])
                    print ('   beta = ', beta[:,jy,ix])
            
                    y0 = empirical_precipvals / beta[0,jy,ix]
                    y1 = empirical_precipvals / beta[1,jy,ix]
                    fitted_CDF0 = stats.gamma.cdf(y0, alpha[0,jy,ix])
                    fitted_CDF1 = stats.gamma.cdf(y1, alpha[1,jy,ix])
                    fitted_CDF = weights[0,jy,ix]*fitted_CDF0 + weights[1,jy,ix]*fitted_CDF1
                    print ('   fitted_CDF = ',fitted_CDF)
                    Dnstat = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
                    print ('   np.abs(fitted_CDF - empirical_CDF) = ', np.abs(fitted_CDF - empirical_CDF))
                    print ('   Dnstat for 2 gamma mixture = ', Dnstat)
                    print ('   weights = ', result_weights_np[:])
                    print ('   result_alpha_beta_np[0,:] = ', result_alpha_beta_np[0,:])
                    print ('   result_alpha_beta_np[1,:] = ', result_alpha_beta_np[1,:])
                    if Dnstat > excessive_threshold:
                        print ('   need to try 3 gamma mixture')
                        result_R = mixtools.gammamixEM(precip_nonzero_R, k=3, maxit=1000) #, 
                            #alpha<-c(result_alpha_beta_np[0,0],result_alpha_beta_np[0,1], 1.0), 
                            #beta<-c(result_alpha_beta_np[1,0],result_alpha_beta_np[1,1], 1.0))
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        #print ('   After mixtools.gammamixEM, jy,ix = ', jy,ix, nxin*nyin, current_time)
                        result_np = np.asarray(result_R, dtype=object)
                        result_weights_np = np.asarray(result_R[1])
                        result_alpha_beta_np = np.asarray(result_R[2])
                        weights[:,jy,ix] = result_weights_np[:]
                        alpha[:,jy,ix] = result_alpha_beta_np[0,:]
                        beta[:,jy,ix] = result_alpha_beta_np[1,:]
                        print ('   weights = ', weights[:,jy,ix])
                        print ('   alpha = ', alpha[:,jy,ix])
                        print ('   beta = ', beta[:,jy,ix])
            
                        y0 = empirical_precipvals / beta[0,jy,ix]
                        y1 = empirical_precipvals / beta[1,jy,ix]
                        y2 = empirical_precipvals / beta[2,jy,ix]
                        fitted_CDF0 = stats.gamma.cdf(y0, alpha[0,jy,ix])
                        fitted_CDF1 = stats.gamma.cdf(y1, alpha[1,jy,ix])
                        fitted_CDF2 = stats.gamma.cdf(y2, alpha[2,jy,ix])
                        fitted_CDF = weights[0,jy,ix]*fitted_CDF0 + weights[1,jy,ix]*fitted_CDF1 + \
                            weights[2,jy,ix]*fitted_CDF2
                        print ('   fitted_CDF = ',fitted_CDF)
                        Dnstat = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
                        print ('   np.abs(fitted_CDF - empirical_CDF) = ', np.abs(fitted_CDF - empirical_CDF))
                        print ('   Dnstat for 3 gamma mixture = ', Dnstat)
            
                except: # bombed off with 2-gamma mixture.  Try three.
                        
                    try:
                        result_R = mixtools.gammamixEM(precip_nonzero_R, k=3, maxit=1000) #, 
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        #print ('   After mixtools.gammamixEM, jy,ix = ', jy,ix, nxin*nyin, current_time)
                        result_np = np.asarray(result_R, dtype=object)
                        result_weights_np = np.asarray(result_R[1])
                        result_alpha_beta_np = np.asarray(result_R[2])
                        weights[:,jy,ix] = result_weights_np[:]
                        alpha[:,jy,ix] = result_alpha_beta_np[0,:]
                        beta[:,jy,ix] = result_alpha_beta_np[1,:]
                        print ('   weights = ', weights[:,jy,ix])
                        print ('   alpha = ', alpha[:,jy,ix])
                        print ('   beta = ', beta[:,jy,ix])
            
                        y0 = empirical_precipvals / beta[0,jy,ix]
                        y1 = empirical_precipvals / beta[1,jy,ix]
                        y2 = empirical_precipvals / beta[2,jy,ix]
                        fitted_CDF0 = stats.gamma.cdf(y0, alpha[0,jy,ix])
                        fitted_CDF1 = stats.gamma.cdf(y1, alpha[1,jy,ix])
                        fitted_CDF2 = stats.gamma.cdf(y2, alpha[2,jy,ix])
                        fitted_CDF = weights[0,jy,ix]*fitted_CDF0 + weights[1,jy,ix]*fitted_CDF1 + \
                            weights[2,jy,ix]*fitted_CDF2
                        print ('   fitted_CDF = ',fitted_CDF)
                        Dnstat = np.max(np.abs(fitted_CDF[i1:] - empirical_CDF[i1:]))
                        print ('   np.abs(fitted_CDF - empirical_CDF) = ', np.abs(fitted_CDF - empirical_CDF))
                        print ('   Dnstat for 3 gamma mixture = ', Dnstat)
                    except: # bombed off with 3-parameter.   Revert to one.
                        weights[0,jy,ix] = 1.0
                        alpha[0,jy,ix] = alpha_hat
                        beta[0,jy,ix] = beta_hat
                        weights[1:,jy,ix] = 0.0
                        alpha[1:,jy,ix] = 1.0
                        beta[1:,jy,ix] = 1.0
                    
            else:
                weights[0,jy,ix] = 1.0
                alpha[0,jy,ix] = alpha_hat
                beta[0,jy,ix] = beta_hat
                weights[1:,jy,ix] = 0.0
                alpha[1:,jy,ix] = 1.0
                beta[1:,jy,ix] = 1.0
            
        else: 
            
            # --- very few positive samples; fit an simple maximum likelihood distribution
            
            pmean = np.mean(precip_ens_1d)
            if pmean == 0.0:  # not long enough sample to find some nonzero precip amounts
                alpha_hat = 1.0
                beta_hat = 100.0
            else:  # use maximum likelihood method of Wilks text to fit non-mixture distribution.
                lnxbar = np.log(pmean)
                meanlnxi = np.mean(np.sum(np.log(precip_ens_1d)))
                D = lnxbar - meanlnxi
                alpha_hat = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
                beta_hat = pmean / alpha_hat
                weights[0,jy,ix] = 1.0
                alpha[0,jy,ix] = alpha_hat
                beta[0,jy,ix] = beta_hat
                weights[1:,jy,ix] = 0.0
                alpha[1:,jy,ix] = 1.0
                beta[1:,jy,ix] = 1.0

# --- save to cPickle file

outfile = master_directory + cmonth+ '_apcp_gamma_parameters_h' + clead + '.cPick'
print ('writing to ', outfile)
ouf = open(outfile, 'wb')
cPickle.dump(weights, ouf)
cPickle.dump(alpha, ouf)
cPickle.dump(beta, ouf)
cPickle.dump(fzero, ouf)
ouf.close()



