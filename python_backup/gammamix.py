##
#Python implementation of mixture of gamma distributions
# Based on the mixtools R package implementation:
# https://github.com/cran/mixtools/blob/master/R/gammamixEM.R
# Notebook at https://github.com/kundajelab/tfmodisco/blob/master/examples/
#             mixture_of_gammas/Mixture%20of%20Gamma%20Distributions.ipynb
##

from __future__ import division, print_function
import numpy as np
from collections import namedtuple
from scipy.stats import gamma
from scipy.optimize import minimize
from scipy.special import digamma
import sys

GammaMixParams = namedtuple("MixParams", ["mix_prop", "alpha", "invbeta", "k"])
GammaMixResult = namedtuple("GammaMixResult", ["params",
                                               "ll", "iteration",
                                               "expected_membership"]) 

#@jit
def gammamix_init(x, mix_prop=None, alpha=None, invbeta=None, k=2):
    """ the original that came with the web code.   not used now.
        didn't like the method of moments estimation. """
    n = len(x)
    if (mix_prop is None):
        mix_prop = np.random.random((k,)) 
        mix_prop = mix_prop/np.sum(mix_prop)
    else:
        k = len(mix_prop)

    if (k==1):
        x_bar = np.array([np.mean(x)])
        x2_bar = np.array([np.mean(np.square(x))])
    else:
        #sort the values
        x_sort = sorted(x) 
        #figure out how many go in each mixing
        #component based on the current mixing
        #parameters
        ind = np.floor(n*np.cumsum(mix_prop)).astype("int")
        #collect the values corresponding to each
        #component to compute the initial alpha and beta
        x_part = []
        x_part.append(x_sort[0:ind[0]])
        for j in range(1,k):
            x_part.append(x_sort[ind[j-1]:ind[j]])
        x_bar = np.array([np.mean(y) for y in x_part])
        x2_bar = np.array([np.mean(np.square(y)) for y in x_part])

    if (alpha is None):
        alpha = np.square(x_bar)/(x2_bar - np.square(x_bar))

    if (invbeta is None):
        invbeta = x_bar/(x2_bar - np.square(x_bar))

    return GammaMixParams(mix_prop=mix_prop,
                          alpha=alpha,
                          invbeta=invbeta, k=k)
                          
#@jit                          
def gammamix_init2(x, mix_prop=None, alpha=None, invbeta=None, k=2):
    """ keeps random aspect of gammamix_init, but uses Thom maximum
    likelihood estimator, per Wilks text. """
    n = len(x)
    if (mix_prop is None):
        mix_prop = 1.0/float(k) + np.random.uniform(low=0.0, high=0.5, size=k)
        mix_prop = mix_prop/np.sum(mix_prop)
        #if k == 2:
        #    mix_prop = np.array([0.8, 0.2])
        #elif k == 3:
        #    mix_prop = np.array([0.6, 0.3, 0.1])
        
        #print ('mix_prop = ', mix_prop)
        alpha = np.zeros((k), dtype=np.float64)
        invbeta = np.zeros((k), dtype=np.float64)
    else:
        k = len(mix_prop)
    mixpropsum = 0.
    ibegin = np.zeros((k), dtype=np.int32)
    iend = np.zeros((k), dtype=np.int32)
    for i in range(k):
        #print ('i = ',i)
        if i == 0:
            mixpropsum = mix_prop[0]
            ibegin[i] = 0
        else:
            mixpropsum = mixpropsum + mix_prop[i]
            ibegin[i] = iend[i-1]
        iend[i] = np.min([int(mixpropsum*float(n)), n])
        pmean = np.mean(x[ibegin[i]:iend[i]])
        lnxbar = np.log(pmean)
        meanlnxi = np.mean(np.log(x[ibegin[i]:iend[i]]))
        D = lnxbar - meanlnxi
        alpha[i] = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
        invbeta[i] = pmean / alpha[i]
    
    return GammaMixParams(mix_prop=mix_prop,
        alpha=alpha, invbeta=invbeta, k=k)   
        
        
def gammamix_init3(x, mix_prop=None, alpha=None, invbeta=None, k=2):
    """ Tom Hamill's home-grown, kludged algorithm that tries to:
    (1) provide arbitrary deterministic sample weights to a 3-parameter mix, and
    (2) use the Thom ML estimator. """
    n = len(x)
    if (mix_prop is None):
        alpha = np.zeros((k), dtype=np.float64)
        invbeta = np.zeros((k), dtype=np.float64)
        weights = np.ones((n,k),dtype=np.float64)
        if k == 2:
            mix_prop = np.array([0.8,0.2])
            weights[0:n//2,0] = 1.75
            weights[n//2:n,0] = 0.25
            weights[0:n//2,1] = 0.25
            weights[n//2:n,1] = 1.75
        elif k == 3:
            mix_prop = np.array([0.6,0.3,0.1])
            weights[0:n//3,0] = 1.75
            weights[n//3:(2*n)//3,0] = 1.0
            weights[(2*n)//3:n,0] = 0.25
            weights[0:n//3,1] = 0.5
            weights[n//3:(2*n)//3,1] = 2.0
            weights[(2*n)//3:n,1] = 0.5
            weights[0:n//3,2] = 0.25
            weights[n//3:(2*n)//3,2] = 1.0
            weights[(2*n)//3:n,2] = 1.75
            
    else:
        k = len(mix_prop)
    mixpropsum = 0.

    for i in range(k):
        if i == 0:
            mixpropsum = mix_prop[0]
        else:
            mixpropsum = mixpropsum + mix_prop[i]
        pmean = np.mean(weights[:,i]*x[:])
        lnxbar = np.log(pmean)
        meanlnxi = np.mean(np.log(weights[:,i]*x[:]))
        D = lnxbar - meanlnxi
        alpha[i] = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
        invbeta[i] = pmean / alpha[i]
    
    return GammaMixParams(mix_prop=mix_prop,
        alpha=alpha, invbeta=invbeta, k=k)           
                               

#@jit
def gamma_component_pdfs(x, theta, k):
    component_pdfs = []
    alpha = theta[0:k]
    invbeta = theta[k:2*k]
    for j in range(k):
        component_pdfs.append(gamma.pdf(x=x, a=alpha[j], scale=invbeta[j])) 
    component_pdfs = np.array(component_pdfs)
    return component_pdfs

#@jit
def log_deriv_gamma_component_pdfs(x, theta, k):
    log_deriv_alpha_component_pdfs = []
    log_deriv_invbeta_component_pdfs = []
    alpha = theta[0:k]
    invbeta = theta[k:2*k]
    for j in range(k):
        log_deriv_invbeta_component_pdfs.append(
            (x/(invbeta[j]**2) - alpha[j]/invbeta[j]))
        log_deriv_alpha_component_pdfs.append(
            (np.log(x) - np.log(invbeta[j]) - digamma(alpha[j])))
    return (np.array(log_deriv_invbeta_component_pdfs),
            np.array(log_deriv_alpha_component_pdfs))

#@jit
def gamma_ll_func_to_optimize(theta, x, expected_membership, mix_prop, k):
    component_pdfs = gamma_component_pdfs(x=x,
                                          theta=theta, k=k)
    if (np.isnan(np.sum(component_pdfs))):
        assert False
    #prevent nan errors for np.log
    component_pdfs = component_pdfs+((component_pdfs == 0)*1e-32)
    #log likelihood
    #if np.min(mix_prop) < 0.001 :
    #    print ('****** gamma_ll_func_to_optimize: mix_prop = ',mix_prop)
    ll =  -np.sum(expected_membership*np.log(
                  mix_prop[:,None]*component_pdfs))
    #log deriv gamma component pdfs
    (log_deriv_invbeta_component_pdfs,
     log_deriv_alpha_component_pdfs) =\
     log_deriv_gamma_component_pdfs(x=x, theta=theta, k=k) 

    log_derivs = np.array(
        list(-np.sum(
             expected_membership
             *log_deriv_alpha_component_pdfs, axis=1))+
        list(-np.sum(
          expected_membership
          *log_deriv_invbeta_component_pdfs, axis=1)))
    
    return ll, log_derivs
                                                          

# -- based on https://github.com/cran/mixtools/blob/master/R/gammamixEM.R

def gammamix_em(x, mix_prop=None, alpha=None, invbeta=None,
                k=2, epsilon=0.001, maxit=1000,
                maxrestarts=20, progress_update=20, verb=False):

    #initialization
    x = np.array(x) 
    if mix_prop is None:
        mix_prop, alpha, invbeta, k = \
            gammamix_init3(x=x, mix_prop=mix_prop, alpha=alpha,\
            invbeta=invbeta, k=k) 
    if verb is True:
        print("   initial vals:",mix_prop, alpha, invbeta, k) 
        sys.stdout.flush()
    theta = np.concatenate([alpha, invbeta],axis=0)

    
    iteration = 0
    mr = 0
    diff = epsilon + 1
    n = len(x)
    
    old_obs_ll = np.sum(np.log(np.sum(
                    mix_prop[:,None]*gamma_component_pdfs(
                        x=x,
                        theta=theta, k=k), axis=0))) 
    ll = [old_obs_ll]

    best_result = None
    best_obs_ll = old_obs_ll

    while ((np.abs(diff) > epsilon) and (iteration < maxit)):
        dens1 = mix_prop[:,None]*gamma_component_pdfs(
                                     x=x,
                                     theta=theta, k=k)
        expected_membership = dens1/np.sum(dens1, axis=0)[None,:] 
        mix_prop_hat = np.mean(expected_membership, axis=1)
        minimization_result = minimize(
            fun=gamma_ll_func_to_optimize,
            x0=theta,
            bounds=[(1e-7,None) for t in theta],
            args=(x, expected_membership, mix_prop, k),
            jac=True) 
        if (minimization_result.success==False):
            print("   Choosing new starting values")
            if (mr==maxrestarts):
                raise RuntimeError("Try a different number of components?") 
            mr += 1 
            mix_prop, alpha, invbeta, k = gammamix_init2(x=x, k=k) 
            theta = np.concatenate([alpha, invbeta],axis=0)
            iteration = 0
            diff = epsilon + 1
            old_obs_ll = np.sum(np.log(np.sum(
                            mix_prop[:,None]*gamma_component_pdfs(
                                x=x,
                                theta=theta, k=k), axis=0))) 
            ll = [old_obs_ll]
        else:
            theta_hat = minimization_result.x 
            alpha_hat = theta_hat[0:k]
            invbeta_hat = theta_hat[k:2*k]


            new_obs_ll = np.sum(np.log(np.sum(
                          mix_prop_hat[:,None]*gamma_component_pdfs(
                            x=x,
                            theta=theta_hat, k=k),axis=0))) 
            diff = new_obs_ll - old_obs_ll
            old_obs_ll = new_obs_ll
            ll.append(old_obs_ll)

            mix_prop = mix_prop_hat
            theta = theta_hat
            alpha = alpha_hat
            invbeta = invbeta_hat
            iteration = iteration + 1

            if (old_obs_ll >= best_obs_ll):
                best_result = GammaMixResult(
                    params=GammaMixParams(mix_prop=mix_prop,
                                          alpha=alpha, invbeta=invbeta, k=k),
                    ll=ll,
                    iteration=iteration,
                    expected_membership=expected_membership)
                best_obs_ll = old_obs_ll
                #if verb:
                #    print("New best!") 
                #    print(GammaMixParams(mix_prop=mix_prop,
                #                          alpha=alpha,
                #                          invbeta=invbeta, k=k))

            if verb:
                if (iteration%progress_update == 0):
                    print("iteration =", iteration,
                          "log-lik diff =", diff,
                          " log-lik =", new_obs_ll) 
                    sys.stdout.flush()

    if (iteration == maxit):
        print("   WARNING! NOT CONVERGENT!")
    print("   Number of iterations=", iteration)

    return best_result



