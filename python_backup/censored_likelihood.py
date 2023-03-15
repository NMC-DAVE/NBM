import numpy as np
import scipy as sp
import math
import os, sys
import matplotlib.pyplot as plt

from numpy import ma
from numpy.random import random_sample
from scipy import stats
from scipy.stats import gamma
from scipy.special import loggamma
from scipy.interpolate import interp1d
from scipy.optimize import minimize


####  Simulate from a 2-component gamma distribution  ###
#
w1 = 0.4
w2 = 0.6
shape1 = 0.2
shape2 = 0.6
scale1 = 3.5
scale2 = 2.2

n = 1000

cmp = np.random.binomial(1,w2,size=n)
n1 = np.sum(cmp==0)
n2 = np.sum(cmp==1)
x = np.zeros(n, dtype=np.float32)
x[cmp==0] = np.random.gamma(shape1,scale1,n1)
x[cmp==1] = np.random.gamma(shape2,scale2,n2)

#
##########################################################



## Fit single component model

pmean = np.mean(x)
lnxbar = np.log(pmean)
meanlnxi = np.mean(np.log(x))
D = lnxbar - meanlnxi
alpha0 = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
invbeta0 = pmean / alpha0


## Fit single component only to upper 20% data, ignoring the censoring issue

c = np.percentile(x,80)

pmeanC1 = np.mean(x[x>c])
lnxbarC1 = np.log(pmeanC1)
meanlnxiC1 = np.mean(np.log(x[x>c]))
D = lnxbarC1 - meanlnxiC1
alphaC1 = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
invbetaC1 = pmeanC1 / alphaC1



## Fit single component only to upper 20% data using censored maximum likelihood

def loglik(par,x_ucs,c,n_usc,n_cs):     # censored, negative log-likelihood function
    T1 = -n_cs * np.log(gamma.cdf(c,par[0],scale=par[1]))
    T2 = np.sum(x_ucs)/par[1] - (par[0]-1)*np.sum(np.log(x_ucs)) + n_ucs*(loggamma(par[0])+par[0]*np.log(par[1]))
    return T1 + T2

bnds = ((0.001,1.5), (0.001,3.0*invbeta0))
pst = [alpha0,invbeta0]                       # start with the parameters fitted to the full data set

x_ucs = x[x>c]
n_ucs = len(x_ucs)
n_cs = len(x) - n_ucs

par_opt = minimize(loglik, pst, args=(x_ucs,c,n_ucs,n_cs), method='L-BFGS-B', bounds=bnds, tol=1e-6).x

alphaC2 = par_opt[0]
invbetaC2 = par_opt[1]

##  Make some plots

prob = np.arange(1,n+1)/(n+1)
xgrd = np.arange(0,15,0.1)

# Plot empirical and fitted CDFs
plt.figure()
plt.scatter(np.sort(x), prob)
plt.plot(xgrd, gamma.cdf(xgrd,alpha0,scale=invbeta0), c='m',label='Single Gamma')
plt.plot(xgrd, gamma.cdf(xgrd,alphaC1,scale=invbetaC1), c='r',label='Single to upper 20%, ignore censor')
plt.plot(xgrd, gamma.cdf(xgrd,alphaC2,scale=invbetaC2), c='b',label='Single to upper 20%, censor')
plt.legend(loc=0)
plt.show()


# Q-Q Plots
plt.figure()
plt.plot(xgrd, xgrd, c='k')
plt.scatter(gamma.ppf(prob,alpha0,scale=invbeta0), np.sort(x), c='m',label='Single Gamma')
plt.scatter(gamma.ppf(prob,alphaC1,scale=invbetaC1), np.sort(x), c='r',label='Single to upper 20%, ignore censor')
plt.scatter(gamma.ppf(prob,alphaC2,scale=invbetaC2), np.sort(x), c='b',label='Single to upper 20%, censor')
plt.legend(loc=0)
plt.show()




