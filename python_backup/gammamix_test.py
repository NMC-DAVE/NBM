"""
gammamix_test.py

"""

import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
from netCDF4 import Dataset
import _pickle as cPickle
import scipy.stats as stats
from gammamix import gammamix_em

infile = 'precipitation_samples.cPick'
inf = open(infile,'rb')
precip_ens_nonzero = cPickle.load(inf)
inf.close()
print ('len = ', len(precip_ens_nonzero))
print ('precip_ens_nonzero[0], [-1] = ', \
    precip_ens_nonzero[0], precip_ens_nonzero[-1])

best_results = gammamix_em(precip_ens_nonzero, \
    mix_prop=None, alpha=None, invbeta=None,
    k=3, epsilon=0.001, maxit=1000,
    maxrestarts=20, verb=False)
    
params = best_results.params
weights = params.mix_prop
alpha = params.alpha
beta = params.invbeta


nz = len(precip_ens_nonzero)
query_these_indices = [ nz//20, nz//10, (3*nz)//20, nz//5, nz//4, (3*nz)//10, \
    (7*nz)//20, (2*nz)//5, (9*nz)//20, nz//2, (11*nz)//20, (3*nz)//5, (13*nz)//20, \
    (7*nz)//10, (3*nz)//4, (4*nz)//5, (17*nz)//20, (35*nz)//40,(9*nz)//10, \
    (37*nz)//40, (19*nz)//20, (39*nz)//40, (79*nz)//80]
print ('   nz, query_these_indices = ', nz, query_these_indices)
nq = len(query_these_indices)
    
empirical_CDF = np.array(query_these_indices, dtype=np.float) / float(nz)     
empirical_precipvals = precip_ens_nonzero[query_these_indices] 

# ---- build the CDF at the precipitation values associated with the quantiles in empirical_CDF

y0 = empirical_precipvals / beta[0]
y1 = empirical_precipvals / beta[1]
y2 = empirical_precipvals / beta[2]
fitted_CDF0 = stats.gamma.cdf(y0, alpha[0])
fitted_CDF1 = stats.gamma.cdf(y1, alpha[1])
fitted_CDF2 = stats.gamma.cdf(y2, alpha[2])
fitted_CDF = weights[0]*fitted_CDF0 + weights[1]*fitted_CDF1 + \
    weights[2]*fitted_CDF2
    
print (' i    precipvals[i], empirical_CDF,    fitted_CDF')
for i in range(nq):
    print (i, empirical_precipvals[i], empirical_CDF[i], fitted_CDF[i])

# ---- build the distribution of fitted precipitation amounts associated with 
#      a particular quantile of the distribution, i.e., the quantile function.


quantile_function = np.zeros(nz, dtype=np.float32)
print ('index   fitted precip   empirical precip')
for iq in range(nq):
    quantile = empirical_CDF[iq] # quantile of interest nz//20, etc.
    idx = np.where(fitted_CDF > quantile) [0] # first instance 
    is_empty = idx.size == 0
    if is_empty == False:
        nidx = len(idx)
        if nidx>0 :
            idx2 = idx[0]
        else:
            idx2 = idx
        if idx2 > 0:
            pbelow = fitted_CDF[idx2-1]
            pabove = fitted_CDF[idx2]            
            frac = (empirical_CDF[iq]- pbelow)/ (pabove-pbelow)
            quantile_function[iq] = (1.-frac)*empirical_precipvals[idx2-1] + \
                frac*empirical_precipvals[idx2]
        else:
            quantile_function[iq] = 0.0 
    print (iq, '   ', quantile_function[iq], '   ', empirical_precipvals[iq])       
    
    