import numpy as np
import scipy.stats as stats
from crps import crps



ens = stats.gamma.rvs(a=1.0, scale=10.0, size=31)
observed = np.zeros((1,1), dtype=np.float64)
observed[0,0] = ens[30]
ensemble = np.zeros((30,1,1), dtype=np.float64)
ensemble[:,0,0] = ens[0:30]
lats = np.zeros((1,1), dtype=np.float64)
lats[0,0] = 30.0

nmembers = 30
ny = 1
nx = 1
thresholds = np.arange(1,50.1,0.1)
nthresh = len(thresholds)


CRPS_domain_average, CRPS_by_gridpoint =  \
    crps (nmembers, ny, nx,  nthresh, \
    thresholds, lats, ensemble, observed)
    