import numpy as np
import scipy.stats as stats

def crps_calc(x,cdf_f,cdf_o):
    n = len(x)
    dx = x[1]-x[0]
    crps = np.sum((1.0/n)*dx*(cdf_f - cdf_o)**2)
    return crps
    

x = np.arange(0.0,1.001,0.01)

cdf_fdeterministic = np.zeros((101), dtype=np.float32)
cdf_fdeterministic[50:] = 1.0  # deterministic forecast at 0.5
print (cdf_fdeterministic[0:-1:5])

cdf_observation = np.zeros((101), dtype=np.float32)
cdf_observation[20:] = 1.0 # obs at 0.2
print (cdf_observation[0:-1:5])

cdf_fprobabilistic = stats.norm.cdf(x, loc=0.5, scale = 0.15)
print (cdf_fprobabilistic[0:-1:5])

crps_deterministic = crps_calc(x,cdf_fdeterministic,cdf_observation)
crps_probabilistic = crps_calc(x,cdf_fprobabilistic,cdf_observation)
print ('CRPS deterministic, probabilistic = ', crps_deterministic, crps_probabilistic)