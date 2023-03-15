import random
import numpy as np
import math
import scipy.stats as stats

# ---- initialize the ensemble

nmembers = 10
niters = 10
alpha = 0.2
mu = 0.0
sigma = 1.0
R = 1.0

xa = np.zeros((nmembers,niters))
xb = np.zeros((nmembers,niters))
#random.seed([-12345])
for imem in range(nmembers):
    xa[imem,0] = random.gauss(mu,sigma)

# ---- now iterate Jeff Anderson's simple 1-d model with nonlinearity, from
# A Non-Gaussian Ensemble Filter Update for Data Assimilation, 2010, MWR to appear

for iter in range(1,niters):
    xb[:,iter] = xa[:,iter-1] + 0.05*(xa[:,iter-1]+alpha*xa[:,iter-1]**2)
    #print 'xa = ',iter-1, xa[:,iter-1]
    #print 'xb = ',iter, xb[:,iter-1]

    obs = random.gauss(0.0, R)

    # calculate first two moments of xb = (xbmean, Pb)

    xb1d = xb[:,iter]
    print xb1d
    xbmean = np.mean(xb1d)
    xbprime = xb1d - xbmean
    stddev = stats.std(xb1d,axis=none)
    Pb = stddev**2                
                                    
    # calculate Kalman gain

    K = Pb / (Pb+R)
    alpha2 = 1.0 / (1.0 + math.sqrt(R/(Pb+R)))

    # ---- perform EnSRF

    xamean = xbmean + K*(obs-xbmean)
    xaprime = xbprime - alpha2*K*xbprime

    # ---- reform analysis ensemble

    xa[:,iter] = xamean + xaprime

print 'Done!'
