def paired_bootstrap(x,y):
    """ given vectors x and y, return a paired bootstrap estimate of the 5th and 95th percentiles
    of the resampled distribution"""
    import numpy as np

    nresa = 100   # 10000 resamplings
    nelts = np.size(x)  # what is the size of the vector x
    sum1 = np.zeros((nresa),dtype=np.float)
    sum2 = np.zeros((nresa),dtype=np.float)
    dist = np.zeros((nresa),dtype=np.float)
    ones = np.ones((nelts),dtype=np.float)
    zeros = np.zeros((nelts),dtype=np.float)
    for i in range(nresa):
        x0 = np.random.rand(nelts)
        iusex = np.where(x0 < 0.5, ones, zeros)
        iusey = 1-iusex
        sum1[i] = np.mean(x*iusey + y*iusex)
        sum2[i] = np.mean(x*iusex + y*iusey)
        dist[i] = sum1[i]-sum2[i]
        #print 'sum1, sum2, dist=',sum1[i],sum2[i],dist[i]

    dsort = np.sort(dist)
    d05 = dsort[np.int(.05*np.float(nresa))]
    d95 = dsort[np.int(.95*np.float(nresa))]
    #print 'd05, d95 = ',d05,d95
    return d05, d95