def det_params_gaussmix(ndates_valid, nlats, nlons, data_input):

    """ try to fit a Gaussian mixture model of two normal distributions to 
         the data at hand, and return the parameters.
    """

    import numpy as np
    import scipy
    import sys
    from sklearn import mixture

    weights = np.zeros((3,nlats,nlons), dtype=np.float64)
    means = np.zeros((3,nlats,nlons), dtype=np.float64)
    stddevs = np.zeros((3,nlats,nlons), dtype=np.float64)
    X = np.zeros((ndates_valid,1), dtype=np.float64)

    # ---- Fit 2 gaussian distribution mixture. Return parameters
    
    print ('determining 3-component Gaussian mixture')
    for ilat in range(nlats):
        print ('ilat = ', ilat)
        for ilon in range(nlons):
            X[:,0] = data_input[:,ilat,ilon]
            clf = mixture.GaussianMixture(n_components=3,\
                covariance_type='spherical',init_params='kmeans', n_init=5)
            #clf = mixture.GaussianMixture(n_components=3,\
            #    covariance_type='spherical',init_params='random',n_init=5)
            clf.fit(X)
            w = clf.weights_
            m = clf.means_
            s = np.sqrt(clf.covariances_)
            if ilat == nlats//2 and ilon == nlons//2:
                print ('overall mean, stddev ',np.mean(X), np.std(X))
                print ('weights = ', w[:])
                print ('means = ', m[:,0])
                print ('stdevs = ',s[:])
            weights[:,ilat,ilon] = w[:]
            means[:,ilat,ilon] = m[:,0]
            stddevs[:,ilat,ilon] = s[:]
    
    return weights, means, stddevs