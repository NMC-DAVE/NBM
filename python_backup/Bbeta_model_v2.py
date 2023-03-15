def Bbeta_model_v2(nlats, nlons, ndates, npts, cyear, clead,\
    warm_or_cold, cpath_Bbeta, bias_corr_3d, already):
    
    """
    form a model for the bias-correction error covariance, 
    using a previously calculated correlation model based on 
    horizontal, vertical distance, land/water.   Covariances
    then formed using beta errors determined grid point by
    grid point from time series of previous decay avg approach
    """

    import numpy as np
    import os, sys
    import _pickle as cPickle

    # --- read correlation of bias to file

    cfile = 'correlation_bias_ERA5grid.cPick'
    inf = open(cfile, 'rb')
    correlation_bias = cPickle.load(inf)
    inf.close()

    # ---- get standard deviation across time dimension
    
    beta_stddev = np.std(bias_corr_3d, axis=0)
    Bbeta = np.zeros((nlats*nlons, nlats*nlons), dtype=np.float32)
    
    ktr1 = 0
    for i1 in range(nlons):
        for j1 in range(nlats):
            ktr2 = 0
            for i2 in range(nlons):
                for j2 in range(nlats):
                    Bbeta[ktr1,ktr2] = beta_stddev[j1,i1]*beta_stddev[j2,i2]*\
                        correlation_bias[j1,i1,j2,i2]
                    ktr2 = ktr2+1
            ktr1 = ktr1+1
        
    return Bbeta