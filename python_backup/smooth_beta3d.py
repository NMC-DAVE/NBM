def smooth_beta3d(ndates, nlats, nlons, beta_3d):
    import numpy as np
    import scipy.ndimage    
    sigma = 3.0
    beta_3d_smoothed = np.zeros((ndates,nlats,nlons), dtype=np.float64)
    for idate in range(ndates):
        b_input = beta_3d[idate,:,:]
        b_output = scipy.ndimage.filters.gaussian_filter(b_input, sigma=sigma, mode='reflect')
        beta_3d_smoothed[idate,:,:] = b_output
    
    return beta_3d_smoothed