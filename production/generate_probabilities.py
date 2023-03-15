def generate_probabilities(nmembers, thresholds, precip_mean, precip_ens, \
    prob, ny, nx, hist_thresholds, closest_hist, dressit, isitunweighted, \
    b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, \
    b1_std_lowrank, b0_mean_midrank, b1_mean_midrank, \
    b0_std_midrank, b1_std_midrank, b0_mean_highrank, \
    b1_mean_highrank, b0_std_highrank, b1_std_highrank) : 
    
    
    """ determine the ensemble probability of exceeding the threshold,
        weighted by closest-member histograms.   Simplify calculations
        and eliminate sums by using CDF of closest-member histograms. 
        Calls a fortran routine for computational efficiency."""
    
    import numpy as np
    from weighted_probability_dressed_f90 import weighted_probability_dressed_f90
    from weighted_probability_f90 import weighted_probability_f90
    from numba import jit
    
    # ---- set weights based on the ensemble-mean amount
    
    nchcats = len(hist_thresholds)+1
    nthresholds = len(thresholds)
    prob = np.zeros((nthresholds, ny, nx), dtype=np.float64)
    zeros = np.zeros((ny, nx), dtype=np.float64)
    ones = np.ones((ny, nx), dtype=np.float64)
    
    if isitunweighted == True:
        prob[:,:,:] = 0.0
        for ithresh, thresh in enumerate(thresholds):
            for imem in range(nmembers):
                onezero = np.where(precip_ens[imem,:,:] >= thresh, ones, zeros)
                prob[ithresh,:,:] = prob[ithresh,:,:] + onezero[:,:]
        prob = prob / float(nmembers)
    else:    
        #print ('   np.shape(precip_ens) = ',np.shape(precip_ens))
        #print ('   nmembers = ', nmembers)
        #print ('   closest_hist.dtype = ', closest_hist.dtype)
        #print ('   np.shape(closest_hist) = ',np.shape(closest_hist)  )
        #print ('   weighted_probability_dressed_f90.__doc__', weighted_probability_dressed_f90.__doc__ )
        #print ('   ny, nx, nmembers, nthresholds, nchcats = ', ny, nx, nmembers, nthresholds, nchcats)
        #print ('   precip_mean.dtype = ', precip_mean.dtype)
        #print ('   np.shape(precip_mean) = ', np.shape(precip_mean))
        #print ('   precip_ens.dtype = ', precip_ens.dtype)
        #print ('   thresholds.dtype = ', thresholds.dtype)
        #print ('   hist_thresholds.dtype = ', hist_thresholds.dtype)
        #print ('   b0_mean_lowrank.dtype = ', b0_mean_lowrank.dtype)
        dress_array = np.array([b0_mean_lowrank, b1_mean_lowrank, \
            b0_std_lowrank, b1_std_lowrank, b0_mean_midrank, \
            b1_mean_midrank, b0_std_midrank, b1_std_midrank, \
            b0_mean_highrank, b1_mean_highrank, b0_std_highrank, \
            b1_std_highrank], dtype=np.float32)
        ndress = 12
        if dressit == True:
            prob = weighted_probability_dressed_f90(precip_mean, precip_ens, \
                thresholds, hist_thresholds, closest_hist, \
                dress_array, ny, nx, nmembers, nthresholds, ndress, nchcats)
                
        else:
            prob = weighted_probability_f90(precip_mean, \
                precip_ens, thresholds, hist_thresholds, \
                closest_hist, ny, nx, nmembers, nthresholds, nchcats)        
                
    a = np.where(prob > 1.0)
    if a: prob[a] = 1.0
    return prob