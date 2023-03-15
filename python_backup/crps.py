def crps (nmembers, ny, nx,  nthresh, thresholds, lats, ensemble, observed):
    """ 
    
    crps:  computes CRPS score from ensemble for one specific case.  
    The routine has been tailored here for precipitation.  Assumes lat-lon
    grid input.
    
    Inputs:
    ------
    nmembers: number of members
    ny: number of grid points in y direction
    nx: number of grid points in x direction
    nthresh
    lats[ny]: latitudes so we can perform cos(latitude) weighting of samples
        to allow for different grid areas
    nthresh: number of thresholds over which to integrate CRPS
    thresholds[nthresh]: values of the precipitation thresholds. 
         Start with nonzero value.
    lats[ny]: in degrees
    ensemble[nmembers, ny, nx]: ensemble precipitation values, assumed mm
    observed[ny, nx]: observed/analyzed precipitation amount    
    
    returns:
    ------- 
    CRPS_domain_average: the area-weighted average CRPS over the whole domain
    CRPS_by_gridpoint: (ny,nx) array of CRPS values for each grid point
    
    """
    
    import numpy as np
    
    # ---- define variables, arrays

    FCDF = np.zeros((nthresh), dtype=np.float64)
    OCDF = np.zeros((nthresh), dtype=np.float64)
    ones = np.ones((nmembers), dtype=np.float64) 
    zeros = np.zeros((nmembers), dtype=np.float64)
    pi = 3.141592654
    CRPS_by_gridpoint = np.zeros((ny,nx), dtype=np.float64)

    # ---- compute CRPS at each grid point + totals needed to compute
    #      area weighted domain average.

    ensemble_sorted = np.sort(ensemble, axis=0)
    #print ('ensemble_sorted[:,0,0] = ', ensemble_sorted[:,0,0])
    #print ('observed[0,0] = ', observed[0,0])
    CRPS_sum = 0.0
    CRPS_sum_coslat = 0.0
    for jy in range(ny):
        coslat = np.cos(2.0*pi*lats[jy]/360.)
        #print ('coslat = ', coslat)
        for ix in range(nx):
            #print ('ithresh, thresh, FCDF, ACDF, CRPS_by_gridpoint = ')
            ensemble_sorted_1d = ensemble_sorted[:,jy,ix]
            for ithresh, thresh in enumerate(thresholds):
                if ithresh == 0:
                    dp = thresholds[0]
                else:
                    dp = thresholds[ithresh] - thresholds[ithresh-1]
                    
                    
                # --- integrate to compute CRPS following Wilks 2011 text,
                #     section 2.5, eq. 8.54a
                    
                a = np.where(thresh < ensemble_sorted_1d, zeros, ones)
                FCDF = np.sum(a) / float(nmembers)
                if thresh < observed[jy,ix]  :
                    ACDF = 0.0
                else:
                    ACDF = 1.0
                CRPS_by_gridpoint[jy,ix] = CRPS_by_gridpoint[jy,ix] + dp*(FCDF-ACDF)**2
                #print (ithresh, thresh, FCDF, ACDF, CRPS_by_gridpoint[jy,ix])
                
            CRPS_sum = CRPS_sum + CRPS_by_gridpoint[jy,ix]*coslat
            CRPS_sum_coslat = CRPS_sum_coslat + coslat
            
    CRPS_domain_average = CRPS_sum / CRPS_sum_coslat 
    
    return CRPS_domain_average, CRPS_by_gridpoint
                    
                