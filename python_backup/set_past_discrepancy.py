def set_past_discrepancy(idate, ilead, ndates, nlats, nlons, \
    analyses_3d, forecast_3d, ones, efold_timescale):
    
    import numpy as np
    
    idate_earliest = idate - 3*efold_timescale
    if idate_earliest < 0: idate_earliest = 0
    numer = np.zeros((nlats,nlons), dtype=np.float64)
    denom = np.zeros((nlats,nlons), dtype=np.float64)
    if idate+1-ilead > 0:
        for jdate in range(idate_earliest, idate+1-ilead):
            weight = np.exp(- (jdate - (idate-ilead))/efold_timescale)
            numer = numer + weight*(analyses_3d[jdate,:,:]-forecast_3d[jdate,:,:])
            denom = denom + weight*ones[:,:]
        past_discrepancy = numer/denom
    else:
        past_discrepancy = numer # effectively zero
        
    return past_discrepancy