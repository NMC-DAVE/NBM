    
def forecast_error_covariance(nlats, nlons, lats, lons, \
    difference_3d, efold, exponenty):
    
    """ 
    forecast_error_covariance.py
    input is 3D array of bias-corrected forecasts (or bias corrections) across US. 
    Output is array of localized forecast-error covariances 
    (or bias-correction covariances), needed in Kalman filter.
    """
    
    import numpy as np
    from datetime import datetime
    
    def set_coslat(lats):
        coslat = np.cos(lats[:,0]*3.1415926/180.)
        return coslat
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('forecast_error_covariance.py.  Current time = ', current_time)
    npts = nlats*nlons
    Bx = np.zeros((npts, npts), dtype=np.float64)
    Bx_localized = np.zeros((npts, npts), dtype=np.float64)
    #Bx_4d = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float64)
    #Bx_localized_4d = np.zeros((nlats, nlons, nlats, nlons), dtype=np.float64)

    coslat = set_coslat(lats)
    
    # --- compute covariance of forecast errors between grid points, and localize
    #     them based on input exponent and horizontal length scale.
    
    ktr1 = 0
    for i1 in range(nlons):
    #for i1 in range(nlons//2, nlons//2+1):
        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        #print ('processing i1 = ',i1,' of ',nlons,'. Current time = ',current_time)
        for j1 in range(nlats):
        #for j1 in range(nlats//2, nlats//2+1):
            x1 = difference_3d[:,j1,i1]
            ktr2 = 0
            for i2 in range(i1,nlons):
                for j2 in range(j1,nlats):
            #for i2 in range(nlons):
            #    for j2 in range(nlats):
                    x2 = difference_3d[:,j2,i2]
                    #hdist = (111./2.) * np.sqrt( ( (i1-i2) * (coslat[j1]+coslat[j2]) / 2. )**2 + \
                    #    (j1-j2)**2)   # 111./2 since 111 km per degree, 1/2 degree grid
                    hdist = (111./2.) * np.sqrt((i1-i2)**2 + (j1-j2)**2)   # 111./2 since 111 km per degree, 1/2 degree grid
                    localizn_factor = np.exp(-(hdist/efold)**exponenty)
                    #print ('j1,i1,j2,i2, hdist, lfactor = ', j1,i1,j2,i2, hdist, localizn_factor)
                    Bx[ktr1,ktr2] = np.cov(x1,x2)[0,1]
                    Bx[ktr2,ktr1] = Bx[ktr1,ktr2] 
                    #Bx_4d[j1,i1,j2,i2] = Bx[ktr1,ktr2] 
                    #Bx_4d[j2,i2,j1,i1] = Bx[ktr1,ktr2] 
                    Bx_localized[ktr1,ktr2] = Bx[ktr1,ktr2]*localizn_factor
                    Bx_localized[ktr2,ktr1] = Bx_localized[ktr1,ktr2]
                    #Bx_localized_4d[j1,i1,j2,i2] = Bx_localized[ktr1,ktr2] 
                    #Bx_localized_4d[j2,i2,j1,i1] = Bx_localized[ktr2,ktr1]
                    ktr2 = ktr2 + 1
            ktr1 = ktr1 + 1
        
    #print (nlats,nlons)   
    #print ('difference_3d[:,nlats//2, nlons//2] = ',difference_3d[:,nlats//2, nlons//2])
    #print ('Bx_4d[nlats//2, nlons//2, nlats//2-10:nlats//2+10, nlons//2] =',\
    #    Bx_4d[nlats//2, nlons//2, nlats//2-10:nlats//2+10, nlons//2])
    #print ('Bx_4d[nlats//2, nlons//2, nlats//2, nlons//2-10:nlons//2+10] =',\
    #    Bx_4d[nlats//2, nlons//2, nlats//2, nlons//2-10:nlons//2+10]) 
    #print ('Bx_localized_4d[nlats//2, nlons//2, nlats//2-10:nlats//2+10, nlons//2] =',\
    #    Bx_localized_4d[nlats//2, nlons//2, nlats//2-10:nlats//2+10, nlons//2])
    #print ('Bx_localized_4d[nlats//2, nlons//2, nlats//2, nlons//2-10:nlons//2+10] =',\
    #    Bx_localized_4d[nlats//2, nlons//2, nlats//2, nlons//2-10:nlons//2+10])
    
    return Bx_localized