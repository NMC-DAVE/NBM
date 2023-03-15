def implement_lowess(forecast, analysis, lsmask, nsamps, \
    nlats, nlons, fcstvals, predicted_anal_vals):
    
    import numpy as np 
    import statsmodels.api as sm
    import sys
    from scipy import interpolate

    lowess = sm.nonparametric.lowess
    
    #print ('forecast = ',forecast[:,1,0])
    #print ('analysis = ',analysis[:,1,0])
    #print ('nlats, nlons = ', nlats, nlons)
    for jy0 in range(nlats):
        for ix0 in range(nlons):
    #for jy in range(nlats//2, nlats//2+1):
    #    for ix in range(nlons//2, nlons//2 + 1):
            #print ('jy0,ix0, lsmask[jy,ix] = ', jy0, ix0, lsmask[jy0,ix0])
            if lsmask[jy0,ix0] == 1:
                a = analysis[:,jy0,ix0]
                f = forecast[:,jy0,ix0]
                fdiff = f[-1] - f[0] 
                #fmin = f[0] - fdiff/10.
                #fmax = f[-1] + fdiff/10.
                fmin = f[0]+0.1
                fmax = f[-1]-0.1
                fincrement = (fmax - fmin) / 25.
			
                # --- predict at these forecast temperature values
			
                f1temps = np.arange(fmin, fmax+0.0001, fincrement)
                #print (jy0, ix0, f1temps)
			
			    # --- implement LOWESS 
			    #     https://www.statsmodels.org/stable/generated/ \
			    #     statsmodels.nonparametric.smoothers_lowess.lowess.html
                #     delta=0.0, 
			
                atemps_predicted = lowess(\
                    a, f,frac=0.2, it=0, delta=0.0, xvals=None,   \
                    is_sorted = False, missing='none',\
                    return_sorted = False)
                #print ('atemps_predicted = ', atemps_predicted)
                #print ('f')
                #f0 = interpolate.interp1d(f1temps,atemps_predicted)
                #print ('np.shape f = ', np.shape(f))
                #print ('np.shape atemps_predicted = ', np.shape(atemps_predicted))
                f0 = interpolate.interp1d(f,atemps_predicted)   
                ainterp = f0(f1temps)
                #if jy0 == 1 and ix0 == 0:
                #    print (lsmask[jy0,ix0])
                #    print ('a = ', a)
                #    print ('f = ', f)
                #    print ('f1temps = ',f1temps)
                #    print ('ainterp = ', ainterp)

                fcstvals[:,jy0,ix0] = f1temps[:]
                predicted_anal_vals[:,jy0,ix0] = ainterp[:]
			
    return fcstvals, predicted_anal_vals 
			