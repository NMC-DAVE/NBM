def compute_closest_member(nmembers_x25, ny, nx, nfour, \
    thresh_low, thresh_mod, thresh_high, ensmean, precip_ens, \
    precip_anal, closest_histogram):
    
    import numpy as np
    
    # --- determine if there is any corrupt ensemble members, and if so, don't 
    #     tally stats for this day.
    
    rm = np.min(precip_ens)
    rma = np.min(precip_anal)
    if rm < 0.0 or rma < 0.0:
        istat = -1
        print ('rm, rma = ', rm, rma)
        print ('identified bad forecast or analysis data, so skip this day.')
    else:
        istat = 0
        for ix in range(nx):
            for jy in range(ny):
            
                # --- determine which member is closest to the analyzed and
                #       how many members have values lower than or equal to analyzed
                
                a = precip_anal[jy,ix]
                iclosest = 1
                rclosest = 9999.
                eclosest = 0.0
                for imem in range(nmembers_x25):
                    e = precip_ens[imem,jy,ix]
                    diff = np.abs(e-a)
                    if diff < rclosest and e > -99:
                        rclosest = diff
                        eclosest = e
                        iclosest = imem

                ibelow = 0
                iequal = 0
                for imem in range(nmembers_x25):

                    e = precip_ens[imem,jy,ix]
                    if imem == iclosest:
                        continue
                    else:
                        if e < eclosest: ibelow = ibelow + 1
                        if e == eclosest: iequal = iequal + 1
                
			        # --- determine the closest_histogram rank 
				
                    if iequal == 0:
                        iclosest = ibelow + 1			
                    else: #with equals, randomly assign rank
                        r = np.random.uniform()*iequal
                        ir = int(r)
                        if ir > iequal: ir = iequal
                        iclosest = ibelow + ir + 1
                
                    em = ensmean[jy,ix]
                    if em < thresh_low:
                        closest_histogram[iclosest,0] = closest_histogram[iclosest,0] + 1
                    elif em >= thresh_low and em < thresh_mod:
                        closest_histogram[iclosest,1] = closest_histogram[iclosest,1] + 1
                    elif em >= thresh_mod and em < thresh_high:
                        closest_histogram[iclosest,2] = closest_histogram[iclosest,2] + 1
                    else:
                        closest_histogram[iclosest,3] = closest_histogram[iclosest,3] + 1

    return closest_histogram, istat
     
