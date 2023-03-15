@jit(nopython=True)
def weighted_probability_dressed(precip_mean, precip_ens, \
    thresholds, thresholds_hist, closest_hist, \
	ny, nx, nmembers, nthresholds, nchcats,\
    b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, b1_std_lowrank,\
	b0_mean_midrank, b1_mean_midrank, b0_std_midrank, b1_std_midrank,\
	b0_mean_highrank, b1_mean_highrank, b0_std_highrank, b1_std_highrank):
	
    
    import scipy.stats as stats
    
	# ---- determine index of point with max ensemble mean precip	
	
    maxindex = precip_mean.argmax()
    jmax, imax = np.unravel_index(maxindex,precip_mean.shape)
	prob = np.zeros((nthresholds, ny, nx), dtype=np.float64)
	prob[:,:,:] = 0.0
	t = thresholds_hist
    
	# ---- process each grid point
	
	for jy in range(ny):
		
		if jy%100 == 0:
            print ('      jy = ',jy,' of ',ny)
		for ix in range(nx)

			# ---- determine the category for the closest-member
			#      histogram CDF to apply for this point based
			#      on quantile-mapped ensemble-mean precipitation
			
			
			p = precip_mean[jy,ix]
			if p == 0.0:
				prob[:,jy,ix] = 0.0 # don't bother sorting and weighting
			else:
				if p <= t[0]:
					icat = 0 # near zero precip
				elif p > t[0] and p < t[1]:
					icat = 1 # light precip
				elif p > t[1] and p < t[2]:
					icat = 2 # light-moderate precip
				elif p > t[2] and p < t[3]:
					icat = 3 # moderate precip
				elif p > t[3] and p < t[4]:
					icat = 4 # moderate-heavy precip
				elif p > t[4] and p < t[5]:
					icat = 5 # heavy precip
				else:
					icat = 6 # extra heavy precip

			
				# ---- Compute final probability as a weighted sum of
				#      kernel densities exceeding the threshold
			
				ens1d[:] = precip_ens[:,jy,ix]
				ens1d = np.sort(ens1d)
                
                for imem in range(nmembers):
					if ens1d[imem] > 0.0:
                        
						if imem == 1:
							pmean = b0_mean_lowrank + ens1d[imem]*b1_mean_lowrank
							pstd = b0_std_lowrank + ens1d[imem]*b1_std_lowrank
						elif imem == nmembers:
							pmean = b0_mean_highrank + ens1d[imem]*b1_mean_highrank
							pstd = b0_std_highrank + ens1d[imem]*b1_std_highrank
						else:	
							pmean = b0_mean_midrank + ens1d[imem]*b1_mean_midrank
							pstd = b0_std_midrank + ens1d[imem]*b1_std_midrank
                            
                        for ithresh in range(nthresholds):
							zscore = (pmean - thresholds[ithresh]) / pstd
							idx = 0
							if zscore < -5.0:
								cum = 0.0
							elif zscore >= -5.0 and zscore <= 5.0:
								cum = stats.norm.cdf(zscore, loc=0, scale=1)
							else
								cum = 1.0
							
							prob[ithresh, jy,ix] = prob[ithresh, jy,ix] + \
								cum * closest_hist[imem,icat]

	print ('  prob[:,jmax,imax] = ', prob[:,jmax,imax])
	print ('  POP prob[jmax,1:nx:10] = ', prob[0,jmax,0:nx:10] )
	
    return prob