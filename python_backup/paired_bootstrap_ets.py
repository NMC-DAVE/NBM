def paired_bootstrap_ets(contab_x,contab_y):
    
    """ given vectors x and y, return a paired bootstrap 
    estimate of the 5th and 95th percentiles
    of the resampled distribution
    
    inputs:
    contab_x: 1st forecast contingency table (raw), dimension [ncasedays,2,2]
    contab_y: 2nd forecast contingency table (QM), dimension [ncasedays,2,2]
    
    """
    
    import numpy as np

    # --------------------
    
    def calc_ets(contab): 
        
        """ calculate the equitable threat score per Wilks """
        aref = (contab[1,1]+ contab[0,1])*(contab[0,1]+contab[1,1])/sum(contab)
        ets = (contab[1,1]-aref)/(contab[1,1]+contab[1,0]+contab[0,1] -aref)
        return ets
        
    def calc_bia(contab): 
        """ calculate the bias per Wilks text"""
        bia = (contab[1,0] + contab[1,1]) / (contab[0,1]+contab[1,1])
        return bia        
        
    # --------------------
        
    nresa = 100   # 10000 resamplings
    nelts, ntwo, ntwoa = np.size(contab_x)   # number cases, 2,2

    ETS1 = np.zeros((nresa), dtype=np.float64) # 1st resampled ETS consistent w. null hypothesis
    ETS2 = np.zeros((nresa), dtype=np.float64) # 2nd resampled ETS consistent w. null hypothesis
    BIA1 = np.zeros((nresa), dtype=np.float64) # 1st resampled BIA consistent w. null hypothesis
    BIA2 = np.zeros((nresa), dtype=np.float64) # 2nd resampled BIA consistent w. null hypothesis
    ETSdiffs = np.zeros((nresa),dtype=np.float) # vector of differences consistent w. null h.
    BIAdiffs = np.zeros((nresa),dtype=np.float) # vector of differences consistent w. null h.
    ones = np.ones((nelts),dtype=np.float) # vector of ones
    zeros = np.zeros((nelts),dtype=np.float) # vector of zeros
    for i in range(nresa):
        
        contab_sum1 = np.zeros((ntwo,ntwo),dtype=np.float)
        contab_sum2 = np.zeros((ntwo,ntwo),dtype=np.float)
        
        x0 = np.random.rand(nelts) # random Uniform [0,1] numbers, dimension nelts
        iusex = np.where(x0 < 0.5, ones, zeros)  # turn this into a list of 1 or 0's based on random #'s
        iusey = 1 - iusex # turn this into an OPPOSITE list of 1 or 0's based on random #'s
        contab_sum1[:,:] = np.mean(contab_x[i,:,:]*iusey + contab_y[i,:,:]*iusex)
        contab_sum2[:,:] = np.mean(contab_x[i,:,:]*iusex + contab_y[i,:,:]*iusey)
        
        # --- new code to calculate ETS and BIA for contingency tables
        
        ETS1[i] = calc_ets(contab_sum1) 
        ETS2[i] = calc_ets(contab_sum2)
        BIA1[i] = calc_bia(contab_sum1)
        BIA2[i] = calc_bia(contab_sum2)
        print ('ETS1,2, BIA1,2 = ', ETS1[i], ETS2[i], BIA1[i], BIA2[i])
        ETSdiffs[i] = ETS1[i]-ETS2[i] # resampled value of ETS diffs. consistent w. null hypothesis
        BIAdiffs[i] = BIA1[i]-BIA2[i] # resampled value of BIA diffs. consistent w. null hypothesis

    # --- sort ETS differences low to high, report back the 5th, 95th percentiles
    
    dsort = np.sort(ETSdiffs)
    d05_ets = dsort[np.int(.05*np.float(nresa))]
    d95_ets = dsort[np.int(.95*np.float(nresa))]
    
    # --- same for bia
    
    dsort = np.sort(BIAdiffs)
    d05_bia = dsort[np.int(.05*np.float(nresa))]
    d95_bia = dsort[np.int(.95*np.float(nresa))]

    return d05_ets, d95_ets, d05_bia, d95_bia