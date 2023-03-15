def quantile_mapping(nlats, nlons, forecast_3d, analyses_3d, \
    beta_3d, lsmask, date_list_forecast, cpath_gefsv12, \
    cpath_era5, clead):
    
    """ quantile mapping.py : From previously determined 
    (estimate_Gaussmix_parameters.py) routine that estimates monthly
    parameters of a 3-mode gaussian mixture model, estimate the
    forecast and analyzed CDFs and apply quantile mapping.  
    Save the corrections in beta_3d """
    
    import numpy as np
    from datetime import datetime
    import sys
    import _pickle as cPickle
    import scipy
    from scipy import stats
        
    # --------------------------------------------------------------  
    
    def populate_CDF_table(avalues, ameans, astds, aweights):
        CDF_table = np.zeros((1001,nlats,nlons), dtype=np.float64)
        Zval_analysis1 = np.zeros((nlats,nlons), dtype=np.float64)
        Zval_analysis2 = np.zeros((nlats,nlons), dtype=np.float64)
        Zval_analysis3 = np.zeros((nlats,nlons), dtype=np.float64)
        for i in range(1001):
            Zval_analysis1[:,:] = (avalues[i] - ameans[0,:,:]) / astds[0,:,:]
            Zval_analysis2[:,:] = (avalues[i] - ameans[1,:,:]) / astds[1,:,:]
            Zval_analysis3[:,:] = (avalues[i] - ameans[2,:,:]) / astds[2,:,:]
            CDF_table[i,:,:] = aweights[0,:,:]*stats.norm.cdf(Zval_analysis1) + \
                aweights[1,:,:]*stats.norm.cdf(Zval_analysis2) + \
                aweights[2,:,:]*stats.norm.cdf(Zval_analysis3)  
        return CDF_table
    
    # --------------------------------------------------------------  
    
    def quantile_map(CDF_forecast, CDF_table, avalues, nlats, nlons):
        
        def find_nearest(vec, value):
            idx = np.abs(vec-value).argmin()
            return idx
        
        for jlat in range(nlats):
            for ilon in range(nlons):
                CDF1d = CDF_table[:,jlat,ilon]
                idx = find_nearest(CDF_forecast[jlat,ilon], CDF1d)
                qmapped_forecast[jlat,ilon] = avalues[idx]
        return qmapped_forecast
    
    # --------------------------------------------------------------   
    
    # ---- load the information for generating CDFs, the climatological distributions 
    
    infile = cpath_gefsv12+'GEFSv12_forecast_Gaussmix2_parameters_f'+clead+'.cPick'
    inf = open(infile, 'rb')
    weights_forecast = cPickle.load(inf)
    means_forecast = cPickle.load(inf)
    stddevs_forecast = cPickle.load(inf)
    inf.close()

    infile = cpath_era5+'ERA5_analyzed_Gaussmix2_parameters_f'+clead+'.cPick'
    inf = open(infile, 'rb')
    weights_analysis = cPickle.load(inf)
    means_analysis = cPickle.load(inf)
    stddevs_analysis = cPickle.load(inf)
    inf.close()
    
    # ---- sequentially loop through dates during the sample,
    #      updating the previous day's bias correction
    #      to the new days fcst vs. obs discrepancy.

    ndates = int(len(date_list_forecast))
    qmapped_forecast = np.zeros((nlats, nlons), dtype=np.float64)
    beta_3d[:,:,:] = 0.0
    avalues = -50.0 + np.arange(1001)/10.  # every 0.1 degrees from -50 to 50.
   
    ifirst = True
    for idate, date in enumerate(date_list_forecast):

        # ---- if the first day of the month, load the relevant month's 
        #      worth of info needed to build CDFs.
        
        cdd = date[6:8]
        if cdd == '01' or ifirst == True:
            ifirst = False
            cmm = date[4:6]
            imonth = int(cmm)-1
            aweights = weights_analysis[imonth,:,:,:]
            ameans = means_analysis[imonth,:,:,:]
            astds = stddevs_analysis[imonth,:,:,:] 
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")   
            CDF_table = populate_CDF_table(avalues, \
                ameans, astds, aweights)
            fweights = weights_forecast[imonth,:,:,:]
            fmeans = means_forecast[imonth,:,:,:]
            fstds = stddevs_forecast[imonth,:,:,:]
        
        # --- determine the cumulative probability for today's forecast relative
        #     to the forecast climatological distribution.
            
        Zval_forecast1 = (forecast_3d[idate,:,:] - fmeans[0,:,:]) / fstds[0,:,:]
        Zval_forecast2 = (forecast_3d[idate,:,:] - fmeans[1,:,:]) / fstds[1,:,:]
        Zval_forecast3 = (forecast_3d[idate,:,:] - fmeans[2,:,:]) / fstds[2,:,:]
        CDF_forecast1 = scipy.stats.norm.cdf(Zval_forecast1, loc=0., scale=1.)
        CDF_forecast2 = scipy.stats.norm.cdf(Zval_forecast2, loc=0., scale=1.)
        CDF_forecast3 = scipy.stats.norm.cdf(Zval_forecast3, loc=0., scale=1.)
        CDF_forecast = fweights[0]*CDF_forecast1 + \
            fweights[1]*CDF_forecast2 + fweights[2]*CDF_forecast3

        # ---- perform quantile mapping 
        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")   
        qmapped_forecast = quantile_map(CDF_forecast, CDF_table, avalues, \
            nlats, nlons)
            
        # ---- calculate the "observation" increment (term in parentheses
        #      in eq. 37 in Dee paper)

        beta_3d[idate,:,:] = lsmask[:,:]* \
            (forecast_3d[idate,:,:] - qmapped_forecast[:,:])
        

    return beta_3d