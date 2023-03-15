def Bxmodel(nlats, nlons, ndates, lats, lons, \
    cyear, clead, warm_or_cold, cpath_Bx, \
    efold, exponenty, forecast_3d, analyses_3d, \
    bias_corr_3d, already):
    
    """
    form a model for the (bias-corrected) forecast-error covariance, 
        including covariance localization.   Save the resulting B_x 
        to file
    """

    import numpy as np
    import os, sys
    import _pickle as cPickle
    #from forecast_error_covariance import forecast_error_covariance
    from forecast_error_covariance_f90 import forecast_error_covariance_f90
    
    npts = int(nlats*nlons)
    Bx_localized = np.zeros((npts, npts), dtype=np.float64)
    
    if already == False:
        
        # ---- bias_correct the forecast before comparing to the analyses

        difference_3d_biascorr = analyses_3d - (forecast_3d - bias_corr_3d)

        # ---- produce estimate of the localized covariance of bias-corrected 
        #      forecast errors between grid points and the localized, inverted 
        #      covariance matrix.

        Bx_localized = \
            forecast_error_covariance_f90(difference_3d_biascorr, efold, \
                exponenty, npts, ndates, nlats, nlons)

        # ---- write the Bx_localized to pickle file.

        outfile = cpath_Bx+'Localized_Bx_'+warm_or_cold+\
            'season_year'+cyear+'_lead='+clead+'_efold'+\
            str(efold)+'.cPick'
        print ('writing to ', outfile)
        ouf = open(outfile,'wb')
        cPickle.dump(Bx_localized, ouf)
        ouf.close()
        
    else:
        # ---- read the Bx_localized from pickle file.

        infile = cpath_Bx+'Localized_Bx_'+warm_or_cold+\
            'season_year'+cyear+'_lead='+clead+'_efold'+\
            str(efold)+'.cPick'
        inf = open(infile,'rb')
        Bx_localized = cPickle.load(inf)
        inf.close()
        
    return Bx_localized

    
