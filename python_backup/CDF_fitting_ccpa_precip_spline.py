def CDF_fitting_ccpa_precip_spline(cmonth, cend_hour):

    """
    CDF_fitting_ccpa_precip_spline.py cmonth cend_hour

    this python script is designed to spline fit an empirical CDF of 
    precipitation. The script is tailored to CCPA precipitation analyses over 
    one of the National Digital Forecast Database domains for the National 
    Blend of Models

    Designed by Tom Hamill, NOAA, with help from Michael Scheuerer, Jan 2021

    """

    import os, sys
    from datetime import datetime
    import numpy as np
    import _pickle as cPickle
    from netCDF4 import Dataset
    import _pickle as cPickle
    import scipy.stats as stats
    import scipy.stats as stats
    from scipy.interpolate import splrep, splev

    # =====================================================================

    def set_domain_boundaries(cdomain):
    
        """ used grib file of 2.5-km blend output grid to determine bounding 
        lat and lon, and from that, the domain bounding indices for the 
        0.25 GEFSv12 reforecast data that will encompass the domain.    
        """
        if cdomain == 'conus': 
            jmin = 93
            jmax = 246
            imin = 368
            imax = 686
        elif cdomain == 'pr':
            jmin = 243
            jmax = 256
            imin = 649
            imax = 667   
        elif cdomain == 'ak':
            jmin = 19
            jmax = 161
            imin = 201
            imax = 967
        else:
            print ('invalid domain.  Exiting.')     
            sys.exit()    
 
        return jmin, jmax, imin, imax
 
    # =====================================================================

    def find_nearest(vec, value):
    
        """ given a vector vec and a particular value, find the index in vec
        that is nearest to value"""
    
        idx = np.abs(vec-value).argmin()
        return idx

    # =====================================================================

    def fraczero_possamps(nsamps, precip_ens):
    
        """
    
        from the vector input sample precip_ens, define the fraction of
        samples with zero precipitation.   For the positive samples, add
        a small random number to deal with the fact that the data was 
        discretized to ~0.1 mm, so that when later creating CDFs we don't 
        have empirical values with lots of tied amounts.  Also, sort the 
        nonzero amounts and return.
    
        """
    
        precip_ens_nonzero = np.delete(precip_ens, \
            np.where(precip_ens <= 0.0))  # censor at 0.0 mm
        precip_ens_nonzero = precip_ens_nonzero + \
            np.random.uniform(low=-0.01,high=0.01,size=len(precip_ens_nonzero))
        precip_ens_nonzero = np.delete(precip_ens_nonzero, \
            np.where(precip_ens_nonzero <= 0.0))  # censor at 0.0 mm
        nz = len(precip_ens_nonzero)    
        
        #precip_ens_nonzero = precip_ens_nonzero + \
        #    np.random.uniform(low=-0.005,high=0.005,size=nz) 
        precip_ens_nonzero = np.sort(precip_ens_nonzero)  
        ntotal = len(precip_ens)
        nzero = ntotal - len(precip_ens_nonzero)
        fraction_zero = float(nzero) / float(ntotal)
    
        return fraction_zero, precip_ens_nonzero, nz
    
    # =====================================================================
    
    # ---- inputs from command line

    #cmonth = sys.argv[1] # '01', '02' etc.
    #cend_hour = sys.argv[2] # 06, 12, 18, 00 -- end hour of 6-h period
    print ('cmonth, cend_hour = ', cmonth, cend_hour)
    imonth = int(cmonth) - 1
    nstride = 1 # do every point
    cdomain = 'conus'

    # ---- set directories, constants

    pflag = True # for print statements
    #master_directory = '/Volumes/Backup Plus/ccpa/'
    master_directory = '/Volumes/NBM/'+cdomain+'_panal/'
    master_directory_out = '/Volumes/NBM/'+cdomain+'_panal/CDF_spline/'
    ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
    ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
    cmonths_early = ['12','01','02','03','04','05','06','07','08','09','10','11']
    cmonths_late =  ['02','03','04','05','06','07','08','09','10','11','12','01']

    # ---- determine the overall number of daily precipitation 
    #      samples across all years for this month and the surrounding
    #      two months

    iearly = int(cmonths_early[imonth])-1
    ilate = int(cmonths_late[imonth])-1

    if imonth != 1:  # not Feb
        nsamps_mid = ndaysomo[imonth]*18
    else:
        nsamps_mid = 4*ndaysomo_leap[imonth] + 14*ndaysomo[imonth]
    
    if iearly != 1:  # not Feb    
        nsamps_early = ndaysomo[iearly]*20
    else:
        nsamps_early = 4*ndaysomo_leap[iearly] + 14*ndaysomo[iearly]
    if ilate != 1:  # not Feb    
        nsamps_late = ndaysomo[ilate]*20
    else:
        nsamps_late = 4*ndaysomo_leap[ilate] + 14*ndaysomo[ilate]
    nsamps = nsamps_mid + nsamps_early + nsamps_late
    print ('nsamps = ', nsamps)

    # ---- read in the previously generated netCDF file with precipitation
    #      for this month and lead time as well as the surrounding
    #      two months.  All dates for this month have
    #      been smushed into one leading index, dimension nsamps,
    #      since the date of the forecast within the month and 
    #      the member number is irrelevant for the distribution 
    #      fitting.
   

    ktr = 0
    for iyear in range(2002,2020):
    
        # --- loop over the month in question and the surrounding 2 months
    
        for cmo in [cmonth, cmonths_early[imonth], cmonths_late[imonth]]:
            imo = int(cmo)-1
            if iyear%4 == 0:
                ndays = ndaysomo_leap[imo]
            else:
                ndays = ndaysomo[imo]
            cyear = str(iyear)    
            infile = master_directory + cyear + cmo + \
                '_ccpa_on_ndfd_grid_6hourly.nc'
            print (infile)
            nc = Dataset(infile)
            yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
            for iday in range(1,ndays+1):
                if iday < 10:
                    cday = '0'+str(iday)
                else:
                    cday = str(iday)
                iyyyymmddhh = int(str(iyear)+cmo+cday+cend_hour)
                idx = np.where(yyyymmddhh_end == iyyyymmddhh)[0]
                precip_in = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
                if iyear == 2002 and iday == 1 and cmo == cmonth:
                    nlats_ndfd, nlons_ndfd = np.shape(precip_in)
                    precip_tseries = np.zeros((nsamps,nlats_ndfd,nlons_ndfd), \
                        dtype=np.float64)
                    missingv = -99.99*np.ones((nlats_ndfd, nlons_ndfd), \
                        dtype=np.float64)
                    lons = nc.variables['lons'][:,:]
                    lats = nc.variables['lats'][:,:]
                precip_in = np.where(precip_in < 500., precip_in, missingv)
                precip_tseries[ktr,:,:] = precip_in[:,:]
                ktr = ktr+1
            nc.close()


    # ---- loop over the grid points and estimate the spline coefficients 

    now = datetime.now()
    begin_time = now.strftime("%H:%M:%S")

    spline_info = np.zeros((17,2,nlats_ndfd,nlons_ndfd), dtype=np.float64) 
    spline_info_inv = np.zeros((17,2,nlats_ndfd,nlons_ndfd), dtype=np.float64)

    usegamma = np.zeros((nlats_ndfd,nlons_ndfd), dtype=np.int32) 
        # flag for whether to use Gamma fit (1)
    indices_to_query = np.zeros((nlats_ndfd,nlons_ndfd,9), dtype=np.float16)
    Dnstat = 0.10*np.ones((nlats_ndfd, nlons_ndfd), dtype=np.float)
    fzero = np.zeros((nlats_ndfd, nlons_ndfd), dtype=np.float)
    cdf_at_indices = np.asarray([ 0.1, 0.25, 0.33333, 0.5, 0.65, 0.8, 0.85, 0.9, 0.95])
        # these are interior knots where to calculate the cum hazard fn.

    print ('******** COMPUTING SPLINE COEFFICIENTS (wetter) or GAMMA PARAMETERS (dry) *********')
    #for jy in range(0,nlats_ndfd,nstride):
    for jy in range(498,499):

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        if jy%10 == 0: print ('   ***** begin time, current time, jy, nlats_ndfd, lat = ',\
            begin_time, current_time, jy, nlats_ndfd, lats[jy,0])
            
        #for ix in range(0,nlons_ndfd, nstride):
        for ix in range (521,522):
        
            # ---- there is a grib round-off error that can give negative
            #      values slightly smaller than teeny precip.  to make sure 
            #      that we don't have either negative values or lots of the 
            #      same tiny values, subtractoff teeny_precip
        
            if pflag == True: print ('******* ',jy,ix)
            precip_ens_1d = precip_tseries[:,jy,ix]
            tp = np.min( [np.abs(np.min(precip_ens_1d)), 0.0] )
            teeny_precip = tp*np.ones(nsamps) # make sure that samples >= 0.0
            precip_ens_1d = precip_ens_1d - teeny_precip[:]
            if pflag == True: print ('min, max precip_ens_1d = ',\
                np.min(precip_ens_1d), np.max(precip_ens_1d))
            # take this grid point's sample, calc fraction_zero, and return
            # the number of nonzero (nz) samples and their sorted values
            fraction_zero, precip_ens_nonzero, nz = \
                fraczero_possamps(nsamps, precip_ens_1d) # return sorted
            if pflag == True and nz > 0 : print ('   precip_ens_nonzero[0:-1] = ',\
                precip_ens_nonzero[0:-1])
            if pflag == True: print ('   nz = ',nz)
            if nz > 40:
            
                # ---- spline fit the CDF to the precipitation values via 
                #      Michael Scheuerer's hazard function (see Fig 3 in
                #      https://doi.org/10.1175/MWR-D-20-0096.1. )
            
                fzero[jy,ix] = fraction_zero 
                empirical_cdf = 1.0/(2.0*nz) + np.arange(nz)/nz
                # ---- set where spline knots are.
                query_these_indices = [ nz//10, nz//4, nz//3, nz//2, (3*nz)//5, \
                    (4*nz)//5, (17*nz)//20, (9*nz)//10, (19*nz)//20]
                indices_to_query[jy,ix,:] = query_these_indices[:]
                empirical_precipvals = precip_ens_nonzero[query_these_indices]
                hazard_function_empirical = -np.log(1.0-empirical_cdf)  
                # ---- compute spline coefficients/  
                spltemp = splrep(precip_ens_nonzero, hazard_function_empirical, \
                    xb=0., task=-1, t = empirical_precipvals)  
                # ---- compute the fitted CDF 
                spline_hazard = splev(precip_ens_nonzero, spltemp)
                spline_cdf = 1.0 - np.exp(-spline_hazard)
            
                # ---- save spline information to numpy array, 
        
                spline_info[:,0,jy,ix] = spltemp[0]
                spline_info[:,1,jy,ix] = spltemp[1]
            
                # ---- in the subsequent quantile mapping, we will want the
                #      analyzed precipitation amount given the quantile.
                #      Accordingly, let's also reverse the data in the spline
                #      and get the spline fits of precipitation amount 
                #      (y) to the cdf (x).
            
                hazard_function_at_indices = -np.log(1.0-cdf_at_indices)
                spltemp_inv = splrep(hazard_function_empirical, precip_ens_nonzero, \
                    xb=0., task=-1, t = cdf_at_indices)
                spline_info_inv[:,0,jy,ix] = spltemp_inv[0]
                spline_info_inv[:,1,jy,ix] = spltemp_inv[1]

                # --- evaluate Dn statistic, goodness of fit.
            
                diff = np.abs(empirical_cdf - spline_cdf)
                Dnstat[jy,ix] = np.max(diff) 
                usegamma[jy,ix] = 0
            
            elif nz < 10:
            
                # --- over the whole training sample, there were not enough 
                #     nonzero precip values recorded.   Flag this as missing data.
            
                spline_info_inv[:,0,jy,ix] = -99.99
                spline_info_inv[:,1,jy,ix] = -99.99
                Dnstat[jy,ix] = -99.99
                fzero[jy,ix] = -99.99
                indices_to_query[jy,ix,:] = -99 
                usegamma[jy,ix] = 0
            
            else:
            
                # ---- too few samples to use splines; fit a single-parameter Gamma 
                #      using Thom (1958) estimator described in Wilks textbook,
                #      Statistical Methods in the Atmospheric Sciences.
            
                empirical_cdf = 1.0/(2.0*nz) + np.arange(nz)/nz
                fzero[jy,ix] = fraction_zero 
                pmean = np.mean(precip_ens_nonzero)
                lnxbar = np.log(pmean)
                if pflag == True: print ('pmean, lnxbar = ', pmean, lnxbar)
                meanlnxi = np.mean(np.log(precip_ens_nonzero))
                if pflag == True: print ('meanlnxi = ', meanlnxi)
                D = lnxbar - meanlnxi
                alpha_hat = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
                beta_hat = pmean / alpha_hat
                if pflag == True: print('alpha_hat, beta_hat = ', alpha_hat, beta_hat)
                indices_to_query[jy,ix,:] = -1 # flag for using Gamma

                spline_info_inv[:,0,jy,ix] = alpha_hat  # smoosh into the spline array
                spline_info_inv[:,1,jy,ix] = beta_hat # smoosh into the spline array
                usegamma[jy,ix] = 1
            
                # --- evaluate Dn statistic, goodness of fit.
            
                y0 = precip_ens_nonzero / beta_hat
                fitted_CDF = stats.gamma.cdf(y0, alpha_hat)
                diff = np.abs(empirical_cdf - fitted_CDF)
                Dnstat[jy,ix] = np.max(diff) 
            
           
           
    sys.exit()
            
    # --- save to cPickle file

    #outfile = master_directory_out + cmonth+'_'+cdomain+\
    #    '_CCPA_spline_info_h' + cend_hour + '.cPick'
    #print ('writing to ', outfile)
    #ouf = open(outfile, 'wb')
    #cPickle.dump(spline_info, ouf)
    #cPickle.dump(spline_info_inv, ouf)
    #cPickle.dump(fzero, ouf)
    #cPickle.dump(usegamma, ouf)
    #cPickle.dump(indices_to_query, ouf)
    #ouf.close()

    outfile = master_directory_out + cmonth+'_'+cdomain+\
        '_CCPA_Dnstat_h' + cend_hour + '.cPick'
    print ('writing to ', outfile)
    ouf = open(outfile, 'wb')
    cPickle.dump(Dnstat, ouf)
    cPickle.dump(lons, ouf)
    cPickle.dump(lats, ouf)
    ouf.close()

    # ---- save to netCDF file

    outfile = master_directory_out + cmonth+'_'+cdomain+\
        '_CCPA_spline_info_h' + cend_hour + 'UTC.nc' 
    print ('writing to ',outfile)
    ncout = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

    xf = ncout.createDimension('xf',nlons_ndfd)
    xvf = ncout.createVariable('xf','i4',('xf',))
    xvf.long_name = "eastward grid point number on NDFD grid"
    xvf.units = "n/a"

    yf = ncout.createDimension('yf',nlats_ndfd)
    yvf = ncout.createVariable('yf','i4',('yf',))
    yvf.long_name = "northward grid point number on NDFD grid"
    yvf.units = "n/a"

    xspd = ncout.createDimension('xspd',17)
    xspdf = ncout.createVariable('xspd','i4',('xspd',))
    xspdf.long_name = "index for spline dimension"
    xspdf.units = "n/a"

    x2 = ncout.createDimension('x2',2)
    x2f = ncout.createVariable('x2','i4',('x2',))
    x2f.long_name = "2nd index for spline dimension"
    x2f.units = "n/a"

    lonsa = ncout.createVariable('lons','f4',('yf','xf',))
    lonsa.long_name = "longitude"
    lonsa.units = "degrees_east"

    latsa = ncout.createVariable('lats','f4',('yf','xf',))
    latsa.long_name = "latitude"
    latsa.units = "degrees_north"

    spline_info_inv_out = ncout.createVariable('spline_info_inv',\
        'f4',('xspd','x2','yf','xf',),
        zlib=True,least_significant_digit=6)
    spline_info_inv_out.units = "n/a"
    spline_info_inv_out.long_name = \
        "Information for computing precipitation from"+\
        "spline inverse (or Gamma CDF for dry points)"
    spline_info_inv_out.missing_value = \
        np.array(-9999.99,dtype=np.float32)

    fzero_out = ncout.createVariable('fzero','f4',('yf','xf',),
        zlib=True,least_significant_digit=6)
    fzero_out.units = "n/a"
    fzero_out.long_name = "fraction_zero"
    fzero_out.missing_value = np.array(-9999.99,dtype=np.float32)

    usegamma_out = ncout.createVariable('usegamma',\
        'i4',('yf','xf',), zlib=True)
    usegamma_out.units = "n/a"
    usegamma_out.long_name = "1 if fit CDF via Gamma distribution, 0 if not"
    usegamma_out.missing_value = np.array(-99,dtype=np.int32)

    # ---- metadata

    ncout.title = "NDFD domain spline inverse coefficients / Gamma parameters "
    ncout.history = "from CDF fitting code by Tom Hamill, PSL"
    ncout.institution =  "psl.noaa.gov"
    ncout.platform = "n/a"
    ncout.references = "n/a"
    
    # ---- copy the outputs to netCDF structures.

    xvf[:] = range(nlons_ndfd)
    yvf[:] = range(nlats_ndfd)
    xspdf[:] = range(17)
    x2f[:] = range(2)
    lonsa[:] = lons[:,:]
    latsa[:] = lats[:,:]
    spline_info_inv_out[:] = spline_info_inv[:,:,:,:]
    fzero_out[:] = fzero[:,:]
    usegamma_out[:] = usegamma[:,:]

    ncout.close()


    istat = 0
    return istat