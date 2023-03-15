def CDF_spline_fitting_forecast_precip_v2(cmonth, clead, cdomain):

    """
    CDF_spline_fitting_forecast_precip_v2.py cmonth clead cdomain
    
    where cmonth = 'Jan', 'Feb' etc
    clead = 

    fit cubic spline to cumulative hazard function of precipitation and save
    for spline parameters for fitting of precipitation CDFs.

    Designed by Tom Hamill, NOAA, with help from Michael Scheuerer, Dec 2020

    """

    import os, sys
    from datetime import datetime
    import numpy as np
    import _pickle as cPickle
    from netCDF4 import Dataset
    import _pickle as cPickle
    import scipy.stats as stats
    from scipy.interpolate import LSQUnivariateSpline, splrep, splev
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    from mpl_toolkits.basemap import Basemap, interp
    from mpl_toolkits.axes_grid1 import make_axes_locatable

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

    def get_surrounding_months(cmonth):

        if cmonth == 'Jan':
            cmonth_early = 'Dec'
            cmonth_late = 'Feb'
        elif cmonth == 'Feb':
            cmonth_early = 'Jan'
            cmonth_late = 'Mar'
        elif cmonth == 'Mar':
            cmonth_early = 'Feb'
            cmonth_late = 'Apr'
        elif cmonth == 'Apr':
            cmonth_early = 'Mar'
            cmonth_late = 'May'
        elif cmonth == 'May':
            cmonth_early = 'Apr'
            cmonth_late = 'Jun'
        elif cmonth == 'Jun':
            cmonth_early = 'May'
            cmonth_late = 'Jul'
        elif cmonth == 'Jul':
            cmonth_early = 'Jun'
            cmonth_late = 'Aug'
        elif cmonth == 'Aug':
            cmonth_early = 'Jul'
            cmonth_late = 'Sep'
        elif cmonth == 'Sep':
            cmonth_early = 'Aug'
            cmonth_late = 'Oct'
        elif cmonth == 'Oct':
            cmonth_early = 'Sep'
            cmonth_late = 'Nov'
        elif cmonth == 'Nov':
            cmonth_early = 'Oct'
            cmonth_late = 'Dec'
        elif cmonth == 'Dec':
            cmonth_early = 'Nov'
            cmonth_late = 'Jan'
        else:
            print ('invalid month')
            sys.exit()
    
        return cmonth_early, cmonth_late

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
        number_zeros = 0
    
        # data discretized, so add random component of this magnitude
        #precip_ens = np.where(precip_ens < 4.0, precip_ens* \
        #    np.random.uniform(low=-0.5,high=1.5,size=len(precip_ens)), precip_ens)
    
        precip_ens_nonzero = np.delete(precip_ens, \
            np.where(precip_ens <= 0.0))  # censor at 0.0 mm
        precip_ens_nonzero = precip_ens_nonzero + \
            np.random.uniform(low=-0.1,high=0.1,size=len(precip_ens_nonzero))
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
    
    # ---- set parameters
        
    cmonth_list = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    cmonthnum = ['01','02','03','03','04','05','06','07','08','09','10','11','12']
    print ('cmonth = ', cmonth)
    print (cmonth_list)
    imonth_index = cmonth_list.index(cmonth)
    jmin, jmax, imin, imax = set_domain_boundaries(cdomain)
    nstride = 1
    cmonth_early, cmonth_late = get_surrounding_months(cmonth)
    pflag = False # for print statements
    master_directory = '/Volumes/NBM/'+cdomain+'_gefsv12/precip/netcdf/'
    master_directory_out = '/Volumes/NBM/'+cdomain+'_gefsv12/CDF_spline/'
    nmembers = 5

    # ---- read in the previously generated netCDF file with precipitation
    #      for this month and lead time.  All members, dates for this 
    #      month have been smushed into one leading index, dimension
    #      nsamps, since the date of the forecast within the month and 
    #      the member number is irrelevant for the distribution fitting.
   
    ncfile = master_directory + cmonth + '_apcp_sfc_h' + clead + '.nc'
    print (ncfile)
    nc = Dataset(ncfile)
    precip_middle = nc.variables['apcp_fcst'][:,jmin:jmax,imin:imax]
    nsamps_middle, ny_gefsv12, nx_gefsv12 = np.shape(precip_middle)
    lons_1d = nc.variables['lons_fcst'][imin:imax]
    lats_1d = nc.variables['lats_fcst'][jmin:jmax]
    nc.close()

    ncfile = master_directory + cmonth_early + '_apcp_sfc_h' + clead + '.nc'
    print (ncfile)
    nc = Dataset(ncfile)
    precip_early = nc.variables['apcp_fcst'][:,jmin:jmax,imin:imax]
    nsamps_early, ny_gefsv12, nx_gefsv12 = np.shape(precip_early)
    nc.close()

    ncfile = master_directory + cmonth_late + '_apcp_sfc_h' + clead + '.nc'
    print (ncfile)
    nc = Dataset(ncfile)
    precip_late = nc.variables['apcp_fcst'][:,jmin:jmax,imin:imax]
    nsamps_late, ny_gefsv12, nx_gefsv12 = np.shape(precip_late)
    nc.close()

    nsamps = nsamps_middle + nsamps_early + nsamps_late

    # ---- more initialization of output storage arrays now that 
    #      we know the array dimensions

    fzero = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.float)
    Dnstat = 0.10*np.ones((ny_gefsv12, nx_gefsv12), dtype=np.float)
   
    # ---- loop over the grid points and estimate the Gamma distributions
    #      for each parameter.  First see if a single Gamma distribution
    #      is appropriate; if not, try a mixture of two.   If that still
    #      doesn't fit well, try a mixture of three.   

    now = datetime.now()
    begin_time = now.strftime("%H:%M:%S")
    plotit = False
 
    spline_info_fcst = np.zeros((ny_gefsv12,nx_gefsv12,2,17), dtype=np.float64) 
    usegamma_fcst = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.int32)
    empirical_precipvals_fcst = np.zeros((9,ny_gefsv12,nx_gefsv12), dtype=np.float16)
    q98 = np.zeros((ny_gefsv12,nx_gefsv12), dtype=np.float32)
                   
    for jy in range(0,ny_gefsv12,nstride):
    #for jy in range(ny_gefsv12//2, ny_gefsv12//2+1):
        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('***** begin time, current time, jy, ny_gefsv12, lat = ',\
            begin_time, current_time, jy, ny_gefsv12, lats_1d[jy])        
        
        for ix in range(0,nx_gefsv12,nstride):
        #for ix in range (nx_gefsv12//2, nx_gefsv12//2+1):
            
            # ---- there is a grib round-off error that can give negative
            #      values slightly smaller than teeny precip.  to make sure 
            #      that we don't have either negative values or lots of the 
            #      same tiny values, so subtract teeny_precip
        
            precip_ens_1d_middle = precip_middle[:,jy,ix]
            precip_ens_1d_early = precip_early[:,jy,ix]
            precip_ens_1d_late = precip_late[:,jy,ix]
            precip_ens_1d = np.concatenate((precip_ens_1d_middle, \
                precip_ens_1d_early, precip_ens_1d_late))
        
            tp = np.min( [np.abs(np.min(precip_ens_1d)), 0.0] )
            teeny_precip = tp*np.ones(nsamps)
            precip_ens_1d = precip_ens_1d - teeny_precip[:]

            fraction_zero, precip_ens_nonzero, nz = \
                fraczero_possamps(nsamps, precip_ens_1d) # return sorted
            fzero[jy,ix] = fraction_zero 
        
            if nz > 50:

                # ---- spline fit with the focus on knots at higher quantiles

                usegamma_fcst[jy,ix] = 0
                query_these_indices = [ nz//10, nz//4, nz//2, (3*nz)//5, \
                    (4*nz)//5, (17*nz)//20, (9*nz)//10, (24*nz)//25, (49*nz)//50]
                empirical_precipvals = precip_ens_nonzero[query_these_indices]
                empirical_precipvals_fcst[:,jy,ix] = empirical_precipvals[:]
                #print ('empirical_precipvals = ', empirical_precipvals)
                q98[jy,ix] = empirical_precipvals[-1]
        
                empirical_cdf = 1.0/(2.0*nz) + np.arange(nz)/nz      
                hazard_function_empirical = -np.log(1.0 - empirical_cdf)
                spltemp = splrep(precip_ens_nonzero, hazard_function_empirical, \
                    xb=0., task=-1, t = empirical_precipvals)   
                spline_hazard = splev(precip_ens_nonzero, spltemp)
                spline_cdf = 1.0 - np.exp(-spline_hazard)       
                diff = np.abs(empirical_cdf - spline_cdf)
        
                # ---- save spline information to numpy array
        
                spline_info_fcst[jy,ix,0,:] = spltemp[0]
                spline_info_fcst[jy,ix,1,:] = spltemp[1]
                Dnstat[jy,ix] = np.max(diff)
        
            elif nz < 10:
            
                # --- over the whole training sample, there were not enough 
                #     nonzero precip values recorded.   Flag this as missing data.
            
                usegamma_fcst[jy,ix] = -1  # flag as insufficient data to do any qmapping
                spline_info_fcst[jy,ix,0,:] = -99.99
                spline_info_fcst[jy,ix,1,:] = -99.99
                Dnstat[jy,ix] = 0.0
                q98[jy,ix] = -99.99

            else:
            
                # ---- too few samples to use splines; fit a single-parameter Gamma 
                #      using Thom (1958) estimator described in Wilks textbook,
                #      Statistical Methods in the Atmospheric Sciences.
                
                usegamma_fcst[jy,ix] = 1  # yes, estimate with Gamma distribution
                empirical_cdf = 1.0/(2.0*nz) + np.arange(nz)/nz
                pmean = np.mean(precip_ens_nonzero)
                lnxbar = np.log(pmean)
                if pflag == True: print ('pmean, lnxbar = ', pmean, lnxbar)
                meanlnxi = np.mean(np.log(precip_ens_nonzero))
                if pflag == True: print ('meanlnxi = ', meanlnxi)
                D = lnxbar - meanlnxi
                alpha_hat = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
                beta_hat = pmean / alpha_hat
                if pflag == True: print('alpha_hat, beta_hat = ', alpha_hat, beta_hat)
                spline_info_fcst[jy,ix,0,:] = alpha_hat  # smoosh into the spline array
                spline_info_fcst[jy,ix,1,:] = beta_hat # smoosh into the spline array
                q98[jy,ix] = precip_ens_nonzero[-1] 
            
                # --- evaluate Dn statistic, goodness of fit.
            
                y0 = precip_ens_nonzero / beta_hat
                fitted_CDF = stats.gamma.cdf(y0, alpha_hat)
                diff = np.abs(empirical_cdf - fitted_CDF)
                Dnstat[jy,ix] = np.max(diff)
           
    # --- save Dnstat goodness of fit to cPickled file

    print ('max Dnstat = ', np.max(Dnstat))
    outfile = master_directory_out + cmonthnum[imonth_index]+'_'+cdomain+\
        '_GEFSv12_Dnstat_h' + clead + '.cPick'
    print ('writing to ', outfile)
    ouf = open(outfile, 'wb')
    cPickle.dump(Dnstat, ouf)
    cPickle.dump(lons_1d, ouf)
    cPickle.dump(lats_1d, ouf)
    ouf.close()

    # ---- save spline information to netCDF file

    #outfile = master_directory_out + cmonthnum[imonth_index]+'_'+cdomain+\
    #    '_GEFSv12_spline_info_h' + clead + '.nc' 
        
    outfile = master_directory_out + cmonth_list[imonth_index]+'_'+cdomain+\
        '_GEFSv12_spline_info_h' + clead + '.nc'
    print ('writing to ',outfile)
    ncout = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

    xf = ncout.createDimension('xf',nx_gefsv12)
    xvf = ncout.createVariable('xf','i4',('xf',))
    xvf.long_name = "eastward grid point number on NDFD grid"
    xvf.units = "n/a"

    yf = ncout.createDimension('yf',ny_gefsv12)
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

    x9 = ncout.createDimension('x9',9)
    x9f = ncout.createVariable('x9','i4',('x9',))
    x9f.long_name = "first dimension of nonzero_indices_of_knots"
    x9f.units = "n/a"

    lonsa = ncout.createVariable('lons','f4',('xf',))
    lonsa.long_name = "longitude"
    lonsa.units = "degrees_east"

    latsa = ncout.createVariable('lats','f4',('yf',))
    latsa.long_name = "latitude"
    latsa.units = "degrees_north"

    spline_info_out = ncout.createVariable('spline_info',\
        'f4',('yf','xf','x2','xspd'),
        zlib=True,least_significant_digit=6)
    spline_info_out.units = "n/a"
    spline_info_out.long_name = \
        "Information for computing precipitation from"+\
        "spline (or Gamma CDF for dry points)"
    spline_info_out.missing_value = \
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
    
    quantile_98 = ncout.createVariable('quantile_98',\
        'f8',('yf','xf',), zlib=True)
    quantile_98.units = "mm"
    quantile_98.long_name = "98th percentile of precipitation CDF"
    quantile_98.missing_value = np.array(-99.99,dtype=np.float32)

    # ---- metadata

    ncout.title = "NDFD domain spline inverse coefficients / Gamma parameters "
    ncout.history = "from CDF fitting code by Tom Hamill, PSL"
    ncout.institution =  "psl.noaa.gov"
    ncout.platform = "n/a"
    ncout.references = "n/a"
    
    # ---- copy the outputs to netCDF structures.

    xvf[:] = range(nx_gefsv12)
    yvf[:] = range(ny_gefsv12)
    xspdf[:] = range(17)
    x2f[:] = range(2)
    x9f[:] = range(9)
    lonsa[:] = lons_1d[:]
    latsa[:] = lats_1d[:]
    spline_info_out[:] = spline_info_fcst[:,:,:,:]
    fzero_out[:] = fzero[:,:]
    usegamma_out[:] = usegamma_fcst[:,:]
    quantile_98[:] = q98[:,:]

    ncout.close()




