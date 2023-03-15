def CDF_spline_fitting_forecast_precip_v4(cmonth, clead, cdomain):

    """
    CDF_spline_fitting_forecast_precip_v4.py cmonth clead cdomain
    
    where cmonth = 'Jan', 'Feb' etc
    clead = '024' or similar (3 digit lead in hours) and 
    cdomain typically is "conus" but permits "pr" "ak" 
    if suitable precipitation analyses can be rustled up for these
    domains.  Not coded for HI or GU, but the domain boundaries could 
    be coded for those, too.

    Purpose: fit cubic spline to cumulative hazard function of precipitation 
    and save these spline parameters for fitting of precipitation CDFs.
    For points with few samples of precipitation, fit Gamma dist. instead.
    
    Designed by Tom Hamill, with help from Michael Scheuerer, Dec 2020 -
    July 2021

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
        
            I've only done this so far for CONUS, PR, AK ....
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

        """ return the surrounding two months, as our spline fit will 
            use data not just from the current month but also those
            surrounding months. 
        """

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
        nonzero amounts (precip_ens_nonzero), count them (nz) and return.
    
        """
        precip_ens_nonzero = np.delete(precip_ens, \
            np.where(precip_ens <= 0.0))  # censor at 0.0 mm
        precip_ens_nonzero = precip_ens_nonzero + \
            np.random.uniform(low=-0.001,high=0.001,size=len(precip_ens_nonzero))
        precip_ens_nonzero = np.delete(precip_ens_nonzero, \
            np.where(precip_ens_nonzero <= 0.0))  # censor at 0.0 mm
        nz = len(precip_ens_nonzero)    
        precip_ens_nonzero = np.sort(precip_ens_nonzero)  
        ntotal = len(precip_ens)
        nzero = ntotal - len(precip_ens_nonzero)
        fraction_zero = float(nzero) / float(ntotal)
        return fraction_zero, precip_ens_nonzero, nz

    # =====================================================================
    
    def define_knot_locations(nz, precip_samples_nonzero):
        
        """ define_knot_locations:  choose the number of knots and the indices 
            in the sorted precipitation samples according to the number of
            positive precipitation amounts.   We define fewer knots for grid 
            points with fewer positive precipitation samples, and we want 
            the precipitation values of the chosen knots to emphasize 
            the upper quantiles of the distribution, as that's where
            we care most about an accurate fit.  This will be described in the
            accompanying Hamill et al. 2022/2023 MWR journal article. 
        """

        # ---- these beta parameters be used to create a beta distribution 
        #      that will emphasize the upper quantiles of the distribution
        
        rp = 3.5
        rq = 1.0
        
        # ---- the number of (interior) knots in the cubic spline will
        #      be set to be no greater than 9.0, and for a sample size
        #      of 100 will be 3.
        
        nknots = min([9,nz//30])
        query_these_indices = []
        cdf_at_indices = []
        
        for iknot in range (1,nknots+1):
            rknot = float(iknot)/(nknots+1)
            xloc = stats.beta.ppf(rknot, rp, rq)
            c = stats.beta.cdf(xloc, rp, rq)
            iloc = int(nz*xloc)
            query_these_indices.append(iloc)
            c = (1./(2.*nz)) + float(iloc)/float(nz)
            cdf_at_indices.append(c)
        
        return nknots, query_these_indices, cdf_at_indices
        
    # =====================================================================
        
    # ---- set general parameters
        
    cmonth_list = ['Jan','Feb','Mar','Apr','May','Jun',\
        'Jul','Aug','Sep','Oct','Nov','Dec']
    cmonthnum = ['01','02','03','03','04','05',\
        '06','07','08','09','10','11','12']
    print ('cmonth = ', cmonth)
    print (cmonth_list)
    imonth_index = cmonth_list.index(cmonth)
    jmin, jmax, imin, imax = set_domain_boundaries(cdomain)
    nstride = 1
    cmonth_early, cmonth_late = get_surrounding_months(cmonth)
    pflag = True # for print statements
    master_directory = '/Volumes/NBM/'+cdomain+'_gefsv12/precip/netcdf/'
    master_directory_out = '/Volumes/NBM/'+cdomain+'_gefsv12/CDF_spline/'
    nmembers = 5 # 5 members in GEFSv12 reforecast

    # ---- read in the previously generated netCDF file with precipitation
    #      for this month and lead time + surrounding months.  
    #      All members, dates for this month have been
    #      smushed into one leading index, dimension
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
    #      we know the array dimensions. This is fraction zero.

    fzero = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.float)
   
    # ---- loop over the grid points and estimate the spline coefficients
    #      for each point.   If very dry, instead estimate Gamma distribution
    #      parameters.   If super dry, do nothing.

    now = datetime.now()
    begin_time = now.strftime("%H:%M:%S")
    plotit = False
 
    spline_info_fcst = np.zeros((ny_gefsv12,nx_gefsv12,2,17), dtype=np.float64) 
    usegamma_fcst = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.int32)
    q99 = np.zeros((ny_gefsv12,nx_gefsv12), dtype=np.float32)
    number_knots = np.zeros((ny_gefsv12,nx_gefsv12), dtype=np.int32)
                   
    for jy in range(0,ny_gefsv12,nstride):
        
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('***** begin time, current time, jy, ny_gefsv12, lat = ',\
            begin_time, current_time, jy, ny_gefsv12, lats_1d[jy])        
        
        for ix in range(0,nx_gefsv12,nstride):
            
            # ---- for this grid point, form a vector of the data 
            #      for this month and the surrounding two months.
            
            precip_ens_1d_middle = precip_middle[:,jy,ix]
            precip_ens_1d_early = precip_early[:,jy,ix]
            precip_ens_1d_late = precip_late[:,jy,ix]
            precip_ens_1d = np.concatenate((precip_ens_1d_middle, \
                precip_ens_1d_early, precip_ens_1d_late))
            
            # ---- there is a grib round-off error that can give negative
            #      values slightly smaller than teeny precip.  to make sure 
            #      that we don't have either negative values or lots of the 
            #      same tiny values, so subtract teeny_precip

            tp = np.min( [np.abs(np.min(precip_ens_1d)), 0.0] )
            teeny_precip = tp*np.ones(nsamps)
            precip_ens_1d = precip_ens_1d - teeny_precip[:]
            
            fraction_zero, precip_ens_nonzero, nz = \
                fraczero_possamps(nsamps, precip_ens_1d) # return sorted
            fzero[jy,ix] = fraction_zero 
        
            if pflag == True: print ('---------  jy, ix = ', jy,ix,' --------- ' )
            if pflag == True: print ('   precip_ens_nonzero[-1] = ', precip_ens_nonzero[-1])
            if pflag == True: print ('   nz, nsamps = ', nz, nsamps)
                
            if nz > 100:

                # ---- spline fit with the focus on knots at higher quantiles

                usegamma_fcst[jy,ix] = 0
                q99[jy,ix] = precip_ens_nonzero[(99*nz)//100]
                empirical_cdf = 1.0/(2.0*nz) + np.arange(nz)/nz      
                hazard_function_empirical = -np.log(1.0 - empirical_cdf)
                
                # ---- set the indices in the sorted precipitation sample 
                #      where spline knots are, the CDF at these values,
                #      the number of knots, the precipitation amount at
                #      chosen percentiles.
                    
                nknots, query_these_indices, cdf_at_indices = \
                    define_knot_locations(nz, precip_ens_nonzero)
                empirical_precipvals = precip_ens_nonzero[query_these_indices]
                spltemp_inv = splrep(precip_ens_nonzero, hazard_function_empirical, \
                    xb=0., task=-1, t = empirical_precipvals)   
                lspline = len(spltemp_inv[0])  # final # knots including @ boundaries
        
                # --- set the number of knots and splines coefficients into arrays to
                #     be dumped to netCDF file.
            
                spline_info_fcst[jy,ix,0,0:lspline] = spltemp_inv[0]
                spline_info_fcst[jy,ix,1,0:lspline] = spltemp_inv[1]
                number_knots[jy,ix] = lspline
        
            elif nz < 10:
            
                # --- over the whole training sample, there were not enough 
                #     nonzero precip values recorded.   Flag this as missing data.
            
                usegamma_fcst[jy,ix] = -1  # flag as insufficient data to do any qmapping
                number_knots[jy,ix] = 0
                spline_info_fcst[jy,ix,0,:] = -99.99
                spline_info_fcst[jy,ix,1,:] = -99.99
                q99[jy,ix] = -99.99

            else:
            
                # ---- too few samples to use splines; fit a single-parameter Gamma 
                #      using Thom (1958) estimator described in Wilks textbook,
                #      Statistical Methods in the Atmospheric Sciences.
                
                usegamma_fcst[jy,ix] = 1  # yes, estimate with Gamma distribution
                number_knots[jy,ix] = 0
                empirical_cdf = 1.0/(2.0*nz) + np.arange(nz)/nz
                pmean = np.mean(precip_ens_nonzero)
                lnxbar = np.log(pmean)
                if pflag == True: print ('  pmean, lnxbar = ', pmean, lnxbar)
                meanlnxi = np.mean(np.log(precip_ens_nonzero))
                if pflag == True: print ('  meanlnxi = ', meanlnxi)
                D = lnxbar - meanlnxi
                alpha_hat = (1.0 + np.sqrt(1.0+4.0*D/3.0)) / (4.0*D)
                beta_hat = pmean / alpha_hat
                if pflag == True: print('   alpha_hat, beta_hat = ', alpha_hat, beta_hat)
                spline_info_fcst[jy,ix,0,:] = alpha_hat  # smoosh into the spline array
                spline_info_fcst[jy,ix,1,:] = beta_hat # smoosh into the spline array
                q99[jy,ix] = precip_ens_nonzero[-1] 
            
                # --- evaluate Dn statistic, goodness of fit.
            
                y0 = precip_ens_nonzero / beta_hat
                fitted_CDF = stats.gamma.cdf(y0, alpha_hat)
           
    # ---- save spline / Gamma distribution information to netCDF file
        
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
    
    quantile_99 = ncout.createVariable('quantile_99',\
        'f8',('yf','xf',), zlib=True)
    quantile_99.units = "mm"
    quantile_99.long_name = "99th percentile of precipitation CDF"
    quantile_99.missing_value = np.array(-99.99,dtype=np.float32)
    
    numberknots_out = ncout.createVariable('number_knots',\
        'i4',('yf','xf',), zlib=True)
    numberknots_out.units = "n/a"
    numberknots_out.long_name = "number of knots when using spline"
    numberknots_out.missing_value = np.array(-99,dtype=np.int32)

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
    quantile_99[:] = q99[:,:]
    numberknots_out[:] = number_knots[:,:]

    ncout.close()




