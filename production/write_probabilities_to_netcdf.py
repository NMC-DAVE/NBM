def write_probabilities_to_netcdf (outfile, ncname, idate, cyyyymmddhh, \
    cyyyymmddhh_fcst,  prob_raw_out, prob_qmapped_out,  prob_qmapped_weighted_out, \
    prob_qmapped_weighted_dressed_out, ensmean_raw_out, ensmean_qmapped_out):    
    
    from netCDF4 import Dataset
    
    """ write the raw and quantile-mapped ensemble probabilities to 
    a netCDF file """
    
    # ---- write probabilities and close

    print ('cyyyymmddhh, cyyyymmddhh_fcst = ', cyyyymmddhh, cyyyymmddhh_fcst)
    print ('   writing probs to netCDF file ',outfile)
    if ncname == 'nc_fullfield':
        nc_fullfield = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        nc_fullfield['yyyymmddhh_init'][idate] = int(cyyyymmddhh)
        nc_fullfield['yyyymmddhh_fcst'][idate] = int(cyyyymmddhh_fcst)
        nc_fullfield['probability_raw'][idate] = prob_raw_out[:,:,:]
        nc_fullfield['probability_qmapped'][idate] = prob_qmapped_out[:,:,:]
        nc_fullfield['probability_qmapped_weighted'][idate] = \
            prob_qmapped_weighted_out[:,:,:]
        nc_fullfield['probability_qmapped_weighted_dressed'][idate] = \
            prob_qmapped_weighted_dressed_out[:,:,:]
        nc_fullfield['ensmean_raw'][idate] = ensmean_raw_out[:,:]
        nc_fullfield['ensmean_qmapped'][idate] = ensmean_qmapped_out[:,:]
        nc_fullfield.close()
    elif ncname == 'nc_thinned':
        nc_thinned = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        nc_thinned['yyyymmddhh_init'][idate] = int(cyyyymmddhh)
        iyyyymmddhh = int(cyyyymmddhh_fcst)
        nc_thinned['yyyymmddhh_fcst'][idate] = iyyyymmddhh
        nc_thinned['probability_raw'][idate] = prob_raw_out[:,:,:]
        nc_thinned['probability_qmapped'][idate] = prob_qmapped_out[:,:,:]
        nc_thinned['probability_qmapped_weighted'][idate] = \
            prob_qmapped_weighted_out[:,:,:]
        nc_thinned['probability_qmapped_weighted_dressed'][idate] = \
            prob_qmapped_weighted_dressed_out[:,:,:]
        nc_thinned['ensmean_raw'][idate] = ensmean_raw_out[:,:]
        nc_thinned['ensmean_qmapped'][idate] = ensmean_qmapped_out[:,:]
        nc_thinned.close()
    else:
        print ('invalid ncname in initalize_netCDF_probfiles = ', ncname)
        sys.exit()
        
    istat = 0
    return istat