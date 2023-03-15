def read_reanalysis_timeseries(cpath_era5, date_list):
    
    """
    read_reanalysis_timeseries.py : read in the ERA-5 reanalysis data desired
    return 3D array with data as well as associated lat/lon
    """

    import os.path as path
    import _pickle as cPickle
    import numpy as np
    import numpy.ma as ma
    import sys

    ndates = len(date_list)
    for idate, date in enumerate(date_list):
 
        #print (idate, date)
        # ---- read the ECMWF ERA5 reanalysis at valid at the forecast date.
 
        cyyyy = date[0:4]
        infile = cpath_era5 + cyyyy + '/t2m_era5_halfdegree_'+date+'.cPick'
        fexist = path.exists(infile)
        #print (infile, fexist)
        if fexist == True:
            inf = open(infile, 'rb')
            analysis = cPickle.load(inf) - 273.16
            analysis = np.flipud(analysis)
            if idate == 0:
                lats = cPickle.load(inf)
                lons = cPickle.load(inf)
                lats = np.flipud(lats)
                nlats, nlons = np.shape(lats)
                #print (nlats, nlons)
                analyses_3d = ma.zeros((ndates,nlats,nlons), dtype=np.float32) 
                #print ('analysis lats[:,0] = ', lats[:,0])
            analyses_3d[idate,:,:] = analysis[:,:]
        else:
            print ('Unable to read ', infile)
            analyses_3d[idate,:,:] = ma.masked
        #print (infile, fexist, idate, date, analysis[55,115])
 
    return analyses_3d, lats, lons