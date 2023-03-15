def read_cloudfcst_timeseries_GEFSv12(cpath_cloud, date_list, clead):
    
    """
    read_cloudfcst_timeseries_GEFSv12.py : read in the cloud cover forecast data desired
    """
    
    import numpy as np
    import numpy.ma as ma
    import pygrib
    import os.path as path
    import sys
    import _pickle as cPickle

    ndates = len(date_list)
    for idate, date in enumerate(date_list):
 
        # ---- read the forecast information for bias corr.
         
        cyear = date[0:4]
        cmmdd = date[4:8]    
        infile = cpath_cloud + cyear + '/'+date+'_lead'+\
            clead+'_cldcover_conus_0.5deg_hour'+clead+'.cPick'
        fexist2 = path.exists(infile)
        if fexist2 == True:
            inf = open(infile,'rb')
            cloud = cPickle.load(inf)
            if cmmdd == '0101':
                nlats, nlons = np.shape(cloud)
                cloud_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            inf.close()
            cloud_3d[idate,:,:] = cloud[:,:]
        else:
            print ('Unable to read ', infile)
            cloud_3d[idate,:,:] = ma.masked
 
    print ('    max, min cloud_3d = ', np.max(cloud_3d), np.min(cloud_3d), np.mean(cloud_3d))
    return cloud_3d
