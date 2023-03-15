    
def read_closest_member_histogram(master_directory_histogram, ctype, \
    cmonth, clead):
    
    from netCDF4 import Dataset
    import numpy as np
    import sys
    
    """ read the closest-member histogram data from netCDF file"""    

    #infile = master_directory_histogram +'closest_member_histogram_'+\
    #    'thinned2_month='+cmonth+'_lead='+clead+'h.nc'
    infile = master_directory_histogram +'closest_member_histogram_'+\
        'thinned_month='+cmonth+'_lead='+clead+'h.nc'
    print ('reading from ',infile)
    
    nc = Dataset(infile,'r')
    
    histogram_thresholds = nc.variables['thresholds'][:]
    nthresh = len(histogram_thresholds)
    closest_member_histogram = nc.variables['closest_hist'][:,:]
    nmemx25, ncats = np.shape(closest_member_histogram)
    
    b0_mean_lowrank = nc.variables['b0_mean_lowrank'][0]
    b1_mean_lowrank = nc.variables['b1_mean_lowrank'][0]
    b0_std_lowrank = nc.variables['b0_std_lowrank'][0]
    b1_std_lowrank = nc.variables['b1_std_lowrank'][0] 
    
    b0_mean_midrank = nc.variables['b0_mean_midrank'][0]
    b1_mean_midrank = nc.variables['b1_mean_midrank'][0]
    b0_std_midrank = nc.variables['b0_std_midrank'][0]
    b1_std_midrank = nc.variables['b1_std_midrank'][0]
    
    b0_mean_highrank = nc.variables['b0_mean_highrank'][0]
    b1_mean_highrank = nc.variables['b1_mean_highrank'][0] 
    b0_std_highrank = nc.variables['b0_std_highrank'][0]  
    b1_std_highrank = nc.variables['b1_std_highrank'][0]
    
    b0_mean_lowrank = np.float32(b0_mean_lowrank)
    b1_mean_lowrank = np.float32(b1_mean_lowrank)
    b0_std_lowrank = np.float32(b0_std_lowrank)
    b1_std_lowrank = np.float32(b1_std_lowrank)
    
    b0_mean_midrank = np.float32(b0_mean_midrank)
    b1_mean_midrank = np.float32(b1_mean_midrank)
    b0_std_midrank = np.float32(b0_std_midrank)
    b1_std_midrank = np.float32(b1_std_midrank)
    
    b0_mean_highrank = np.float32(b0_mean_highrank)
    b1_mean_highrank = np.float32(b1_mean_highrank)
    b0_std_highrank = np.float32(b0_std_highrank)
    b1_std_highrank = np.float32(b1_std_highrank)
    
    nc.close()
    
    #print ('histogram_thresholds = ', histogram_thresholds)
    #print ('closest_member_histogram[0::10,0] = ', closest_member_histogram[0::10,0])
    #print ('closest_member_histogram[0::10,1] = ', closest_member_histogram[0::10,1])
    #print ('closest_member_histogram[0::10,2] = ', closest_member_histogram[0::10,2])
    #print ('closest_member_histogram[0::10,3] = ', closest_member_histogram[0::10,3])
    #print ('closest_member_histogram[0::10,4] = ', closest_member_histogram[0::10,4])
    #print ('closest_member_histogram[0::10,5] = ', closest_member_histogram[0::10,5])
    #print ('closest_member_histogram[0::10,6] = ', closest_member_histogram[0::10,6])
    #print ('nmemx25, ncats = ', nmemx25, ncats)
    #print ('b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, b1_std_lowrank = ',\
    #    b0_mean_lowrank, b1_mean_lowrank, b0_std_lowrank, b1_std_lowrank)
    #print ('b0_mean_midrank, b1_mean_midrank, b0_std_midrank, b1_std_midrank = ',\
    #    b0_mean_midrank, b1_mean_midrank, b0_std_midrank, b1_std_midrank)
    #print ('b0_mean_highrank, b1_mean_highrank, b0_std_highrank, b1_std_highrank = ',\
    #    b0_mean_highrank, b1_mean_highrank, b0_std_highrank, b1_std_highrank)
    
    return histogram_thresholds, closest_member_histogram, \
        nmemx25, ncats, nthresh, b0_mean_lowrank, \
        b1_mean_lowrank, b0_std_lowrank, b1_std_lowrank, \
        b0_mean_midrank, b1_mean_midrank, b0_std_midrank, \
        b1_std_midrank, b0_mean_highrank, b1_mean_highrank, \
        b0_std_highrank, b1_std_highrank
    