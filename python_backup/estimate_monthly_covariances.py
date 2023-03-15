"""
estimate_monthly_covariances.py

"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
import numpy as np
import numpy.ma as ma
import scipy.signal as signal
import scipy.stats as stats
import scipy
from astropy.convolution import convolve
from numba import jit
from forecast_error_covariance_twod_fourd_v2_f90 import forecast_error_covariance_twod_fourd_v2_f90
from reformat_gain_to_4d_f90 import reformat_gain_to_4d_f90
from estimate_analysis_error_covariance_f90 import estimate_analysis_error_covariance_f90

# --------------------------------------------------------------   

def form_diagonal_matrix(npts, vary):
    R = vary*np.identity(npts, dtype=np.float64)
    return R

# --------------------------------------------------------------

@jit
def calculate_cov(x1, x2, n):
    x1mean = np.sum(x1) / n
    x2mean = np.sum(x2) / n
    cov = np.sum( (x1-x1mean) * (x2-x2mean) ) / (n-1)
    return cov

# ---- various initialization

clead = sys.argv[1]
iskip = int(clead)//24
cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_beta = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_random = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_gain = '/Volumes/Backup Plus/gefsv12/t2m/'
efold_bias = 1500.
efold_random = 600.
efold_analysis = 100.
exponenty = 2.0
anal_err_var = 1.0


# ---- load the bias correction time series

bias_file = 'bias_correction_2000_2018_convolved_lead'+clead+'.cPick'
inf = open(bias_file, 'rb')
bias_3d = cPickle.load(inf)
date_list_anal = cPickle.load(inf)
inf.close()

# ---- load the random data time series

random_file = 'random_error_2000_2018_lead'+clead+'.cPick'
inf = open(random_file, 'rb')
random_3d = cPickle.load(inf)
inf.close()

ndates, nlats, nlons = np.shape(random_3d)
npts = nlats*nlons
print ('forming R ... ')
R = form_diagonal_matrix(npts, anal_err_var)
    
# ---- loop thru months

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

for imonth, cmonth in enumerate(cmonths):
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('**** PROCESSING month = ', cmonth, current_time)
    
    # ---- loop through dates and process day by day

    datevalid_indices = np.zeros(ndates,dtype=np.int32)
    for idate, date in enumerate(date_list_anal):
        thismonth = int(str(date)[4:6]) - 1

        if imonth == 0:
            if thismonth == imonth or thismonth == 11 or \
                thismonth == imonth+1: datevalid_indices[idate] = 1
        elif imonth == 11:
            if thismonth == imonth or thismonth == 0 or \
                thismonth == imonth-1: datevalid_indices[idate] = 1
        else:
            if thismonth == imonth or thismonth == imonth+1 or \
                thismonth == imonth-1: datevalid_indices[idate] = 1
                
    ndates_valid = np.sum(datevalid_indices)
    bias_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)
    random_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)

    ktr = 0
    for idate, date in enumerate(date_list_anal):
        if datevalid_indices[idate] == 1:
            bias_validdates[ktr,:,:] = bias_3d[idate,:,:]
            random_validdates[ktr,:,:] = random_3d[idate,:,:]
            #print (idate,date,bias_validdates[ktr,nlats/6,nlons/8], bias_validdates[ktr,nlats/3,55])
            ktr = ktr+1
            
    # ---- compute covariances in 2-D and 4-D.   Localize

    Bbeta_localized_2D = np.zeros((npts,npts), dtype=np.float64)
    Bbeta_localized_4D = np.zeros((nlats,nlons,nlats,nlons), dtype=np.float64)
    Brandom_localized_2D = np.zeros((npts,npts), dtype=np.float64)
    Brandom_localized_4D = np.zeros((nlats,nlons,nlats,nlons), dtype=np.float64)


    #print (forecast_error_covariance_twod_fourd_v2_f90.__doc__)
    print ('   computing bias covariances ... ', nlats, nlons)
    Bbeta_localized_2D, Bbeta_localized_4D = \
        forecast_error_covariance_twod_fourd_v2_f90(bias_validdates, \
        efold_bias, npts, ndates_valid, nlats, nlons)
        
    print ('   computing random covariances ... ', nlats, nlons)
    Brandom_localized_2D, Brandom_localized_4D = \
        forecast_error_covariance_twod_fourd_v2_f90(random_validdates, \
        efold_random, npts, ndates_valid, nlats, nlons)

    # ---- write the Bbeta_localized 2D to pickle file.

    outfile = cpath_beta+'Localized_Bbeta_2D_'+cmonth+\
        '_lead='+clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(Bbeta_localized_2D, ouf)
    ouf.close()
    
    # ---- write the Bbeta_localized 4D to pickle file.

    outfile = cpath_beta+'Localized_Bbeta_4D_'+cmonth+\
        '_lead='+clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(Bbeta_localized_4D, ouf)
    ouf.close()
    
    # ---- write the Brandom_localized 2D to pickle file.

    outfile = cpath_random+'Localized_Brandom_2D_'+cmonth+\
        '_lead='+clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(Brandom_localized_2D, ouf)
    ouf.close()
    
    # ---- write the Brandom_localized 4D to pickle file.

    outfile = cpath_random+'Localized_Brandom_4D_'+cmonth+\
        '_lead='+clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(Brandom_localized_4D, ouf)
    ouf.close()
    
    # ---- compute Kalman gain for bias and for forecast, and write to file
    
    Bbeta_plus_Bx_plus_R = R + Brandom_localized_2D + Bbeta_localized_2D
    #Bbeta_plus_Bx_plus_R = scipy.linalg.pinvh(Bbeta_plus_Bx_plus_R, cond=0.000001)
    Bbeta_plus_Bx_plus_R_inv = np.linalg.inv(Bbeta_plus_Bx_plus_R)
    #Bbeta_localized_2D = scipy.linalg.pinvh(Bbeta_localized_2D, cond=0.000001)
    print ('   inverse calculated')
    
    Kalman_gain_beta = np.matmul(Bbeta_localized_2D, Bbeta_plus_Bx_plus_R_inv)
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    Kalman_gain_beta_4d = reformat_gain_to_4d_f90(Kalman_gain_beta, nlats, nlons)
    
    gain_outfile = cpath_gain + 'GEFSv12_KFgain_'+cmonth+'_lead'+clead+'.cPick'
    print ('   writing Kalman gain to ', gain_outfile)
    ouf = open(gain_outfile, 'wb')
    cPickle.dump(Kalman_gain_beta_4d, ouf)
    ouf.close()
    print ('   done writing')




    
    

    
    
    