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
from astropy.convolution import convolve
from numba import jit
from forecast_error_covariance_twod_fourd_f90 import forecast_error_covariance_twod_fourd_f90
from reformat_gain_to_4d_f90 import reformat_gain_to_4d_f90
from estimate_analysis_error_covariance_f90 import estimate_analysis_error_covariance_f90

# --------------------------------------------------------------   

def form_diagonal_matrix(npts, vary):
    B = vary*np.identity(npts, dtype=np.float64)
    return B

# --------------------------------------------------------------

@jit
def calculate_cov(x1, x2, n):
    x1mean = np.sum(x1) / n
    x2mean = np.sum(x2) / n
    cov = np.sum( (x1-x1mean) * (x2-x2mean) ) / (n-1)
    return cov

# ---- various initialization

#print (forecast_error_covariance_twod_fourd_f90.__doc__)

clead = sys.argv[1]
iskip = int(clead)//24
cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_beta = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_random = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_gain = '/Volumes/Backup Plus/gefsv12/t2m/'
efold = 500.
efold_analysis = 200.
exponenty = 2.0
anal_err_var = 1.0

# ---- load the bias correction time series

diff_file = 'differences_2000_2018_lead'+clead+'.cPick'
inf = open(diff_file, 'rb')
differences_3d = cPickle.load(inf)
date_list_anal = cPickle.load(inf)
inf.close()

ndates, nlats, nlons = np.shape(differences_3d)
npts = nlats*nlons
#R = form_diagonal_matrix(npts, anal_err_var)
print ('forming R ... ')
R = np.zeros((npts, npts), dtype=np.float64)
#R = estimate_analysis_error_covariance_f90(npts, anal_err_var, efold_analysis)
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
    diff_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float32)

    ktr = 0
    for idate, date in enumerate(date_list_anal):
        if datevalid_indices[idate] == 1:
            diff_validdates[ktr,:,:] = differences_3d[idate,:,:]
            ktr = ktr+1
            
    # ---- compute covariances in 2-D and 4-D.   Localize

    B_localized_2D = np.zeros((npts,npts), dtype=np.float64)
    B_localized_4D = np.zeros((nlats,nlons,nlats,nlons), dtype=np.float64)
    print ('   computing B ... ', nlats, nlons)
    
    B_localized_2D, B_localized_4D = \
        forecast_error_covariance_twod_fourd_f90(diff_validdates, efold, \
    	exponenty, npts, ndates_valid, nlats, nlons)

    # ---- write the Bbeta_localized 2D to pickle file.

    outfile = cpath_beta+'Localized_B_2D_'+cmonth+\
        '_lead='+clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(B_localized_2D, ouf)
    ouf.close()
    
    # ---- write the Bbeta_localized 4D to pickle file.

    outfile = cpath_beta+'Localized_B_4D_'+cmonth+\
        '_lead='+clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(B_localized_4D, ouf)
    ouf.close()
    
    # ---- compute Kalman gain and write to file
    
    B_plus_R = R + B_localized_2D
    B_plus_R_inv = np.linalg.inv(B_plus_R)
    Kalman_gain = np.matmul(B_localized_2D, B_plus_R_inv)
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    Kalman_gain_4d = reformat_gain_to_4d_f90(Kalman_gain, nlats, nlons)
    
    gain_outfile = cpath_gain + 'GEFSv12_KFgain_together_'+cmonth+'_lead'+clead+'.cPick'
    print ('   writing Kalman gain to ', gain_outfile)
    ouf = open(gain_outfile, 'wb')
    cPickle.dump(Kalman_gain_4d, ouf)
    ouf.close()
    print ('   done writing')




    
    

    
    
    