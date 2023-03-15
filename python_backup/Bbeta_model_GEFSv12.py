
"""
Bbeta_model_GEFSv12.py
thin gridded bias time series data down to just the desired 3 months.  
form a model for the bias-correction error covariance, 
including covariance localization.   Save the resulting B_beta 
to file and return localized covariance.  If already, just read from file.
"""

#from fec_2D_4D_f90 import fec_2D_4D_f90a
import numpy as np
import os, sys
import _pickle as cPickle
from datetime import datetime

# ----------------------------------------------------------

def calculate_cov(x1, x2, n):   
    x1mean = np.sum(x1) / n
    x2mean = np.sum(x2) / n
    cov = np.sum( (x1-x1mean) * (x2-x2mean) ) / (n-1)
    return cov
    
# ----------------------------------------------------------

cseason = sys.argv[1]
clead = sys.argv[2]
cpath_beta = '/Volumes/Backup Plus/gefsv12/t2m/beta/'
efold = 800.0
exponenty = 2.0
already = False

if already == False:
     
    # ---- read in lat/lons
       
    infile = '/Volumes/Backup Plus/gefsv12/t2m/gefsv12_latlon_subset.cPick'
    inf = open(infile,'rb')
    latsf = cPickle.load(inf)
    lonsf = cPickle.load(inf)
    nlats, nlons = np.shape(latsf)
    npts = nlats*nlons
    inf.close()
        
    # ---- produce estimate of the localized covariance of bias-corrected 
    #      forecast errors between grid points and the localized, inverted 
    #      covariance matrix.

    if cseason == 'JFM': 
        mmddhh_begin = 10100
        mmddhh_end = 33100
    elif cseason == 'AMJ': 
        mmddhh_begin = 40100
        mmddhh_end = 63000
    elif cseason == 'JAS': 
        mmddhh_begin = 70100
        mmddhh_end = 93000
    elif cseason == 'OND': 
        mmddhh_begin = 100100
        mmddhh_end = 123100
            
    # ---- load the bias correction file

    bias_file = 'bias_correction_decayavg_lead'+clead+'.cPick'
    inf = open(bias_file, 'rb')
    bias_3d = cPickle.load(inf)
    date_list_anal = cPickle.load(inf)
    ndates = len(date_list_anal)
    inf.close()     

    # ---- make list of valid dates
        
    date_list_good = []
    index_good = []
    for idate, date in enumerate(date_list_anal):
        immddhh = int(date[2:8])
        if immddhh >= mmddhh_begin and immddhh <= mmddhh_end:
            date_list_good.append(date)
            index_good.append(idate)
    ndatesub = len(date_list_good)
                
    bias_3d_season = np.zeros((ndatesub, nlats, nlons), dtype=np.float32)
    for idate,date in enumerate(date_list_good):
        bias_3d_season[idate,:,:] = bias_3d[index_good[idate],:,:]
        
    # ---- compute covariances in 2-D space, not 4-D.   Localize
        
    Bbeta_localized = np.zeros((npts,npts), dtype=np.float64)
    Bbeta_localized_4D = np.zeros((nlats,nlons,nlats,nlons), dtype=np.float64)
    #Bbeta_localized, Bbeta_localized_4D = \
    #    fec_2D_4D_f90a(bias_3d_season, \
    #    efold, exponenty, npts, ndatesub, nlats, nlons)         
    
    ktr1 = 0
    for i1 in range(nlons):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('processing ',i1, current_time)
        for j1 in range(nlats):
            x1 = bias_3d_season[:,j1,i1]
            ktr2 = 0
            for i2 in range(nlons):
                for j2 in range(nlats):
                    x2 = bias_3d_season[:,j2,i2]
                    hdist = (111./2.) * np.sqrt( np.float(i1-i2)**2 + np.float(j1-j2)**2)
                    localizn_factor = np.exp(-(hdist/efold)**exponenty)
                    covv = calculate_cov(x1,x2,ndatesub)
                    Bbeta_localized[ktr1,ktr2] = covv*localizn_factor
                    Bbeta_localized[ktr2,ktr1] = Bbeta_localized[ktr1,ktr2]
                    Bbeta_localized_4D[j1,i1,j2,i2] = Bbeta_localized[ktr1,ktr2]
                    Bbeta_localized_4D[j2,i2,j1,i1] = Bbeta_localized[ktr1,ktr2]
                    ktr2 = ktr2 + 1
            ktr1 = ktr1 + 1   
                
    # ---- write the Bx_localized and Bx_localized_inverse to pickle file.

    outfile = cpath_beta+'Localized_Bbeta_'+cseason+\
            '_lead='+clead+'_'+str(efold)+'.cPick'
    print ('writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(Bbeta_localized_4D, ouf)
    cPickle.dump(Bbeta_localized, ouf)
    ouf.close()
else:
    infile = cpath_beta+'Localized_Bbeta_'+cseason+\
        '_lead='+clead+'_'+str(efold)+'.cPick'
    inf = open(infile,'rb')
    Bbeta_localized_4D = cPickle.load(inf)
    Bbeta_localized = cPickle.load(inf)
    inf.close()
        
print ('Done')