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
#import scipy.stats.linregress as linregress

# --------------------------------------------------------------   

def form_diagonal_matrix(npts, vary):
    R = vary*np.identity(npts, dtype=np.float64)
    return R

# --------------------------------------------------------------

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

date_end = dateshift('2018123100',-2*int(clead))
date_list_anal = daterange('2000010100',date_end,24)

ndates = len(date_list_anal)
date_list_forecast = []
for i in range(ndates):
    date_list_forecast.append(dateshift(date_list_anal[i], int(clead)))
    
date_list_anal_ver = []
for i in range(ndates):
    date_list_anal_ver.append(dateshift(date_list_forecast[i], int(clead)))

# ---- loop through dates and process day by day

ndatektr = 0
ndatektr_yearly = np.zeros((18), dtype=np.int32)
for idate, date in enumerate(date_list_anal):
    
    # ---- read reanalysis appropriate at the time of this forecast for bias corr.
    
    rem = idate%30
    if rem == 0: print ('processing date = ', idate, date)
    cyear = date[0:4]
    cmm = date[4:6]
    cmmdd = date[4:8]
    imm = int(cmm)
    iyear = int(cyear)-2000
    datef = date_list_forecast[idate]
    cyearf = datef[0:4]
    
    infile = cpath_era5 +cyearf+'/t2m_era5_halfdegree_'+datef+'.cPick'
    fexist1 = os.path.exists(infile)
    if fexist1 == True:
        inf = open(infile, 'rb')
        analysis = cPickle.load(inf) - 273.16
        analysis = np.flipud(analysis)
        if idate == 0: 
            lats = cPickle.load(inf)
            lons = cPickle.load(inf)
            nlats, nlons = np.shape(lats)
            lats = np.flipud(lats)
            analysis_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            forecast_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            slope = np.zeros((12,nlats,nlons), dtype=np.float64)
            intercept = np.zeros((12,nlats,nlons), dtype=np.float64)
    else:
        print ('1. did not find file ', infile)
            
    # ---- read the forecast information for bias corr.
            
    cyear = date[0:4]    
    cpath_forecast = cpath_gefsv12+cyear+'/'
    infile = cpath_forecast + date + '_lead'+clead+'_conus_0.5deg_hour'+clead+'.cPick'
    fexist2 = os.path.exists(infile)
    if fexist2 == True:
        inf = open(infile,'rb')
        forecast = cPickle.load(inf)
        inf.close()
        if cmmdd == '0101':
            infile = '/Volumes/Backup Plus/gefsv12/t2m/gefsv12_latlon_subset.cPick'
            inf = open(infile,'rb')
            latsf = cPickle.load(inf)
            lonsf = cPickle.load(inf)
            inf.close()
    else:
        print ('2. did not find file ', infile)

    analysis_3d[idate,:,:] = analysis[:,:] 
    forecast_3d[idate,:,:] = forecast[:,:] 
    
    
# ---- print forecast/analysis

print ('sample date anal fcst ')
for idate in range(100):
    print (idate, analysis_3d[idate,nlats//2, nlons//2])
    
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
    forecast_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)
    analysis_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)

    ktr = 0
    for idate, date in enumerate(date_list_anal):
        if datevalid_indices[idate] == 1:
            analysis_validdates[ktr,:,:] = analysis_3d[idate,:,:]
            forecast_validdates[ktr,:,:] = forecast_3d[idate,:,:]
            ktr = ktr+1
    
    print ('valid idate, analysis, fcst ')
    for idate in range(ktr):
        print (idate, analysis_validdates[ktr,nlats//2, nlons//2], forecast_validdates[ktr,nlats//2, nlons//2])      
    # ---- perform linear regression analysis for each grid point
    
    for ix in range(nlons):
        for jy in range(nlats):
            x = forecast_validdates[:,jy,ix]
            y = analysis_validdates[:,jy,ix]
            slope_out, intercept_out, r_value, p_value, std_err = scipy.stats.linregress(x,y)
            slope[imonth,jy,ix] = slope_out
            intercept[imonth,jy,ix] = intercept_out
            if ix == nlons//2 and jy == nlats//2
                print ('slope, intercerpt = ', slope_out[imonth, jy, ix], intercept[imonth,jy,ix])
    
    # ---- write the slope and intercept for this month to file.

    print ('   mean slope = ',np.mean(slope[imonth,:,:]))
    print ('   mean intercept = ',np.mean(intercept[imonth,:,:]))
    
    outfile = cpath_beta+'MOS_slope_intercept_'+cmonth+\
        '_lead='+clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(slope[imonth,:,:], ouf)
    cPickle.dump(intercept[imonth,:,:], ouf)
    ouf.close()
    print ('   done writing')
    
    
    
   




    
    

    
    
    