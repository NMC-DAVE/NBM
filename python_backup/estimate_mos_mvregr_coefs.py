"""
estimate_mos_mvregr_coefs.py

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
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
import math
import scipy
import statsmodels.api as sm
from decayavg_biascorrection2 import decayavg_biascorrection2
import matplotlib.pyplot as plt
from quantile_mapping import quantile_mapping

# --------------------------------------------------------------  

def set_alpha(clead):
    if clead == '12':
        alpha = 0.16
    elif clead == '24':
        alpha = 0.14  
    elif clead == '36':
        alpha = 0.08
    elif clead == '48':
        alpha = 0.06 
    elif clead == '60':
        alpha = 0.06 
    elif clead == '72':
        alpha = 0.06 
    elif clead == '84':
        alpha = 0.06
    elif clead == '96':
        alpha = 0.04 
    elif clead == '108':
        alpha = 0.04 
    elif clead == '120':
        alpha = 0.04  
    elif clead == '132':
        alpha = 0.04
    return alpha
    
# -------------------------------------------------------------- 

def read_lsmask():

    infilename = 'lsmask_0p5.cPick'
    print (infilename)
    inf = open(infilename, 'rb')
    lsmask = cPickle.load(inf)
    inf.close()
    return lsmask
  
# --------------------------------------------------------------   

clead = sys.argv[1]
alpha = set_alpha(clead)
ilead = int(clead)
efold_timescale = 8.0
iskip = int(clead) // 24
if clead == '12' or clead == '36' or clead == '60' or \
    clead == '84' or clead == '108'  or clead == '132':
    iskip = iskip+1
lsmask = read_lsmask()
already_qmap = False

# ---- various initialization

cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'

date_end = dateshift('2018123100',-2*int(clead))
date_list_anal = daterange('2000010100',date_end,24)
ndates = len(date_list_anal)

date_list_forecast = []
for i in range(ndates):
    date_list_forecast.append(dateshift(date_list_anal[i], int(clead)))
    
date_list_anal_ver = []
for i in range(ndates):
    date_list_anal_ver.append(dateshift(date_list_forecast[i], int(clead)))

# ---- read estimated 2000-2018 climatology to file

if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
    infile = cpath_era5 + 'ERA5_temperature_climatology_12UTC.cPick'
else:
    infile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'

print ('reading from ', infile)
inf = open(infile,'rb')
climo_temps_estimated = cPickle.load(inf)
inf.close()

# ---- loop through dates and process day by day

ndatektr = 0
ndatektr_yearly = np.zeros((18), dtype=np.int32)
for idate, date in enumerate(date_list_anal):
    
    # ---- read reanalysis appropriate at the time of this forecast for bias corr.
    
    rem = idate%30
    if rem == 0: print ('processing date = ', idate, date)

    datef = date_list_forecast[idate]
    cyearf = datef[0:4]
    iyear = int(cyearf)-2000
    cmm = datef[4:6]
    imm = int(cmm)
    idd = int(datef[6:8])
    iyear_full = int(cyearf)
    julday = dayofyear(iyear_full, imm, idd) - 1
    if julday > 364: julday = 364
    
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
            differences_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            analyses_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            forecast_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            analyses_3d_raw = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            forecast_3d_raw = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            qmap_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            ones = np.ones((nlats, nlons), dtype=np.float64)
            print ('arrays established')
    else:
        print ('1. did not find file ', infile)
            
    # ---- read the t2m forecast information for bias corr.
            
    cyear = date[0:4]    
    cpath_forecast = cpath_gefsv12+cyear+'/'
    infile = cpath_forecast + date + '_lead'+clead+'_conus_0.5deg_hour'+clead+'.cPick'
    fexist2 = os.path.exists(infile)
    if fexist2 == True:
        inf = open(infile,'rb')
        forecast = cPickle.load(inf)
        inf.close()
    else:
        print ('2. did not find file ', infile)       
    
    # ---- subtract off climatology, save info to 3D arrays.

    #analyses_3d[idate,:,:] = analysis[:,:] - climo_temps_estimated[julday,:,:]
    #forecast_3d[idate,:,:] = forecast[:,:] - climo_temps_estimated[julday,:,:]
    analyses_3d[idate,:,:] = analysis[:,:] 
    forecast_3d[idate,:,:] = forecast[:,:] 
    analyses_3d_raw[idate,:,:] = analysis[:,:] 
    forecast_3d_raw[idate,:,:] = forecast[:,:]

# ---- get 2000-2018 decaying average bias correction
    
beta_3d = np.zeros((ndates, nlats, nlons), dtype=np.float64)
npts = ndates
beta_3d = decayavg_biascorrection2(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, lsmask, \
    date_list_forecast, alpha)    
    
# ---- get 2000-2018 quantile mapping

if already_qmap == False:      
    print ('performing quantile mapping')
    qmap_3d = quantile_mapping(nlats, nlons, \
        forecast_3d_raw, analyses_3d_raw, qmap_3d, lsmask, \
        date_list_forecast, cpath_gefsv12, cpath_era5, clead)
    
    outfile = 'qmap_f'+clead+'.cPick'
    ouf = open(outfile,'wb')
    cPickle.dump(qmap_3d, ouf)
    ouf.close()
else:
    infile = 'qmap_f'+clead+'.cPick'
    inf = open(infile,'rb')
    qmap_3d = cPickle.load(inf)
    inf.close()
    
# ---- loop thru months

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
regr_coefs = np.zeros((12,4,nlats,nlons), dtype=np.float64)
#regr_coefs = np.zeros((12,3,nlats,nlons), dtype=np.float64)
tvalues_monthly = np.zeros((12,4),dtype=np.float64)
#tvalues_monthly = np.zeros((12,3),dtype=np.float64)

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
    forecast_validdates_raw = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)
    analysis_validdates_raw = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)
    beta_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)
    qmap_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float64)

    ktr = 0
    for idate, date in enumerate(date_list_anal):
        if datevalid_indices[idate] == 1:
            analysis_validdates[ktr,:,:] = analyses_3d[idate,:,:]
            forecast_validdates[ktr,:,:] = forecast_3d[idate,:,:]
            qmap_validdates[ktr,:,:] = qmap_3d[idate,:,:] 
            #if cmonth == "Jun":
            #    print ('id,d,qm,f,qm3d = ',idate, date, qmap_validdates[ktr,59,72], \
            #        forecast_3d_raw[idate,59,72], qmap_3d[idate,59,72])
            
            if idate - iskip >= 0:
                beta_validdates[ktr,:,:] = beta_3d[idate-iskip,:,:]
            else:
                beta_validdates[ktr,:,:] = 0.0
            ktr = ktr+1
        
    # ---- perform multiple linear regression analysis for each grid point

    ktrsoil = 0
    for jy in range(nlats):
        for ix in range(nlons):
    #for jy in [59]:
    #    for ix in [72]:
            
            if lsmask[jy,ix] == 1:
                
                X = np.zeros((ktr,4), dtype=np.float64)
                #X = np.zeros((ktr,3), dtype=np.float64)
                X[:,0] = 1.0
                X[:,1] = forecast_validdates[:,jy,ix]
                X[:,2] = beta_validdates[:,jy,ix]
                X[:,3] = qmap_validdates[:,jy,ix]
            
                y = np.zeros((ktr), dtype=np.float64)
                y[:] = analysis_validdates[:,jy,ix]

                model = sm.OLS(y, X)
                results = model.fit()
                tvalues = results.tvalues
            
                tvalues_monthly[imonth,:] = tvalues_monthly[imonth,:] + \
                    np.abs(tvalues[:])
                regr_coefs[imonth,:,jy,ix] = results.params
                
                #print ('min, max analysis_validdates = ', np.min(analysis_validdates[:,59,72]), np.max(analysis_validdates[:,59,72]))
                #print ('min, max forecast_validdates = ', np.min(forecast_validdates[:,59,72]), np.max(forecast_validdates[:,59,72]))
                #print ('min, max beta_validdates = ', np.min(beta_validdates[:,59,72]), np.max(beta_validdates[:,59,72]))
                #print ('min, max qmap_validdates = ', np.min(qmap_validdates[:,59,72]), np.max(qmap_validdates[:,59,72]))
                #print ('min, max y = ', np.min(y), np.max(y))  
                #print ('regr_coefs= ', regr_coefs[imonth,:,jy,ix] )
                
                
            else: 
                
                regr_coefs[imonth,:,jy,ix] = 0.0      
                
    # --- get average t values and report
    
    tvalues_monthly[imonth,:] = tvalues_monthly[imonth,:] / np.sum(lsmask)
    print ('    avg tvalue[const] = ',tvalues_monthly[imonth,0])
    print ('    avg tvalue[t2m] = ',tvalues_monthly[imonth,1])
    print ('    avg tvalue[beta] = ',tvalues_monthly[imonth,2])
    print ('    avg tvalue[qmap] = ',tvalues_monthly[imonth,3])

    outfile = cpath_gefsv12+'MOS_multiple_regr_'+cmonth+\
       '_lead='+clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(regr_coefs[imonth,:,:,:], ouf)
    ouf.close()
    
sys.exit()   
    
