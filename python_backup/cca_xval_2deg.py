"""
cca.py: user inputs the month of choice and the forecast lead time. program
then loads in the HUC precipitation and the ERA-5 analysis data,
and it munges the data to make it suitable for Canonical Correlation Analysis.
It calls scipy.stats library routine to perform the CCA and saves output.
"""

from netCDF4 import Dataset
import numpy as np
from dateutils import daterange, datetohrs, dateshift
from datetime import datetime
import sys
import os
from os import path
import numpy.ma as ma
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap
import scipy.stats as stats
from sklearn.cross_decomposition import PLSRegression

# ==================================================================
    
def load_era5(cyear, cmonth, cvar):
    """ read the ERA5 reanalysis data at 4 vertical levels, 
    NHem, averaged to 2 degrees """
    infile = 'cca/'+cyear+cmonth+'_'+cvar+'_era5_2deg.cPick'
    #print (infile)
    inf = open(infile,'rb')
    input_data = cPickle.load(inf)
    ntimes, nlevels, ny, nx = np.shape(input_data)
    lon = cPickle.load(inf)
    lat = cPickle.load(inf)
    yyyymmddhh = cPickle.load(inf)
    inf.close()
    return yyyymmddhh, input_data, ntimes, nlevels, ny, nx, lon, lat

# ==================================================================
    
def load_era5_monolevel(cyear, cmonth, cvar):
    
    """ read the ERA5 reanalysis data at single levels, 
    NHem, averaged to 2 degrees """
    infile = 'cca/'+cyear+cmonth+'_'+cvar+'_era5_2deg.cPick'
    #print (infile)
    inf = open(infile,'rb')
    input_data = cPickle.load(inf)
    ntimes, ny, nx = np.shape(input_data)
    lon = cPickle.load(inf)
    lat = cPickle.load(inf)
    yyyymmddhh = cPickle.load(inf)
    inf.close()
    return yyyymmddhh, input_data, ntimes,  ny, nx, lon, lat
    
# ==================================================================
    
def convert_to_yyyymmddhh(precipDates):
    
    # --- convert Matt's HUC date array to yyyymmddhh format, making the 
    #     assumption that the 7am local time is close enough to 12Z.
    
    npdates, nymd = np.shape(precipDates)
    yyyymmddhh = []
    for i in range(npdates):
        
        yyyy = str(precipDates[i,0])
        imm = precipDates[i,1]
        if imm < 10:
            cmm = '0'+str(imm)
        else:
            cmm = str(imm)
        idd = precipDates[i,2]
        if idd < 10:
            cdd = '0'+str(idd)
        else:
            cdd = str(idd)
        yyyymmddhh.append(int(yyyy+cmm+cdd+'12'))
    return yyyymmddhh
    
# ==================================================================

def compute_n_sampledays(cmonth):
    """ compute the number of days of samples """    

    imonth = int(cmonth)
    if imonth == 1 or imonth == 3 or imonth == 5 or imonth == 7 \
    or imonth == 8 or imonth == 10 or imonth == 12:
        ndays_total = 31*nyears
    elif imonth == 4 or imonth == 6 or imonth == 9 or imonth == 11:
        ndays_total = 30*nyears 
    elif imonth == 2:
        ndays_total = 9*29 + (nyears-9)*28 
    return ndays_total
    
# ==================================================================

def read_HUC_data():
            
    """ read in the HUC data provided by Matt Switanek """

    f1 = open("cca/PRISM_pr_hucs_19810101-20180930.pickle",'rb')
    precipHUCs = cPickle.load(f1) #the daily precip accum in mm (days,hucs)
    ndayhucs, nhucs = np.shape(precipHUCs)
    precipDates = cPickle.load(f1) #the dates
    hucLats = cPickle.load(f1) #centroid lats of hucs
    hucLons = cPickle.load(f1) #centroid lons of hucs
    hucShapes = cPickle.load(f1) #embedded lists of huc boundaries
    hucDivision4 = cPickle.load(f1) #the div 4 numeric codes of the hucs
    f1.close()  
    return precipHUCs, ndayhucs, nhucs, precipDates, hucLats,\
         hucLons, hucShapes, hucDivision4      
        
# ==================================================================        
        
def set_before_after(cmonth):
    if cmonth == '01':
        cmonth_before = '12'
        cmonth_after = '02'
    elif cmonth == '02':
        cmonth_before = '01'
        cmonth_after = '03'
    elif cmonth == '03':
        cmonth_before = '02'
        cmonth_after = '04'
    elif cmonth == '04':
        cmonth_before = '03'
        cmonth_after = '05'
    elif cmonth == '05':
        cmonth_before = '04'
        cmonth_after = '06'
    elif cmonth == '06':
        cmonth_before = '05'
        cmonth_after = '07'
    elif cmonth == '07':
        cmonth_before = '06'
        cmonth_after = '08'
    elif cmonth == '08':
        cmonth_before = '07'
        cmonth_after = '09'
    elif cmonth == '09':
        cmonth_before = '08'
        cmonth_after = '10'
    elif cmonth == '10':
        cmonth_before = '09'
        cmonth_after = '11'
    elif cmonth == '11':
        cmonth_before = '10'
        cmonth_after = '12'
    elif cmonth == '12':
        cmonth_before = '11'
        cmonth_after = '01'
    else:
        print ('invalid month! ', cmonth)
        sys.exit()        

    return cmonth_before, cmonth_after 

# ==================================================================

def control_load_era5(cyears, cmonth, cmonth_before, cmonth_after, ndays_total):
    """ ERA5 was upscaled on Tom Hamill's mac24 in cca directory from
    Cathy Smith data store of ERA5 data, copied to /Public/thamill, 
    then ftped to Tom's home computer """

    ktr = 0
    for iyear, cyear in enumerate(cyears):
        
        for cm in [cmonth, cmonth_before, cmonth_after]:
    
            yyyymmddhh, temp, ntimes, nlevels, ny, nx, lon, lat = \
                load_era5(cyear, cm, 'air')
            yyyymmddhh, shum, ntimes, nlevels, ny, nx, lon, lat = \
                load_era5(cyear, cm, 'shum')    
            yyyymmddhh, uwnd, ntimes, nlevels, ny, nx, lon, lat = \
                load_era5(cyear, cm, 'uwnd')    
            yyyymmddhh, vwnd, ntimes, nlevels, ny, nx, lon, lat = \
                load_era5(cyear, cm, 'vwnd') 
            yyyymmddhh, mslp, ntimes, ny, nx, lon, lat = \
                load_era5_monolevel(cyear, cm, 'prmsl')
            yyyymmddhh, pwat, ntimes, ny, nx, lon, lat = \
                load_era5_monolevel(cyear, cm, 'pr_wtr')    
                
            zero_store = np.zeros((ntimes, nlevels, ny, nx), dtype=np.float32)    
            if iyear == 0 and cm == cmonth:
                temp_store = np.zeros((ndays_total,4,ny//2, nx), dtype=np.float64) 
                shum_store = np.zeros((ndays_total,4,ny//2, nx), dtype=np.float64)
                uwnd_store = np.zeros((ndays_total,4,ny//2, nx), dtype=np.float64)
                vwnd_store = np.zeros((ndays_total,4,ny//2, nx), dtype=np.float64)
                mslp_store = np.zeros((ndays_total, ny//2, nx), dtype=np.float64)
                pwat_store = np.zeros((ndays_total, ny//2, nx), dtype=np.float64)
                yyyymmddhh_store = np.zeros((ndays_total), dtype=np.int32)
                
            shum = np.where(shum < 0.0, zero_store, shum)
            temp_store[ktr:ktr+ntimes,:,:,:] = temp[:,:,0:ny//2,:]
            shum_store[ktr:ktr+ntimes,:,:,:] = shum[:,:,0:ny//2,:]**0.3333
            uwnd_store[ktr:ktr+ntimes,:,:,:] = uwnd[:,:,0:ny//2,:]
            vwnd_store[ktr:ktr+ntimes,:,:,:] = vwnd[:,:,0:ny//2,:]
            mslp_store[ktr:ktr+ntimes,:,:]   = mslp[:,0:ny//2,:]
            pwat_store[ktr:ktr+ntimes,:,:]   = pwat[:,0:ny//2,:]
            yyyymmddhh_store[ktr:ktr+ntimes] = yyyymmddhh[:]
            ktr = ktr + ntimes 
    
    return zero_store, temp_store, shum_store, uwnd_store, \
        vwnd_store, mslp_store, pwat_store, yyyymmddhh_store, \
        ktr, nx, ny//2        

# ==================================================================

def standard_normal(data_store,ktr):
    """ computed standardized anomalies, return in same data structure """
    
    ndims = np.ndim(data_store)    
    data_stddev = np.std(data_store, axis=0)
    data_mean = np.mean(data_store, axis=0)
    for idate in range(ktr):
        if ndims == 4:
            data_store[idate,:,:,:] = (data_store[idate,:,:,:] - data_mean[:,:,:]) / data_stddev[:,:,:]
        else:
            data_store[idate,:,:] = (data_store[idate,:,:] - data_mean[:,:]) / data_stddev[:,:]
    return data_store
    
# ==================================================================
       
def munge_predictors_predictand_to_2darray(nx, ny, ktr, nhucs, \
    hucnumber, temp_store, shum_store, uwnd_store, vwnd_store, \
    mslp_store, pwat_store, precipHUCs, yyyymmddhh_hucs, \
    yyyymmddhh_store, ifhour):       

    # --- rearrange data; first index is time, second index is 1-d composite vector
    #     of the standardized temperature, specific humidity, and u- and v-wind 
    #     components

    ngrid = nx*ny*4
    X = np.zeros((ktr, ngrid*4 + 2*nx*ny), dtype=np.float64) # 4 as there are four fields
    #Y = np.zeros((ktr, 1), dtype=np.float64)
    Y = np.zeros((ktr, nhucs), dtype=np.float64)
    for idate in range(ktr):
        temp_1d = np.reshape(temp_store[idate,:,:,:], ngrid)
        shum_1d = np.reshape(shum_store[idate,:,:,:], ngrid)
        uwnd_1d = np.reshape(uwnd_store[idate,:,:,:], ngrid)
        vwnd_1d = np.reshape(vwnd_store[idate,:,:,:], ngrid)
        mslp_1d = np.reshape(mslp_store[idate,:,:], ngrid//4)
        pwat_1d = np.reshape(pwat_store[idate,:,:], ngrid//4)
        X[idate,0:ngrid] = temp_1d[:]
        X[idate,ngrid:2*ngrid] = shum_1d[:]
        X[idate,2*ngrid:3*ngrid] = uwnd_1d[:]
        X[idate,3*ngrid:4*ngrid] = vwnd_1d[:]
        X[idate,4*ngrid:4*ngrid+nx*ny] = mslp_1d[:]
        X[idate,4*ngrid+nx*ny:4*ngrid+2*nx*ny] = pwat_1d[:]
                
        # --- for an n.5 -day forecast, pluck the HUC data offset by +n.5 days.
        #     got to get this data individually for the given month of each year.
        #     yyyymmddhh_store has the vector of initial condition dates.
        #     yyyymmddhh_hucs has the vector of HUC verification period (end)
    
        fcst_date = int(dateshift(str(yyyymmddhh_store[idate]), ifhour))
        timeindex = yyyymmddhh_hucs.index(fcst_date)
        Y[idate,:] = precipHUCs[timeindex,:]**0.3333
        #Y[idate,0] = precipHUCs[timeindex,hucnumber]**0.3333
        
    return X, Y

# ==================================================================

def separate_train_validation (X, Y, ixval, nxval, ktr):

    """separate data into training and validation parts
    for cross validation """
    
    ixstart = (ktr*ixval) // nxval
    ixend = ixstart + ktr//nxval
    ival = [*range(ixstart, ixend)]
    nval = len(ival)
    itrain = [*range(ktr)]
    del itrain[ixstart:ixend]
    ntrain = len(itrain)
    Xtrain = X[itrain,:]
    Ytrain = Y[itrain,:]
    Xval = X[ival,:]
    Yval = Y[ival,:]
    return Xtrain, Xval, Ytrain, Yval, ntrain, \
        nval, itrain, ival

# ==================================================================
        
# --- get inputs from command line
#cmonth = sys.argv[1] # 01 to 12
#clead = sys.argv[2] # forecast lead time in days, e.g., 2.5   Use half days so the 00 UTC
cyears = \
    ['1981', '1982', '1983', '1984','1985', '1986', \
    '1987', '1988', '1989', '1990', '1991', '1992', \
    '1993', '1994', '1995', '1996', '1997', '1998', \
    '1999', '2000', '2001', '2002', '2003', '2004', \
    '2005', '2006', '2007', '2008', '2009', '2010', \
    '2011', '2012', '2013', '2014', '2015', '2016'] # , '2017']
nyears = len(cyears)
   
#for cmonth in ['01','02','03','04','05','06','07','08','09','10','11']:

# --- perform cross validation over 6 groups of 6 years.

nxval = 6
hucnumber = 112
for cmonth in ['03']:
    cmonth_before, cmonth_after = set_before_after(cmonth)
    nc = compute_n_sampledays(cmonth) 
    na = compute_n_sampledays(cmonth_after) 
    nb = compute_n_sampledays(cmonth_before)
    ndays_total = na+nb+nc
    for clead in ['1.5']:
        ifhour = int(float(clead)*24)
        
        # ---- read analyzed precipitation in HUC4 boundaries. 
        
        precipHUCs, ndayhucs, nhucs, precipDates, hucLats,\
            hucLons, hucShapes, hucDivision4 = read_HUC_data()
            
        # ---- convert the precipitation dates into yyyymmddhh format used w. ERA5 data

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('Converting HUC precipitation dates.  Current time is ', current_time)
        yyyymmddhh_hucs = convert_to_yyyymmddhh(precipDates)

        # ---- for the chosen month, load the ERA5 analysis data over the multiple 
        #      years

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ("Loading ERA5 data.  Current time = ", current_time)
        zero_store, temp_store, shum_store, uwnd_store, vwnd_store, mslp_store, \
            pwat_store, yyyymmddhh_store, ktr, nx, ny = \
            control_load_era5(cyears, cmonth, cmonth_before, cmonth_after, ndays_total)
    
        # ---- reshape arrays into predictors (X) and predictand (Y).   Y is
        #      exponentiated (**0.3333)

        temp_store = standard_normal(temp_store, ktr)
        shum_store = standard_normal(shum_store, ktr)
        uwnd_store = standard_normal(uwnd_store, ktr)
        vwnd_store = standard_normal(vwnd_store, ktr)
        mslp_store = standard_normal(mslp_store, ktr)
        pwat_store = standard_normal(pwat_store, ktr)

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ("Calling munge_predictors_predictand_to_2darray .  Current time = ", current_time)
        X, Y = munge_predictors_predictand_to_2darray \
            (nx, ny, ktr, nhucs, hucnumber, temp_store, \
            shum_store, uwnd_store, vwnd_store, mslp_store, \
            pwat_store, precipHUCs, yyyymmddhh_hucs, \
            yyyymmddhh_store, ifhour)    

        # --- convert to a standard normal deviate. apply standard deviation to Y, HUC data  
    
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ("Computing mean and std dev.  Current time = ", current_time)
        Ystd = np.std(Y,axis=0)
        Ymean = np.mean(Y,axis=0)
        for idate in range(ktr):
            Y[idate,:] = (Y[idate,:]-Ymean[:]) / Ystd[:]
        Ypred_full = np.copy(Y)  

        for ixval in range(nxval): 

            print ('**** performing cross validation ',ixval,' of ',nxval)
            
            # --- separate the data into training and validation parts. itrain
            #     and ival contain the associated indices of training, validation data
            #     within the X, Y arrays
            
            print ('calling separate train_validation, ktr = ', ktr)
            print ('np.shape(X) = ', np.shape(X))
            print ('np.shape(Y) = ', np.shape(Y))
            Xtrain, Xval, Ytrain, Yval, ntrain, nval, itrain, ival = \
                separate_train_validation (X, Y, ixval, nxval, ktr)

            # --- perform the partial-least-squares regression on the
            #     training data.

            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print ("Starting PLSR.  Current time = ", current_time)
            plsr = PLSRegression(n_components=40, scale=False) # , max_iter=150)
            plsr.fit(Xtrain,Ytrain)
            Ypred = plsr.predict(Xval)
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print ("Finished PLSR.  Current time = ", current_time)
            
            # ---- reconstitute the full vector of forecast data, including
            #      back transformation to original space.  Back out the 
            #      indexing for cross_validation 
              
            Ypred_full[ival,:] = (Ypred[:,:]*Ystd[:] + Ymean[:])**3.0
            
        # ---- save cross-validated predicted precipitation forecasts to file
        
        for i in range(len(yyyymmddhh_store)):
            if yyyymmddhh_store[i] == 2003031700: 
                print (i-2,yyyymmddhh_store[i-2], Ypred_full[i-2,hucnumber])
                print (i-1,yyyymmddhh_store[i-1], Ypred_full[i-1,hucnumber])
                print (i,yyyymmddhh_store[i], Ypred_full[i,hucnumber])
                print (i+1,yyyymmddhh_store[i+1], Ypred_full[i+1,hucnumber])
                print (i+2,yyyymmddhh_store[i+2], Ypred_full[i+2,hucnumber])

        outfile = 'cca/PLSR_regression_data_month='+cmonth+'_lead='+clead+'days.cPick'
        print ('writing to ', outfile)
        ouf = open(outfile, 'wb')
        cPickle.dump(Ypred_full, ouf)
        cPickle.dump(yyyymmddhh_store, ouf)
        ouf.close()
      