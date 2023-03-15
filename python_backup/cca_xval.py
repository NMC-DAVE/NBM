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
    NHem, averaged to 5 degrees """
    infile = 'cca/'+cyear+cmonth+'_'+cvar+'_era5_5deg.cPick'
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
    
def convert_to_yyyymmddhh(precipDates):
    
    # --- convert Matt's HUC date array to yyyymmddhh format, making the 
    #     assumption that the 7am local time is close enough to 12Z.
    
    npdates, nymd = np.shape(precipDates)
    yyyymmddhh = []
    #print ('npdates, nymd = ', npdates, nymd)
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
        #print (precipDates[i,0], precipDates[i,1], \
        #precipDates[i,2], int(yyyy+cmm+cdd+'12'))
        #if i == 1000: sys.exit()
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

def control_load_era5(cyears, cmonth):
    """ ERA5 was upscaled on Tom Hamill's mac24 in cca directory from
    Cathy Smith data store of ERA5 data, copied to /Public/thamill, 
    then ftped to Tom's home computer """

    ktr = 0
    for iyear, cyear in enumerate(cyears):
    
        yyyymmddhh, temp, ntimes, nlevels, ny, nx, lon, lat = \
            load_era5(cyear, cmonth, 'air')
        yyyymmddhh, shum, ntimes, nlevels, ny, nx, lon, lat = \
            load_era5(cyear, cmonth, 'shum')    
        yyyymmddhh, uwnd, ntimes, nlevels, ny, nx, lon, lat = \
            load_era5(cyear, cmonth, 'uwnd')    
        yyyymmddhh, vwnd, ntimes, nlevels, ny, nx, lon, lat = \
            load_era5(cyear, cmonth, 'vwnd') 
                
        if iyear == 0:
            zero_store = np.zeros((ntimes, nlevels, ny, nx), dtype=np.float32)
            temp_store = np.zeros((ndays_total,4,ny//2, nx), dtype=np.float64) 
            shum_store = np.zeros((ndays_total,4,ny//2, nx), dtype=np.float64)
            uwnd_store = np.zeros((ndays_total,4,ny//2, nx), dtype=np.float64)
            vwnd_store = np.zeros((ndays_total,4,ny//2, nx), dtype=np.float64)
            yyyymmddhh_store = np.zeros((ndays_total), dtype=np.int32)
                
        shum = np.where(shum < 0.0, zero_store, shum)
        temp_store[ktr:ktr+ntimes,:,:,:] = temp[:,:,0:ny//2,:]
        shum_store[ktr:ktr+ntimes,:,:,:] = shum[:,:,0:ny//2,:]**0.3333
        uwnd_store[ktr:ktr+ntimes,:,:,:] = uwnd[:,:,0:ny//2,:]
        vwnd_store[ktr:ktr+ntimes,:,:,:] = vwnd[:,:,0:ny//2,:]
        yyyymmddhh_store[ktr:ktr+ntimes] = yyyymmddhh[:]
        ktr = ktr + ntimes 
    return zero_store, temp_store, shum_store, uwnd_store, \
        vwnd_store, yyyymmddhh_store, ktr, nx, ny//2,        

# ==================================================================

def standard_normal(data_store,ktr):
    """ computed standardized anomalies, return in same data structure """
    
    data_stddev = np.std(data_store, axis=0)
    data_mean = np.mean(data_store, axis=0)
    for idate in range(ktr):
        data_store[idate,:,:,:] = (data_store[idate,:,:,:] - data_mean[:,:,:]) / data_stddev[:,:,:]
    return data_store
    
# ==================================================================
       
def munge_predictors_predictand_to_2darray(nx, ny, ktr, nhucs, temp_store, \
    shum_store, uwnd_store, vwnd_store, precipHUCs, yyyymmddhh_hucs, \
    yyyymmddhh_store, ifhour):       

    # --- rearrange data; first index is time, second index is 1-d composite vector
    #     of the standardized temperature, specific humidity, and u- and v-wind 
    #     components

    ngrid = nx*ny*4
    X = np.zeros((ktr, ngrid*4), dtype=np.float64) # 4 as there are four fields
    Y = np.zeros((ktr, nhucs), dtype=np.float64)
    print ('ktr = ', ktr)
    for idate in range(ktr):
        temp_1d = np.reshape(temp_store[idate,:,:,:], ngrid)
        shum_1d = np.reshape(shum_store[idate,:,:,:], ngrid)
        uwnd_1d = np.reshape(uwnd_store[idate,:,:,:], ngrid)
        vwnd_1d = np.reshape(vwnd_store[idate,:,:,:], ngrid)
        X[idate,0:ngrid] = temp_1d[:]
        X[idate,ngrid:2*ngrid] = shum_1d[:]
        X[idate,2*ngrid:3*ngrid] = uwnd_1d[:]
        X[idate,3*ngrid:4*ngrid] = vwnd_1d[:]
                
        # --- for an n.5 -day forecast, pluck the HUC data offset by +n.5 days.
        #     got to get this data individually for the given month of each year.
        #     yyyymmddhh_store has the vector of initial condition dates.
        #     yyyymmddhh_hucs has the vector of HUC verification period (end)
    
        #print ('str(yyyymmddhh_store[idate] = ', str(yyyymmddhh_store[idate]))
        fcst_date = int(dateshift(str(yyyymmddhh_store[idate]), ifhour))
        timeindex = yyyymmddhh_hucs.index(fcst_date)
        Y[idate,:] = precipHUCs[timeindex,:]**0.3333
        
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

    #ival = range(ixval,ktr,nxval)
    #itrain = []
    #for ix in range(nxval):
    #    if itrain != ixval:
    #        itrain = itrain + range(ix,ktr,nxval)
    #ntrain = len(itrain)
    
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
for cmonth in ['03']:
    ndays_total = compute_n_sampledays(cmonth)
    for clead in ['1.5']:
        ifhour = int(float(clead)*24)
        for ixval in range(nxval): 
            
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

            zero_store, temp_store, shum_store, uwnd_store, vwnd_store, \
                yyyymmddhh_store, ktr, nx, ny = control_load_era5(cyears, cmonth)
    
            # --- convert to a standard normal deviate.

            temp_store = standard_normal(temp_store, ktr)
            shum_store = standard_normal(shum_store, ktr)
            uwnd_store = standard_normal(uwnd_store, ktr)
            vwnd_store = standard_normal(vwnd_store, ktr)

            X, Y = munge_predictors_predictand_to_2darray \
                (nx, ny, ktr, nhucs, temp_store, \
                shum_store, uwnd_store, vwnd_store, precipHUCs, \
                yyyymmddhh_hucs, yyyymmddhh_store, ifhour)      

            # --- apply standard deviation to Y, HUC data  
    
            Ystd = np.std(Y,axis=0)
            Ymean = np.mean(Y,axis=0)
            for idate in range(ktr):
                Y[idate,:] = (Y[idate,:]-Ymean[:]) / Ystd[:]

            # --- separate the data into training and validation parts
            
            Xtrain, Xval, Ytrain, Yval, ntrain, nval, itrain, ival = \
                separate_train_validation (X, Y, ixval, nxval, ktr)

            # --- perform the partial-least-squares regression on the
            #     training data.

            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print ("Starting PLSR.  Current time = ", current_time)
            plsr = PLSRegression(n_components=160, scale=False)
            plsr.fit(Xtrain,Ytrain)
            Ypred = plsr.predict(Xval)
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print ("Finished PLSR.  Current time = ", current_time)
            
            # ---- reconstitute the full vector of forecast data, including
            #      back transformation to original space.

            if ixval == 0: Ypred_full = np.copy(Y)
            print ('np.shape(Ypred_full[ival,:]) = ', np.shape(Ypred_full[ival,:]))
            print ('np.shape(Ypred[:,:]) = ', np.shape(Ypred[:,:]))
            print ('np.shape(Ystd[:]) = ', np.shape(Ystd[:]))
            print ('np.shape(Ymean[:]) = ', np.shape(Ymean[:]))
            
            Ypred_full[ival,:] = (Ypred[:,:]*Ystd[:] + Ymean[:])**3.0
            
        # ---- save cross-validated predicted precipitation forecasts to file

        outfile = 'cca/PLSR_regression_data_month='+cmonth+'_lead='+clead+'days.cPick'
        print ('writing to ', outfile)
        ouf = open(outfile, 'wb')
        cPickle.dump(Ypred_full, ouf)
        cPickle.dump(yyyymmddhh_store, ouf)
        ouf.close()
      