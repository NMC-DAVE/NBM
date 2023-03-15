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

rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'
rcParams['legend.fontsize']='xx-small'
rcParams['legend.fancybox']=True


#def load_latlon():
    
#    infile = 'cca/199901_air_era5.nc'
#    nc = Dataset(infile)
#    lat = nc.variables['lat'][:] 
#    lon = nc.variables['lon'][:]
#    nc.close()
    
def load_era5(cyear, cmonth, cvar):
    
    infile = 'cca/'+cyear+cmonth+'_'+cvar+'_era5.cPick'
    #print (infile)
    inf = open(infile,'rb')
    input_data = cPickle.load(inf)
    ntimes, nlevels, ny, nx = np.shape(input_data)
    lon = cPickle.load(inf)
    lat = cPickle.load(inf)
    yyyymmddhh = cPickle.load(inf)
    inf.close()
    return yyyymmddhh, input_data, ntimes, nlevels, ny, nx, lon, lat
    
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
        #print (precipDates[i,0], precipDates[i,1], precipDates[i,2], int(yyyy+cmm+cdd+'12'))
        #if i == 1000: sys.exit()
    return yyyymmddhh
        
# --- get inputs from command line

#cmonth = sys.argv[1] # 01 to 12
#clead = sys.argv[2] # forecast lead time in days, e.g., 2.5   Use half days so the 00 UTC
cyears = ['1981', '1982', '1983', '1984','1985', '1986', '1987', '1988', '1989', '1990', \
    '1991', '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', \
    '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011']
nyears = len(cyears)
   
#for cmonth in ['01','02','03','04','05','06','07','08','09','10','11']:
for cmonth in ['03']:
    for clead in ['1.5']:

        # ---- ECMWF initial dates line up with the 12 UTC HUC dates.

        ifhour = int(float(clead)*24)
        print ('forecast lead in hours: ', ifhour)
        imonth = int(cmonth)
        if imonth == 1 or imonth == 3 or imonth == 5 or imonth == 7 \
        or imonth == 8 or imonth == 10 or imonth == 12:
            ndays_total = 31*nyears
        elif imonth == 4 or imonth == 6 or imonth == 9 or imonth == 11:
            ndays_total = 30*nyears 
        elif imonth == 2:
            ndays_total = 7*29 + (nyears-7)*28

        # --- read in the HUC data provided by Matt Switanek

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('Reading HUC data. Current time is ', current_time)
        f1 = open("PRISM_pr_hucs_19810101-20180930.pickle",'rb')
        precipHUCs = cPickle.load(f1) #the daily precip accumulations in mm (days,hucs)
        ndayhucs, nhucs = np.shape(precipHUCs)
        precipDates = cPickle.load(f1) #the dates
        hucLats = cPickle.load(f1) #centroid lats of hucs
        hucLons = cPickle.load(f1) #centroid lons of hucs
        hucShapes = cPickle.load(f1) #embedded lists of huc boundaries
        hucDivision4 = cPickle.load(f1) #the division 4 numeric codes of the hucs
        #print ('np.shape(precipHUCs) = ', np.shape(precipHUCs))
        #print ('precipDates[0:-1:100] = ', precipDates[0:-1:100]) # [2012    3   19]
        #print ('hucLons = ', hucLons) # negative for west
        #print ('hucLats = ', hucLats)
        #print ('hucShapes.__doc__ = ',hucShapes.__doc__ )
        #print ('hucDivision4 = ', hucDivision4)
        f1.close()
        #for i in range(len(hucLons)):
        #    print (i,hucLons[i],hucLats[i])
        #sys.exit()
        

        # ---- convert the precipitation dates into yyyymmddhh format

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('Converting HUC precipitation dates.  Current time is ', current_time)
        yyyymmddhh_hucs = convert_to_yyyymmddhh(precipDates)
        #print ('HUC dates = ', yyyymmddhh_hucs[0:100])
        #sys.exit()

        # ---- for the chosen month, load the ERA5 analysis data over the multiple 
        #      years

        temp_store = np.zeros((ndays_total,4,45,180), dtype=np.float64) # 4 levels, 45 lats, 180 lons
        shum_store = np.zeros((ndays_total,4,45,180), dtype=np.float64)
        uwnd_store = np.zeros((ndays_total,4,45,180), dtype=np.float64)
        vwnd_store = np.zeros((ndays_total,4,45,180), dtype=np.float64)
        yyyymmddhh_store = np.zeros((ndays_total), dtype=np.int32)

        ktr = 0
        for iyear, cyear in enumerate(cyears):
    
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            #print ("Loading ERA-5 reanalysis data from disk for year ",cyear,". Current time = ", current_time)
            yyyymmddhh, temp, ntimes, nlevels, ny, nx, lon, lat = load_era5(cyear, cmonth, 'air')
            yyyymmddhh, shum, ntimes, nlevels, ny, nx, lon, lat = load_era5(cyear, cmonth, 'shum')    
            yyyymmddhh, uwnd, ntimes, nlevels, ny, nx, lon, lat = load_era5(cyear, cmonth, 'uwnd')    
            yyyymmddhh, vwnd, ntimes, nlevels, ny, nx, lon, lat = load_era5(cyear, cmonth, 'vwnd') 
            zero_store = np.zeros((ntimes, nlevels, ny, nx), dtype=np.float32)
            shum = np.where(shum < 0.0, zero_store, shum)
            #print (yyyymmddhh) 
            #print ('ntimes = ', ntimes)
            temp_store[ktr:ktr+ntimes,:,:,:] = temp[:,:,0:45,:]
            shum_store[ktr:ktr+ntimes,:,:,:] = shum[:,:,0:45,:]**0.3333
            uwnd_store[ktr:ktr+ntimes,:,:,:] = uwnd[:,:,0:45,:]
            vwnd_store[ktr:ktr+ntimes,:,:,:] = vwnd[:,:,0:45,:]
            yyyymmddhh_store[ktr:ktr+ntimes] = yyyymmddhh[:]
            ktr = ktr + ntimes
    
        # --- determine the standard deviation across all time samples, and normalize by this.

        #print ('yyyymmddhh_store from ERA5 = ',yyyymmddhh_store[0:100])
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('Determining standard deviation and dividing by this.  Current time =  ', current_time)
        temp_stddev = np.std(temp_store, axis=0)
        shum_stddev = np.std(shum_store, axis=0)
        uwnd_stddev = np.std(uwnd_store, axis=0)
        vwnd_stddev = np.std(vwnd_store, axis=0)
        temp_mean = np.mean(temp_store, axis=0)
        shum_mean = np.mean(shum_store, axis=0)
        uwnd_mean = np.mean(uwnd_store, axis=0)
        vwnd_mean = np.mean(vwnd_store, axis=0)
        print ('max, min temp_mean = ', np.max(temp_mean), np.min(temp_mean))
        print ('max, min temp_stddev = ', np.max(temp_stddev), np.min(temp_stddev))

        for idate in range(ktr):
            temp_store[idate,:,:,:] = (temp_store[idate,:,:,:] - temp_mean[:,:,:]) / temp_stddev[:,:,:]
            shum_store[idate,:,:,:] = (shum_store[idate,:,:,:] - shum_mean[:,:,:]) / shum_stddev[:,:,:]
            uwnd_store[idate,:,:,:] = (uwnd_store[idate,:,:,:] - uwnd_mean[:,:,:]) / uwnd_stddev[:,:,:]
            vwnd_store[idate,:,:,:] = (vwnd_store[idate,:,:,:] - vwnd_mean[:,:,:]) / vwnd_stddev[:,:,:]
        print ('max, min temp_store after normalization = ', np.max(temp_store), np.min(temp_store))
        print ('max, min shum_store after normalization = ', np.max(shum_store), np.min(shum_store))
        print ('max, min uwnd_store after normalization = ', np.max(uwnd_store), np.min(uwnd_store))
        print ('max, min vwnd_store after normalization = ', np.max(vwnd_store), np.min(vwnd_store))
    
        # --- rearrange data; first index is time, second index is 1-d composite vector
        #     of the standardized temperature, specific humidity, and u- and v-wind 
        #     components

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ('Munging data into X, Y arrays. Current time = ', current_time)
        ngrid = nx*45*4
        X = np.zeros((ktr, ngrid*4), dtype=np.float64) # 4 as there are four fields
        Y = np.zeros((ktr, nhucs), dtype=np.float64)
        #Y = np.zeros((ktr, 1), dtype=np.float64)
        print ('ktr = ', ktr)
        for idate in range(ktr):
    
            if idate%30 == 0:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                #print ('Processing idate = ',idate,' of ',ktr,'.  Current time = ', current_time)
            temp_1d = np.reshape(temp_store[idate,:,:,:], ngrid)
            shum_1d = np.reshape(shum_store[idate,:,:,:], ngrid)
            uwnd_1d = np.reshape(uwnd_store[idate,:,:,:], ngrid)
            vwnd_1d = np.reshape(vwnd_store[idate,:,:,:], ngrid)
            #print ('idate, max, min temp_1d = ',idate,np.max(temp_1d), np.min(temp_1d))
            #print ('idate, max, min shum_1d = ',idate,np.max(shum_1d), np.min(shum_1d))
            #print ('idate, max, min uwnd_1d = ',idate,np.max(uwnd_1d), np.min(uwnd_1d))
            #print ('idate, max, min vwnd_1d = ',idate,np.max(vwnd_1d), np.min(vwnd_1d))

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
            if idate%30 == 0:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                #print ('Finding HUC data for date = ', fcst_date,'. Current time is ', current_time)
            timeindex = yyyymmddhh_hucs.index(fcst_date)
            Y[idate,:] = precipHUCs[timeindex,:]**0.3333
  
        # --- apply standard deviation to Y, HUC data  
    
        print ('max, min X = ', np.max(X), np.min(X))
        print ('max, min Y after power transform ', np.max(Y), np.min(Y))
        Ystd = np.std(Y,axis=0)
        Ymean = np.mean(Y,axis=0)
        for idate in range(ktr):
            Y[idate,:] = (Y[idate,:]-Ymean[:]) / Ystd[:]
        print ('max, min Y after normalization', np.max(Y), np.min(Y))

        # --- perform the canonical correlation analysis

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ("Performing CCA.  Current time = ", current_time)
        plsr = PLSRegression(n_components=160, scale=False)
        plsr.fit(X,Y)
        Ypred = plsr.predict(X)
        print ('max, min Y = ', np.max(Y), np.min(Y))
        print ('max, min Ypred = ', np.max(Ypred), np.min(Ypred))
        for i in range(ktr):
            print (i,Ypred[i,163], Y[i,163])
        X_c, Y_c = plsr.transform(X, Y)
        print ('shape X_c, Y_c = ', np.shape(X_c), np.shape(Y_c))
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print ("Done!  Current time = ", current_time)

        # ---- save to file

        outfile = 'cca/cca_data_month='+cmonth+'_lead='+clead+'days.cPick'
        print ('writing to ', outfile)
        ouf = open(outfile, 'wb')
        cPickle.dump(Ystd, ouf)
        cPickle.dump(Ymean, ouf)
        cPickle.dump(X_c, ouf)
        cPickle.dump(Y_c, ouf)
        cPickle.dump(Ypred, ouf)
        cPickle.dump(yyyymmddhh_store, ouf)
        ouf.close()




      