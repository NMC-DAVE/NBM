"""
determine_raw_forecast_climo.py

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

# --------------------------------------------------------------   

# ---- various initialization

clead = sys.argv[1]
cvariable = '2t'
#cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_forecast = '/Volumes/Backup Plus/gefsv12/t2m/'
date_list_anal = daterange('2000010100','2018123100',24)
#date_list_anal = daterange('2000010112','2018123112',24)
ndates = len(date_list_anal)
knots = [0.1*3.14159*2., 0.2*3.14159*2.,0.3*3.14159*2., \
    0.4*3.14159*2.,0.5*3.14159*2.,0.6*3.14159*2,\
    0.7*3.14159*2., 0.8*3.14159*2., 0.9*3.14159*2.]

# ---- loop through dates and process day by day

ndatektr = 0
ndatektr_yearly = np.zeros((18), dtype=np.int32)
radians_doy = np.zeros((ndates), dtype=np.float32)
for idate, date in enumerate(date_list_anal):
    
    # ---- read reanalysis appropriate at the time of 
    #      this forecast for bias corr.
    
    rem = idate%30
    if rem == 0: print ('processing date = ', idate, date)
    cyear = date[0:4]
    cmm = date[4:6]
    cmmdd = date[4:8]
    imm = int(cmm)
    idd = int(date[6:8])
    iyear = int(cyear)-2000
    iyear_full = int(cyear)
    julday = dayofyear(iyear_full, imm, idd)
    radians_doy[idate] = np.where(iyear_full%4 == 0,  \
        2.*math.pi*float(julday)/365., \
        2.*math.pi*float(julday)/364.) +\
        (iyear_full-2000)*0.00001 # extra term to aid in sorting
    
    #infile = cpath_era5 +cyearf+'/t2m_era5_halfdegree_'+date+'.cPick'
    infile = cpath_forecast + cyear + '/'+date+'_lead'+\
        clead+'_conus_0.5deg_hour'+clead+'.cPick'
    fexist1 = os.path.exists(infile)
    #print (infile, fexist1)
    if fexist1 == True:
        inf = open(infile, 'rb')
        forecast = cPickle.load(inf) 
        if idate == 0: 
            nlats, nlons = np.shape(forecast)
            forecast_3d = ma.zeros((ndates,nlats,nlons), dtype=np.float64)
            climo_temps_estimated = np.zeros((365,nlats,nlons))
            climo_temps_stddev = np.zeros((365,nlats,nlons))
        forecast_3d[idate,:,:] = forecast[:,:]
    else:
        forecast_3d[idate,:,:].mask = True
        print ('missing data for idate, date = ',idate, date)
        
        
# ---- estimate mean temperature for each julian day with cubic spline

x = np.arange(0.,2.*math.pi,2.*math.pi/365.)
print ('shape of x = ', len(x))
for jy in range(nlats):
    for ix in range(nlons):
        temps = forecast_3d[:,jy,ix]
        rads = radians_doy[:]
        indices = np.argsort(rads)
        temps_sorted = temps[indices]
        rads_sorted = rads[indices]
        spltemp = splrep(rads_sorted, temps_sorted, \
            xb=0., xe=2.*math.pi, k=3, task=-1, per=1, t=knots)
        climo_temps_estimated[:,jy,ix] = splev(x, spltemp)            
            
# ---- for every Julian day and grid point, determine the 
#      standard deviation around the estimated mean.  Use +/- 30 days

for idate, date in enumerate(date_list_anal):
    cyear = date[0:4]
    cmm = date[4:6]
    cmmdd = date[4:8]
    imm = int(cmm)
    idd = int(date[6:8])
    iyear = int(cyear)-2000
    iyear_full = int(cyear)
    cyearf = date[0:4]
    julday = dayofyear(iyear_full, imm, idd) - 1
    if julday > 364: julday = 364
    forecast_3d[idate,:,:] = forecast_3d[idate,:,:] - \
        climo_temps_estimated[julday,:,:]

iuse = np.zeros(ndates, dtype=np.int32)
sumx = np.zeros((nlats, nlons), dtype = np.float64)
sumx2 = np.zeros((nlats, nlons), dtype = np.float64)
for iday in range(365):
    iuse = np.zeros((ndates,nlats,nlons), dtype=np.int32)
    for idate, date in enumerate(date_list_anal):
        cyear = date[0:4]
        cmm = date[4:6]
        cmmdd = date[4:8]
        imm = int(cmm)
        idd = int(date[6:8])
        iyear = int(cyear)-2000
        iyear_full = int(cyear)
        cyearf = date[0:4]
        julday = dayofyear(iyear_full, imm, idd) 
        if julday > 364: julday = 364
        imin = np.min([np.abs(iday - julday), \
            np.abs(iday - julday + 365), np.abs(iday - julday -365)])
        if imin < 30:
            iuse[idate,:,:] = 1

    sumx = np.sum(forecast_3d*iuse, axis=0) 
    sumx2 = np.sum(forecast_3d*forecast_3d*iuse, axis=0)
    nd = np.sum(iuse[:,0,0])
    climo_temps_stddev[iday,:,:] = np.sqrt((sumx2 - sumx**2/nd)/(nd-1))
    print (iday, climo_temps_stddev[iday,nlats//2,nlons//2])
     
# ---- save estimated 2000-2018 climatology to file

#outfile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'
outfile = cpath_forecast + 'GEFSv12_temperature_climatology_lead='+\
    clead+'h.cPick'
print ('writing to ', outfile)
ouf = open(outfile,'wb')
cPickle.dump(climo_temps_estimated, ouf)
cPickle.dump(climo_temps_stddev, ouf)
ouf.close()
                