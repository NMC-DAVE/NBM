"""
determine_era5_climatology.py

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

cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
#date_list_anal = daterange('2000010100','2018123100',24)
date_list_anal = daterange('2000010112','2018123112',24)
ndates = len(date_list_anal)
knots = [0.1*3.14159*2., 0.2*3.14159*2.,0.3*3.14159*2., 0.4*3.14159*2.,\
    0.5*3.14159*2.,0.6*3.14159*2,0.7*3.14159*2., 0.8*3.14159*2., 0.9*3.14159*2.]

# ---- loop through dates and process day by day

ndatektr = 0
ndatektr_yearly = np.zeros((18), dtype=np.int32)
radians_doy = np.zeros((ndates), dtype=np.float32)
for idate, date in enumerate(date_list_anal):
    
    # ---- read reanalysis appropriate at the time of this forecast for bias corr.
    
    rem = idate%30
    if rem == 0: print ('processing date = ', idate, date)
    cyear = date[0:4]
    cmm = date[4:6]
    cmmdd = date[4:8]
    imm = int(cmm)
    idd = int(date[6:8])
    iyear = int(cyear)-2000
    iyear_full = int(cyear)
    cyearf = date[0:4]
    julday = dayofyear(iyear_full, imm, idd)
    radians_doy[idate] = np.where(iyear_full%4 == 0,  \
        2.*math.pi*float(julday)/365., \
        2.*math.pi*float(julday)/364.) +\
        (iyear_full-2000)*0.00001 # extra term to aid in sorting
    
    infile = cpath_era5 +cyearf+'/t2m_era5_halfdegree_'+date+'.cPick'
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
            analyses_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            climo_temps_estimated = np.zeros((365,nlats,nlons))
            climo_temps_stddev = np.zeros((365,nlats,nlons))
    analyses_3d[idate,:,:] = analysis[:,:]

x = np.arange(0.,2.*math.pi,2.*math.pi/365.)
print ('shape of x = ', len(x))
for jy in range(nlats):
    for ix in range(nlons):
        temps = analyses_3d[:,jy,ix]
        rads = radians_doy[:]
        indices = np.argsort(rads)
        temps_sorted = temps[indices]
        rads_sorted = rads[indices]
        spltemp = splrep(rads_sorted, temps_sorted, \
            xb=0., xe=2.*math.pi, k=3, task=-1, per=1, t=knots)
        climo_temps_estimated[:,jy,ix] = splev(x, spltemp)
        #if jy == nlats//2 and ix == nlons//2:
        #    print ('climo_temps sample = ',climo_temps_estimated[0:-1:5,jy,ix])
            
            
# ---- for every day and grid point, determine the standard deviation around
#      the estimated mean.

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
    analyses_3d[idate,:,:] = analyses_3d[idate,:,:] - climo_temps_estimated[julday,:,:]

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
        imin = np.min([np.abs(iday - julday), np.abs(iday - julday + 365), np.abs(iday - julday -365)])
        if imin < 30:
            iuse[idate,:,:] = 1
        #print (np.abs(iday - julday), np.abs(iday - julday + 365), np.abs(iday - julday -365))
    #print (iday, date, julday, iuse[0:40,0,0], iuse[-40:-1,0,0])
    #sys.exit()

    sumx = np.sum(analyses_3d*iuse, axis=0) 
    sumx2 = np.sum(analyses_3d*analyses_3d*iuse, axis=0)
    nd = np.sum(iuse[:,0,0])
    climo_temps_stddev[iday,:,:] = np.sqrt((sumx2 - sumx**2/nd)/(nd-1))
    print (iday, climo_temps_stddev[iday,nlats//2,nlons//2])
     
# ---- save estimated 2000-2018 climatology to file

#outfile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'
outfile = cpath_era5 + 'ERA5_temperature_climatology_12UTC.cPick'
print ('writing to ', outfile)
ouf = open(outfile,'wb')
cPickle.dump(climo_temps_estimated, ouf)
cPickle.dump(lats, ouf)
cPickle.dump(lons, ouf)
cPickle.dump(climo_temps_stddev, ouf)
ouf.close()
                