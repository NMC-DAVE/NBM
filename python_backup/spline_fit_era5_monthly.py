""" spline_fit_era5_monthly.py
"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import _pickle as cPickle
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
import math


cyears = ['2000', '2001', '2002', '2003', '2004', 
    '2005', '2006', '2007', '2008', '2009',\
    '2010', '2011', '2012', '2013', '2014',\
    '2015', '2016', '2017', '2018', '2019']

cmmddhhs = ['010100', '020100', '030100', '040100', '050100',  '060100', \
    '070100',  '080100', '090100',  '100100', '110100',  '120100']

# ---- read in 20 years of monthly ERA5 t2m at 00 UTC, and average this.

for imonth, cmmddhh in enumerate(cmmddhhs):
    for cyear in cyears:

        infilename = 'ecmwf/t2m_climo_era5_halfdegree_'+cyear+cmmddhh+'.cPick'
        print (infilename)
        inf = open(infilename, 'rb')
        temp_halfdegree = cPickle.load(inf)
        if cyear == cyears[0] and imonth == 0:
            ny, nx = np.shape(temp_halfdegree)
            temp_mean = np.zeros((12,ny,nx), dtype=np.float32)
            latsa_halfdegree = cPickle.load(inf)
            lonsa_halfdegree = cPickle.load(inf)
            inf.close()
        print (temp_halfdegree[ny//2, nx//2])
        temp_mean[imonth,:,:] = temp_mean[imonth,:,:] + temp_halfdegree[:,:]
    temp_mean[imonth,:,:] = temp_mean[imonth,:,:] / len(cyears)
    
# ---- spline fit to the day of the year.
    
knots = [0.1*3.14159*2., 0.2*3.14159*2.,0.3*3.14159*2., 0.4*3.14159*2.,\
    0.5*3.14159*2.,0.6*3.14159*2,0.7*3.14159*2., 0.8*3.14159*2., 0.9*3.14159*2.]

# ---- will repeat the ERA5 time series 3x over for spline fit; get the associated
#      radian related to day of year.
    
radians = np.zeros(36, dtype=np.float32)
rktr = 0
for iyear in range(3):
    for imonth in range(12):
        fracyear = (float(imonth)+0.5) / 12.
        radians[rktr] = 2.*3.1415926*(iyear-1.+fracyear)
        rktr = rktr+1

# ---- loop over months, days, and spline fit

temp_mean_dayofyear = np.zeros((365, ny, nx), dtype=np.float32)
x = np.arange(0.,2.*math.pi,2.*math.pi/365.)
for ix in range(nx):
    print ('processing ix = ', ix)
    for jy in range(ny):
        tm = temp_mean[:,jy,ix] 
        temp_mean_3x = np.tile(tm,3)

        spltemp = splrep(radians, temp_mean_3x, \
            xb=0., xe=2.*math.pi, k=3, task=-1, per=1, t=knots)
        temp_mean_dayofyear[:,jy,ix] = splev(x, spltemp)
        if ix == nx//2 and jy == ny//2:
            print (radians)
            print (temp_mean_3x)
            print ('temp_mean_dayofyear[0:365:5] = ', temp_mean_dayofyear[0:365:5, jy, ix])
    
# ---- save to cPickle file.
        
outfilename = 'ecmwf/t2m_climo_daily_era5_halfdegree.cPick'
print (outfilename)
ouf = open(outfilename, 'wb')
cPickle.dump(temp_mean_dayofyear, ouf)
cPickle.dump(latsa_halfdegree, ouf)
cPickle.dump(lonsa_halfdegree, ouf)
ouf.close()    

