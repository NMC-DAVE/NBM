""" create_conserved_timeseries.py """
from dateutils import daterange, datetohrs, dayofyear
import sys
import os
import os.path
import numpy as np
from os import path
import numpy.ma as ma
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset

def set_Lwater(temp): # return the Latent heat of evaporation, f(temp)
    Lwater = 2500.8 - 2.36*temp + 0.016*temp**2 - 0.00006*temp**3
    return Lwater

date_list_1999 = daterange('2000010106','2003123106', 24)
date_list_2003 = daterange('2004010106','2007123106', 24)
date_list_2007 = daterange('2008010106','2011123106', 24)
date_list_2011 = daterange('2012010106','2015123106', 24)
date_list_2015 = daterange('2016010106','2019123106', 24)
date_list = date_list_1999 + date_list_2003 + \
    date_list_2007 + date_list_2011 + date_list_2015
ndates = len(date_list)
decimalyear = ma.zeros((ndates), dtype=np.float64)
total_precip = ma.zeros((ndates), dtype=np.float64)
surface_pressure = ma.zeros((ndates), dtype=np.float64)
precipitation_rate = ma.zeros((ndates), dtype=np.float64)
evaporation_rate = ma.zeros((ndates), dtype=np.float64)
temp_2m = ma.zeros((ndates), dtype=np.float64)
hourssince1CE_2000 = datetohrs('2000010100')

for idate, date in enumerate(date_list):
    hourssince2000  = datetohrs(date_list[idate]) - hourssince1CE_2000
    decimalyear[idate] = 2000. + hourssince2000 /(24.*365.25)
    if decimalyear[idate] < 2004:
        cpath = '/Volumes//Backup Plus/bfg/1999/'
    elif decimalyear[idate] >= 2004 and decimalyear[idate] < 2008:
        cpath = '/Volumes//Backup Plus/bfg/2003/'
    elif decimalyear[idate] >= 2008 and decimalyear[idate] < 2012:
        cpath = '/Volumes//Backup Plus/bfg/2007/'
    elif decimalyear[idate] >= 2012 and decimalyear[idate] < 2016:
        cpath = '/Volumes//Backup Plus/bfg/2011/'  
    else:      
        cpath = '/Volumes//Backup Plus/bfg/2015/' 
        
    infile = cpath + 'bfg_'+date+'_fhr00_control2.nc4'
    fexist = path.exists(infile)
    #print (infile, fexist)
    if fexist == True:       
        
        try:
            nc = Dataset(infile)
            presssfc = nc.variables['pressfc'][0,:,:]
            prate_avesfc = nc.variables['prate_avesfc'][0,:,:]
            lhtfl_avesfc = nc.variables['lhtfl_avesfc'][0,:,:]
            precipitable_water_pressure = nc.variables['pwatacol'][0,:,:] *9.81
            
            temp_2m = nc.variables['tmp2m'][0,:,:] - 273.16
            if idate == 0:
                lon = nc.variables['lon'][:]
                lat = nc.variables['lat'][:]
                lon2d, lat2d = np.meshgrid(lon,lat)
                coslat = np.cos(lat2d*3.1415926/180.)
                coslatsum = np.sum(coslat)
                
            #print ('mean precipitable_water_pressure = ', \
            #    np.sum(precipitable_water_pressure*coslat)/coslatsum)
            # ---- from the 2-meter temperature, estimate the latent heat of evaporation
            #      Lwater calculated from https://en.wikipedia.org/wiki/Latent_heat
            #      as a cubic function of temperature in degrees C.

            Lwater = set_Lwater(temp_2m)
            #Lwater = 2500.*np.ones((ny,nx), dtype=np.float32) # a common approximation

            # Evaporation rate can be calculated from the latent heat flux divided by
            # the latent heat of evaporation of water (Lwater).
            # E = (latent heat flux)/Lwater , where latent heat flux
            # read in from netCDF file.
            #
            # Do the units work out?
            #
            # Numerator's units: evaporative heat flux units: W/m**2
            # 1 W = 1 J/s = 1 Nm/s = 1 kg*m**2/s**3 , so evaporative heat flux
            # units are (kg m**2/s**3)*(1/m**2) = kg/s**3
            #
            # Denominator's units: Lwater as calculated there has units of
            # J/gm = (kg*m**2/s**2)/ gm.  Hence multiply by 0.001 kg/gm
            # to get expressed in m**2/s**2 .
            #
            # So, the final units are (kg/s**3) / (m**2/s**2) =
            #   kg/s**3 * s**2/m**2 = kg/(m**2 s)
            #
            # which coincides with the precipitation rate units in the grib table
            # https://www.nco.ncep.noaa.gov/pmb/docs/on388/table2.html
            #

            evap = lhtfl_avesfc / (Lwater*1000.)
            #total_precip[idate] = np.sum(tprcpsfc*coslat) / coslatsum
            surface_pressure[idate] = (np.sum(presssfc*coslat) - \
                np.sum(precipitable_water_pressure*coslat)) / (coslatsum*100.)
            precipitation_rate[idate] = 2.0 * np.sum(prate_avesfc*coslat) / coslatsum
            evaporation_rate[idate] = np.sum(evap*coslat) / coslatsum
            nc.close() 
        except (RuntimeError, OSError) as e:
            surface_pressure[idate] = ma.masked
            precipitation_rate[idate] = ma.masked
            evaporation_rate[idate] = ma.masked
        
    else:
        #total_precip[idate] = ma.masked
        surface_pressure[idate] = ma.masked
        precipitation_rate[idate] = ma.masked
        evaporation_rate[idate] = ma.masked
    print (idate,  precipitation_rate[idate], \
        evaporation_rate[idate], surface_pressure[idate])

# ---- save to cPickle file.

outfile = 'precip_budget_and_pressure_reanalysis_timeseries.cPick'
ouf = open(outfile, 'wb')
cPickle.dump(decimalyear, ouf)
#cPickle.dump(total_precip, ouf)
cPickle.dump(surface_pressure, ouf)
cPickle.dump(precipitation_rate, ouf)
cPickle.dump(evaporation_rate, ouf)
ouf.close()
    
    