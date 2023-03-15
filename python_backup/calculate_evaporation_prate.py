"""
calculate_evaporation_prate.py

read in global time series of analyzed precipitation rate,
latent heat flux, and surface temperature.   From LHF and
temperature, calculate the evaporative flux.   Get average
over the globe and write out time series.
""" 


from netCDF4 import Dataset
import numpy as np
from dateutils import daterange
import sys
import os
import os.path
from os import path
import numpy.ma as ma
import _pickle as cPickle

def set_Lwater(temp): # return the Latent heat of evaporation, f(temp)
    Lwater = 2500.8 - 2.36*temp + 0.016*temp**2 - 0.00006*temp**3
    return Lwater

# ----- develop a table of the latent heat of condensation, -25 to 40 by 0.1

cstreams = ['1999','2003','2007','2011','2015']
date_list_1999 = daterange('2000010100','2003123100',24)
date_list_2003 = daterange('2004010100','2007123100',24)
date_list_2007 = daterange('2008010100','2011123100',24)
date_list_2011 = daterange('2012010100','2015123100',24)
date_list_2015 = daterange('2016010100','2019123100',24)
ndates_total = len(date_list_1999) + len(date_list_2003) + \
    len(date_list_2007) + len(date_list_2011) + len(date_list_2015)
prate_timeseries = np.zeros((ndates_total), dtype=np.float64)
evap_timeseries  = np.zeros((ndates_total), dtype=np.float64)
print (ndates_total)

# --- read in sample lat/lon indices; set up cos(latitude) grid

infile = '/Users/Tom/python/gefsv12/1999/bfg_2003123100_fhr00_control2.nc4'
nc = Dataset(infile)
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]
nlons = len(lon)
nlats = len(lat)
cosfac_field = np.zeros((nlats,nlons), dtype=np.float64)
nc.close()
for i in range(nlats):
    cosfac_field[i,:] = np.cos(lat[i]*3.1415926/180.)

# --- loop over streams

ktr = 0
for istream, cstream in enumerate(cstreams):

    print (istream, cstream )
    # ---- determine the dates

    if cstream == '1999':
        date_list = date_list_1999
    elif cstream == '2003':
        date_list = date_list_2003
    elif cstream == '2007':
        date_list = date_list_2007
    elif cstream == '2011':
        date_list = date_list_2011
    else:
        date_list = date_list_2015
                    
    print ('***** processing stream = ', cstream)
    # ---- loop thru each date in the date_list
    
    for idate, date in enumerate(date_list):
        
        # --- read in from netCDF file
        
        infile1 = '/Users/Tom/python/gefsv12/'+cstream+'/bfg_'+date+'_fhr00_control2.nc4'
        fexist = path.exists(infile1)
        if fexist == True:
            
            nc = Dataset(infile1)
            #print (infile1)
            latent = nc.variables['lhtfl_avesfc'][0,:,:]
            prate = nc.variables['prate_avesfc'][0,:,:]
            temperature = nc.variables['tmp2m'][0,:,:] - 273.15 # convert to deg C
            ny, nx = np.shape(temperature)
            #print ('max, min temperature = ', np.max(temperature), np.min(temperature))
        
            # ---- from the 2-meter temperature, estimate the latent heat of evaporation
            #      Lwater calculated from https://en.wikipedia.org/wiki/Latent_heat
            #      as a cubic function of temperature in degrees C.
        
            Lwater = set_Lwater(temperature)
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
        
            evap = latent / (Lwater*1000.)
            nc.close()
            evap_timeseries[ktr] = np.sum(evap*cosfac_field) / np.sum(cosfac_field)
            #if istream == 0:
            prate_timeseries[ktr] = 2.0 * np.sum(prate*cosfac_field) / np.sum(cosfac_field) 
            #else:
            #    prate_timeseries[ktr] = np.sum(prate*cosfac_field) / np.sum(cosfac_field) 
            print (idate, ktr, date, evap_timeseries[ktr], prate_timeseries[ktr], \
                np.sum(Lwater*cosfac_field) / np.sum(cosfac_field)) 
                
            if evap_timeseries[ktr] > 0.1: evap_timeseries[idate-1] # bad data point filter
            if prate_timeseries[ktr] > 0.1: prate_timeseries[idate-1] 
                
        else:
            
            evap_timeseries[ktr] = -99.99
            prate_timeseries[ktr] = -99.99
        
        ktr = ktr+1
        
# ---- save to file.

outfile = 'evap_prate.cPick'
print ('writing to ', outfile)
ouf = open(outfile,'wb')
cPickle.dump(evap_timeseries, ouf)
cPickle.dump(prate_timeseries, ouf)
ouf.close()