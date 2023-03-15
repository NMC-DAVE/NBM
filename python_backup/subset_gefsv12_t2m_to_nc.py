"""
subset_gefsv12_t2m_to_nc.py

for input lead time, develop a netCDF file of the time series of 
t2m reforecasts for this lead time.

coded by: Tom Hamill, Oct 2021, tom.hamill@noaa.gov 

"""

import os, sys
from datetime import datetime # Jeff Whitaker's datetime library
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from dateutils import daterange, dateshift
import pygrib # grib-reading routine
import scipy.stats as stats

# =====================================================================        
        
def read_gribdata(gribfilename, endStep):
    istat = -1
    fexist_grib = False
    fexist_grib = os.path.exists(gribfilename)
    print (gribfilename, endStep)
    t2m_input = []
    lats_full = []
    lons_full = []
    if fexist_grib:
        try:
            fcstfile = pygrib.open(gribfilename)
            grb = fcstfile.select(shortName='2t',endStep=endStep)[0]
            t2m_input = grb.values
            lats_full, lons_full = grb.latlons()
            istat = 0
            fcstfile.close()
        except IOError:
            print ('   IOError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1
        except NameError:
            print ('   NameError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1    
        except ValueError:
            print ('   ValueError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1
        except RuntimeError:
            print ('   RuntimeError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1
        except UnboundLocalError:
            print ('   UnboundLocalError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1
    return istat, t2m_input, lats_full, lons_full
       
# =====================================================================    

def get_domain_subset_of_gefs(cyyyymmddhh, master_directory_in, \
    cmem, clead):

    """ read in the global forecast.  Subset to CONUS. """

    nib = 886  # precomputed for set of lat-lons that surrounds
    nie = nib+318 # the NDFD CONUS 2.5-km Lambert conformal domain.
    njb = 131 
    nje = njb+153 
    ny = 153
    nx = 318
    
    # ---- read in forecast grid covering the whole globe.
    
    infile = master_directory_in + cmem + '/tmp_2m_' + cyyyymmddhh + \
        '_'+cmem+'.grib2' 
    print ('   reading from ', infile)
    endStep = int(clead)
    istat, t2m_input, lats_full, lons_full = \
        read_gribdata(infile, endStep)
    
    # ---- subset for CONUS.
    
    if istat == 0:
        t2m = t2m_input[njb:nje,nib:nie]
        lons_1D = lons_full[0,nib:nie]-360.
        lats_1D = lats_full[njb:nje,0]
    else:
        t2m = -99.99*np.ones((ny,nx), dtype=np.float32)
        lats_1D = -99.99*np.ones((ny), dtype=np.float32)
        lons_1D = -99.99*np.ones((nx), dtype=np.float32)

    return t2m, lons_1D, lats_1D, ny, nx
    
# =====================================================================    

def initialize_netCDF(master_directory_ncout, clead, cyear, \
    nx, ny, nmembers, lons_in, lats_in):
                    
    """ initialize the output netCDF file """

    outfile = master_directory_ncout + cyear +'_t2m_conus_'+clead+'h.nc'
    print ('initializing netCDF file ',outfile)
    ncout = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
    
    # --- initialize dimensions, variable names
        
    xf = ncout.createDimension('xf',nx)
    xvf = ncout.createVariable('xf','f4',('xf',))
    xvf.long_name = "eastward grid point number on 1/4-degree lat-lon grid"
    xvf.units = "n/a"

    yf = ncout.createDimension('yf',ny)
    yvf = ncout.createVariable('yf','f4',('yf',))
    yvf.long_name = "northward grid point number on 1/4-degree lat-lon grid"
    yvf.units = "n/a"
    
    enssize = ncout.createDimension('enssize',nmembers)
    memberv = ncout.createVariable('memberv','f4',('enssize',))
    memberv.long_name = "member number (control = 0)"
    memberv.units = "n/a"

    sample = ncout.createDimension('sample',None)
    samplev = ncout.createVariable('samplev','i4',('sample',))
    samplev.units = "n/a"
    
    lons_out = ncout.createVariable('lons','f4',('xf',))
    lons_out.long_name = "longitude"
    lons_out.units = "degrees_east"

    lats_out = ncout.createVariable('lats','f4',('yf',))
    lats_out.long_name = "latitude"
    lats_out.units = "degrees_north"

    yyyymmddhh_init = ncout.createVariable('yyyymmddhh_init','i4',('sample',))
    yyyymmddhh_init.longname = "Initial condition date/time in yyyymmddhh format"

    yyyymmddhh_fcst = ncout.createVariable('yyyymmddhh_fcst','i4',('sample',))
    yyyymmddhh_fcst.longname = "Forecast valid date/time in yyyymmddhh format"

    # --- declare the single-level variable information on lat-lon grid

    t2m = ncout.createVariable('t2m','f4',\
        ('sample','enssize','yf','xf',),zlib=True,least_significant_digit=3)
    t2m.units = "n/a"
    t2m.long_name = '2-meter temperature in deg C.'
    t2m.valid_range = [-50., 45.]
    t2m.missing_value = np.array(-99.99,dtype=np.float32)
    
    # ---- metadata

    ncout.title = 'GEFSv12 reforecast t2m forecast data '+\
        ' for this lead on 0.25-deg grid surrounding the CONUS'
    ncout.history = "GEFSv12 implemented at NCEP/EMC Sep 2020"
    ncout.institution =  "NCEP/EMC and PSL"
    ncout.platform = ""
    ncout.references = ""

    # ---- initialize

    xvf[:] = np.arange(nx)
    yvf[:] = np.arange(ny)
    memberv[:] = range(nmembers)
    print ('max, min lons_in = ', np.max(lons_in), np.min(lons_in))
    lons_out[:] = lons_in[:]
    lats_out[:] = lats_in[:]
    
    return ncout

# =====================================================================    
    
def write_t2m_to_netcdf (isamp, cyyyymmddhh, \
    cyyymmddhh_fcst, t2m_out, ncout):    
    
    """ write the t2m record to netCDF file """

    ncout['yyyymmddhh_init'][isamp] = int(cyyyymmddhh)
    ncout['yyyymmddhh_fcst'][isamp] = int(cyyyymmddhh_fcst)
    ncout['t2m'][isamp] = t2m_out[:,:,:] - 273.16
    istat = 0
    return istat
    
# =====================================================================    

clead = sys.argv[1]  # 3 digits

master_directory_in = '/Volumes/NBM/gefsv12/t2m/'
master_directory_ncout = '/Volumes/NBM/gefsv12/t2m/conus_netCDF/'

isamp = 0
cmembers = ['c00','p01','p02','p03','p04']
nmembers = len(cmembers)

#for iyear in range(2006,2007):
for iyear in range(2000,2019):
    cyear = str(iyear)
    date_begin = cyear+'010100'
    date_end = cyear+'123100'
    #date_begin = cyear+'040100'
    #date_end = cyear+'040200'
    date_list = daterange(date_begin, date_end,24)
    for idate, cyyyymmddhh_init in enumerate(date_list):
        cyyyymmddhh_fcst = dateshift(cyyyymmddhh_init, int(clead))
    
        for imem, cmem in enumerate(cmembers):
        
            # --- read in grib file t2m data for this date and lead time
        
            t2m, lons_1D, lats_1D, ny, nx =  \
                get_domain_subset_of_gefs(cyyyymmddhh_init, \
                master_directory_in, cmem, clead)

            # --- first time through?  Initialize output netCDF file.
        
            if idate == 0 and imem == 0:
                t2m_store = np.zeros((nmembers,ny,nx), dtype=np.float32)
                ncout = initialize_netCDF(master_directory_ncout, clead, \
                    cyear, nx, ny, nmembers, lons_1D, lats_1D)
                
            t2m_store[imem,:,:] = t2m[:,:]
        
        # --- write this netCDF record
        
        istat = write_t2m_to_netcdf (idate, cyyyymmddhh_init, \
            cyyyymmddhh_fcst, t2m_store, ncout)
    
    ncout.close() 

        