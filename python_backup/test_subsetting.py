"""
test_subsetting.py cyyyymmddhh clead 

"""
import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pygrib
from quantile_mapping_gamma_mixture_v2_f90 import \
    quantile_mapping_gamma_mixture_v2_f90
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
from qmapping_spline import qmapping_spline
import _pickle as cPickle
import scipy.stats as stats


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'



# =====================================================================

# --- read grib data on a single level

def read_gribdata(gribfilename, endStep):
    istat = -1
    fexist_grib = False
    fexist_grib = os.path.exists(gribfilename)
    #print (gribfilename, endStep)
    if fexist_grib:
        try:
            fcstfile = pygrib.open(gribfilename)
            grb = fcstfile.select(shortName='tp',endStep=endStep)[0]
            precip_realtime = grb.values
            lats_full, lons_full = grb.latlons()
            istat = 0
            fcstfile.close()
        except IOError:
            print ('   IOError in read_gribdata reading ', \
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
    return istat, precip_realtime, lats_full, lons_full
    
# =====================================================================
    
def get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead):

    # ---- get the desired 2021 GEFSv12 forecast as grib file downloaded
    #      from NOMADS server. Then subset to big domain encompassing
    #      all NDFD domains, and subset again to processing domain of
    #      interest.

    nib = 886 # ~lon 220E
    nie = nib+318 # up to ~lon 310E
    njb = 131 # lat ~ 80.5
    nje = njb+153 # down to lat ~ -30.5
    
    # ---- get the whole globe.
    
    input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
    infile = input_directory + cyyyymmddhh + \
        '_ge'+cmem+'.t00z.pgrb2s.0p25.f' + clead 
    print ('  reading from ', infile)
    endStep = int(clead)
    #print ('endStep = ', endStep)
    istat, precip_realtime, lats_full, lons_full = \
        read_gribdata(infile, endStep)
    
    # ---- subset for the big grid encompassing all NBM domains.
    #      The eastern boundary crosses Greenwich meridian, so
    #      fiddle to make sure longitudes in ascending order.
    
    precip_realtime = precip_realtime[njb:nje,nib:nie]
    lons_1D_realtime = lons_full[0,nib:nie]-360.
    lats_1D_realtime = lats_full[njb:nje,0]
    print ('njb,nje, nib,nie =  ', njb,nje, nib,nie)
    print ('min, max lons_1D_realtime = ',\
        np.min(lons_1D_realtime), np.max(lons_1D_realtime))
    print ('min, max lats_1D_realtime = ',\
        np.min(lats_1D_realtime), np.max(lats_1D_realtime))
    return precip_realtime, lons_1D_realtime, lats_1D_realtime
                
# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
  
master_directory_forecast_spline = '/Volumes/NBM/'+cdomain+'_gefsv12/CDF_spline/'
   
# ---- load spline information for forecast netCDF file.  
#      Should be same lat/lon as real-time, but check

infile = master_directory_forecast_spline + \
    ccmonth + '_' + cdomain + \
    '_GEFSv12_spline_info_h' + clead + '.nc' 
nc = Dataset(infile)    
lons_spline_gefsv12_1d = nc.variables['lons'][:]
lats_spline_gefsv12_1d = nc.variables['lats'][:]
print ('min, max lons_spline_gefsv12_1d = ', \
    np.min(lons_spline_gefsv12_1d), \
    np.max(lons_spline_gefsv12_1d) )
print ('min, max lats_spline_gefsv12_1d = ', \
    np.min(lats_spline_gefsv12_1d), \
    np.max(lats_spline_gefsv12_1d) )
ny_gefsv12 = len(lats_spline_gefsv12_1d)
nx_gefsv12 = len(lons_spline_gefsv12_1d)
print ('ny_gefsv12, nx_gefsv12 = ', ny_gefsv12, nx_gefsv12)
nc.close()

cmem = 'c00'
precip_realtime, lons_1d_realtime, lats_1d_realtime = \
    get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead) 
nx_gefsv12_v2 = len(lons_1d_realtime)
ny_gefsv12_v2 = len(lats_1d_realtime)
print ('ny_gefsv12_v2, nx_gefsv12_v2 = ', ny_gefsv12_v2, nx_gefsv12_v2)



