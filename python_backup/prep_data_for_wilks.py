import os, sys
import numpy as np
import numpy.ma as ma
from dateutils import daterange, dateshift
import pygrib # grib-reading routine
from write_to_f77_unformatted_f90 import write_to_f77_unformatted_f90

# =====================================================================

def find_nearest(vec, value):

    """ given a vector vec and a particular value, find the index in vec
    that is nearest to value"""

    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================
    
def get_domain_subset_of_gefs(find_near, rlon, rlat, \
    cyyyymmddhh, cmem, clead, ilon, ilat):

    """ read in the global forecast.  Subset to CONUS. """
    
    # ---- read in forecast grid covering the whole globe.
    
    cycle = cyyyymmddhh[8:10]
    input_directory = '/Projects/Hydro_Forcings/EMC/retro_forecasts/'        
    infile = input_directory + 'apcp_'+\
        cyyyymmddhh+'_'+cmem+'_.f'+clead+'.grib2'
        
    print ('   reading from ', infile)
    endStep = int(clead)
    istat, precip_realtime, lats_full, lons_full = \
        read_gribdata(infile, endStep)
    
    if find_nearest == True:
        lons_1d = lons_full[0,:]
        lats_1d = lats_full[:,0]
        ilon = find_nearest(lons_1d, rlon)   
        jlat = find_nearest(lats_1d, rlat)
        print ('nearest longitude point to ',rlon,' is ',ilon, lons_1d[ilon])
        print ('nearest latitude point to ',rlat,' is ',jlat, lats_1d[jlat])
        find_nearest = False
    
    # ---- subset for +/- 10 grid points around chosen center point
    
    precip_realtime = precip_realtime[jlat-10:jlat+11, ilon-10:ilon+11]
    lons_1D_realtime = lons_1d[ilon-10:ilon+11]-360.
    lats_1D_realtime = lats_1d[jlat-10:jlat+11]

    return precip_realtime, lons_1D_realtime, lats_1D_realtime, find_near, \
        ilon, jlat

# =====================================================================

def read_gribdata(gribfilename, endStep):
    
    """ read grib data"""
    
    istat = -1
    fexist_grib = False
    fexist_grib = os.path.exists(gribfilename)
    print ('   reading ',gribfilename, fexist_grib)
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


cyyyymm = sys.argv[1] # year/month to process from command line
clead = sys.argv[2] # 024, etc. 3 digits.

ndaysomo = [31,28,31,  30,31,30,  31,31,30, 31,30,31] 
cmembers = ['c00','p01', 'p02','p03','p04','p05','p06','p07','p08','p09','p10',\
    'p11', 'p12','p13','p14','p15','p16','p17','p18','p19','p20',\
    'p21', 'p22','p23','p24','p25','p26','p27','p28','p29','p30']
nmembers = len(cmembers)

imm = int(cyyyymm[4:6])-1
ndays = ndaysomo[imm]
cyyyymmddhh_begin = cyyyymm+'0100'
cyyyymmddhh_end = cyyyymm+str(ndays)+'00'
date_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end, 24)
ilon = 0
jlat = 0
find_near = True
rlon = -76.5 # Ithaca NY
rlat = 42.5


for idate, cyyyymmddhh in enumerate(date_list):  
    for imem, cmem in enumerate(cmembers):
        
        precip_realtime, lons_1D_realtime, lats_1D_realtime, \
            find_near, ilon, jlat = get_domain_subset_of_gefs( \
            find_near, rlon, rlat, cyyyymmddhh, cmem, clead, \
            ilon, jlat)
            
        if imem == 0:
            ny = len(lats_1d_realtime)
            nx = len(lons_1d_realtime)
            precip_realtime_ens = np.zeros((nmembers,ny,nx), dtype=np.float32)
        
        precip_realtime_ens[imem,:,:] = precip_realtime[:,:]
        
    # ---- write this set of ensemble members to a fortran unformatted file
    #      after flipping i and j indices
    
    outfile = '/data/thamill/wilks/'+clead+'/apcp_'+cyyyymmddhh+'_f'+clead+'.dat'
    istat = write_to_f77_unformatted_f90(outfile, precip_realtime_ens, \
        lons_1d_realtime, lats_1d_realtime, nmembers, ny, nx)
