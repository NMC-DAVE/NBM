"""

python ccpa_to_netCDF_halfdegree.py 

For chosen month (01 to 12), extract merged CCPA and MSWEP
and upscale to 0.5 degrees and save to a new netCDF file.

Tom Hamill, Jan 2022

"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import scipy.stats as stats
import pygrib
from netCDF4 import Dataset
from dateutils import hrs_since_day1CE_todate, \
    dateto_hrs_since_day1CE, hrstodate, datetohrs, dateshift
from mpl_toolkits.basemap import Basemap, interp
from upscale_precip_f90 import upscale_precip_f90

# ========================================================================

# ---- define the nearest CCPA grid to each 1/2-degree point

def define_nearest_ndfd_point(lons_halfdegree, lats_halfdegree, \
    ny_half, nx_half, lons_ccpa, lats_ccpa, ny_ccpa, nx_ccpa) :
    
    
    infile = 'nearest_NDFD_to_halfdegree.cPick'
    fexist = os.path.exists(infile)

    if fexist == True:
        inf = open(infile,'rb')
        jnear = cPickle.load(inf)
        inear = cPickle.load(inf)
        #print ('jnear[:, nx_half//2] = ',jnear[:, nx_half//2])
        #print ('inear[ny_half//2,:] = ',inear[ny_half//2,:])
        #sys.exit()
        dmin = cPickle.load(inf)
        inf.close()
    else:
        print ('2. did not find file ', infile)
    
        jnear = np.zeros((ny_half, nx_half), dtype=int)
        inear = np.zeros((ny_half, nx_half), dtype=int)
        dmin = np.zeros((ny_half, nx_half), dtype=np.float32)
        for ilon, rlon in enumerate(lons_halfdegree):
        #for ilon, rlon in zip([137],lons_halfdegree[137:138]):
            #print ('processing ilon, rlon = ',ilon, rlon)
            distx = np.abs(lons_ccpa - rlon)
            for jlat, rlat in enumerate(lats_halfdegree):
            #for jlat, rlat in zip([39],lats_halfdegree[39:40]):
                #print ('processing jlat, rlat = ',jlat, rlat)
                #print ('rlon, rlat = ', rlon, rlat)
                disty = np.abs(lats_ccpa - rlat)
                dist = distx**2 + disty**2
                #print ('np.shape(dist)', np.shape(dist))
                #print ('dist(ny_ccpa//2, 0:nx_ccpa:10) = ', dist[ny_ccpa//2, 0:nx_ccpa:10])
                j,i = np.unravel_index(dist.argmin(), dist.shape)
                jnear[jlat,ilon] = j 
                inear[jlat,ilon] = i
                #print ('unraveled indices j,i, ny_ccpa, nx_ccpa = ', j,i, ny_ccpa, nx_ccpa)
                dmin[jlat,ilon] = dist[j,i]
                #print ('dmin[jlat,ilon] = ', dmin[jlat,ilon])
                #print ('jlat, ilon, lons_halfdegree[ilon],'+\
                #    ' lats_halfdegree[jlat], lons_ccpa[j,i], lats_ccpa[j,i] = ',\
                #    jlat,ilon,lons_halfdegree[ilon], lats_halfdegree[jlat], \
                #    lons_ccpa[j,i], lats_ccpa[j,i])
                #print ('dist[j-10:j+10,i-10:nx_ccpa] = ',dist[j-10:j+10,i-10:nx_ccpa] )
                #print (j,i,ny_ccpa, nx_ccpa)
        #sys.exit()
        
        outfile = 'nearest_NDFD_to_halfdegree.cPick'
        ouf = open(outfile, 'wb')
        cPickle.dump(jnear, ouf)
        cPickle.dump(inear, ouf)
        cPickle.dump(dmin, ouf)
        ouf.close()
        
    return jnear, inear, dmin

# ========================================================================

# ---- get averaging box min and max

def define_averaging_box(jnear, inear, dmin, ny_half, nx_half, ny_ccpa, nx_ccpa):
    jboxmin = -99*np.ones((ny_half, nx_half), dtype=int)
    jboxmax = -99*np.ones((ny_half, nx_half), dtype=int)
    iboxmin = -99*np.ones((ny_half, nx_half), dtype=int)
    iboxmax = -99*np.ones((ny_half, nx_half), dtype=int)
    for j in range(ny_half):
        jminus = np.max([j-1,0])
        jplus = np.min([ny_half-1,j+1])
        for i in range(nx_half):
            iminus = np.max([i-1,0])
            iplus = np.min([nx_half-1,i+1])
            if dmin[j,i] < 0.001:
            
                iboxmin[j,i] = max([(inear[j,iminus]+inear[j,i])//2,0])
                iboxmax[j,i] = min([(inear[j,iplus]+inear[j,i])//2,nx_ccpa-1])
                jboxmin[j,i] = max([(jnear[jplus,i]+jnear[j,i])//2,0])
                jboxmax[j,i] = min([(jnear[jminus,i]+jnear[j,i])//2,ny_ccpa-1])
                
                if iboxmin[j,i] < 0 or iboxmax[j,i] > nx_ccpa-1 or \
                jboxmin[j,i] < 0 or jboxmax[j,i] > ny_ccpa-1:
                    iboxmin[j,i] = -99
                    iboxmax[j,i] = -99
                    jboxmin[j,i] = -99
                    jboxmax[j,i] = -99
            
    #print ('jboxmin[:,nx_half//2] = ',jboxmin[:,nx_half//2])
    #print ('jboxmax[:,nx_half//2] = ',jboxmax[:,nx_half//2])
    #print ('iboxmin[ny_half//2,:] = ',iboxmin[ny_half//2,:])
    #print ('iboxmax[ny_half//2,:] = ',iboxmax[ny_half//2,:])  
    #sys.exit()              

    return jboxmin, jboxmax, iboxmin, iboxmax
                

# ========================================================================

def initialize_netCDF_halfdegree(outfile, lats_halfdegree, \
    lons_halfdegree, nx_half, ny_half, ndates):

    # ---- open netCDF output file and deal with all the variable definition and
    #      such.
    
    ncout = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

    xf = ncout.createDimension('xf',nx_half)
    xvf = ncout.createVariable('xf','i4',('xf',))
    xvf.long_name = "eastward grid point number on 1/2-degree grid"
    xvf.units = "n/a"

    yf = ncout.createDimension('yf',ny_half)
    yvf = ncout.createVariable('yf','i4',('yf',))
    yvf.long_name = "northward grid point number on 1/2-degree grid"
    yvf.units = "n/a"

    sample = ncout.createDimension('sample',None)
    samplev = ncout.createVariable('samplev','i4',('sample',))
    samplev.units = "index to time dimension, that's all"

    lonsa = ncout.createVariable('lons','f4',('xf',))
    lonsa.long_name = "longitude"
    lonsa.units = "degrees_east"

    latsa = ncout.createVariable('lats','f4',('yf',))
    latsa.long_name = "latitude"
    latsa.units = "degrees_north"

    yyyymmddhh_begin = ncout.createVariable('yyyymmddhh_begin','i4',('sample',))
    yyyymmddhh_begin.longname = \
        "Precip accumulation period beginning in yyyymmddhh format"

    #print (yyyymmddhh_begin)
    yyyymmddhh_end = ncout.createVariable('yyyymmddhh_end','i4',('sample',))
    yyyymmddhh_end.longname = \
        "Precip accumulation period ending in yyyymmddhh format"
    #print (yyyymmddhh_end)

    # --- declare the single-level variable information on lat-lon grid

    apcp_analysis = ncout.createVariable('apcp_analysis','f4',('sample','yf','xf',),
        zlib=True,least_significant_digit=4)
    apcp_analysis.units = "mm"
    apcp_analysis.long_name = \
        "Upscaled 6-h accumulated combined CCPA/MSWEP analysis on CONUS 1/2-degree grid"
    apcp_analysis.valid_range = [0.,1000.]
    apcp_analysis.missing_value = np.array(-99.99,dtype=np.float32)

    # ---- initialize

    xvf[:] = range(nx_half)
    yvf[:] = range(ny_half)
    lonsa[:] = lons_halfdegree[:]
    latsa[:] = lats_halfdegree[:]

    # ---- metadata

    ncout.title = "Combined CCPA/MSWEP 1/2-degree precip analysis upscaled from NDFD grid synthesis"
    ncout.history = "Interpolated CCPA provided by Yan Luo, NCEP/EMC"
    ncout.institution =  "NCEP/EMC"
    ncout.platform = "Precipitation analysis"
    ncout.references = "n/a"
    
    istat = 0
    ncout.close()
    return istat
    
    
# ========================================================================
    
# ---- define 1/2-degree lats and lons

lats_halfdegree = np.arange(57., 19.49, -0.5)
lons_halfdegree = np.arange(-138.0, -59.4, 0.5)
nx_half = len(lons_halfdegree)
ny_half = len(lats_halfdegree)

# ---- get the month and end time from the commmand line.  The first 00
#      hour analysis of the month will need to access the data from
#      the previous month.

ccpa_directory = '/Volumes/NBM/conus_panal/'

#for iyear in range(2002,2020):
for iyear in range(2020,2021):
    print ('processing year = ', iyear)
    #for imonth, cmonth in enumerate(\
    #    ['01','02','03','04','05','06','07','08','09','10','11','12']):
    for imonth, cmonth in enumerate(['01']):
    #for imonth, cmonth in enumerate(['01','02']):
        
        print ('  processing month = ',imonth)
        cyear = str(iyear)
        infile = ccpa_directory + cyear + cmonth + '_ccpa_on_ndfd_grid_6hourly.nc'
        #print ('  reading ', infile)
        ncin = Dataset(infile)

        apcp_anal = ncin.variables['apcp_anal'][:,:,:]
        yyyymmddhh_begin_in = ncin.variables['yyyymmddhh_begin'][:]
        yyyymmddhh_end_in = ncin.variables['yyyymmddhh_end'][:]
        ndates = len(yyyymmddhh_begin_in)
        
        outfile = ccpa_directory + cyear + cmonth + '_ccpa_on_halfdegree_grid_6hourly.nc'
        #print ('  initializing netCDF file.')
        istat = initialize_netCDF_halfdegree(outfile, lats_halfdegree, \
            lons_halfdegree, nx_half, ny_half, ndates)
        
        # ---- first pass through, define averaging box
        
        if iyear == 2020 and imonth == 0:
            lons_ccpa = ncin.variables['lons'][:,:]
            lats_ccpa = ncin.variables['lats'][:,:]
            ny_ccpa, nx_ccpa = np.shape(lats_ccpa)
            conusmask = ncin.variables['conusmask'][:,:]
            jnear, inear, dmin = define_nearest_ndfd_point(\
                lons_halfdegree, lats_halfdegree, ny_half, \
                nx_half, lons_ccpa, lats_ccpa, ny_ccpa, nx_ccpa)
            jboxmin, jboxmax, iboxmin, iboxmax = \
                define_averaging_box(jnear, inear, dmin, \
                ny_half, nx_half, ny_ccpa, nx_ccpa)

        ncin.close()
        
        apcp_anal_upscaled_today = -99.99*np.ones((ny_half,nx_half), dtype=np.float32)

        ktr = 0
        print ('  writing to ',outfile)
        ncout = Dataset(outfile,'a',format='NETCDF4_CLASSIC')
        for idate in range(ndates):
            
            #print ('processing date = ', idate, ' of ', ndates)
            
            # ---- upscale the data
            
            apcp_anal_today = apcp_anal[idate,:,:]
            #print ('idate, apcp_anal_today[jnear[39,137], inear[39,137]]= ', \
            #    idate, apcp_anal_today[jnear[39,137], inear[39,137]])
            apcp_anal_upscaled_today = upscale_precip_f90(apcp_anal_today,jboxmin, \
                jboxmax, iboxmin, iboxmax, ny_half, nx_half, ny_ccpa, nx_ccpa)
            #print ('  idate, max apcp_anal_today, upscaled = ', \
            #    idate, np.max(apcp_anal_today), np.max(apcp_anal_upscaled_today))
            #print ('apcp_anal_upscaled_today[39,137] = ', \
            #    apcp_anal_upscaled_today[39,137])
            #sys.exit()

            # ---- write netCDF record to file.

            ncout['yyyymmddhh_begin'][ktr] = yyyymmddhh_begin_in[ktr]
            ncout['yyyymmddhh_end'][ktr] = yyyymmddhh_end_in[ktr]
            ncout['apcp_analysis'][ktr] = apcp_anal_upscaled_today[:,:]
            ktr = ktr + 1

        print ('  closing ', outfile)
        ncout.close()






