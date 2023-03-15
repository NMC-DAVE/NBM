"""
python thin_and_upscale_ccpa.py cyyyymm

In order to produce and validate the thinned (every 10th grid pt)
and upscaled GEFSv12 forecasts, raw and quantile mapped, 
we will need thinned and upscaled the NDFD grid merged 
CCPA/MSWEP precipitation analyses. This routine creates these.

Tom Hamill

"""

import os, sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
from PIL import Image

# =====================================================================

def initialize_netCDF_outfile(outfile, ny, nx, \
    lats_in, lons_in, conusmask_in):
        
    """ initialize the netCDF files for writing probability output"""

    nc = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
    print ('initializing ', outfile)
        
    xf = nc.createDimension('xf',nx)
    xvf = nc.createVariable('xf','f4',('xf',))
    xvf.long_name = "eastward grid point number"
    xvf.units = "n/a"

    yf = nc.createDimension('yf',ny)
    yvf = nc.createVariable('yf','f4',('yf',))
    yvf.long_name = "northward grid point number"
    yvf.units = "n/a"
    
    sample = nc.createDimension('sample',None)
    samplev = nc.createVariable('samplev','i4',('sample',))
    samplev.units = "n/a"

    lons_out = nc.createVariable('lons','f4',('yf','xf',))
    lons_out.long_name = "longitude"
    lons_out.units = "degrees_east"

    lats_out = nc.createVariable('lats','f4',('yf','xf',))
    lats_out.long_name = "latitude"
    lats_out.units = "degrees_north"
    
    conusmask_out = nc.createVariable('conusmask','f4',('yf','xf',))
    conusmask_out.long_name = "land water fraction"
    conusmask_out.units = "n/a"

    yyyymmddhh_begin = nc.createVariable('yyyymmddhh_begin','i4',('sample',))
    yyyymmddhh_begin.longname = \
        "Precip analysis period beginning date/time in yyyymmddhh format"

    yyyymmddhh_end = nc.createVariable('yyyymmddhh_end','i4',('sample',))
    yyyymmddhh_end.longname = \
        "Precip analysis period ending date/time in yyyymmddhh format"

    # --- declare the single-level variable information on lat-lon grid

    apcp_anal = nc.createVariable('apcp_anal','f4',\
        ('sample','yf','xf',),zlib=True,least_significant_digit=3)
    apcp_anal.units = "n/a"
    apcp_anal.long_name = 'analyzed precipitation amount (mm)'
    apcp_anal.valid_range = [0.,300.0]
    apcp_anal.missing_value = np.array(-99.99,dtype=np.float32)
    
    # ---- metadata

    nc.title = 'NDFD CONUS precipitation analyses from CCPA/MSWEP combo, upscaled or thinned'
    nc.history = " "
    nc.institution =  "NCEP/EMC and PSL primarily, MSWEP from Princeton"
    nc.platform = ""
    nc.references = ""

    # ---- initialize

    xvf[:] = np.arange(nx)
    yvf[:] = np.arange(ny)
    lons_out[:] = lons_in[:,:]
    lats_out[:] = lats_in[:,:]
    conusmask_out[:] = conusmask_in[:]
    istat = 0
    
    return istat, nc
    
# =====================================================================

# ---- get the month and end time from the commmand line.  The first 00
#      hour analysis of the month will need to access the data from
#      the previous month.

cyyyymm = sys.argv[1] # 202001 etc
nstride = 10  # thin to every 10th grid point 
master_directory = '/Volumes/NBM/conus_panal/'

infile = master_directory + cyyyymm + \
    '_ccpa_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile)
yyyymmddhh_begin_in = nc.variables['yyyymmddhh_begin'][:]
yyyymmddhh_end_in = nc.variables['yyyymmddhh_end'][:]
conusmask_in = nc.variables['conusmask'][:,:]
precip_in = np.squeeze(nc.variables['apcp_anal'][:,:,:])
ndates, ny_ndfd, nx_ndfd = np.shape(precip_in)
lons_in = nc.variables['lons'][:,:]
lats_in = nc.variables['lats'][:,:]
nc.close()

# ---- upscale the precipitation analyses as per 
#      quantile_map_ensemble_bymonth_andlead_v4.py

for idate in range(ndates):
    
    # ---- initialize stuff first time through
    
    if idate == 0:
        
        im = Image.fromarray(lons_in)
        imu = im.resize(size=(ny_ndfd//10, nx_ndfd//10),resample=Image.BOX)
        lons_upscaled = np.transpose(np.asarray(imu))
        
        im = Image.fromarray(lats_in)
        imu = im.resize(size=(ny_ndfd//10, nx_ndfd//10),resample=Image.BOX)
        lats_upscaled = np.transpose(np.asarray(imu))
        
        ny_upscaled, nx_upscaled = np.shape(lons_upscaled)
        
        im = Image.fromarray(conusmask_in)
        imu = im.resize(size=(ny_ndfd//10, nx_ndfd//10),resample=Image.BOX)
        conusmask_upscaled = np.transpose(np.asarray(imu))
            
        precip_upscaled = np.zeros((ny_upscaled, nx_upscaled), \
            dtype=np.float32)

        outfile_upscaled = infile = master_directory + cyyyymm + \
            '_ccpa_on_ndfd_grid_6hourly_upscaled.nc'
        istat, nc = initialize_netCDF_outfile(outfile_upscaled, \
            ny_upscaled, nx_upscaled, lats_upscaled, \
            lons_upscaled, conusmask_upscaled)

    # --- upscale each day's precipitation
    
    p = precip_in[idate,:,:]
    im = Image.fromarray(p)
    imu = im.resize(size=(ny_ndfd//10, nx_ndfd//10),resample=Image.BOX)
    precip_upscaled[:,:] = np.transpose(np.asarray(imu))
    
    # --- write precipitation analysis record and date stuff
    
    nc['yyyymmddhh_begin'][idate] = yyyymmddhh_begin_in[idate]
    nc['yyyymmddhh_end'][idate] = yyyymmddhh_end_in[idate]
    nc['apcp_anal'][idate] = precip_upscaled[:,:]
    

nc.close()
print ('writing to ', outfile_upscaled,' completed.')

# ---- thin the precipitation analyses 
    
for idate in range(ndates):
    
    if idate == 0: 
        lons_thinned = lons_in[5::nstride, 5::nstride]
        lats_thinned = lats_in[5::nstride, 5::nstride] 
        conusmask_thinned = conusmask_in[5::nstride, 5::nstride] 
        ny_thinned, nx_thinned = np.shape(lons_thinned)          
        precip_thinned = np.zeros((ny_thinned, nx_thinned), \
            dtype=np.float32)
        outfile_thinned = master_directory + cyyyymm + \
            '_ccpa_on_ndfd_grid_6hourly_thinned.nc'
        istat, nc = initialize_netCDF_outfile(outfile_thinned, \
            ny_thinned, nx_thinned, lats_thinned, \
            lons_thinned, conusmask_thinned)

    # --- thin the precip array     
        
    p = precip_in[idate,:,:]
    precip_thinned = p[5::nstride, 5::nstride]

    # --- write precipitation analysis record and date stuff
    
    nc['yyyymmddhh_begin'][idate] = yyyymmddhh_begin_in[idate]
    nc['yyyymmddhh_end'][idate] = yyyymmddhh_end_in[idate]
    nc['apcp_anal'][idate] = precip_thinned[:,:]
    

nc.close()
print ('writing to ', outfile_thinned,' completed.')

