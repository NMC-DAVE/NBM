
import xarray as xr
import matplotlib.pyplot as plt
from dateutils import daterange
from datetime import datetime
import dask.array as da
from rechunker import rechunk
import sys, os
import zarr as zarr

# ---- develop a string list of the dates to read in.

date_list = daterange('2000010100','2000013100',24) 
    # Jeff Whitaker function, creates a list of dates
    # including the end date, spanned by 24 h in this case
ndates= len(date_list)

# ---- specify directory names for input data, for working, for
#      output zarr chunks.

data_directory = '/Volumes/NBM/gefsv12/t2m/c00/'
data_directory_zarr_temporary = '/Volumes/NBM/gefsv12/t2m/c00/zarr_temporary/'
data_directory_zarr_chunks = '/Volumes/NBM/gefsv12/t2m/c00/zarr_chunks/'

# ---- note to self.  To open multiple files simultaneously 
#      in parallel using Dask delayed, use open_mfdataset():
#      xr.open_mfdataset('my/files/*.nc', parallel=True).
#      Then this loop can go away.


# --- loop over dates, read in grib files and concatenate data 
#     structures.

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('before opening grib files.  time = ', current_time)
for idate, date in enumerate(date_list):
    infile = data_directory + 'tmp_2m_'+date+'_c00.grib2'
    print (infile)
    ds = xr.open_dataset(infile, engine='cfgrib')    
    if idate == 0:
        ds_concatenated = ds
    else:
        ds_concatenated = xr.concat((ds_concatenated, ds), dim='time', \
            data_vars='minimal', coords='different', compat='equals', \
            positions=None, join='outer', combine_attrs='override')
print (ds_concatenated)

# --- note to self.   In addition to rechunking, we may want to transpose 
#     array coordinates.  see Xarray DataArray.transpose(*dims, 
#     transpose_coords=True, missing_dims='raise')[source]


# --- perform a cleanup in case we've used temporary or output 
#     directories before.

cmd = 'rm -rf /Volumes/NBM/gefsv12/t2m/c00/zarr_temporary/*'
istat = os.system(cmd)
cmd = 'rm -rf /Volumes/NBM/gefsv12/t2m/c00/zarr_chunks/*'
istat = os.system(cmd)

# --- write concatenated data to a temporary zarr file

tempfile = data_directory_zarr_temporary + 'tempfile.zarr'
print ('writing to zarr file ',tempfile)
ds_concatenated.to_zarr(tempfile)

# --- read back in.   Seems pointless.   My guess is that
#     ds_concatenated is xarray structure, ds_concatenated_zarr
#     is a zarr structure, and we need the data in this
#     structure to perform the following rechunking.

ds_concatenated_zarr = zarr.open(tempfile)
ds_concatenated_zarr.t2m.encoding = {}
print (ds_concatenated_zarr)
print (ds_concatenated_zarr.tree())

# --- extract the actual t2m data from the data structure

t2m_array = ds_concatenated_zarr['t2m']
print ('t2m_array.info',t2m_array.info)

# ---- again, clean up the temporary working directory

cmd = 'rm -rf /Volumes/NBM/gefsv12/t2m/c00/zarr_temporary/*'
istat = os.system(cmd)

# ---- in this example, I want every distinct lead time (step)
#      to be a separate chunk, and I want the globe split up
#      into tiles 180x180 grid points.   But I want all forecast
#      initial dates (time) to be in the same chunk.  The 

#target_chunks_dict = {'time':ndates, 'step': 1, 'latitude':180, 'longitude':180}
intermediate = data_directory_zarr_temporary + 'intermediate.zarr'
target = data_directory_zarr_chunks + 'tmp_2m_chunked.zarr'

#t2m_array_plan = rechunk(t2m_array, \
#    target_chunks = target_chunks_dict, max_mem=512000000,\
#    target_store=target, temp_store=intermediate)
t2m_array_plan = rechunk(t2m_array, \
    target_chunks = (ndates,1,180,180), max_mem=512000000,\
    target_store=target, temp_store=intermediate)
    
result = t2m_array_plan.execute()
print ('result.chunks = ',result.chunks)
    
# ---- how long did that take?

now = datetime.now()
current_time = now.strftime("%H:%M:%S")

print ('ending.  time = ', current_time)
