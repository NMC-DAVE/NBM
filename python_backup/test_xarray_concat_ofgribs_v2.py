
import xarray as xr
import matplotlib.pyplot as plt
from dateutils import daterange
from datetime import datetime
import dask.array as da

date_list = daterange('2000010100','2000013000',24) 
    # Jeff Whitaker function, creates a list of dates
    # including the end date, spanned by 24 h in this case
data_directory = '/Volumes/NBM/gefsv12/t2m/c00/'
data_directory_zarr_temporary = '/Volumes/NBM/gefsv12/t2m/c00/zarr_temporary/'
data_directory_zarr_chunks = '/Volumes/NBM/gefsv12/t2m/c00/zarr_chunks/'

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('before opening grib files.  time = ', current_time)


#To open multiple files simultaneously in parallel using Dask delayed, use open_mfdataset():
#xr.open_mfdataset('my/files/*.nc', parallel=True)

for idate, date in enumerate(date_list):
    infile = data_directory + 'tmp_2m_'+date+'_c00.grib2'
    print (infile)
    #ds = xr.load_dataset(infile, engine='cfgrib')
    ds = xr.open_dataset(infile, engine='cfgrib')    
    if idate == 0:
        ds_concatenated = ds
    else:
        ds_concatenated = xr.concat((ds_concatenated, ds), dim='time', \
            data_vars='minimal', coords='different', compat='equals', \
            positions=None, join='outer', combine_attrs='override')
            
            
# --- convers from xarray data structure to dask data frame

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('after grib, convert to dask data frame.  time = ', current_time)
#ds_dask = ds_concatenated.to_dask_dataframe()

# now following guidance at https://discourse.pangeo.io/t/
# best-practices-to-go-from-1000s-of-netcdf-files-to-analyses-on-a-hpc-cluster/588/11
# Tom Augspurger May 20

# --- split only pass, splitting the global array into 180x180 grids. 

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('split only pass.  time = ', current_time)
ds_dask2 = ds_concatenated.rechunk(chunks={"time":1, "step": 1, "latitude":180, "longitude":180})

# ---- write to temporary file(s)

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('writing to temporary file.  time = ', current_time)
outfile = data_directory_zarr_temporary +'temp.zarr'
da_dask2.to_zarr(outfile,overwrite=True)

# ---- read output of step 2

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('reading output of step 2.  time = ', current_time)
da_dask3.from_zarr(outfile)

# ---- merge pass.  Chunk all times together.

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('merging all times.  time = ', current_time)
da_dask4.rechunk(chunks={"time":-1})

# ---- rewrite to the final file

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('writing final file.  time = ', current_time)
outfile = data_directory_zarr_chunks +'tmp_2m_chunked.zarr'
da_dask4.to_zarr(outfile, overwrite=True)

# ---- how long did that take?

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('ending.  time = ', current_time)


print ('da_dask4 = ',da_dask4) 
print ('da_dask4.dims = ',da_dask4.dims)
print ('da_dask4.coords = ',da_dask4.coords)
print ('da_dask4.attrs = ',da_dask4.attrs)

#ds_concatenated.to_zarr(store=data_directory_zarr_chunks, \
#    chunk_store=data_directory_zarr_chunks, mode='w')
#ds_concatenated.to_zarr(store=data_directory)





#I'd like to see a time breakdown from each of these three steps

#%time ds = xr.open_mfdataset('your files *')

#%time delayed_obj = ds.to_zarr(store=store, ..., compute=False)

#%time delayed_obj.compute()


print ('done writing')
