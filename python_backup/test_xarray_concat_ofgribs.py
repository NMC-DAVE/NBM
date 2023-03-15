
import xarray as xr
import matplotlib.pyplot as plt
from dateutils import daterange
from datetime import datetime

date_list = daterange('2000010100','2000013000',24) 
    # Jeff Whitaker function, creates a list of dates
    # including the end date, spanned by 24 h in this case
data_directory = '/Volumes/NBM/gefsv12/t2m/c00/'
data_directory_zarr_chunks = '/Volumes/NBM/gefsv12/t2m/c00/zarr_chunks/'

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('time = ', current_time)

for idate, date in enumerate(date_list):
    infile = data_directory + 'tmp_2m_'+date+'_c00.grib2'
    #ds = xr.load_dataset(infile, engine='cfgrib')
    ds = xr.open_dataset(infile, engine='cfgrib', \
        chunks={"step": 1, "latitude":180, "longitude":180})    
    if idate == 0:
        ds_concatenated = ds
    else:
        ds_concatenated = xr.concat((ds_concatenated, ds), dim='time', \
            data_vars='minimal', coords='different', compat='equals', \
            positions=None, join='outer', combine_attrs='override')

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('time = ', current_time)
print ('ds = ',ds) 
print ('ds_concatenated.dims = ',ds_concatenated.dims)
print ('ds_concatenated.coords = ',ds_concatenated.coords)
print ('ds_concatenated.attrs = ',ds_concatenated.attrs)

ds_concatenated.to_zarr(store=data_directory_zarr_chunks, \
    chunk_store=data_directory_zarr_chunks, mode='w')
#ds_concatenated.to_zarr(store=data_directory)

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('time = ', current_time)


#I'd like to see a time breakdown from each of these three steps

#%time ds = xr.open_mfdataset('your files *')

#%time delayed_obj = ds.to_zarr(store=store, ..., compute=False)

#%time delayed_obj.compute()


print ('done writing')
