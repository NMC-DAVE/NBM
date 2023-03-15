from netCDF4 import Dataset

def initialize(outfile):
    ncout = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
    sample = ncout.createDimension('sample',None)
    samplev = ncout.createVariable('samplev','i4',('sample',))
    probability_raw = ncout.createVariable('probability_raw','f4',\
        ('sample',),zlib=True,least_significant_digit=6)
    istat = 0
    return ncout

def write_data(isamp, pout, ncout):
    ncout['samplev'][isamp] = isamp
    ncout['probability_raw'][isamp] = pout
    istat = 0
    return istat
    
outfile1 = 'test1.nc'
ncout1 = initialize(outfile1)

outfile2 = 'test2.nc'
ncout2 = initialize(outfile2)

isamp = 0
pout = 0.3
istat1 = write_data(isamp, pout, ncout1)

pout = 0.7
istat2 = write_data(isamp, pout, ncout2)

ncout1.close()
ncout2.close()
