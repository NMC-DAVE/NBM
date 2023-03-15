from netCDF4 import Dataset

def test_nc(rootgrp):
    print(rootgrp.data_model)
    istat = 0
    return istat

rootgrp1 = Dataset("test1.nc", "w", format="NETCDF4")
rootgrp2 = Dataset("test2.nc", "w", format="NETCDF4")

istat1 = test_nc(rootgrp1)
istat2 = test_nc(rootgrp2)