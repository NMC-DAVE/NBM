""" 
control_reforecast_2netcdf_240h.py

control the execution of the production of netCDF files for 
precipitation from grib files """

from reforecast_2netcdf_240h import reforecast_2netcdf_240h
from dateutils import daterange, dateshift, dayofyear, splitdate
import sys

cdayend_noleap = ['31','28','31','30','31','30',  '31','31','30','31','30','31']
cdayend_leap = ['31','29','31','30','31','30',  '31','31','30','31','30','31']
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#outfiledir = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'
outfiledir = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'

date_list = []
for imonth in range(1,13):  #range(2,3):
    for iyear in range(2000, 2020):
        cyear = str(iyear)
        if imonth < 10:
            cmonth = '0'+str(imonth)
        else:
            cmonth = str(imonth)
        cdaystart = '01'
        if iyear%4 == 0:
            cdayend = cdayend_leap[imonth-1]
        else:
            cdayend = cdayend_noleap[imonth-1]
        cyyyymmddhh_start = cyear+cmonth+'0100'
        cyyyymmddhh_end = cyear+cmonth+cdayend+'00'
        date_list_thisyear = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)
        date_list.extend(date_list_thisyear)
    
    outfilename = outfiledir + cmonths[imonth-1] +'_apcp_mean_10day.nc'
    print ('   start, end dates: ',date_list[0], date_list[-1])
    print ('   writing to ', outfilename)
    istat = reforecast_2netcdf_240h(date_list, outfilename)
    
