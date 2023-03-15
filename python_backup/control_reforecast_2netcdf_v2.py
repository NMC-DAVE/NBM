""" 
control_reforecast_2netcdf_v2.py

control the execution of the production of netCDF files for 
precipitation from grib files """

from reforecast_2netcdf_v2 import reforecast_2netcdf_v2
from dateutils import daterange, dateshift, dayofyear, splitdate
import sys

clead = sys.argv[1] # Enter lead time in hours.   Enter 06 for 6
   # Do every 3 hours between 03 and 240.
ilead = int(clead)
cdayend_noleap = ['31','28','31','30','31','30',  '31','31','30','31','30','31']
cdayend_leap = ['31','29','31','30','31','30',  '31','31','30','31','30','31']
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
outfiledir = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'

date_list = []
for imonth in range(1,2):

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
        print ('iyear, date_list_thisyear ', iyear, date_list_thisyear )
        date_list.extend(date_list_thisyear)
    
    outfilename = outfiledir + cmonths + '+'_apcp_h'+clead+'.nc'
    print ('   start, end dates: ',date_list[0], date_list[-1])
    print ('   writing to ', outfilename)
    istat = reforecast_2netcdf_v2(date_list, ilead, outfilename)