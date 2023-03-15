"""
wget_python.py -- a simple script to extract 20 years of GEFSv12
    reforecast data of precipitation, all members
"""
from dateutils import daterange # utility of Jeff Whitaker's
import os, sys

for iyear in range (2000, 2020):
    cyear = str(iyear)
    cyyyymmddhh_begin = cyear+'120100'
    cyyyymmddhh_end = cyear+'123100'
    
    #if iyear%4 == 0:
    #    cyyyymmddhh_begin = cyear+'020100'
    #    cyyyymmddhh_end = cyear+'022900'
    #else:
    #    cyyyymmddhh_begin = cyear+'020100'
    #    cyyyymmddhh_end = cyear+'022800'
        
    date_list = daterange(cyyyymmddhh_begin, cyyyymmddhh_end,24)
    for idate, date in enumerate(date_list):
        cyear = date[0:4]
        for imem, cmem in enumerate(['c00','p01','p02','p03','p04']):
    
            url = 'https://noaa-gefs-retrospective.s3.amazonaws.com/GEFSv12/reforecast/'+\
                cyear + '/'+date+ '/'+cmem+'/Days:1-10/apcp_sfc_'+date+'_'+cmem+'.grib2'
            print (url)
            cmd = 'wget '+url
            #print (cmd)
            try:
                istat = os.system(cmd)
            except:
                print ('could not find data for ', date)
        