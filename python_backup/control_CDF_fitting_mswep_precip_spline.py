# control_CDF_fitting_mswep_precip_spline.py 
import os, sys
from CDF_fitting_mswep_precip_spline_v2 import CDF_fitting_mswep_precip_spline_v2

#cmonths = ['01','02','03','04','05','06','07','08','09','10','11','12']
cend_hours = ['00','06','12','18']

#cmonths = ['02']
#cend_hours = ['06']

cmonths = ['03']
#cend_hours = ['18']

for cmonth in cmonths:
    for cend_hour in cend_hours:
        #try:
        istat = CDF_fitting_mswep_precip_spline_v2(cmonth, cend_hour)
        #except:
        #    print ('Didnt work! ', cmonth, cend_hour)