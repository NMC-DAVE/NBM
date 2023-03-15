# test_control_CDF_spline_fitting_forecast_precip.py 
import os, sys
from CDF_spline_fitting_forecast_precip_v4 import CDF_spline_fitting_forecast_precip_v4

#cmonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
#    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] 
    
cmonths = ['Mar'] 
cleads = ['024']
cdomain = 'conus'

for cmonth in cmonths:
    for clead in cleads:
        #try:
        print (cmonth, clead)
        istat = CDF_spline_fitting_forecast_precip_v4(\
            cmonth, clead, cdomain)
        #except:
        #    print ('Didnt work! ', cmonth, clead, cdomain)
