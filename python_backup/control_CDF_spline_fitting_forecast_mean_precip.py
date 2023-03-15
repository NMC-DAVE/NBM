# control_CDF_spline_fitting_forecast_mean_precip.py 

import os, sys
from CDF_spline_fitting_forecast_precip_mean import CDF_spline_fitting_forecast_precip_mean

#cmonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
#    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
cmonths = ['Apr', 'May', 'Jun', \
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
cleads = ['06','12','18','24','30','36', \
    '42','48','54','60','66', \
    '72','78','84','90','96', \
    '102','108','114','120']
    
#cmonths = ['Jan']
#cleads = ['06']
        
cdomain = 'conus'

for cmonth in cmonths:
    for clead in cleads:
        #try:
        print (cmonth, clead)
        istat = CDF_spline_fitting_forecast_precip_mean(\
            cmonth, clead, cdomain)
        #except:
        #    print ('Didnt work! ', cmonth, clead, cdomain)
