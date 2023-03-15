""" control_CDF_spline_fitting_forecast_precip_mean_AMJ.py  

    control the spline fitting of the mean for the months April, May, June """

import os, sys
from CDF_spline_fitting_forecast_precip_mean_v4 import CDF_spline_fitting_forecast_precip_mean_v4

cmonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] 

cleads = ['006','012','018','024','030', '036','042','048','054','060', \
          '066','072','078','084','090', '096','102','108','114','120', \
          '126','132','138','144','150', '156','162','168']
        
cdomain = 'conus'

for cmonth in cmonths:
    for clead in cleads:

        print (cmonth, clead)
        istat = CDF_spline_fitting_forecast_precip_mean_v4(\
            cmonth, clead, cdomain)
