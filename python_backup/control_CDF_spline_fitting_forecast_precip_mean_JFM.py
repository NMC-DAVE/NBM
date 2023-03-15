# control_CDF_spline_fitting_forecast_precip.py 
import os, sys
from CDF_spline_fitting_forecast_precip_mean_v4 import CDF_spline_fitting_forecast_precip_mean_v4

#cmonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
#    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'] 
    
cmonths = ['Jan', 'Feb', 'Mar'] 

cleads = ['006','012','018','024','030', '036','042','048','054','060', \
          '066','072','078','084','090', '096','102','108','114','120', \
          '126','132','138','144','150', '156','162','168']
    
#cmonths = ['Jan']
#cleads = ['024']
#cleads = ['018']
        
cdomain = 'conus'

for cmonth in cmonths:
    for clead in cleads:
        #try:
        print (cmonth, clead)
        istat = CDF_spline_fitting_forecast_precip_mean_v4(\
            cmonth, clead, cdomain)
        #except:
        #    print ('Didnt work! ', cmonth, clead, cdomain)
