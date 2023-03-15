""" control_precip_forecast_2netcdf.py """
import os

cmonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

cleads = ['006','012','018','024','030', '036','042','048','054','060', \
          '066','072','078','084','090', '096','102','108','114','120', \
          '126','132','138','144','150', '156','162','168','174','180', \
          '186','192','198','204','210', '216','222','228','234','240']
          
for cmonth in cmonths:
    for clead in cleads:
        cmd = 'python precip_forecast_2netcdf.py '+clead+' '+ cmonth
        istat = os.system(cmd)