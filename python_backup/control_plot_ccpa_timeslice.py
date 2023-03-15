""" control_plot_ccpa_timeslice.py """
import os, sys
from dateutils import daterange

date_list = daterange('2019010200','2019033100',24)

for idate, date in enumerate(date_list):

    cmd = 'python plot_ccpa_timeslice.py '+ date + ' 00'
    try:
        istat = os.system(cmd)
    except:
        print ('this didnt work!')
        
    cmd = 'python plot_ccpa_timeslice.py '+ date + ' 06'
    try:
        istat = os.system(cmd)
    except:
        print ('this didnt work!')
        
    cmd = 'python plot_ccpa_timeslice.py '+ date + ' 12'
    try:
        istat = os.system(cmd)
    except:
        print ('this didnt work!')
        
    cmd = 'python plot_ccpa_timeslice.py '+ date + ' 18'
    try:
        istat = os.system(cmd)
    except:
        print ('this didnt work!')
    