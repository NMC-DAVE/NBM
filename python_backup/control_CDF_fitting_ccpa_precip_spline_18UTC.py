# control_CDF_fitting_ccpa_precip_spline_18UTC.py 
import os, sys
from CDF_fitting_ccpa_precip_spline_flexiknot \
    import CDF_fitting_ccpa_precip_spline_flexiknot

cmonths = ['01','02','03','04','05','06','07','08','09','10','11','12']
cend_hours = ['18']

for cmonth in cmonths:
    for cend_hour in cend_hours:
        istat = CDF_fitting_ccpa_precip_spline_flexiknot(cmonth, cend_hour)
