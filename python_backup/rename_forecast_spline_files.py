# rename_forecast_spline_files.py

import os, sys

cmonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
cleads = ['006','012','018','024','030','036', \
    '042','048','054','060','066', \
    '072','078','084','090','096', \
    '102','108','114','120','126', \
    '132','138','144','150','156', \
    '162','168','174','180','186', \
    '192','198','204','210','216', \
    '222','228','234','240']
    
cmonthnum = ['01','02','03','04','05','06','07','08','09','10','11','12']
    
#cmonths = ['Jan']
#cleads = ['006']
        
cdomain = 'conus'
master_directory_out = '/Volumes/NBM/'+cdomain+'_gefsv12/CDF_spline/'

for imonth in range(12):
    for clead in cleads:
        
        current_outfile = master_directory_out + cmonthnum[imonth]+'_'+cdomain+\
            '_GEFSv12_spline_info_h' + clead + '.nc'
        new_outfile = master_directory_out + cmonths[imonth]+'_'+cdomain+\
            '_GEFSv12_spline_info_h' + clead + '.nc'
        cmd = 'mv '+current_outfile+' '+new_outfile
        print (cmd)   
        
        try:
            istat = os.system(cmd)
        except:
            print ('could not find data for ', date)