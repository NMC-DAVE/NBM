""" control_brier_reliability_all.py """
import os, sys

types = ['probability_raw', 'probability_qmapped', \
    'probability_qmapped_weighted', \
    'probability_qmapped_weighted_dressed']

#cleads = ['012','024','036','048','060','072','084','096','108','120',\
#    '132','144','156','168','180','192','204','216','228','240'] # '036',
#cleads = ['024','072','120']
cleads = ['108','120','132','240']
#cleads = ['120']
    
#cyyyymms = ['201712']
#cyyyymms = ['201712','201801','201802','201803','201804',\
#    '201805','201806','201807','201808','201809',\
#    '201810','201811','201812','201901','201902',\
#    '201903','201904','201905','201906','201907',\
#    '201908','201909','201910','201911']
cyyyymms = ['201712','201801','201802','201803',\
    '201810','201811','201812','201901','201902',\
    '201903','201910','201911']
    
#cyyyymms = ['201810']

for cyyyymm in cyyyymms:        
    for clead in cleads:  
        for ctype in types:

            #cmd = 'python calculate_brier_reliability.py '+ \
            #    cyyyymm + ' ' + clead + ' ' + ctype
            cmd = 'python calculate_brier_reliability_any.py '+ \
                cyyyymm + ' ' + clead + ' ' + ctype
            print (cmd)

            try:
                istat = os.system(cmd)
            except:
                print ('this didnt work!')
                        
