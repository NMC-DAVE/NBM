""" control_brier_reliability.py """
import os, sys

#types = ['thinned', 'upscaled']
types = ['thinned']

#cleads = ['006','012','018','024','030','036', \
#    '042','048','054','060','066','072','078',\
#    '084','090','096','102','108','114','120',\
#    '126','132','138','144','150','156','162',\
#    '168','174','180','186','192','198','204',\
#    '210','216','222','228','234','240']
#cleads = ['024','072','120']
cleads = ['024']
    
#cyyyymms = ['201810']
cyyyymms = ['201712','201801','201802','201803','201804',\
    '201805','201806','201807','201808','201809',\
    '201810','201811','201812','201901','201902',\
    '201903','201904','201905','201906','201907',\
    '201908','201909','201910','201911']

    
for ctype in types:
    for clead in cleads:
        for cyyyymm in cyyyymms:
            #cmd = 'python calculate_brier_reliability.py '+ \
            #    cyyyymm + ' ' + clead + ' ' + ctype
            cmd = 'python calculate_brier_reliability_weighted.py '+ \
                cyyyymm + ' ' + clead + ' ' + ctype
            print (cmd)

            try:
                istat = os.system(cmd)
            except:
                print ('this didnt work!')
                        
