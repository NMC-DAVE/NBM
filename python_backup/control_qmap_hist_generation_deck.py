""" control_qmap_hist_generation_deck.py """
import os, sys

cseason = sys.argv[1] # 'cool', 'warm'
clead = sys.argv[2] # '024' etc, 3-digit number
    
if cseason == 'cool':
    cyyyymms = ['201712','201801','201802','201803',\
        '201810','201811','201812','201901','201902',\
        '201903','201910','201911']
else:
    cyyyymms = ['201804','201805','201806','201807',\
        '201808','201809','201904','201905','201906',\
        '201907','201908','201909']

outfile = 'qmap_chist_'+cseason+'_'+clead+'.deck'
print ('writing to ', outfile)
ouf = open(outfile,'w')
ctext = '#!/bin/bash'    
print (ctext, file=ouf)

for cyyyymm in cyyyymms:   
    cyyyymmddhh = cyyyymm + '0100' 
    cmd = 'python qmap_reforecast_for_histogram_stats_thinned.py '+ \
        cyyyymmddhh + ' ' + clead
    print (cmd,file=ouf)
ouf.close()
                        
