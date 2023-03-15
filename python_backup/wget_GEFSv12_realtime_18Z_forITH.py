"""
wget_python_apcp_multilead_18Z.py 
"""
from dateutils import daterange # utility of Jeff Whitaker's
import os, sys

cyyyymmdd = sys.argv[1]
cleadb = sys.argv[2]
cleade = sys.argv[3]

cyear = cyyyymmdd[0:4]

#cmembers = ['c00','p01', 'p02','p03','p04','p05','p06','p07','p08','p09','p10',\
#    'p11', 'p12','p13','p14','p15','p16','p17','p18','p19','p20',\
#    'p21', 'p22','p23','p24','p25','p26','p27','p28','p29','p30']
    
cmembers = ['c00','p01']
output_directory = '/Volumes/NBM/ITH/'
#cmembers = ['c00']

for ilead in range(int(cleadb),int(cleade)+1,3):
    if ilead < 10:
        clead = '00'+str(ilead)
    elif ilead >= 10 and ilead < 100: 
        clead = '0'+str(ilead)
    else:
        clead = str(ilead)
    print ('------- processing clead = ',clead)
    for cmem in cmembers:
        print (' ------ cmem = ', cmem)
        output_file_precip = output_directory + cyyyymmdd+'18_precip_ge'+ cmem+\
            '.t18z.pgrb2s.0p25.f'+clead
        output_file_precip_ITH = output_directory + cyyyymmdd+'18_precip_ITH_ge'+ cmem+\
            '.t18z.pgrb2s.0p25.f'+clead
        output_file_t2m = output_directory + cyyyymmdd+'18_t2m_ge'+ cmem+\
            '.t18z.pgrb2s.0p25.f'+clead
        output_file_t2m_ITH = output_directory + cyyyymmdd+'18_t2m_ITH_ge'+ cmem+\
            '.t18z.pgrb2s.0p25.f'+clead
        url = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/gens/prod/gefs.'+\
            cyyyymmdd+'/18/atmos/pgrb2sp25/ge'+ cmem+\
            '.t18z.pgrb2s.0p25.f'+clead
        print (url)
        inv_url = url+'.idx'
        print (inv_url)
        grep_text = 'grep -F ":APCP:" '
        cmd = './get_inv.pl ' + inv_url + ' | ' + grep_text + ' | ./get_grib.pl ' + url+ ' '+output_file_precip
        istat = os.system(cmd)
        print (cmd)
        grep_text = 'grep -F ":TMP:" '
        cmd = './get_inv.pl ' + inv_url + ' | ' + grep_text + ' | ./get_grib.pl ' + url+ ' '+output_file_t2m
        istat = os.system(cmd)
        print (cmd)
    
        # --- extract subset over ITH and delete original file
        
        cmd = 'wgrib2 '+output_file_precip+' -small_grib 282.5:284.5 41.5:43.5 '+output_file_precip_ITH
        print (cmd)
        istat = os.system(cmd)
        cmd = 'rm '+output_file_precip
        print (cmd)
        #istat = os.system(cmd)
        
        cmd = 'wgrib2 '+output_file_t2m+' -small_grib 282.5:284.5 41.5:43.5 '+output_file_t2m_ITH
        print (cmd)
        istat = os.system(cmd)
        cmd = 'rm '+output_file_t2m
        print (cmd)
        #istat = os.system(cmd)
        
        #try:
        #    istat = os.system(cmd)
        #except:
        #    print ('could not find data for ', cyyyymmdd, clead)
        