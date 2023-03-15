import sys
import csv

filename = sys.argv[1]
print ('reading ', filename)
f = open(filename, 'r')
ktr = 0
sat_dict = {ktr:{'sat':'XXX', \
                'inst':'XXX', \
                'channel':'-99',\
                'date_start':'19990101',\
                'date_end':'19990101'}}
ktr = 1
try:
    for row in f:
        sat = row.split()  
        nchannels = len(sat) - 7
        for i in range(7,7+nchannels): 
            sat_dict.update({ktr:{'sat':sat[0], 'inst':sat[5],'channel':sat[i],\
                'date_start':sat[1],'date_end':sat[3]}})
            outfile = 'satfiles/'+sat[0]+'_inst='+sat[5]+'_channel='+sat[i]+'.txt'
            print (outfile)
            f = open(outfile,'a')
            stringo = sat[1]+" "+sat[3]+"\n"
            f.write(stringo)
            f.close()
            ktr = ktr+1
            
finally:
    f.close()
