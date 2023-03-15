#!/usr/bin/env python
import pygrib, sys
if len(sys.argv) < 2:
    print ('grib_list <grib_filename> to get long listing')
    print ('grib_list <grib_filename> -s to get short listing')
    sys.exit(1)
fname = sys.argv[1]
short = False
if len(sys.argv) > 2 and sys.argv[2] == '-s':
    short = True
grbs = pygrib.open(fname)
if short:
    for grb in grbs:
        print (grb)
else:
    for grb in grbs:
        print ('------message %d------' %grb.messagenumber)
        for k in grb.keys():
            if k.startswith('mars'): continue
            if k == 'values' or k == 'codedValues': continue
            if grb.is_missing(k):
                print (k,'= MISSING')
            else:
                try:
                    grb[k]
                    print (k,'=',grb[k])
                except:
                    print (k,'= NOT FOUND')
grbs.close()