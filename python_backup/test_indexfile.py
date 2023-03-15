
import numpy as np
import sys
import pygrib
import os

gribfilename = 'test/apcp_sfc_2000010100_c00.grib2'
idxfilename = 'test/apcp_sfc_2000010100_c00.pyidx'

validityDate = 20000111
validityTime = 600
fexist_grib = False
fexist_grib = os.path.exists(gribfilename)
fexist_idx  = False
fexist_idx  = os.path.exists(idxfilename)
print ('fexist_grib, fexist_idx  = ', fexist_grib, fexist_idx)
if fexist_grib and fexist_idx:
    try:
        #fcstfile = pygrib.open(gribfilename)
        fcstfile = pygrib.index(idxfilename)
        #print ('   q:',validityDate, validityTime)
        grb = fcstfile.select(shortName='tp',\
            validityDate=validityDate, \
            validityTime=validityTime)[0]
        gv = grb.values
        istat = 0
        fcstfile.close()
        ny, nx = np.shape(gv)
        print ('ny, nx = ', ny, nx)
        print ('sample data  = ', gv[ny//2,0:nx:10])
    except IOError:
        print ('   IOError in read_gribdata reading ', \
            gribfilename, validityDate, validityTime)
    except ValueError:
        print ('   ValueError in read_gribdata reading ', \
            gribfilename, validityDate, validityTime)
    except RuntimeError:
        print ('   RuntimeError in read_gribdata reading ', \
            gribfilename, validityDate, validityTime)

