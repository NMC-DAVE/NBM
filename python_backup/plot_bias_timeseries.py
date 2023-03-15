"""
plot_bias_timeseries.py

"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
# --------------------------------------------------------------   


def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# ---- various initialization

clead = sys.argv[1]
cmonth = sys.argv[2] 
clon = sys.argv[3]
clat = sys.argv[4]
rlon = float(clon)
rlat = float(clat)
iskip = int(clead)//24
cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_beta = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_random = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_gain = '/Volumes/Backup Plus/gefsv12/t2m/'

if cmonth == 'Jan':
    cmonth3 = 'DJF'
elif cmonth == 'Feb':
    cmonth3 = 'JFM'
elif cmonth == 'Mar':
    cmonth3 = 'FMA'
elif cmonth == 'Apr':
    cmonth3 = 'MAM'
elif cmonth == 'May':
    cmonth3 = 'AMJ'
elif cmonth == 'Jun':
    cmonth3 = 'MJJ'
elif cmonth == 'Jul':
    cmonth3 = 'JJA'
elif cmonth == 'Aug':
    cmonth3 = 'JAS'
elif cmonth == 'Sep':
    cmonth3 = 'ASO'
elif cmonth == 'Oct':
    cmonth3 = 'SON'
elif cmonth == 'Nov':
    cmonth3 = 'OND'
elif cmonth == 'Dec':
    cmonth3 = 'NDJ'


# - get lat/lon indice

infile = cpath_era5 +'2001/t2m_era5_halfdegree_2001010100.cPick'
inf = open(infile, 'rb')
analysis = cPickle.load(inf) - 273.16
analysis = np.flipud(analysis)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
nlats, nlons = np.shape(lats)
lats = np.flipud(lats)
lons_1d = lons[0,:]
lats_1d = lats[:,0]
imin = find_nearest(lons_1d, rlon)
jmin = find_nearest(lats_1d, rlat)
inf.close()

# ---- load the bias correction time series

bias_file = 'bias_correction_2000_2018_convolved_lead'+clead+'.cPick'
inf = open(bias_file, 'rb')
bias_3d = cPickle.load(inf)
print ('jmin, imin = ', jmin, imin)
#for i in range(20):
#    ioff = 365*i
#    print (bias_3d[0+ioff:90+ioff:5,jmin,imin])
#sys.exit()
date_list_anal = cPickle.load(inf)
inf.close()

# ---- load the random data time series

random_file = 'random_error_2000_2018_lead'+clead+'.cPick'
inf = open(random_file, 'rb')
random_3d = cPickle.load(inf)
inf.close()

ndates, nlats, nlons = np.shape(random_3d)
npts = nlats*nlons
    
# ---- loop thru months

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']


for imonth, cmonth in enumerate(cmonths):
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print ('**** PROCESSING month = ', cmonth, current_time)
    
    # ---- loop through dates and process day by day

    datevalid_indices = np.zeros(ndates,dtype=np.int32)
    
    for idate, date in enumerate(date_list_anal):
        thismonth = int(str(date)[4:6]) - 1

        if imonth == 0:
            if thismonth == imonth or thismonth == 11 or \
                thismonth == imonth+1: datevalid_indices[idate] = 1
        elif imonth == 11:
            if thismonth == imonth or thismonth == 0 or \
                thismonth == imonth-1: datevalid_indices[idate] = 1
        else:
            if thismonth == imonth or thismonth == imonth+1 or \
                thismonth == imonth-1: datevalid_indices[idate] = 1
                
    ndates_valid = np.sum(datevalid_indices)
    bias_byyear = np.zeros((19,90), dtype=np.float32)



    ktrdatev = 0
    iyyyy_old = 1999
    for idate, date in enumerate(date_list_anal):
        if datevalid_indices[idate] == 1:
            d = date_list_anal[idate]
            iyyyy = int(d[0:4])
            print (idate, date, d, iyyyy, iyyyy_old)
            if iyyyy != iyyyy_old:
                iyyyy_old = iyyyy
                iy = iyyyy - 2000
                ktr = 0
            if ktr < 90: bias_byyear[iy,ktr] = bias_3d[idate,jmin,imin]
            ktr = ktr+1
            
    # ---- Q-Q plot of forecast data

    f = plt.figure(figsize=(6.5,5.))

    ax = f.add_axes([.15,.15,.79,.75])
    ax.set_title('Bias time series'+cmonth3+',\n lon = '+clon+' lat = '+clat+' lead = '+clead+'h',fontsize=13)
    for iyear in range(19):
        ax.plot(range(90), bias_byyear[iyear,:],'k-',lw=2)
    
    ax.set_ylim(-6,6)
    ax.grid(True,lw=0.25,color='LightGray')
    ax.set_xlim(0,90)
    ax.set_xlabel('Day number in the 3-month series',fontsize=11)
    ax.set_ylabel('Bias estimate (deg C)',fontsize=11)

    figname = 'bias_tseries_'+cmonth3+'_lon='+clon+'_lat='+clat+'_lead='+clead+'.pdf'
    plt.savefig(figname,dpi=400)
    print ('Plot done', figname)
    plt.close()
    

    
    
    