"""

python plot_closest_member_histogram.py cmonth clead 
    where cmonth = 01,12
    clead = 3 digits, e.g., 024
    ctype = 'thinned' 'upscaled'

This script will (a) plot the closest-member histogram for the
desired month (integer digit), lead time in hours, and type
(whether data was saved for "thinned" output or "upscaled.")
It not only plots the raw closets-member histogram, but it also
performs a Savitzky-Golay smoothing of the interior values,
and it saves those interior values to a netCDF file for later
use in weighting the sorted, quantile-mapped ensemble.

This code is adapted to the GEFSv12 pre-production parallel
runs and their quantile-mapped output.   These data were
saved for 1 Dec 2017 to 30 Nov 2019.

Coded by Tom Hamill, 4 Oct 2021

"""

import numpy as np
from dateutils import daterange
import os, sys
import _pickle as cPickle
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.optimize import curve_fit
from matplotlib import rcParams
import scipy.signal as signal
import scipy.stats as stats
from matplotlib import rcParams

# ==================================================================    
    
# --- read in the desired month (01 - 12)

cmonth = sys.argv[1] # 01-12
clead = sys.argv[2] # 3 digits, e.g., 024
ctype = 'thinned2' #sys.argv[3] # thinned or upscaled
makeplot = True

imonth = int(cmonth)
master_directory_histogram_out = '/Volumes/NBM/conus_gefsv12/'+ctype+'/chist/'
master_directory_histogram_in = '/Volumes/NBM/conus_gefsv12/'+ctype+'/'


if imonth == 1:
    date_list1 = daterange('2017120100','2018022800', 24)
    date_list2 = daterange('2018120100','2019022800', 24)
    date_list3 = []
    cmonthy = 'Dec-Jan-Feb'
elif imonth == 2: 
    date_list1 = daterange('2018010100','2018033100', 24)
    date_list2 = daterange('2019010100','2019033100', 24)
    date_list3 = []
    cmonthy = 'Jan-Feb-Mar'
elif imonth == 3: 
    date_list1 = daterange('2018020100','2018043000', 24)
    date_list2 = daterange('2019020100','2019043000', 24)
    date_list3 = []
    cmonthy = 'Feb-Mar-Apr'
elif imonth == 4: 
    date_list1 = daterange('2018030100','2018053100', 24)
    date_list2 = daterange('2019030100','2019053100', 24)
    date_list3 = []
    cmonthy = 'Mar-Apr-May'
elif imonth == 5: 
    date_list1 = daterange('2018040100','2018063000', 24)
    date_list2 = daterange('2019040100','2019063000', 24)
    date_list3 = []
    cmonthy = 'Apr-May-Jun'
elif imonth == 6: 
    date_list1 = daterange('2018050100','2018073100', 24)
    date_list2 = daterange('2019050100','2019073100', 24)
    date_list3 = []
    cmonthy = 'May-Jun-Jul'
elif imonth == 7: 
    date_list1 = daterange('2018060100','2018083100', 24)
    date_list2 = daterange('2019060100','2019083100', 24)
    date_list3 = []
    cmonthy = 'Jun-Jul-Aug'
elif imonth == 8: 
    date_list1 = daterange('2018070100','2018093000', 24)
    date_list2 = daterange('2019070100','2019093000', 24)
    date_list3 = []
    cmonthy = 'Jul-Aug-Sep'
elif imonth == 9: 
    date_list1 = daterange('2018080100','2018100100', 24)
    date_list2 = daterange('2019080100','2019100100', 24)
    date_list3 = []
    cmonthy = 'Aug-Sep-Oct'
elif imonth == 10: 
    date_list1 = daterange('2018090100','2018113000', 24)
    date_list2 = daterange('2019090100','2019113000', 24)
    date_list3 = []
    cmonthy = 'Sep-Oct-Nov'
elif imonth == 11: 
    date_list1 = daterange('2017120100','2017123100', 24)
    date_list2 = daterange('2018100100','2019123100', 24)
    date_list3 = daterange('2019100100','2019113000', 24)
    cmonthy = 'Oct-Nov-Dec'
elif imonth == 12: 
    date_list1 = daterange('2017120100','2018013100', 24)
    date_list2 = daterange('2018110100','2019013100', 24)
    date_list3 = daterange('2019110100','2019113000', 24)
    cmonthy = 'Nov-Dec-Jan'


date_list = date_list1 + date_list2 + date_list3

for idate, date in enumerate(date_list):
    
    # --- read in c-m histogram for day of interest.
    
    infile = master_directory_histogram_in +\
        'closest_histogram'+date+'_'+clead+'.cPick'
    print (infile)
    fexist = False
    fexist = os.path.exists(infile)
    if fexist:
        inf = open(infile, 'rb')
        closest_histogram_input = cPickle.load(inf)
        print (np.shape(closest_histogram_input))
        inf.close()
    
    # --- set up arrays if this is the first time through
    
    if idate == 0:
        nmembers, ncats = np.shape(closest_histogram_input)
        closest_histogram = np.zeros((nmembers, ncats), dtype=np.float64)
    
    # --- add to running total
    
    if fexist:
        closest_histogram = closest_histogram + closest_histogram_input

#print ('np.shape(closest_histogram) = ', np.shape(closest_histogram))
closest_histogram_raw = np.copy(closest_histogram)
closest_histogram_cdf = np.copy(closest_histogram)
for i in range(7):
    closest_histogram[:,i] = closest_histogram[:,i]/np.sum(closest_histogram[:,i])
    for j in range(nmembers):
        closest_histogram_cdf[j,i] = np.sum(closest_histogram[0:j,i])

# ---- now make closest histogram plots

if makeplot == True:
    fig = plt.figure(figsize=(9.,6.))
    plt.suptitle('Closest-member histograms for '+cmonthy+\
        ', lead = +'+clead+' h',fontsize=18)
    rcParams['legend.fontsize']='small'
    for i in range(1,7):
    
        if i == 1:
            ctitle = '(a) 0.01 mm '+r'$\leq \bar{p}$ < 0.1 mm '
            axlocn = [0.09,0.55,0.22,0.32]
        elif i == 2:
            ctitle = '(b) 0.1 mm '+r'$\leq \bar{p}$ < 0.5 mm '
            axlocn = [0.42,0.55,0.22,0.32]
        elif i == 3:
            ctitle = '(c) 0.5 mm '+r'$\leq \bar{p}$ < 2 mm'
            axlocn = [0.75,0.55,0.22,0.32]
        elif i == 4:
            ctitle = '(d) 2 mm '+r'$\leq \bar{p} < $ 6 mm'
            axlocn = [0.09,0.09,0.22,0.32]
        elif i == 5:
            ctitle = '(e) 6 mm '+r'$\leq \bar{p} < $ 15 mm'
            axlocn = [0.42,0.09,0.22,0.32]
        else:
            ctitle = '(f) '+r'$\bar{p} \geq$ 15 mm'
            axlocn = [0.75,0.09,0.22,0.32]
        
        a1 = fig.add_axes(axlocn)
        a1.set_title(ctitle,fontsize=11)
        a1.set_xlabel('Fraction between lowest & highest rank',fontsize=8)
        a1.set_ylabel('Fraction of analyzed closest\nto this sorted member',fontsize=10)
        a1.set_xlim(0,1)
        a1.set_ylim(0.0001,0.5)
        a1.set_yscale('log')
        a1.grid(color='Gray',lw=0.2,linestyle='--')

        # --- apply Savitzky-Golay smoother to deal with sampling variability for
        #     internal ranks.
    
        csavgol = np.copy(closest_histogram[:,i])
        csavgol = signal.savgol_filter(csavgol, 9, 1, mode='mirror')
        closest_histogram[17:-17,i] = csavgol[17:-17]
    
        ctext_low_N = "{:.3f}".format(closest_histogram[0,i]/ \
            np.sum(closest_histogram[:,i]))
        ctext_high_N = "{:.3f}".format(closest_histogram[-1,i]/ \
            np.sum(closest_histogram[:,i]))    

        a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
            closest_histogram_raw[:,i]/np.sum(closest_histogram_raw[:,i]),'-',\
            color='Red',linewidth=1.5,label='Raw')
    
        a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
            2.*closest_histogram[:,i]/np.sum(closest_histogram[:,i]) ,'-',\
            color='RoyalBlue',label=r'Smoothed x 2',linewidth=2.5)
            
        a1.text(0.03,float(ctext_low_N),r'$\leftarrow\ $'+ctext_low_N)
        a1.text(0.67,float(ctext_high_N),ctext_high_N+r'$\rightarrow\ $')        
        
        a1.legend(loc=9)

    plot_title = 'closest_member_histogram_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)

# ---- rescale the closest-member histograms as weights summing to 1.0
#      before saving the data.

closest_histogram[:,0] = 1.0/float(nmembers)
for icat in range(1,ncats):
    closest_histogram[:,icat] = \
        closest_histogram[:,icat] / np.sum(closest_histogram[:,icat])

# ---- make plots of the closest-member histogram CDFs

makeplot = True
if makeplot == True:
    fig = plt.figure(figsize=(6.,6.5))

    ctitle = 'Closest-member histogram CDFs for\n'+cmonthy+\
        ', lead = +'+clead+' h'
    axlocn = [0.11,0.11,0.85,0.79]

    a1 = fig.add_axes(axlocn)
    a1.set_title(ctitle,fontsize=17)
    a1.set_xlabel('Fraction between lowest and highest rank',fontsize=14)
    a1.set_ylabel('Cumulative probability',fontsize=14)
    a1.set_xlim(0,1)
    a1.set_ylim(0,1)
    a1.grid(color='Gray',lw=0.2,linestyle='--')

    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,0],color='Red',linewidth=2,label='< 0.01')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,1],color='RoyalBlue',linewidth=2,label='0.01 to 0.1')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,2],color='LimeGreen',linewidth=2,label='0.1 to 0.5')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,3],color='Gray',linewidth=2,label='0.5 to 2.0')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,4],color='Black',linewidth=2,label='2.0 to 6.0')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,5],color='Purple',linewidth=2,label='6.0 to 15.0')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,6],color='Orange',linewidth=2,label='> 15.0')
    
    a1.legend(loc=0)

    plot_title = 'closest_member_histogram_CDF_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)


# ===================================================================
# --- save closest-member histogram to netCDF file
# ===================================================================

outfile = master_directory_histogram_out +\
    'closest_member_histogram_'+ctype+'_month='+cmonth+'_lead='+clead+'h.nc'
print ('writing to ',outfile)
nc = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

# --- initialize dimensions, variable names

#print ('ncats = ', ncats)
ncatsd = nc.createDimension('ncatsd',ncats)
ncatsv = nc.createVariable('ncatsv','i4',('ncatsd',))
ncatsv.long_name = "category numbers of histogram"
ncatsv.units = "n/a"

nm1 = ncats-1
#print ('nm1 = ', nm1)
ncats_minus = nc.createDimension('ncats_minus',nm1)
ncatsv_minus = nc.createVariable('ncatsv_minus','i4',('ncats_minus',))
ncatsv_minus.long_name = "between-categories thresholds"
ncatsv_minus.units = "n/a"

nmem = nc.createDimension('nmem',nmembers)
nmemv = nc.createVariable('nmemv','i4',('nmem'))
nmemv.long_name = "member numbers (x25 here, with stencil)"
nmemv.units = "n/a"

# --- declare the thresholds between categories

thresholds = nc.createVariable('thresholds','f4',('ncats_minus'),\
    zlib=True,least_significant_digit=6)
thresholds.units = 'mm'
thresholds.long_name = 'Thresholds between closest-member histogram categories'
thresholds.valid_range = [0.,25.]

# --- declare the closest-member histogram

closest_hist = nc.createVariable('closest_hist','f8',\
    ('nmem','ncatsd'),zlib=True,least_significant_digit=9)
closest_hist.units = 'fractional weight'
closest_hist.long_name = 'Closest-member histogram for '+\
    '25-fold ensemble and various precipitation categories'
closest_hist.valid_range = [0.,1.]


# ---- metadata

nc.title = 'closest_member histogram for this month, forecast lead '+\
    'using 5x5 stencil of points and GEFSv12 ensemble (31 members)'
nc.history = 'Tom Hamill/Diana Stovern, munging the GEFSv12 data'
nc.institution =  "NOAA Physical Sciences Lab"

# ---- initialize

ncatsv = range(ncats)
ncatsv_minus = range(ncats-1)
nmemv = range(nmembers)
thresholds[:] = [0.01, 0.1, 0.5, 2.0, 6.0, 15.0]
closest_hist[:] = closest_histogram[:,:]
print ('closing netCDF file')
nc.close()









