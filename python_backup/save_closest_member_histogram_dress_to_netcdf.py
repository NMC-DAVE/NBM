"""

python save_closest_member_histogram_dress_to_netcdf.py cmonth clead 
    where cmonth = 01,12
    clead = 3 digits, e.g., 024
    

This script will generate closest-member histogram and 
fitted closest-member dressing statistics for the
desired month (integer digit), lead time in hours.
It not only plots the raw closest-member histogram, but it also
performs a Savitzky-Golay smoothing of the interior values,
and it saves those interior values to a netCDF file for later
use in weighting the sorted, quantile-mapped ensemble.

This code is adapted to the GEFSv12 pre-production parallel
runs and their quantile-mapped output.   These data were
saved for 1 Dec 2017 to 30 Nov 2019.

Coded by Tom Hamill, 17 Dec 2021

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
from sklearn.linear_model import LinearRegression

# ==================================================================    
    
def populate_datelists(imonth):
    """ for the integer month, 1-12, get a list of the retro dates and 
    a character string to use in possible plots """
    
    if imonth == 1:
        date_list1 = daterange('2017120100','2018022800', 24)
        date_list2 = daterange('2018120100','2019022800', 24)
        date_list3 = []
        cmonthy = 'January'
    elif imonth == 2: 
        date_list1 = daterange('2018010100','2018033100', 24)
        date_list2 = daterange('2019010100','2019033100', 24)
        date_list3 = []
        cmonthy = 'February'
    elif imonth == 3: 
        date_list1 = daterange('2018020100','2018043000', 24)
        date_list2 = daterange('2019020100','2019043000', 24)
        date_list3 = []
        cmonthy = 'March'
    elif imonth == 4: 
        date_list1 = daterange('2018030100','2018053100', 24)
        date_list2 = daterange('2019030100','2019053100', 24)
        date_list3 = []
        cmonthy = 'April'
    elif imonth == 5: 
        date_list1 = daterange('2018040100','2018063000', 24)
        date_list2 = daterange('2019040100','2019063000', 24)
        date_list3 = []
        cmonthy = 'May'
    elif imonth == 6: 
        date_list1 = daterange('2018050100','2018073100', 24)
        date_list2 = daterange('2019050100','2019073100', 24)
        date_list3 = []
        cmonthy = 'June'
    elif imonth == 7: 
        date_list1 = daterange('2018060100','2018083100', 24)
        date_list2 = daterange('2019060100','2019083100', 24)
        date_list3 = []
        cmonthy = 'July'
    elif imonth == 8: 
        date_list1 = daterange('2018070100','2018093000', 24)
        date_list2 = daterange('2019070100','2019093000', 24)
        date_list3 = []
        cmonthy = 'August'
    elif imonth == 9: 
        date_list1 = daterange('2018080100','2018100100', 24)
        date_list2 = daterange('2019080100','2019100100', 24)
        date_list3 = []
        cmonthy = 'September'
    elif imonth == 10: 
        date_list1 = daterange('2018090100','2018113000', 24)
        date_list2 = daterange('2019090100','2019113000', 24)
        date_list3 = []
        cmonthy = 'October'
    elif imonth == 11: 
        date_list1 = daterange('2017120100','2017123100', 24)
        date_list2 = daterange('2018100100','2019123100', 24)
        date_list3 = daterange('2019100100','2019113000', 24)
        cmonthy = 'November'
    elif imonth == 12: 
        date_list1 = daterange('2017120100','2018013100', 24)
        date_list2 = daterange('2018110100','2019013100', 24)
        date_list3 = daterange('2019110100','2019113000', 24)
        cmonthy = 'December'
    date_list = date_list1 + date_list2 + date_list3
    return date_list, cmonthy

# ==================================================================    

def read_closest_member_stats_and_dress_stats(date_list, \
    master_directory_histogram_in, clead):
    
    """ from the daily file output of closest-member statistics
        and dressing statistics, tally up information to determine
        them over multiple days and months. """

    for idate, date in enumerate(date_list):
  
        # --- read in c-m histogram and stats to determine
        #     dresssing mean and std dev for day of interest.
    
        cyyyymm = date[0:6]
        infile = master_directory_histogram_in + cyyyymm + '/' +\
            'closest_histogram'+date+'_'+clead+'.cPick' 
        fexist = False
        fexist = os.path.exists(infile)
        print (infile, fexist)
        if fexist:
            inf = open(infile, 'rb')
            closest_histogram_input = cPickle.load(inf)
            sumxi_low = cPickle.load(inf)
            sumxi2_low = cPickle.load(inf)
            nsamps_low = cPickle.load(inf)
            sumxi_mid = cPickle.load(inf)
            sumxi2_mid = cPickle.load(inf)
            nsamps_mid = cPickle.load(inf)        
            sumxi_high = cPickle.load(inf)
            sumxi2_high = cPickle.load(inf)
            nsamps_high = cPickle.load(inf)
            #print (np.shape(closest_histogram_input))
            inf.close()
    
            # --- set up arrays if this is the first time through
    
            if idate == 0:
                nmembers, ncats = np.shape(closest_histogram_input)
                closest_histogram = np.copy(closest_histogram_input) 
                    #np.zeros((nmembers, ncats), dtype=np.float64)
                sumxi_low_allcases = np.copy(sumxi_low)
                sumxi2_low_allcases = np.copy(sumxi2_low)
                nsamps_low_allcases = np.copy(nsamps_low)
                sumxi_mid_allcases = np.copy(sumxi_mid)
                sumxi2_mid_allcases = np.copy(sumxi2_mid)
                nsamps_mid_allcases = np.copy(nsamps_mid)
                sumxi_high_allcases = np.copy(sumxi_high)
                sumxi2_high_allcases = np.copy(sumxi2_high)
                nsamps_high_allcases = np.copy(nsamps_high)
            else:
                closest_histogram = closest_histogram + closest_histogram_input
                sumxi_low_allcases = sumxi_low_allcases + sumxi_low
                sumxi2_low_allcases = sumxi2_low_allcases + sumxi2_low
                nsamps_low_allcases = nsamps_low_allcases + nsamps_low
                sumxi_mid_allcases = sumxi_mid_allcases + sumxi_mid
                sumxi2_mid_allcases = sumxi2_mid_allcases + sumxi2_mid
                nsamps_mid_allcases = nsamps_mid_allcases + nsamps_mid
                sumxi_high_allcases = sumxi_high_allcases + sumxi_high
                sumxi2_high_allcases = sumxi2_high_allcases + sumxi2_high
                nsamps_high_allcases = nsamps_high_allcases + nsamps_high

    return closest_histogram, sumxi_low_allcases, \
        sumxi2_low_allcases, nsamps_low_allcases, \
        sumxi_mid_allcases, sumxi2_mid_allcases, \
        nsamps_mid_allcases, sumxi_high_allcases, \
        sumxi2_high_allcases, nsamps_high_allcases
        
        
# ==================================================================    
    
def calc_mean_stddev(sumxi, sumxi2, nsamps):
    ns = len(nsamps)
    dressmean = np.zeros((ns), dtype=np.float64)
    dress_std = np.zeros((ns), dtype=np.float64)
    a = np.where(nsamps > 1)
    dressmean[a] = sumxi[a] / nsamps[a]
    dress_std[a] = np.sqrt( (sumxi2[a] - \
        sumxi[a]**2/nsamps[a]) / (nsamps[a]-1))
    return dressmean, dress_std
    
# ==================================================================

def weighted_linear_regression(precip_amt, dressmean, dress_std, nsamps):

    lns = len(nsamps)
    weights = 0.00000001*np.ones((len(nsamps)), dtype=np.float64)
    a = np.where(nsamps > 0)
    weights[a] = nsamps[a]
    weights = weights / weights.max()
    regr = LinearRegression()
    X = np.zeros((lns,1), dtype=np.float64)
    X[:,0] = precip_amt
    
    # ---- fit coefficients to the mean

    regr.fit(X, dressmean, sample_weight=weights)
    b0_mean = float(regr.intercept_)
    b1_mean = float(regr.coef_)
    
    # ---- fit coefficients to the spread
    regr.fit(X, dress_std, sample_weight=weights)
    b0_std = float(regr.intercept_)
    b1_std = float(regr.coef_)

    return b0_mean, b1_mean, b0_std, b1_std
        
        
    
# ==================================================================   
# MAIN PROGRAM STARTS HERE
# ==================================================================    

# --- read in the desired month (01 - 12)

cmonth = sys.argv[1] # 01-12
clead = sys.argv[2] # 3 digits, e.g., 024
ctype = '2' #sys.argv[3] #  on grid with every other point
ctype_out = ''
makeplot = True
precip_amt = np.arange(251)/10.

imonth = int(cmonth)
#master_directory_histogram_out = '/Volumes/NBM/conus_gefsv12/'+ctype+'/chist/'
#master_directory_histogram_in = '/Volumes/NBM/conus_gefsv12/'+ctype+'/'

master_directory_histogram_in = '/Volumes/NBM/chist/'
master_directory_histogram_out = '/Volumes/NBM/chist/'

# ---- get a list of retro run dates to read in, and a character string
#      describing the 3-month period centered on the month of interest

date_list, cmonthy = populate_datelists(imonth)

# ---- using the list of dates, return the closest-member histograms and
#      information that can be used to determine dressing statistics.
  
closest_histogram, sumxi_low_allcases, sumxi2_low_allcases, \
    nsamps_low_allcases, sumxi_mid_allcases, sumxi2_mid_allcases, \
    nsamps_mid_allcases, sumxi_high_allcases, sumxi2_high_allcases, \
    nsamps_high_allcases = read_closest_member_stats_and_dress_stats(\
    date_list, master_directory_histogram_in, clead)

# ---- determine the CDF of the closest-member histograms.

nmembers, ncats = np.shape(closest_histogram)

closest_histogram[:,0] = 99999

closest_histogram_raw = np.copy(closest_histogram.astype(np.float64))
closest_histogram_smoothed = np.copy(closest_histogram.astype(np.float64))
closest_histogram_cdf = np.copy(closest_histogram.astype(np.float64))

for i in range(7):
    closest_histogram_raw[:,i] = closest_histogram_raw[:,i] / \
        np.sum(closest_histogram_raw[:,i])
    closest_histogram_smoothed[:,i] = closest_histogram_smoothed[:,i] / \
        np.sum(closest_histogram_smoothed[:,i])
    for j in range(nmembers):
        closest_histogram_cdf[j,i] = np.sum(closest_histogram_raw[0:j,i])
    for i in range(1,7):
        csavgol = np.copy(closest_histogram_raw[:,i])
        csavgol2 = signal.savgol_filter(csavgol, 9, 1, mode='mirror')
        closest_histogram_smoothed[17:-17,i] = csavgol2[17:-17]          

# ---- calculate regression estimates of the mean and standard deviation
#      of dressing statistics as a function of the best-member precipitation
#      amount.

dressmean_smallest, dress_std_smallest = \
    calc_mean_stddev(sumxi_low_allcases, \
    sumxi2_low_allcases, nsamps_low_allcases)
b0_mean_smallest, b1_mean_smallest,\
    b0_std_smallest, b1_std_smallest = \
    weighted_linear_regression(\
    precip_amt, dressmean_smallest, \
    dress_std_smallest, nsamps_low_allcases)

dressmean_middle, dress_std_middle= \
    calc_mean_stddev(sumxi_mid_allcases, \
    sumxi2_mid_allcases, nsamps_mid_allcases)
b0_mean_mid, b1_mean_mid,\
    b0_std_mid, b1_std_mid = \
    weighted_linear_regression(\
    precip_amt, dressmean_middle, \
    dress_std_middle, nsamps_mid_allcases)

dressmean_highest, dress_std_highest = \
    calc_mean_stddev(sumxi_high_allcases, \
    sumxi2_high_allcases, nsamps_high_allcases)
b0_mean_largest, b1_mean_largest,\
    b0_std_largest, b1_std_largest = \
    weighted_linear_regression(\
    precip_amt, dressmean_highest, \
    dress_std_highest, nsamps_high_allcases)

# ===================================================================
# --- save closest-member histogram and dressing stats to netCDF file
# ===================================================================

outfile = master_directory_histogram_out +\
    'closest_member_histogram_'+ctype_out+\
    '_month='+cmonth+'_lead='+clead+'h.nc'
print ('writing to ',outfile)
nc = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

# --- initialize dimensions, variable names

ncatsd = nc.createDimension('ncatsd',ncats)
ncatsv = nc.createVariable('ncatsv','i4',('ncatsd',))
ncatsv.long_name = "category numbers of histogram"
ncatsv.units = "n/a"

nm1 = ncats-1
ncats_minus = nc.createDimension('ncats_minus',nm1)
ncatsv_minus = nc.createVariable('ncatsv_minus','i4',('ncats_minus',))
ncatsv_minus.long_name = "between-categories thresholds"
ncatsv_minus.units = "n/a"

nmem = nc.createDimension('nmem',nmembers)
nmemv = nc.createVariable('nmemv','i4',('nmem'))
nmemv.long_name = "member numbers (x25 here, with stencil)"
nmemv.units = "n/a"

one = nc.createDimension('one',1)
onev = nc.createVariable('onev','i4',('one'))
onev.long_name = "1-dimension for constant"
onev.units = "n/a"

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

# --- declare dressing statistics

b0_mean_lowrank = nc.createVariable('b0_mean_lowrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b0_mean_lowrank.units = 'mm'
b0_mean_lowrank.long_name = 'intercept regression coefficient for mean and '+\
    'for smallest sorted member dressing statistics.'
b0_mean_lowrank.valid_range = [-10.,10.]

b1_mean_lowrank = nc.createVariable('b1_mean_lowrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b1_mean_lowrank.units = 'mm'
b1_mean_lowrank.long_name = 'slope regression coefficient for mean and '+\
    'for smallest sorted member dressing statistics.'
b1_mean_lowrank.valid_range = [-10.,10.]

b0_std_lowrank = nc.createVariable('b0_std_lowrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b0_std_lowrank.units = 'mm'
b0_std_lowrank.long_name = 'intercept regression coefficient for std dev and '+\
    'for smallest sorted member dressing statistics.'
b0_std_lowrank.valid_range = [-10.,10.]

b1_std_lowrank = nc.createVariable('b1_std_lowrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b1_std_lowrank.units = 'mm'
b1_std_lowrank.long_name = 'slope regression coefficient for std dev and '+\
    'for smallest sorted member dressing statistics.'
b1_std_lowrank.valid_range = [-10.,10.]


b0_mean_midrank = nc.createVariable('b0_mean_midrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b0_mean_midrank.units = 'mm'
b0_mean_midrank.long_name = 'intercept regression coefficient for mean and '+\
    'for intermediate sorted member dressing statistics.'
b0_mean_midrank.valid_range = [-10.,10.]

b1_mean_midrank = nc.createVariable('b1_mean_midrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b1_mean_midrank.units = 'mm'
b1_mean_midrank.long_name = 'slope regression coefficient for mean and '+\
    'for intermediate sorted member dressing statistics.'
b1_mean_midrank.valid_range = [-10.,10.]

b0_std_midrank = nc.createVariable('b0_std_midrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b0_std_midrank.units = 'mm'
b0_std_midrank.long_name = 'intercept regression coefficient for std dev and '+\
    'for intermediate sorted member dressing statistics.'
b0_std_midrank.valid_range = [-10.,10.]

b1_std_midrank = nc.createVariable('b1_std_midrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b1_std_midrank.units = 'mm'
b1_std_midrank.long_name = 'slope regression coefficient for std dev and '+\
    'for intermediate sorted member dressing statistics.'
b1_std_midrank.valid_range = [-10.,10.]


b0_mean_highrank = nc.createVariable('b0_mean_highrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b0_mean_highrank.units = 'mm'
b0_mean_highrank.long_name = 'intercept regression coefficient for mean and '+\
    'for highest sorted member dressing statistics.'
b0_mean_highrank.valid_range = [-10.,10.]

b1_mean_highrank = nc.createVariable('b1_mean_highrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b1_mean_highrank.units = 'mm'
b1_mean_highrank.long_name = 'slope regression coefficient for mean and '+\
    'for highest sorted member dressing statistics.'
b1_mean_highrank.valid_range = [-10.,10.]

b0_std_highrank = nc.createVariable('b0_std_highrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b0_std_highrank.units = 'mm'
b0_std_highrank.long_name = 'intercept regression coefficient for std dev and '+\
    'for highest sorted member dressing statistics.'
b0_std_highrank.valid_range = [-10.,10.]

b1_std_highrank = nc.createVariable('b1_std_highrank','f8',\
    ('one'),zlib=True,least_significant_digit=6)
b1_std_highrank.units = 'mm'
b1_std_highrank.long_name = 'slope regression coefficient for std dev and '+\
    'for highest sorted member dressing statistics.'
b1_std_highrank.valid_range = [-10.,10.]


# ---- metadata

nc.title = 'closest_member histogram and dressing stats for this month,  '+\
    'forecast lead using 5x5 stencil of points and GEFSv12 ensemble (31 members)'
nc.history = 'Tom Hamill/Diana Stovern, munging the GEFSv12 data'
nc.institution =  "NOAA Physical Sciences Lab"

# ---- initialize and copy data to records

ncatsv = range(ncats)
ncatsv_minus = range(ncats-1)
nmemv = range(nmembers)
thresholds[:] = [0.01, 0.1, 0.5, 2.0, 6.0, 15.0]
closest_hist[:] = closest_histogram_raw[:,:]

b0_mean_lowrank[:] = np.max([b0_mean_smallest, 0.0001])
b1_mean_lowrank[:] = np.max([b1_mean_smallest, 0.0001])
b0_std_lowrank[:] = np.max([b0_std_smallest, 0.0001])
b1_std_lowrank[:] = np.max([b1_std_smallest, 0.0001])

b0_mean_midrank[:] = np.max([b0_mean_mid, 0.0001])
b1_mean_midrank[:] = np.max([b1_mean_mid, 0.0001])
b0_std_midrank[:] = np.max([b0_std_mid, 0.0001])
b1_std_midrank[:] = np.max([b1_std_mid, 0.0001])

b0_mean_highrank[:] = np.max([b0_mean_largest, 0.0001])
b1_mean_highrank[:] = np.max([b1_mean_largest, 0.0001])
b0_std_highrank[:] = np.max([b0_std_largest, 0.0001])
b1_std_highrank[:] = np.max([b1_std_largest, 0.0001])

print('b0_mean_smallest = ',b0_mean_smallest)
print('b1_mean_smallest = ',b1_mean_smallest)
print('b0_std_smallest = ',b0_std_smallest)
print('b1_std_smallest = ',b1_std_smallest)

print('b0_mean_mid = ',b0_mean_mid)
print('b1_mean_mid = ',b1_mean_mid)
print('b0_std_mid = ',b0_std_mid)
print('b1_std_mid = ',b1_std_mid)

print('b0_mean_largest = ',b0_mean_largest)
print('b1_mean_largest = ',b1_mean_largest)
print('b0_std_largest = ',b0_std_largest)
print('b1_std_largest = ',b1_std_largest)

nc.close()

# ===================================================================
# ---- now make closest histogram plots
# ===================================================================

if makeplot == True:
    fig = plt.figure(figsize=(9.,6.))
    plt.suptitle('Closest-member histograms for '+cmonthy+\
        ', lead = +'+clead+' h',fontsize=18)
    rcParams['legend.fontsize']='small'
    for i in range(1,7):
    
        if i == 1:
            ctitle = '(a) 0.01 mm '+r'$\leq \overline{\tilde{x}}^f_i$ < 0.1 mm '
            axlocn = [0.09,0.55,0.22,0.32]
        elif i == 2:
            ctitle = '(b) 0.1 mm '+r'$\leq \overline{\tilde{x}}^f_i$ < 0.5 mm '
            axlocn = [0.42,0.55,0.22,0.32]
        elif i == 3:
            ctitle = '(c) 0.5 mm '+r'$\leq \overline{\tilde{x}}^f_i$ < 2 mm'
            axlocn = [0.75,0.55,0.22,0.32]
        elif i == 4:
            ctitle = '(d) 2 mm '+r'$\leq \overline{\tilde{x}}^f_i < $ 6 mm'
            axlocn = [0.09,0.09,0.22,0.32]
        elif i == 5:
            ctitle = '(e) 6 mm '+r'$\leq \overline{\tilde{x}}^f_i < $ 15 mm'
            axlocn = [0.42,0.09,0.22,0.32]
        else:
            ctitle = '(f) '+r'$\overline{\tilde{x}}^f_i \geq$ 15 mm'
            axlocn = [0.75,0.09,0.22,0.32]
        
        a1 = fig.add_axes(axlocn)
        a1.set_title(ctitle,fontsize=11)
        a1.set_xlabel('Fraction between lowest & highest rank',fontsize=8)
        a1.set_ylabel('Fraction of analyzed closest\nto this sorted member',fontsize=8)
        a1.set_xlim(0,1)
        a1.set_ylim(0.0001,0.1)
        a1.set_yscale('log')
        a1.grid(color='Gray',lw=0.2,linestyle='--')

        # --- apply Savitzky-Golay smoother to deal with sampling variability for
        #     internal ranks.
    
        ctext_low_N = "{:.3f}".format(closest_histogram_raw[0,i])
        ctext_high_N = "{:.3f}".format(closest_histogram_raw[-1,i])    

        a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
            closest_histogram_raw[:,i],'-',\
            color='Red',linewidth=1.) # ,label='Raw'
    
        #a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        #    2.*closest_histogram_smoothed[:,i] ,'-',\
        #    color='RoyalBlue',label=r'Smoothed x 2',linewidth=2.)
            
        a1.text(0.03,float(ctext_low_N),r'$\leftarrow\ $'+ctext_low_N)
        a1.text(0.67,float(ctext_high_N),ctext_high_N+r'$\rightarrow\ $')        
        
        #a1.legend(loc=9)

    plot_title = 'closest_member_histogram_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)

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
        closest_histogram_cdf[:,0],color='Red',linewidth=1.2,label='< 0.01')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,1],color='RoyalBlue',linewidth=1.2,label='0.01 to 0.1')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,2],color='LimeGreen',linewidth=1.2,linestyle='-.',label='0.1 to 0.5')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,3],color='Gray',linewidth=1.2,linestyle='--',label='0.5 to 2.0')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,4],color='Black',linewidth=1.2,label='2.0 to 6.0')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,5],color='Purple',linewidth=1.2,linestyle='--',label='6.0 to 15.0')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_histogram_cdf[:,6],color='Orange',linewidth=1.2,linestyle='--',label='> 15.0')
    
    a1.legend(loc=0)

    plot_title = 'closest_member_histogram_CDF_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    

# ---- make plots of the closest-member histogram exceedance function

closest_hist_exceedance = np.copy(closest_histogram_raw)
for icat in range(7):
    for imem in range(nmembers):
        closest_hist_exceedance[imem,icat] = \
            np.sum(closest_histogram_raw[nmembers-imem:nmembers,icat])

makeplot = True
if makeplot == True:
    fig = plt.figure(figsize=(6.,6.5))
    rcParams['legend.fontsize']='large'
    ctitle = 'Exceedance distribution functions\nfor '+cmonthy+\
        ', lead = +'+clead+' h'
    axlocn = [0.11,0.11,0.85,0.79]

    a1 = fig.add_axes(axlocn)
    a1.set_title(ctitle,fontsize=17)
    a1.set_xlabel('Fraction of members exceeding the threshold',fontsize=14)
    a1.set_ylabel('Exceedance probability',fontsize=14)
    a1.set_xlim(0,1)
    a1.set_ylim(0,1)
    a1.grid(color='Gray',lw=0.2,linestyle='--')

    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_hist_exceedance[:,0],color='Red',linewidth=1.2,\
        label=r'$\overline{\tilde{x}}^f_i$ < 0.01 mm')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_hist_exceedance[:,1],color='RoyalBlue',\
        linewidth=1.2,label=r'0.01 mm $\leq \overline{\tilde{x}}^f_i$ < 0.1 mm ')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_hist_exceedance[:,2],color='LimeGreen',\
        linewidth=1.2,linestyle='-.',label=r'0.1 mm $\leq \overline{\tilde{x}}^f_i$ < 0.5 mm')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_hist_exceedance[:,3],color='Gray',\
        linewidth=1.2,linestyle='--',label=r'0.5 mm $\leq \overline{\tilde{x}}^f_i$ < 2 mm')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_hist_exceedance[:,4],color='Black',\
        linewidth=1.2,label=r'2 mm $\leq \overline{\tilde{x}}^f_i $<  6 mm')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_hist_exceedance[:,5],color='Orchid',\
        linewidth=1.2,linestyle='--',label=r'6 mm $\leq \overline{\tilde{x}}^f_i $<  15 mm')
    a1.plot(np.real(range(1,nmembers+1))/float(nmembers), \
        closest_hist_exceedance[:,6],color='Orange',\
        linewidth=1.2,linestyle='--',label=r'$\overline{\tilde{x}}^f_i \geq$ 15 mm')
        

    a1.legend(loc=0)

    plot_title = 'closest_member_histogram_EDF_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)    


# ---- now make dressing statistics plots

if makeplot == True:
    #fig = plt.figure(figsize=(6.,9.))
    fig = plt.figure(figsize=(3.4,9.))
    #plt.suptitle('Dressing statistics for '+cmonthy+', lead = +'+clead+' h',fontsize=16)
    plt.suptitle('Dressing statistics for\n'+cmonthy+', lead = +'+clead+' h',fontsize=16)
    for i in range(3):
    
        if i == 0:
            ctitle = r'(a) Lowest member'
            #axlocn = [0.12,0.69,0.78,0.22]
            axlocn = [0.19,0.69,0.63,0.2]
            dressmean = dressmean_smallest
            dress_std = dress_std_smallest
            xtop = 4.
            meanpredict = b0_mean_smallest + b1_mean_smallest*precip_amt
            stdpredict = b0_std_smallest + b1_std_smallest*precip_amt
            nsamps = nsamps_low_allcases
            ymax = 10000000
            xmax = 25
        elif i == 1:
            ctitle = r'(b) Intermediate members'
            axlocn = [0.19,0.38,0.63,0.2]
            dressmean = dressmean_middle
            dress_std = dress_std_middle
            xtop = 25.
            meanpredict = b0_mean_mid + b1_mean_mid*precip_amt
            stdpredict = b0_std_mid + b1_std_mid*precip_amt
            nsamps = nsamps_mid_allcases
            ymax = 10000000
            xmax = 25
        else:
            ctitle = r'(c) Highest member'
            axlocn = [0.19,0.07,0.63,0.2]
            dressmean = dressmean_highest
            dress_std = dress_std_highest
            xtop = 50.
            meanpredict = b0_mean_largest + b1_mean_largest*precip_amt
            stdpredict = b0_std_largest + b1_std_largest*precip_amt
            nsamps = nsamps_high_allcases
            ymax = 10000000
            xmax = 25.
        
        a1 = fig.add_axes(axlocn)
        a1.set_title(ctitle,fontsize=11)
        a1.set_xlabel('Precipitation amount of\nclosest dressed member (mm)',fontsize=9)
        a1.set_ylabel('Dressing distribution\nmean and spread (mm)',fontsize=9)
        a1.set_xlim(0,25.)
        a1.set_ylim(0.001, xtop)
        a1.grid(color='Gray',lw=0.2,linestyle='--',axis='x')

        rcParams['legend.fontsize']='small'
        a1.plot(precip_amt,dressmean,'-',\
            color='Red',linewidth=1.5,label='Empirical mean')
        a1.plot(precip_amt,dress_std,'-',\
            color='RoyalBlue',linewidth=1.5,label=r'Empirical $\sigma$')
        a1.plot(precip_amt,meanpredict,'-',\
            color='Red',linewidth=1.5,linestyle=':',label='Regressed mean')
        a1.plot(precip_amt,stdpredict,'-',\
            color='RoyalBlue',linewidth=1.5,linestyle=':', label=r'Regressed $\sigma$')
        a1.plot([0.,0.],[0.0, 0.-0], linewidth=1.25, label='Sample size',color='Gray')
        a1.plot(precip_amt, precip_amt,color='LightGray',linewidth=0.5)
        if i == 0: a1.legend(loc=4)
        
        
        a2 = a1.twinx()  # instantiate a second axes that shares the same x-axis
        a2.set_ylabel('Number of samples')  # we already handled the x-label with ax1
        a2.set_yscale('log')
        a2.plot(precip_amt,nsamps, color='Gray', lw='1.25',label='Sample size')
        a2.set_ylim(1,ymax)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
    
    plot_title = 'dress_stats_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)



    fig = plt.figure(figsize=(9., 3.3))
    plt.suptitle('Dressing statistics for '+cmonthy+', lead = +'+clead+' h',fontsize=16)
    for i in range(3):
    
        if i == 0:
            ctitle = r'(a) Lowest member'
            #axlocn = [0.12,0.69,0.78,0.22]
            axlocn = [0.04,0.18,0.22,0.62]
            dressmean = dressmean_smallest
            dress_std = dress_std_smallest
            xtop = 4.
            meanpredict = b0_mean_smallest + b1_mean_smallest*precip_amt
            stdpredict = b0_std_smallest + b1_std_smallest*precip_amt
            nsamps = nsamps_low_allcases
            ymax = 10000000
            xmax = 25
        elif i == 1:
            ctitle = r'(b) Intermediate members'
            axlocn = [0.38,0.18,0.22,0.62]
            dressmean = dressmean_middle
            dress_std = dress_std_middle
            xtop = 25.
            meanpredict = b0_mean_mid + b1_mean_mid*precip_amt
            stdpredict = b0_std_mid + b1_std_mid*precip_amt
            nsamps = nsamps_mid_allcases
            ymax = 10000000
            xmax = 25
        else:
            ctitle = r'(c) Highest member'
            axlocn = [0.72,0.18,0.22,0.62]
            dressmean = dressmean_highest
            dress_std = dress_std_highest
            xtop = 50.
            meanpredict = b0_mean_largest + b1_mean_largest*precip_amt
            stdpredict = b0_std_largest + b1_std_largest*precip_amt
            nsamps = nsamps_high_allcases
            ymax = 10000000
            xmax = 25.
        
        a1 = fig.add_axes(axlocn)
        a1.set_title(ctitle,fontsize=11)
        a1.set_xlabel('Precipitation amount of\nclosest dressed member (mm)',fontsize=7)
        a1.set_ylabel('Dressing distribution mean and spread (mm)',fontsize=7)
        a1.set_xlim(0,25.)
        a1.set_ylim(0.001, xtop)
        a1.grid(color='Gray',lw=0.2,linestyle='--',axis='x')

        rcParams['legend.fontsize']='x-small'
        a1.plot(precip_amt,dressmean,'-',\
            color='Red',linewidth=1.5,label='Empirical mean')
        a1.plot(precip_amt,dress_std,'-',\
            color='RoyalBlue',linewidth=1.5,label=r'Empirical $\sigma$')
        a1.plot(precip_amt,meanpredict,'-',\
            color='Red',linewidth=1.5,linestyle=':',label='Regressed mean')
        a1.plot(precip_amt,stdpredict,'-',\
            color='RoyalBlue',linewidth=1.5,linestyle=':', label=r'Regressed $\sigma$')
        a1.plot([0.,0.],[0.0, 0.-0], linewidth=1.25, label='Sample size',color='Gray')
        a1.plot(precip_amt, precip_amt,color='LightGray',linewidth=0.5)
        if i == 0: a1.legend(loc=0)
        
        
        a2 = a1.twinx()  # instantiate a second axes that shares the same x-axis
        a2.set_ylabel('Number of samples',fontsize=7,color='Gray')  # we already handled the x-label with ax1
        a2.set_yscale('log')
        a2.plot(precip_amt,nsamps, color='Gray', lw='1.25',label='Sample size')
        a2.set_ylim(1,ymax)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
    
    plot_title = 'dress_stats_horiz_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)