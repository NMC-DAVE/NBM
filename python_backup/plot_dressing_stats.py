"""

python plot_dressing_stats.py cmonth clead
    where cmonth = 01,12
    clead = 3 digits, e.g., 024

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
import numpy.ma as ma
from sklearn.linear_model import LinearRegression

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

    print ('nsamps = ', nsamps)
    print ('len(nsamps) = ', len(nsamps))
    lns = len(nsamps)
    weights = 0.00000001*np.ones((len(nsamps)), dtype=np.float64)
    a = np.where(nsamps > 0)
    weights[a] = nsamps[a]
    weights = weights / weights.max()
    regr = LinearRegression()
    X = np.zeros((lns,1), dtype=np.float64)
    X[:,0] = precip_amt
    
    # ---- fit coefficients to the mean
    #print ('weights = ', weights)
    #print ('precip_amt = ', precip_amt)
    #print ('dressmean = ', dressmean)
    #print ('np.shape(weights) = ', np.shape(weights))
    #print ('np.shape(precip_amt) = ', np.shape(precip_amt))
    #print ('np.shape(dressmean) = ', np.shape(dressmean))
    regr.fit(X, dressmean, sample_weight=weights)
    #print ('')
    b0_mean = regr.intercept_
    b1_mean = regr.coef_
    print ('b0_mean, b1_mean = ', b0_mean, b1_mean)
    # ---- fit coefficients to the spread
    regr.fit(X, dress_std, sample_weight=weights)
    b0_std = regr.intercept_
    b1_std = regr.coef_

    return b0_mean, b1_mean, b0_std, b1_std

# ================================================================== 

# --- read in the desired month (01 - 12)

cmonth = sys.argv[1] # 01-12
clead = sys.argv[2] # 3 digits, e.g., 024
ctype = 'thinned2' # sys.argv[3] # thinned or upscaled
makeplot = True

imonth = int(cmonth)
master_directory_histogram = '/Volumes/NBM/conus_gefsv12/'+ctype+'/'

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


#date_list = date_list1 + date_list2 + date_list3
date_list = date_list1 
ifirst = True
for idate, date in enumerate(date_list):
    
    # --- read in c-m histogram for day of interest.
    
    infile = master_directory_histogram +\
        'closest_histogram'+date+'_'+clead+'.cPick'
    print (infile)
    fexist = False
    fexist = os.path.exists(infile)
    if fexist:
        file_size = os.path.getsize(infile)
        #print (file_size)
    
    if fexist and file_size > 28000:
        inf = open(infile, 'rb')
        print (infile)
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
        print (sumxi_high)
        #print (sumxi_high)
        print (nsamps_high)        
        inf.close()
    
        # --- set up arrays if this is the first time through
    
        if ifirst == True:
            ifirst = False
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
            sumxi_low_allcases = sumxi_low_allcases + sumxi_low
            sumxi2_low_allcases = sumxi2_low_allcases + sumxi2_low
            nsamps_low_allcases = nsamps_low_allcases + nsamps_low
            sumxi_mid_allcases = sumxi_mid_allcases + sumxi_mid
            sumxi2_mid_allcases = sumxi2_mid_allcases + sumxi2_mid
            nsamps_mid_allcases = nsamps_mid_allcases + nsamps_mid
            sumxi_high_allcases = sumxi_high_allcases + sumxi_high
            sumxi2_high_allcases = sumxi2_high_allcases + sumxi2_high
            nsamps_high_allcases = nsamps_high_allcases + nsamps_high

precip_amt = np.arange(251)/10.

# ---- now make closest histogram plots

if makeplot == True:
    fig = plt.figure(figsize=(6.,9.))
    plt.suptitle('Dressing statistics\n'+cmonthy+\
        ', lead = +'+clead+' h',fontsize=18)
    for i in range(3):
    
        if i == 0:
            ctitle = r'(a) Smallest member'
            axlocn = [0.17,0.66,0.78,0.2]
            dressmean, dress_std = \
                calc_mean_stddev(sumxi_low_allcases, \
                sumxi2_low_allcases, nsamps_low_allcases)
            xtop = 3.
            b0_mean_smallest, b1_mean_smallest,\
                b0_std_smallest, b1_std_smallest = \
                weighted_linear_regression(\
                precip_amt, dressmean, dress_std, nsamps_low_allcases)
            #print ('b0, b1 mean smallest = ', b0_mean_smallest, b1_mean_smallest )
            meanpredict = b0_mean_smallest + b1_mean_smallest*precip_amt
            stdpredict = b0_std_smallest + b1_std_smallest*precip_amt
        elif i == 1:
            ctitle = r'(b) Middle-ranked members'
            axlocn = [0.17,0.36,0.78,0.2]
            dressmean, dress_std = \
                calc_mean_stddev(sumxi_mid_allcases, \
                sumxi2_mid_allcases, nsamps_mid_allcases)
            xtop = 3.
            b0_mean_mid, b1_mean_mid,\
                b0_std_mid, b1_std_mid = \
                weighted_linear_regression(\
                precip_amt, dressmean, dress_std, nsamps_mid_allcases)
            meanpredict = b0_mean_mid + b1_mean_mid*precip_amt
            stdpredict = b0_std_mid + b1_std_mid*precip_amt
        else:
            ctitle = r'(c) Largest member'
            axlocn = [0.17,0.06,0.78,0.2]
            dressmean, dress_std = \
                calc_mean_stddev(sumxi_high_allcases, \
                sumxi2_high_allcases, nsamps_high_allcases)
            xtop = 50.
            b0_mean_largest, b1_mean_largest,\
                b0_std_largest, b1_std_largest = \
                weighted_linear_regression(\
                precip_amt, dressmean, dress_std, nsamps_high_allcases)
            meanpredict = b0_mean_largest + b1_mean_largest*precip_amt
            stdpredict = b0_std_largest + b1_std_largest*precip_amt

        #print (np.shape(dressmean))
        #print (np.shape(dress_std))
        
        a1 = fig.add_axes(axlocn)
        a1.set_title(ctitle,fontsize=14)
        a1.set_xlabel('Precipitation amount of closest dressed member',fontsize=12)
        a1.set_ylabel('Mean and spread',fontsize=11)
        a1.set_xlim(0,25.)
        a1.set_ylim(0.001, xtop)
        #a1.set_yscale('log')
        a1.grid(color='Gray',lw=0.2,linestyle='--')

        a1.plot(precip_amt,dressmean,'-',\
            color='Red',linewidth=1.5,label='Mean')
        a1.plot(precip_amt,dress_std,'-',\
            color='RoyalBlue',linewidth=1.5,label='Std Dev')
        a1.plot(precip_amt,meanpredict,'-',\
            color='Red',linewidth=.5)
        #if i == 0:
        #   print ('meanpredict = ', meanpredict)
        #    sys.exit()
        a1.plot(precip_amt,stdpredict,'-',\
            color='RoyalBlue',linewidth=.5)
        a1.plot(precip_amt, precip_amt,color='LightGray',linewidth=0.5)
        a1.legend(loc=9)

    plot_title = 'dress_stats_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)


