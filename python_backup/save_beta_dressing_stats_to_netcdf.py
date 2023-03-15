"""

python save_beta_dressing_stats_to_netcdf.py cmonth clead 
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
from matplotlib import rcParams
from sklearn.linear_model import LinearRegression
#import statsmodels
import statsmodels.api as sm
lowess = sm.nonparametric.lowess

# ==================================================================    
    
def populate_datelists(imonth):
    """ for the integer month, 1-12, get a list of the retro dates and 
    a character string to use in possible plots """
    
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
    return date_list, cmonthy

# ==================================================================    

def read_beta_stats(date_list, \
    master_directory_histogram_in, clead):
    
    """ from the daily file output of closest-member statistics
        and dressing statistics, tally up information to determine
        them over multiple days and months. """

    ndates = len(date_list)
    n251 = 251
    for idate, date in enumerate(date_list):
  
        # --- read in c-m histogram and stats to determine
        #     dresssing mean and std dev for day of interest.
    
        infile = master_directory_histogram_in + \
            'beta_stats_for_closest_histogram' + \
            date + '_' + clead + '.cPick'
        fexist = False
        fexist = os.path.exists(infile)
        print (infile, fexist)
        if fexist:
        
            inf = open(infile, 'rb')
            sum_fracrank = cPickle.load(inf)
            sum_fracrank_squared = cPickle.load(inf)
            nsamps_fracrank = cPickle.load(inf)            
            inf.close()
            
            # --- set up arrays if this is the first time through
    
            if idate == 0:
                n251 = len(sum_fracrank)
                sum_fracrank_allcases = np.copy(sum_fracrank)
                sum_fracrank_squared_allcases = np.copy(sum_fracrank_squared)
                nsamps_fracrank_allcases = np.copy(nsamps_fracrank)
                sum_fracrank_bycase = np.zeros((ndates,n251), dtype=np.float64)
                sum_fracrank_squared_bycase = np.zeros((ndates,n251), dtype=np.float64)
                nsamps_fracrank_bycase = np.zeros((ndates,n251), dtype=int)
            else:
                sum_fracrank_allcases = sum_fracrank_allcases + sum_fracrank
                sum_fracrank_squared_allcases = \
                    sum_fracrank_squared_allcases + sum_fracrank_squared
                nsamps_fracrank_allcases = \
                    nsamps_fracrank_allcases + nsamps_fracrank
                    
            sum_fracrank_bycase[idate,:] = sum_fracrank[:]
            sum_fracrank_squared_bycase[idate,:] = sum_fracrank_squared[:]
            nsamps_fracrank_bycase[idate,:] = nsamps_fracrank[:]

    return n251, ndates, sum_fracrank_allcases, sum_fracrank_squared_allcases, \
        nsamps_fracrank_allcases, sum_fracrank_bycase,\
        sum_fracrank_squared_bycase, nsamps_fracrank_bycase
        
        
# ==================================================================    
    
def calc_mean_stddev(sum_fracrank, sum_fracrank_squared, nsamps):
    #print ('in calc_mean_stddev')
    ns = len(nsamps)
    beta_mean = -99.99*np.ones((ns), dtype=np.float64)
    beta_std = -99.99*np.ones((ns), dtype=np.float64)
    a = np.where(nsamps > 100)
    beta_mean[a] = sum_fracrank[a] / nsamps[a]
    beta_std[a] = np.sqrt( (sum_fracrank_squared[a] - \
        sum_fracrank[a]**2 / nsamps[a]) / (nsamps[a]-1) )
    #print ('before loop')
    #for isamp in range(251):
    #    print (isamp, beta_mean[isamp], beta_std[isamp])
    #sys.exit()
    return beta_mean, beta_std
            
# ==================================================================     
    
def calc_alpha_beta(beta_mean, beta_std):
    ns = len(beta_mean)
    a = np.where(beta_mean > 0.0)
    alpha = -99.99*np.ones((ns), dtype=np.float64)
    beta = -99.99*np.ones((ns), dtype=np.float64)
    alpha[a] = beta_mean[a]**2 * (1.-beta_mean[a]) / beta_std[a]**2 - beta_mean[a]
    beta[a] = alpha[a]*(1.-beta_mean[a])/beta_mean[a]
    return alpha, beta

# ==================================================================    

# ==================================================================   
# MAIN PROGRAM STARTS HERE
# ==================================================================    

# --- read in the desired month (01 - 12)

cmonth = sys.argv[1] # 01-12
clead = sys.argv[2] # 3 digits, e.g., 024
ctype = 'thinned2' #sys.argv[3] # thinned on grid with every other point
ctype_out = 'thinned'
makeplot = True
precip_amt = np.arange(251)/5.

imonth = int(cmonth)
#master_directory_histogram_out = '/Volumes/NBM/conus_gefsv12/'+ctype+'/chist/'
#master_directory_histogram_in = '/Volumes/NBM/conus_gefsv12/'+ctype+'/'

master_directory_histogram_in = '/Volumes/NBM/conus_gefsv12/thinned2/'
master_directory_histogram_out = '/Volumes/NBM/chist/'

# ---- get a list of retro run dates to read in, and a character string
#      describing the 3-month period centered on the month of interest

date_list, cmonthy = populate_datelists(imonth)
ndates = len(date_list)
#print (date_list)

# ---- using the list of dates, return the info needed to calculate
#      mean and std dev to set beta distribution parameters.


n251, ndates, sum_fracrank_allcases, sum_fracrank_squared_allcases, \
    nsamps_fracrank_allcases, sum_fracrank_bycase,\
    sum_fracrank_squared_bycase, nsamps_fracrank_bycase = \
    read_beta_stats(date_list, master_directory_histogram_in, clead)
      
# ---- calculate regression estimates of the mean and standard deviation
#      of dressing statistics as a function of the best-member precipitation
#      amount.

mean_for_beta, std_for_beta = calc_mean_stddev(\
    sum_fracrank_allcases, sum_fracrank_squared_allcases, \
    nsamps_fracrank_allcases)
    
# --- calculate beta distribution parameters where there are sufficient
#     sample sizes

alpha_hat, beta_hat = calc_alpha_beta(mean_for_beta, std_for_beta)
print ('i   precip   beta_mean beta_std nsamps, std_for_beta,  alpha, beta')
for i in range(251):
    print (i,precip_amt[i], nsamps_fracrank_allcases[i], mean_for_beta[i],  \
        std_for_beta[i], alpha_hat[i], beta_hat[i])    
    
# ---- implement 100x resampling with replacement to determine uncertainty
#      in beta parameters

nresa = 100
mean_for_beta_resamp = np.zeros((nresa,n251), dtype=np.float64)
std_for_beta_resamp= np.zeros((nresa,n251), dtype=np.float64)
mean_for_beta_resamp = np.zeros((n251), dtype=np.float64)
std_for_beta_resamp = np.zeros((n251), dtype=np.float64)
alpha_hat_resamp = np.zeros((nresa,n251), dtype=np.float64)
beta_hat_resamp = np.zeros((nresa,n251), dtype=np.float64)
for iresa in range(nresa):
    
    idates = np.random.randint(low=0, high=ndates, size=ndates)
    sum_fracrank_resamp = np.sum(sum_fracrank_bycase[idates,:],axis=0)
    sum_fracrank_squared_resamp = np.sum(sum_fracrank_squared_bycase[idates,:],axis=0)
    nsamps_fracrank_resamp = np.sum(nsamps_fracrank_bycase[idates,:],axis=0)
    mean_for_beta_resamp, std_for_beta = calc_mean_stddev(\
        sum_fracrank_resamp, sum_fracrank_squared_resamp, \
        nsamps_fracrank_resamp)
    alpha_hat_resamp[iresa,:], beta_hat_resamp[iresa,:] = \
        calc_alpha_beta(mean_for_beta, std_for_beta)
        
# --- apply LOWESS regression        
        
alpha_fitted = lowess(alpha_hat[1:], \
    precip_amt[1:], frac=0.15, it=3, delta=2.0, xvals=None, \
    is_sorted=False, missing='drop', return_sorted=False)
        
beta_fitted = lowess(beta_hat[1:], \
    precip_amt[1:], frac=0.15, it=3, delta=2.0, xvals=None, \
    is_sorted=False, missing='drop', return_sorted=False)
alpha_fitted = np.insert(alpha_fitted,0,alpha_hat[0])
beta_fitted = np.insert(beta_fitted,0,alpha_hat[0])

# --- get 5,25,75,95 of resampled distribution parameters

alpha_hat_05 = np.zeros(len(alpha_hat), dtype=np.float64)
alpha_hat_25 = np.zeros(len(alpha_hat), dtype=np.float64)
alpha_hat_75 = np.zeros(len(alpha_hat), dtype=np.float64)
alpha_hat_95 = np.zeros(len(alpha_hat), dtype=np.float64)
beta_hat_05 = np.zeros(len(beta_hat), dtype=np.float64)
beta_hat_25 = np.zeros(len(beta_hat), dtype=np.float64)
beta_hat_75 = np.zeros(len(beta_hat), dtype=np.float64)
beta_hat_95 = np.zeros(len(beta_hat), dtype=np.float64)
for i251 in range(n251):
    alpha_sorted = np.sort(alpha_hat_resamp[:,i251])
    alpha_hat_05[i251] = alpha_sorted[4]
    alpha_hat_25[i251] = alpha_sorted[24]
    alpha_hat_75[i251] = alpha_sorted[74]
    alpha_hat_95[i251] = alpha_sorted[94]
    beta_sorted = np.sort(beta_hat_resamp[:,i251])
    beta_hat_05[i251] = beta_sorted[4]
    beta_hat_25[i251] = beta_sorted[24]
    beta_hat_75[i251] = beta_sorted[74]
    beta_hat_95[i251] = beta_sorted[94]

# ===================================================================
# --- save closest-member histogram and dressing stats to netCDF file
# ===================================================================

outfile = master_directory_histogram_out +\
    'beta_statistics_'+ctype_out+\
    '_month='+cmonth+'_lead='+clead+'h.nc'
print ('writing to ',outfile)
nc = Dataset(outfile,'w',format='NETCDF4_CLASSIC')

# --- initialize dimensions, variable names

n251d = nc.createDimension('n251d',n251)
n251v = nc.createVariable('n251d','i4',('n251d',))
n251v.long_name = "indices into the alpha, beta array"
n251v.units = "n/a"

# --- declare alpha, beta 

precip_amt = nc.createVariable('precip_amt','f8',\
    ('n251d'),zlib=True,least_significant_digit=6)
precip_amt.units = 'mm'
precip_amt.long_name = 'ensemble-mean precipitation amount'
precip_amt.valid_range = [0.0,25.0]

alpha = nc.createVariable('alpha','f8',\
    ('n251d'),zlib=True,least_significant_digit=6)
alpha.units = 'mm'
alpha.long_name = 'fitted alpha coefficient for beta distribution'
alpha.valid_range = [0.0,25.0]

beta = nc.createVariable('beta','f8',\
    ('n251d'),zlib=True,least_significant_digit=6)
beta.units = 'mm'
beta.long_name = 'fitted alpha coefficient for beta distribution'
beta.valid_range = [0.0,25.0]

# ---- metadata

nc.title = 'beta distribution coefficients used to estimate closest-member histogram'
nc.history = 'Tom Hamill, from the GEFSv12 data'
nc.institution =  "NOAA Physical Sciences Lab"

# ---- initialize and copy data to records

n251v = range(n251)
precip_amt = np.arange(n251)/5.
alpha[:] = alpha_hat[:]
beta[:] = beta_hat[:]

nc.close()

# ===================================================================
# ---- now make closest histogram plots estimated from beta distribution
# ===================================================================

beta_pdf = np.zeros((251,1001), dtype=np.float64)
beta_cdf = np.zeros((251,1001), dtype=np.float64)
beta_pdf_lowess = np.zeros((251,1001), dtype=np.float64)
beta_cdf_lowess = np.zeros((251,1001), dtype=np.float64)
chist = np.zeros((251,1000), dtype=np.float64)
chist_lowess = np.zeros((251,1000), dtype=np.float64)
x = np.arange((1001), dtype=np.float64) / 1000.
print (x[0::10])
x1000 = 0.0005 + np.arange((1000), dtype=np.float64) / 1000.

if makeplot == True:
    fig = plt.figure(figsize=(9.,6.))
    plt.suptitle(r'$\beta$ closest-member histograms for '+cmonthy+\
        ', lead = +'+clead+' h',fontsize=18)
    rcParams['legend.fontsize']='xx-small'
    
    for idx in range(251):
        beta_pdf[idx,:] = stats.beta.pdf(x,alpha_hat[idx], beta_hat[idx])
        beta_cdf[idx,:] = stats.beta.cdf(x,alpha_hat[idx], beta_hat[idx])
        beta_pdf_lowess[idx,:] = stats.beta.pdf(x,alpha_fitted[idx], beta_fitted[idx])
        beta_cdf_lowess[idx,:] = stats.beta.cdf(x,alpha_fitted[idx], beta_fitted[idx])
        for j in range(1000):
            chist[idx,j] = beta_cdf[idx,j+1] - beta_cdf[idx,j]
            chist_lowess[idx,j] = beta_cdf_lowess[idx,j+1] - beta_cdf_lowess[idx,j]
    
    for i in range(1,7):
    
        if i == 1:
            ctitle = '(a) 10 mm '
            idx = 50
            axlocn = [0.09,0.55,0.22,0.32]
            print ('10 mm: alpha, beta = ',alpha_hat[idx], beta_hat[idx])
        elif i == 2:
            ctitle = '(b) 20 mm '
            idx = 100
            axlocn = [0.42,0.55,0.22,0.32]
            print ('20 mm: alpha, beta = ',alpha_hat[idx], beta_hat[idx])
        elif i == 3:
            ctitle = '(c) 26 mm'
            idx = 130
            axlocn = [0.75,0.55,0.22,0.32]
            print ('26  mm: alpha, beta = ',alpha_hat[idx], beta_hat[idx])
        elif i == 4:
            ctitle = '(d) 34 mm'
            idx = 170
            axlocn = [0.09,0.09,0.22,0.32]
            print ('34 mm: alpha, beta = ',alpha_hat[idx], beta_hat[idx])
        elif i == 5:
            ctitle = '(e) 42 mm'
            idx = 210
            axlocn = [0.42,0.09,0.22,0.32]
            print ('42 mm: alpha, beta = ',alpha_hat[idx], beta_hat[idx])
        else:
            ctitle = '(f) 50 mm'
            idx = 250
            axlocn = [0.75,0.09,0.22,0.32]
            print ('50 mm: alpha, beta = ',alpha_hat[idx], beta_hat[idx])


        #print (idx, j, chist[idx,0::33])
        
        a1 = fig.add_axes(axlocn)
        a1.set_title(ctitle,fontsize=11)
        a1.set_xlabel('Fraction between lowest & highest rank',fontsize=8)
        a1.set_ylabel('Estimated closest-member histogram\nfraction from beta distribution',fontsize=8)
        a1.set_xlim(0,1)
        a1.set_ylim(0.0001,0.1)
        a1.set_yscale('log')
        a1.grid(color='Gray',lw=0.2,linestyle='--')
        a1.plot(x1000,chist[idx,:],color='Red',lw=1.2,label=r'closest-member histogram from fitted $\beta$')
        a1.plot(x1000,chist_lowess[idx,:],color='RoyalBlue',lw=1.2,label='Lowess smoothed parameters')
        if i==2: a1.legend(loc=0)

    plot_title = 'beta_closest_member_histogram_'+ctype+'_'+cmonthy+'_'+clead+'.png'
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

    a1.plot(x, beta_cdf[1,:],color='Red',linewidth=1.2,label='0.2 mm')
    a1.plot(x, beta_cdf[10,:],color='RoyalBlue',linewidth=1.2,label='2 mm')
    a1.plot(x, beta_cdf[50,:],color='LimeGreen',linewidth=1.2,linestyle='-.',label='10 mm')
    a1.plot(x, beta_cdf[100,:],color='Gray',linewidth=1.2,linestyle='--',label='20 mm')
    a1.plot(x, beta_cdf[170,:],color='Black',linewidth=1.2,label='34 mm')
    a1.plot(x, beta_cdf[250,:],color='Purple',linewidth=1.2,linestyle='--',label='50 mm')
    a1.legend(loc=0)

    plot_title = 'beta_closest_member_histogram_CDF_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    
    
    fig = plt.figure(figsize=(6.,6.5))

    for idx in range(251):
        beta_cdf[idx,:] = stats.beta.cdf(x,alpha_hat[idx], beta_hat[idx])
        for j in range(1000):
            chist[idx,j] = beta_cdf[idx,j+1] - beta_cdf[idx,j]

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

    a1.plot(x, beta_cdf[50,:],color='Black',linewidth=1.2,label='10 mm')
    a1.plot(x, beta_cdf[100,:],color='Red',linewidth=1.2,label='20 mm')
    a1.plot(x, beta_cdf[130,:],color='RoyalBlue',linewidth=1.2,label='26 mm')
    a1.plot(x, beta_cdf[170,:],color='LimeGreen',linewidth=1.2,linestyle='-.',label='34 mm')
    a1.plot(x, beta_cdf[210,:],color='Gray',linewidth=1.2,linestyle='--',label='42 mm')
    a1.plot(x, beta_cdf[250,:],color='Purple',linewidth=1.2,linestyle='--',label='50 mm')
 
    a1.legend(loc=0)

    plot_title = 'beta_closest_member_histogram_part2_CDF_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    
    
    fig = plt.figure(figsize=(6.,6.5))
    ctitle = 'LOWESS closest-member histogram CDFs for\n'+cmonthy+\
        ', lead = +'+clead+' h'
    axlocn = [0.11,0.11,0.85,0.79]

    a1 = fig.add_axes(axlocn)
    a1.set_title(ctitle,fontsize=17)
    a1.set_xlabel('Fraction between lowest and highest rank',fontsize=14)
    a1.set_ylabel('Cumulative probability',fontsize=14)
    a1.set_xlim(0,1)
    a1.set_ylim(0,1)
    a1.grid(color='Gray',lw=0.2,linestyle='--')

    a1.plot(x, beta_cdf_lowess[50,:],color='Black',linewidth=1.2,label='10 mm')
    a1.plot(x, beta_cdf_lowess[100,:],color='Red',linewidth=1.2,label='20 mm')
    a1.plot(x, beta_cdf_lowess[130,:],color='RoyalBlue',linewidth=1.2,label='26 mm')
    a1.plot(x, beta_cdf_lowess[170,:],color='LimeGreen',linewidth=1.2,linestyle='-.',label='34 mm')
    a1.plot(x, beta_cdf_lowess[210,:],color='Gray',linewidth=1.2,linestyle='--',label='42 mm')
    a1.plot(x, beta_cdf_lowess[250,:],color='Purple',linewidth=1.2,linestyle='--',label='50 mm')
 
    a1.legend(loc=0)

    plot_title = 'beta_closest_member_histogram_part2_lowess_CDF_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    
    
    
    fig = plt.figure(figsize=(6.,8.))


    ctitle = r'Beta dist. parameter $\hat{\alpha}$ of closest-ranked member,'+'\n'+cmonthy+\
        ', lead = +'+clead+' h'
    axlocn = [0.11,0.57,0.85,0.36]

    a1 = fig.add_axes(axlocn)
    a1.set_title(ctitle,fontsize=13)
    a1.set_xlabel('Ensemble-mean precipitation amount (mm)',fontsize=11)
    a1.set_ylabel(r'$\hat{\alpha}$',fontsize=11)
    a1.set_xlim(0,50)
    a1.set_ylim(0,1.5)
    a1.grid(color='Gray',lw=0.2,linestyle='--')
    #for iresa in range(nresa):
    #    a1.plot(precip_amt,alpha_hat_resamp[iresa,:],color='Red',linewidth=0.05 )
    #    a1.plot(precip_amt,beta_hat_resamp[iresa,:],color='RoyalBlue',linewidth=0.05 )
    a1.fill_between(precip_amt, alpha_hat_95, alpha_hat_05, color='LightPink',zorder=0,alpha=0.6) 
    a1.fill_between(precip_amt, alpha_hat_75, alpha_hat_25, color='LightCoral',zorder=1,alpha=0.6)      
    a1.plot(precip_amt, alpha_hat,color='Red',linewidth=1.5,label=r'$\hat{\alpha}$',zorder=2)
    a1.plot(precip_amt, alpha_fitted,color='Red',linestyle=':',\
        linewidth=1.5,label=r'LOWESS $\hat{\alpha}$')
    a1.legend(loc=0)
    
    
    ctitle = r'Beta dist. parameter $\hat{\beta}$ of closest-ranked member,'+'\n'+cmonthy+\
        ', lead = +'+clead+' h'    
    axlocn = [0.11,0.07,0.85,0.36]

    a1 = fig.add_axes(axlocn)
    a1.set_title(ctitle,fontsize=13)
    a1.set_xlabel('Ensemble-mean precipitation amount (mm)',fontsize=11)
    a1.set_ylabel(r'$\hat{\beta}$',fontsize=11)
    a1.set_xlim(0,50)
    a1.set_ylim(0,1.5)
    a1.grid(color='Gray',lw=0.2,linestyle='--')
    
    #for iresa in range(nresa):
    #    a1.plot(precip_amt,beta_hat_resamp[iresa,:],color='RoyalBlue',linewidth=0.05 )
            
    a1.fill_between(precip_amt, beta_hat_95, beta_hat_05, color='#cceeff',zorder=0,alpha=0.6) 
    a1.fill_between(precip_amt, beta_hat_75, beta_hat_25, color='#66ccff',zorder=1,alpha=0.6) 
    a1.plot(precip_amt, beta_hat,color='#0088cc',linewidth=1.5,label=r'$\hat{\beta}$',zorder=2)
    a1.plot(precip_amt, beta_fitted,color='#0088cc',linestyle=':',\
        linewidth=1.5,label=r'LOWESS $\hat{\beta}$')
    a1.legend(loc=0)
    
    plot_title = 'beta_parameters_'+ctype+'_'+cmonthy+'_'+clead+'.png'
    fig.savefig(plot_title, dpi=400)
    print ('saving plot to file = ',plot_title)






# ---- make plots of the closest-member histogram exceedance function

closest_hist_exceedance = np.copy(chist_lowess)
for idx in range(251):
    for j in range(1000):
        closest_hist_exceedance[idx,j] = np.sum(chist_lowess[idx,1000-j:1000])


fig = plt.figure(figsize=(6.,6.5))
rcParams['legend.fontsize']='large'
ctitle = 'Beta-distribution exceedance distribution functions\nfor '+cmonthy+\
    ', lead = +'+clead+' h'
axlocn = [0.11,0.11,0.85,0.79]

a1 = fig.add_axes(axlocn)
a1.set_title(ctitle,fontsize=15)
a1.set_xlabel('Fraction of members exceeding the threshold',fontsize=14)
a1.set_ylabel('Exceedance probability',fontsize=14)
a1.set_xlim(0,1)
a1.set_ylim(0,1)
a1.grid(color='Gray',lw=0.2,linestyle='--')

a1.plot(np.real(range(1,1001))/float(1000), \
    closest_hist_exceedance[50,:],color='Red',linewidth=1.2,\
    label=r'$\overline{\tilde{x}}^f_i \simeq$ 10 mm')
a1.plot(np.real(range(1,1001))/float(1000), \
    closest_hist_exceedance[100,:],color='RoyalBlue',\
    linewidth=1.2,label=r'$\overline{\tilde{x}}^f_i \simeq$ 20 mm')
a1.plot(np.real(range(1,1001))/float(1000), \
    closest_hist_exceedance[130,:],color='LimeGreen',\
    linewidth=1.2,linestyle='-.',label=r'$\overline{\tilde{x}}^f_i \simeq$ 26 mm')
a1.plot(np.real(range(1,1001))/float(1000), \
    closest_hist_exceedance[170,:],color='Gray',\
    linewidth=1.2,linestyle='--',label=r'$\overline{\tilde{x}}^f_i \simeq$ 34 mm')
a1.plot(np.real(range(1,1001))/float(1000), \
    closest_hist_exceedance[210,:],color='Black',\
    linewidth=1.2,label=r'$\overline{\tilde{x}}^f_i \simeq $ 42 mm')
a1.plot(np.real(range(1,1001))/float(1000), \
    closest_hist_exceedance[250,:],color='Orchid',\
    linewidth=1.2,linestyle='--',label=r'$\overline{\tilde{x}}^f_i \simeq$ 50 mm')

a1.legend(loc=0)

plot_title = 'closest_member_histogram_beta_EDF_'+ctype+'_'+cmonthy+'_'+clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
