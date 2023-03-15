# python define_important_casedates.py warm/cold ncases cexclude
import _pickle as cPickle
import sys, os
import numpy as np
import math
from matplotlib.path import Path
from shapely.geometry import Polygon

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import shapefile

# ---- input the season and total number of reforecast cases from command line

cseason = sys.argv[1] # warm or cool
ctotal_ncases = sys.argv[2] 
cexclude_hucs = sys.argv[3] # if 1, then exclude HUCs 9,13,16

# --- per feedback from Sunghee, provide ability to save data with or
#     without excluding HUCs

if cexclude_hucs == '1':
    huc_use = [1,1,1,1,1, 1,1,1,0,1, 1,1,0,1,1, 0,1,1,1]
    cexclude = '_no9_13_16'
else:
    huc_use = [1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1, 1,1,1,1]
    cexclude = '_allhucs'
    
total_ncases_out = int(ctotal_ncases)
print ('total_ncases_out = ', total_ncases_out)
print ('cexclude_hucs = ',cexclude_hucs)
print ('cexclude = ', cexclude)
print ('np.sum(huc_use) = ', np.sum(huc_use))

# ---- set the parameters for how to divide up cases, some selected by the
#      average mean precip over the HUC; others by the precip max at the top
#      20 grid points within the HUC, and the rest by CONUS-averaged mean precip

nhucs = 19 # 18 HUC-2's + whole conus [#19]
fraction_by_each_huc_bymean = \
    0.3 / float(np.sum(huc_use)-1) 
    # the fraction of cases by HUC mean 
fraction_by_each_huc_bymax20 = \
    0.6 / float(np.sum(huc_use)-1) 
    # the fraction of cases excl. CONUS
    # by HUC max at top 20 points within that HUC
fraction_conuswide = 1.0 - \
    fraction_by_each_huc_bymean*(np.sum(huc_use)-1) -  \
    fraction_by_each_huc_bymax20*(np.sum(huc_use)-1) \
    # the fraction of cases by CONUS mean
print ('fractions bymean bymax20 conus = ', \
    fraction_by_each_huc_bymean, fraction_by_each_huc_bymax20,\
    fraction_conuswide)
cases_for_each_huc_bymean = int(total_ncases_out*fraction_by_each_huc_bymean)
cases_for_each_huc_bymax20 = int(total_ncases_out*fraction_by_each_huc_bymax20)
cases_for_conus = total_ncases_out - (nhucs-1)*cases_for_each_huc_bymean - \
    (nhucs-1)*cases_for_each_huc_bymax20
print ('cases_for_each_huc_bymean, cases_for_each_huc_bymax20, cases_for_conus = ',\
     cases_for_each_huc_bymean, cases_for_each_huc_bymax20, cases_for_conus)
print ('total number of cases = ',\
    (np.sum(huc_use)-1)*cases_for_each_huc_bymean + \
    (np.sum(huc_use)-1)*cases_for_each_huc_bymax20 + \
    cases_for_conus)

# ---- do some preprocessing, loading in the data and getting the indices associated
#      with sorted data.  The last HUC is for CONUS.

for ihuc in range(19):
    
    # --- read the yyyymmddhh of the initial date, the mean precip in this HUC over the 
    #     240 h period, and the mean of the precip at the grid points with 20 largest 
    #     values.   
    
    infile = cseason + '_precip_stats_huc2number'+str(ihuc+1)+'.cPick'
    print ('reading from ', infile)
    inf = open(infile, 'rb')
    yyyymmddhh = cPickle.load(inf)
    mean_precip = cPickle.load(inf)
    max20_precip = cPickle.load(inf)
    max20_meanlon = cPickle.load(inf)
    max20_meanlat = cPickle.load(inf)    
    inf.close()
    
    # ---- Apply argsort to get the associated indices of sorted low to high values.
    #      Initialize arrays also.
        
    argsort_mean = np.argsort(mean_precip)
    argsort_max20 = np.argsort(max20_precip)
    if ihuc == 0: # first time through
        ndates = len(yyyymmddhh)
        argsort_mean_allhucs = np.zeros((19,ndates), dtype=np.int32)
        argsort_max20_allhucs = np.zeros((19,ndates), dtype=np.int32)
        weighting_mean_allhucs = np.zeros((19,ndates), dtype=np.float64)
        weighting_max20_allhucs = np.zeros((19,ndates), dtype=np.float64)
        lon_max20_allhucs = np.zeros((19,ndates), dtype=np.float64)
        lat_max20_allhucs = np.zeros((19,ndates), dtype=np.float64)
    argsort_mean_allhucs[ihuc,:] = argsort_mean[:]
    argsort_max20_allhucs[ihuc,:] = argsort_max20[:]
    
    # ---- Develop an initial weighting estimate based on the precipitation
    #      value divided by the maximum precipitation value
    
    weighting_mean_allhucs[ihuc,:] = mean_precip[:] / mean_precip[argsort_mean[-1]]
    weighting_max20_allhucs[ihuc,:] = max20_precip[:] / max20_precip[argsort_max20[-1]]
    lon_max20_allhucs[ihuc,:] = max20_meanlon[:]
    lat_max20_allhucs[ihuc,:] = max20_meanlat[:]

# ---- develop a set of cases across US based on ***ensemble mean*** in each real HUC.
#      Process the full CONUS (stored in the last HUC index) later.

print ('******** beginning loop thru cases')
casedates = np.zeros((ndates), dtype=np.int32)
which_huc = np.zeros((ndates), dtype=np.int16)
meanlon = -99.99*np.ones((ndates), dtype=np.float32) # used to ID center lon when selected by max20
meanlat = -99.99*np.ones((ndates), dtype=np.float32) # used to ID center lat when selected by max20
casenum = 0
for icase in range(cases_for_each_huc_bymean):
    print ('processing cases by mean', icase, cases_for_each_huc_bymean)
    for ihuc in range(18): # 18, so excl. CONUS
        
        if huc_use[ihuc] == 1:

            # ---- find the date with the maximum weight for mean precip.  Adjust
            #      weights of nearby cases downward somewhat so that they are less
            #      likely to be chosen so that we don't choose too many cases
            #      that are clustered around a few chosen dates.
        
            weights_mean = weighting_mean_allhucs[ihuc,:]
            #print ('len(weights_mean) = ', len(weights_mean))
            indices = np.argsort(weights_mean)
        
            if indices[-1] != 0:
                iminus = indices[-1]-1
            else:
                iminus = indices[-1]
            if indices[-1] != ndates-1:
                iplus = indices[-1]+1
            else:
                iplus = indices[-1]
        
            if indices[-1] > 1:
                iminus2 = indices[-1]-2
            else:
                iminus2 = indices[-1]
            if indices[-1] < ndates-2:
                iplus2 = indices[-1]+2
            else:
                iplus2 = indices[-1]
        
            casedates[indices[-1]] = 1 # choose the date with max weight
            which_huc[indices[-1]] = ihuc+1
        
            if ihuc == 17: print ('casedate, huc, precip, wt = ',\
                yyyymmddhh[indices[-1]], which_huc[indices[-1]], \
                mean_precip[indices[-1]], weights_mean[indices[-1]])
            weighting_mean_allhucs[:,indices[-1]] = 0.0 
                # zero out this date's weight so not chosen again
            weighting_mean_allhucs[:,iminus] = \
                weighting_mean_allhucs[:,iminus]*0.7 #de-emph the surr. dates
            weighting_mean_allhucs[:,iplus] = \
                weighting_mean_allhucs[:,iplus]*0.7
            weighting_max20_allhucs[:,indices[-1]] = 0.0 
                # zero out this date's weight so not chosen again
            weighting_max20_allhucs[:,iminus] = \
                weighting_max20_allhucs[:,iminus]*0.7 #de-emph the surr. dates
            weighting_max20_allhucs[:,iplus] = \
                weighting_max20_allhucs[:,iplus]*0.7
            weighting_max20_allhucs[:,iminus2] = \
                weighting_max20_allhucs[:,iminus2]*0.85 #de-emph the surr. dates
            weighting_max20_allhucs[:,iplus2] = \
                weighting_max20_allhucs[:,iplus2]*0.85
            casenum = casenum + 1
    
# ---- develop a set of cases across US based on ***ensemble max20 *** in each real HUC.
#      Process the full CONUS (stored in the last HUC index) later.

print ('******** beginning loop thru cases')
for icase in range(cases_for_each_huc_bymax20):
    print ('processing cases by max20', icase, cases_for_each_huc_bymax20)
    for ihuc in range(18):
        
        if huc_use[ihuc] == 1:
            
            # ---- find the date with the maximum weight for max20 precip.  Adjust
            #      weights of nearby cases downward somewhat so that they are less
            #      likely to be chosen.
    
            weights = weighting_max20_allhucs[ihuc,:]
            indices = np.argsort(weights)
    
            if indices[-1] != 0:
                iminus = indices[-1]-1
            else:
                iminus = indices[-1]
            if indices[-1] != ndates-1:
                iplus = indices[-1]+1
            else:
                iplus = indices[-1]
    
            if indices[-1] > 1:
                iminus2 = indices[-1]-2
            else:
                iminus2 = indices[-1]
            if indices[-1] < ndates-2:
                iplus2 = indices[-1]+2
            else:
                iplus2 = indices[-1]
            
            casedates[indices[-1]] = 1
            which_huc[indices[-1]] = ihuc+1
            meanlon[indices[-1]] = lon_max20_allhucs[ihuc,indices[-1]]
            meanlat[indices[-1]] = lat_max20_allhucs[ihuc,indices[-1]]
        
            if ihuc == 17: print ('casedate, huc, wt = ',yyyymmddhh[indices[-1]], \
                which_huc[indices[-1]], max20_precip[indices[-1]], weights[indices[-1]])
            weighting_mean_allhucs[:,indices[-1]] = 0.0 
                # zero out this date's weight so not chosen again
            weighting_mean_allhucs[:,iminus] = \
                weighting_mean_allhucs[:,iminus]*0.7 #de-emph the surr. dates
            weighting_mean_allhucs[:,iplus] = \
                weighting_mean_allhucs[:,iplus]*0.7
            weighting_max20_allhucs[:,indices[-1]] = 0.0 
                # zero out this date's weight so not chosen again
            weighting_max20_allhucs[:,iminus] = \
                weighting_max20_allhucs[:,iminus]*0.7 #de-emph the surr. dates
            weighting_max20_allhucs[:,iplus] = \
                weighting_max20_allhucs[:,iplus]*0.7
            weighting_max20_allhucs[:,iminus2] = \
                weighting_max20_allhucs[:,iminus2]*0.85 #de-emph the surr. dates
            weighting_max20_allhucs[:,iplus2] = \
                weighting_max20_allhucs[:,iplus2]*0.85
            casenum = casenum + 1
        

# ---- use the remaining dates to select cases with CONUS-wide impact that
#      haven't been selected already.
    
print ('******** cases_for_conus',cases_for_conus)
for icase in range(cases_for_conus):
    
    # ---- find the date with the maximum weight for mean precip.  Adjust
    #      weights of nearby cases downward somewhat so that they are less
    #      likely to be chosen so that we don't choose too many cases
    #      that are clustered around a few chosen dates.
        
    weights_mean = weighting_mean_allhucs[18,:]
    indices = np.argsort(weights_mean)
    
    if indices[-1] != 0:
        iminus = indices[-1]-1
    else:
        iminus = indices[-1]
        
    if indices[-1] != ndates-1:
        iplus = indices[-1]+1
    else:
        iplus = indices[-1]
    
    if indices[-1] > 1:
        iminus2 = indices[-1]-2
    else:
        iminus2 = indices[-1]
    if indices[-1] < ndates-2:
        iplus2 = indices[-1]+2
    else:
        iplus2 = indices[-1]
    
    casedates[indices[-1]] = 1 # choose the date with max weight
    which_huc[indices[-1]] = 19
    print ('casedate, huc, precip, wt = ',yyyymmddhh[indices[-1]], \
        which_huc[indices[-1]], mean_precip[indices[-1]], weights_mean[indices[-1]])
    weighting_mean_allhucs[:,indices[-1]] = 0.0 # zero out this date's weight so not chosen again
    weighting_mean_allhucs[:,iminus] = weighting_mean_allhucs[:,iminus]*0.7 #de-emph the surr. dates
    weighting_mean_allhucs[:,iplus] = weighting_mean_allhucs[:,iplus]*0.7 
    weighting_mean_allhucs[:,iminus2] = weighting_mean_allhucs[:,iminus2]*0.85 #de-emph the surr. dates
    weighting_mean_allhucs[:,iplus2] = weighting_mean_allhucs[:,iplus2]*0.85     

# ---- save list of chosen cases.  Also save lon/lat of mean of max20 grid points

outfile = 'case_list_'+cseason+'season_ncases'+ctotal_ncases+cexclude+'.txt'
print ('writing case dates to ', outfile)
ouf = open(outfile,'w')
print ('casedates[0:30] = ', casedates[0:30])
print ('yyyymmddhh[0:30] = ', yyyymmddhh[0:30])
print ('which_huc[0:30] = ', which_huc[0:30])
for idate in range(ndates):
    if casedates[idate] == 1:
        print(yyyymmddhh[idate], which_huc[idate], meanlon[idate], meanlat[idate], file=ouf)
ouf.close()