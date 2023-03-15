"""
using ECMWF forecast and ERA5 verification data, implement Kalman-filter
type regression bias correction approach to estimate forecast bias.   Grid covers
CONUS, approximately.  Coded by Tom Hamill, May-Jun 2020.
tom.hamill@noaa.gov
"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from os import path
import statsmodels.api as sm
from scipy.optimize import curve_fit
import scipy.stats.mstats as stats
import matplotlib.pyplot as plt

# =====================================================================

def gamma_exponential(data_array, rho_horiz, vdconst, gamma):
    term = np.sqrt(data_array[0,:]**2/rho_horiz**2 + \
        data_array[1,:]**2/vdconst**2)
    expterm_gamma = np.exp(-term**(gamma)) 
    print ('min, max expterm_gamma = ', np.min(expterm_gamma), np.max(expterm_gamma))     
    return expterm_gamma   
    
# =====================================================================
    
def estimate_variance(forecast_deviation):
    season_variance = np.var(forecast_deviation, axis=0)
    return season_variance

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
ilead = int(clead)
ndstart = ilead // 24
datadir_reanl = '/Users/Tom/python/ecmwf/'
datadir = '/Volumes/Backup Plus/ecmwf/'
makeplots = True

cvariable = '2t'
tally_statistics = False
#datestart = dateshift('2018110100',ilead)
dateend = dateshift('2018123100',-ilead)

date_list_anal_warmseason = daterange('2018040100','2018093000',24) # initial time of the current forecast
date_list_anal_coolseason = daterange('2018010100','2018033100',24)  + \
   daterange('2018100100',dateend,24) # initial time of the current forecast


ndates_warmseason = len(date_list_anal_warmseason)
ndates_coolseason = len(date_list_anal_coolseason)
date_list_fcst_warmseason = []
date_list_fcst_coolseason = []
for idate in range(ndates_warmseason):
    date_list_fcst_warmseason.append(dateshift(date_list_anal_warmseason[idate],ilead)) # initial times of fcst
for idate in range(ndates_coolseason):
    date_list_fcst_coolseason.append(dateshift(date_list_anal_coolseason[idate],ilead)) # initial times of fcst

# ---- read in the terrain elevation

terrain_elevation = np.load('terrain_halfdeg.npy')

# ---- loop over warm-season dates and update bias estimates

for idate, datea in enumerate(date_list_anal_warmseason):
    
    datef = date_list_fcst_warmseason[idate]

    # ---- read the ECMWF ERA5 reanalysis at valid at the forecast date.
    
    infile = datadir_reanl + 't2m_era5_halfdegree_'+datef+'.cPick'
    #print (infile)
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats)
        print (nlats, nlons)
        analyzed_yearly_warmseason = ma.zeros((ndates_warmseason,nlats,nlons), dtype=np.float32) 
        forecast_yearly_warmseason = ma.zeros((ndates_warmseason,nlats,nlons), dtype=np.float32)    
    inf.close()
    analyzed_yearly_warmseason[idate,:,:] = analysis[:,:] 
    
    # ---- read the control forecast at this lead time and initial date

    infile = datadir + cvariable+'_'+datea+'_f'+clead+'.grib2'  
    fexist = path.exists(infile)
    print (infile, fexist)
    if fexist == True:
        #print (infile)
        grbfile = pygrib.open(infile) 
        grb = grbfile.select()[0] 
        forecast = grb.values
        grbfile.close()
        forecast_yearly_warmseason[idate,:,:] = forecast[:,:] 
    else:
        forecast_yearly_warmseason[idate,:,:] = ma.masked
  
# ---- loop over cool-season dates 
      
for idate, datea in enumerate(date_list_anal_coolseason):
    
    datef = date_list_fcst_coolseason[idate]

    # ---- read the ECMWF ERA5 reanalysis at valid at the forecast date.
    
    infile = datadir_reanl + 't2m_era5_halfdegree_'+datef+'.cPick'
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats)
        print (nlats, nlons)
        analyzed_yearly_coolseason = ma.zeros((ndates_warmseason,nlats,nlons), dtype=np.float32) 
        forecast_yearly_coolseason = ma.zeros((ndates_warmseason,nlats,nlons), dtype=np.float32)    
    inf.close()
    analyzed_yearly_coolseason[idate,:,:] = analysis[:,:] 
    
    # ---- read the control forecast at this lead time and initial date

    infile = datadir + cvariable+'_'+datea+'_f'+clead+'.grib2'  
    fexist = path.exists(infile)
    print (infile, fexist)
    if fexist == True:
        #print (infile)
        grbfile = pygrib.open(infile) 
        grb = grbfile.select()[0] 
        forecast = grb.values
        grbfile.close()
        forecast_yearly_coolseason[idate,:,:] = forecast[:,:] 
    else:
        forecast_yearly_coolseason[idate,:,:] = ma.masked
    
# ---- for each grid point, first produce a linear regression in order 
#      to get a sense of the standard error of the regression coefficients.
    
forecast_deviation_warmseason = forecast_yearly_warmseason - analyzed_yearly_warmseason
forecast_deviation_coolseason = forecast_yearly_coolseason - analyzed_yearly_coolseason

# ---- estimate warm-season variance and cold-season variance.

warm_season_variance = estimate_variance(forecast_deviation_warmseason)
cool_season_variance = estimate_variance(forecast_deviation_coolseason)

# ---- populate arrays of time series correlation, horiz and vert distances, and then 

kmperdeg = 111.
kmperdeg_d2 = 111./2.
corr_list_warmseason = []
corr_list_coolseason = []
hdist_list = []
vdiff_list = []
for ix1 in range(nlons):
#for ix1 in [nlons//2]:
    print ('processing ix1 ',ix1, nlons)
    imin = max(0,ix1-10)
    imax = min(ix1+10,nlons)
    for jy1 in range(nlats):
    #for jy1 in [nlats//2]:
        fd1c = forecast_deviation_coolseason[:,jy1,ix1]
        fd1w = forecast_deviation_warmseason[:,jy1,ix1]
        #print ('fd1c = ', fd1c[:])
        #print ('fd1w = ', fd1w[:])
        if ix1 == 0: 
            coslat1 = np.cos(lats[jy1,0]*3.14159/180.)
        jmin = max(0,jy1-10)
        jmax = min(jy1+10,nlats)    
        for ix2 in range(imin, imax):
            idiff = np.abs(ix1-ix2)
            for jy2 in range(jmin, jmax):
                fd2c = forecast_deviation_coolseason[:,jy2,ix2]
                fd2w = forecast_deviation_warmseason[:,jy2,ix2]
                jdiff = np.abs(jy1-jy2)
                if ix2 == 0: 
                    coslat2 = np.cos(lats[jy2,0]*3.14159/180.)
                cosavg = (coslat1+coslat2)/2.
                cc = np.corrcoef(fd1c,fd2c)
                cw = np.corrcoef(fd1w,fd2w)
                hdist = np.sqrt((idiff*kmperdeg_d2)**2 + (cosavg*jdiff*kmperdeg_d2)**2)
                r = np.random.uniform(low=0.95,high=1.05)
                hdist = hdist*r
                vdiff = np.abs(terrain_elevation[jy1,ix1] - terrain_elevation[jy2,ix2])
                r = np.random.uniform(low=0.95,high=1.05)
                vdiff = vdiff*r
                #print ('ix2,jy2, hdist, vdiff = ', ix2,jy2, hdist, vdiff)
                #print ('corr cool warm = ', cc[0,1], cw[0,1])
                corr_list_warmseason.append(cw[0,1])
                corr_list_coolseason.append(cc[0,1])
                hdist_list.append(hdist)
                vdiff_list.append(vdiff)

#sys.exit()

# ---- move the data of interest to arrays and process.

corr_arr_warmseason = np.squeeze(np.array(corr_list_warmseason))
corr_arr_coolseason = np.squeeze(np.array(corr_list_coolseason))
print ('shape corr_arr_warmseason corr_arr_coolseason', \
    np.shape(corr_arr_warmseason), np.shape(corr_arr_coolseason))
hdist_arr = np.squeeze(np.array(hdist_list))
vdiff_arr = np.squeeze(np.array(vdiff_list))
ones = np.ones(len(corr_arr_warmseason), dtype=np.float32)


sample_var_coolseason = 1.0 - corr_arr_coolseason**2
sample_var_warmseason = 1.0 - corr_arr_warmseason**2
sample_var_coolseason = np.where(corr_arr_coolseason < 0.25, 1.0 - 0.0625*ones, sample_var_coolseason)
sample_var_warmseason = np.where(corr_arr_warmseason < 0.25, 1.0 - 0.0625*ones, sample_var_warmseason)
sample_var_coolseason = np.where(sample_var_warmseason < 0.05, 0.05*ones, sample_var_coolseason)
sample_var_warmseason = np.where(sample_var_warmseason < 0.05, 0.05*ones, sample_var_warmseason)
data_array = np.zeros((2,len(hdist_list)),dtype=np.float32)
data_array = np.zeros((2,len(hdist_list)),dtype=np.float32)
data_array[0,:] = hdist_arr[:]
data_array[1,:] = vdiff_arr[:]
print ('min, max sample_var_coolseason = ', \
    np.min(sample_var_coolseason), np.max(sample_var_coolseason))
print ('min, max sample_var_warmseason = ', \
    np.min(sample_var_warmseason), np.max(sample_var_warmseason))
print ('min, max corr_arr_coolseason', \
    np.min(corr_arr_coolseason), np.max(corr_arr_coolseason))
print ('min, max corr_arr_warmseason', \
    np.min(corr_arr_warmseason), np.max(corr_arr_warmseason))
print ('number of masked warm, cold corr ', \
    ma.count_masked(corr_arr_warmseason), \
    ma.count_masked(corr_arr_coolseason))   
        
# ------ estimate covariance parameters warm season
        
rho_horiz = 200.  # set to conservative values determined from Jul data.
vdconst = 200.
gamma = 1.8  
            
print ('processing warm season ')
popt, pcov = curve_fit(gamma_exponential, \
    data_array, corr_arr_warmseason, \
    p0=(rho_horiz, vdconst, gamma), \
    check_finite=True,  method='trf', diff_step=0.000001, \
    #sigma = sample_var_warmseason, absolute_sigma=False, \
    sigma = ones, absolute_sigma=False, \
    bounds = ([0.0,0.0,1.0], [2500.,5000.,2.0]))    
    
rho_horiz_warmseason = popt[0]
vdconst_warmseason = popt[1]
gamma_warmseason = popt[2]

# ------ estimate covariance parameters cool season

rho_horiz = 200.  # set to conservative values determined from Jul data.
vdconst = 200.
gamma = 1.8  
            
print ('processing cool season')
popt, pcov = curve_fit(gamma_exponential, \
    data_array, corr_arr_coolseason, \
    p0=(rho_horiz, vdconst, gamma), \
    check_finite=True,  method='trf', diff_step=0.000001, \
    #sigma = sample_var_coolseason, absolute_sigma=False, \
    sigma = ones, absolute_sigma=False, \
    bounds = ([0.0,0.0,1.0], [2500.,5000.,2.0]))    
    
rho_horiz_coolseason = popt[0]
vdconst_coolseason = popt[1]
gamma_coolseason = popt[2]
    
#rho_horiz_sprd = np.sqrt(pcov[0,0])
#vdconst_sprd = np.sqrt(pcov[1,1])
#gamma_sprd = np.sqrt(pcov[2,2])
        
print ('warm season: rho_horiz, vdconst, gamma = ', \
    rho_horiz_warmseason,  vdconst_warmseason, \
    gamma_warmseason)
print ('cool season: rho_horiz, vdconst, gamma = ', \
    rho_horiz_coolseason,  vdconst_coolseason, \
    gamma_coolseason)

# ----   make plot of horizontal covariance functions, 
#        now for gamma-exponential with the correct data
        
if makeplots == True:

    plot_title = 'horiz_error_corr_lead'+clead+'.pdf'
    fig = plt.figure(figsize=(6.5, 9.))
    dist = np.arange(0,1001,1)

    a1 = fig.add_axes([.14,.56,.82,.4])
    title = '(a) '+clead+'-h cool-season forecast horizontal error correlation'
    a1.set_title(title,fontsize=12)
    term = dist/rho_horiz_coolseason
    corrfn = np.exp(-term**gamma_coolseason) # gamma-exponential
    a1.scatter(hdist_arr[0:-1:50], corr_arr_coolseason[0:-1:50], marker='o', s=0.1, color='Gray',zorder=8)
    a1.plot(range(1001),corrfn,'-',color='Red',linewidth=3,zorder=11)
    a1.plot([0,1000],[0,0],'-',color='Black',linewidth=1,zorder=10)
    a1.set_xlim(0,1000)
    a1.set_ylim(-0.4,1.02)
    a1.set_xlabel('Euclidean distance (km)',fontsize=13)
    a1.set_ylabel('Correlation',fontsize=13)
    a1.grid(color='LightGray',linewidth=0.3)
    
    a1 = fig.add_axes([.14,.06,.82,.4])
    title = '(b) '+clead+'-h warm-season forecast horizontal error correlation'
    a1.set_title(title,fontsize=12)
    term = dist/rho_horiz_warmseason
    corrfn = np.exp(-term**gamma_warmseason) # gamma-exponential
    a1.scatter(hdist_arr[0:-1:50], corr_arr_warmseason[0:-1:50], marker='o', s=0.1, color='Gray',zorder=8)
    a1.plot(range(1001),corrfn,'-',color='Red',linewidth=3,zorder=11)
    a1.plot([0,1000],[0,0],'-',color='Black',linewidth=1,zorder=10)
    a1.set_xlim(0,1000)
    a1.set_ylim(-0.4,1.02)
    a1.set_xlabel('Euclidean distance (km)',fontsize=13)
    a1.set_ylabel('Correlation',fontsize=13)
    a1.grid(color='LightGray',linewidth=0.3)
    
    print ('saving plot to file = ',plot_title)
    plt.savefig(plot_title)
    plt.close()
    print ('Plot done')
        
# ---- save data assimilation parameter estimates to file
        
outfile = 'gamma_exponential_horiz_vert_lead'+clead+'.pydat'
ouf = open(outfile,'wb')
cPickle.dump(rho_horiz_warmseason, ouf)
cPickle.dump(vdconst_warmseason, ouf)
cPickle.dump(gamma_warmseason, ouf)
cPickle.dump(warm_season_variance, ouf) 
cPickle.dump(rho_horiz_coolseason, ouf)
cPickle.dump(vdconst_coolseason, ouf)
cPickle.dump(gamma_coolseason, ouf)
cPickle.dump(cool_season_variance, ouf)      
ouf.close()

print ("\a")
        
