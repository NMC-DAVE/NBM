"""
generate_error_stats_decayavg.py 

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
import matplotlib.pyplot as plt
from os import path
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'
rcParams["contour.negative_linestyle"]='solid'


# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
calpha = sys.argv[2]
makeplots = False

ilead = int(clead)
datadir_reanl = '/Users/Tom/python/ecmwf/'
datadir = '/Users/Tom/python/ecmwf/'
cvariable = '2t'
date_anal_start = dateshift('2019010100',ilead)
date_anal_end = dateshift('2019123100',-ilead)
date_list_anal = daterange(date_anal_start,date_anal_end,24)
#date_list_anal = daterange(datestart,'2019013100',24)
ndates = len(date_list_anal)
rmse_kf = ma.zeros((ndates), dtype=np.float32)
rmse_decayavg = ma.zeros((ndates), dtype=np.float32)
rmse_raw = ma.zeros((ndates), dtype=np.float32)
date_list_fcst = []
date_list_bcorr = []
frac2019 = np.zeros((ndates), dtype=np.float32)
for idate in range(ndates):
    date_list_fcst.append(dateshift(date_list_anal[idate],ilead)) # initial times of fcst
    
    
# ---- read in the time-dependent ERA5 climatology of t2m

infilename = datadir_reanl+'t2m_climo_daily_era5_halfdegree.cPick'
print (infilename)
inf = open(infilename, 'rb')
climo_yearly = cPickle.load(inf)
print ('shape climo_yearly = ', np.shape(climo_yearly))
latsa_halfdegree = cPickle.load(inf)
lonsa_halfdegree = cPickle.load(inf)
nlats, nlons = np.shape(latsa_halfdegree)
inf.close()    

for idate, datea in enumerate(date_list_anal):
    
    datef = date_list_fcst[idate]
    if datea == '2019010100': dstart = idate
    #print ('------ processing analysis, forecast dates = ', datea, datef)
    yyyy,mm,dd,hh = splitdate(datea) 
    doy = dayofyear(yyyy,mm,dd)
    fracyear = doy/365.

    # ---- read the ECMWF ERA5 reanalysis at the date of the forecast.
    
    infile = datadir_reanl + 't2m_era5_halfdegree_'+datef+'.cPick'
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats)
        npts = nlats*nlons 
    inf.close()
    
    # ---- read the control forecast at this lead time and initial date
 
    infile = datadir + cvariable+'_'+datea+'_f'+clead+'.grib2' 
    fexist = path.exists(infile)
    if fexist:
         
        grbfile = pygrib.open(infile) 
        grb = grbfile.select()[0] 
        forecast = grb.values
        grbfile.close()
    
        # ---- read the decaying average bias corrections estimate for this date.
    
        infilename = datadir + 'bias_decayavg_alpha'+calpha+'_'+datea+'_f'+clead+'.cPick'
        print (infilename)
        inf = open(infilename, 'rb')
        bias_decayavg = cPickle.load(inf)
        print ('max, min bias_decayavg = ', ma.max(bias_decayavg), ma.min(bias_decayavg))
        inf.close()
    
        # ---- read in the Kalman filter regression coefficients
    
        infilename = datadir + 'bias_KFregression_'+datea+'_f'+clead+'.cPick'
        inf = open(infilename, 'rb')
        regr_slope = cPickle.load(inf)
        regr_intercept = cPickle.load(inf)
        print ('max, min regr_slope = ', ma.max(regr_slope), ma.min(regr_slope))
        print ('max, min regr_intercept = ', ma.max(regr_intercept), ma.min(regr_intercept))
        inf.close
        forecast_dev = forecast[:,:] - climo_yearly[doy,:,:]
        forecast_kf = forecast_dev*regr_slope + regr_intercept + climo_yearly[doy,:,:]
    
        # ---- change the forecast to deviation from climatology, and pre
    
        rmse_kf[idate] = ma.sqrt(ma.sum((forecast_kf-analysis)**2)/(npts-1.))
        rmse_raw[idate] = ma.sqrt(ma.sum((forecast-analysis)**2)/(npts-1.))
        rmse_decayavg[idate] = ma.sqrt(ma.sum(((forecast-bias_decayavg)-analysis)**2)/(npts-1.))
    else:
        rmse_kf[idate] = ma.masked
        rmse_raw[idate] = ma.masked
        rmse_decayavg[idate] = ma.masked
        

rmse_kf_mean = ma.mean(rmse_kf)
rmse_raw_mean = ma.mean(rmse_raw)
rmse_decay_mean = ma.mean(rmse_decayavg)

print ('alpha, rmse raw decay kf = ',calpha, rmse_raw_mean, rmse_decay_mean, rmse_kf_mean )
    
    
    
    
            

