"""

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

def initialize(ilead, alpha):
    
    # ---- initialize various
    
    R = 0.25**2  # 1.0**2 # observation-error variance
    B = (0.1 + (ilead/96.))**2 # (0.5 + (ilead/48.))**2 # estimated forecast random-error variance, which grows with lead time
    Bbeta = np.zeros((5,5), dtype=np.float32) # simplified error covariance for seas'ly depdt coefficients
    Bbeta[0,0] = alpha *0.25 
    Bbeta[1,1] = alpha *0.25
    Bbeta[2,2] = alpha *0.25
    Bbeta[3,3] = alpha *0.25
    Bbeta[4,4] = alpha *0.25
    
    KF_betahat = np.zeros((5), dtype=np.float32)
   
    return R, B, Bbeta, KF_betahat

# =====================================================================

def cosfac_sinfac (date):
    
    # ---- compute cos, sin of Julian day/365.
    
    yyyy,mm,dd,hh = splitdate(date) 
    doy = dayofyear(yyyy,mm,dd)
    fracyear = doy/365.
    fac = 2.*3.14159*(np.real(doy)/365.)
    cosfac = np.cos(fac)
    sinfac = np.sin(fac)
    fac = 4.*3.14159*(np.real(doy)/365.)
    cos2fac = np.cos(fac)
    sin2fac = np.sin(fac)
    return cosfac, sinfac, cos2fac, sin2fac, fracyear
    
# =====================================================================

def decayavg_bias(alpha, obs, forecast, bias):
    
    # ---- compute the bog-standard decaying average bias correction estimate
       
    bias = (1-alpha)*bias + alpha*(forecast-obs)
    return bias

# =====================================================================

def seasonalKFbias(cosfac, sinfac, cos2fac, sin2fac, Bbeta, B, R, \
    KF_betahat, obsy, fcsty, biasy):
        
    # ---- estimate the Kalman gain for the bias correction.
    
    L = np.array([1.0, sinfac, cosfac, sin2fac, cos2fac])
    BbetaLT = np.matmul(Bbeta[:,:], np.transpose(L))
    LBbetaLT = np.matmul(L,BbetaLT)
    LBbetaLT_plus_B_plus_R = LBbetaLT + B + R
    LBbetaLT_plus_B_plus_R_inv = 1.0 / LBbetaLT_plus_B_plus_R
    Kfgain_beta = BbetaLT * LBbetaLT_plus_B_plus_R_inv
    print ('np.shape(Kfgain_beta) = ', np.shape(Kfgain_beta))
    print ('obsy, fcsty, biasy = ', obsy, fcsty, biasy)

    # ---- update bias estimate with new data

    for i in range(5):
        KF_betahat[i] = KF_betahat[i] - \
            Kfgain_beta[i]*(obsy - (fcsty - biasy))      
    biasy = L[0]*KF_betahat[0] + \
        L[1]*KF_betahat[1] + L[2]*KF_betahat[2] + \
        L[3]*KF_betahat[3] + L[4]*KF_betahat[4]
    return KF_betahat, biasy

# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
clonb = sys.argv[2]
clatb = sys.argv[3]
clone = sys.argv[4]
clate = sys.argv[5]
rlonb = float(clonb)
rlatb = float(clatb)
rlone = float(clone)
rlate = float(clate)
alpha = 0.02
calpha = str(alpha)

ilead = int(clead)
datadir = '/Users/Tom/python/ecmwf/'
cvariable = '2t'
datestart = dateshift('2018110100',ilead)
date_list_anal = daterange(datestart,'2019123100',24)
ndates = len(date_list_anal)
date_list_fcst = []
for idate in range(ndates):
    date_list_fcst.append(dateshift(date_list_anal[idate],-ilead)) # initial times of fcst

forecast_box = np.zeros((ndates), dtype=np.float32) 
analysis_box = np.zeros((ndates), dtype=np.float32) 
bias_decayavg = np.zeros((ndates), dtype=np.float32)
bias_seasonalKF = np.zeros((ndates), dtype=np.float32)
frac2019 = np.zeros((ndates), dtype=np.float32)  
    
# ---- call initialization routine

R, B, Bbeta, KF_betahat = initialize(ilead, alpha)

for idate, datea in enumerate(date_list_anal):
    
    datef = date_list_fcst[idate]
    if datea == '2019010100': dstart = idate
    print ('------ processing analysis, forecast dates = ', datea, datef)

    # ---- read the ECMWF ERA5 reanalysis at this analysis date.
    
    infile = datadir + 't2m_era5_halfdegree_'+datea+'.cPick'
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats) 
        imin = np.argmin(np.abs(lons[0,:]-rlonb))
        imax = np.argmin(np.abs(lons[0,:]-rlone))
        jmin = np.argmin(np.abs(lats[:,0]-rlate))
        jmax = np.argmin(np.abs(lats[:,0]-rlatb))
        print ('imin, imax, jmin, jmax = ', imin, imax, jmin, jmax)
        #sys.exit()
    inf.close()
    
    # ---- read the ECMWF control forecast at this lead time and initial date
 
    infile = datadir + cvariable+'_'+datef+'_f'+clead+'.grib2'  
    grbfile = pygrib.open(infile) 
    grb = grbfile.select()[0] 
    forecast = grb.values
    grbfile.close()
    
    # ---- read the ERA5 analysis valid at this date.
    
    infilename = datadir+'t2m_era5_halfdegree_'+datea+'.cPick'
    inf = open(infilename, 'rb')
    obs = cPickle.load(inf)
    inf.close()    
    
    forecast_box[idate] = np.mean(forecast[jmin:jmax,imin:imax])
    analysis_box[idate] = np.mean(analysis[jmin:jmax,imin:imax])
    
    cosfac, sinfac, cos2fac, sin2fac, fracyear = cosfac_sinfac (datea)
    if int(datea[0:4]) < 2019:
        frac2019[idate] = fracyear-1.0
    else:
        frac2019[idate] = fracyear
    
    # ---- produce estimate of standard decaying-average bias correction

    bias = bias_decayavg[idate-1]
    print ('alpha, obs, analysis_box[idate], forecast_box[idate] = ', alpha, \
        analysis_box[idate], forecast_box[idate], bias)
    bias_decayavg[idate] = decayavg_bias(alpha, analysis_box[idate], \
        forecast_box[idate], bias)
    
    # ---- produce estimate of Kalman filter bias correction with seasonal variability.
    
    bias = bias_seasonalKF[idate-1]
    KF_betahat, bias_seasonalKF[idate] = seasonalKFbias(cosfac, sinfac, \
        cos2fac, sin2fac, Bbeta, B, R, KF_betahat, \
        analysis_box[idate], forecast_box[idate], bias)
    
# --- make plot    

print ('frac2019[0:] = ',frac2019[0:])
print ('frac2019[dstart:] = ',frac2019[dstart:])
    
fig = plt.figure(figsize=(9.,5.6))

axloc = [0.08,0.11,0.68,0.82]
a1 = fig.add_axes(axloc)
a1.set_title('Mean forecasts and observations, lead = '+clead+' h',fontsize=16)
a1.plot(frac2019[dstart:],forecast_box[dstart:]-analysis_box[dstart:],'.',color='Black',linewidth=0.3)
#a1.plot(np.arange(ndates)/365.,analysis_box,'.',color='Red',linewidth=0.3)
a1.plot(frac2019[dstart:],bias_decayavg[dstart:], 'k-',lw=2,label=r'Decaying average, $\alpha$='+calpha)
a1.plot(frac2019[dstart:],bias_seasonalKF[dstart:], 'r-',lw=2,\
    label='Kalman filter permitting seasonal\nand subseasonal bias dependence')

a1.plot([0,1],[0,0],'k-',lw=1)
a1.set_xlabel('Fraction of calendar year', fontsize=14)
a1.set_ylabel('Temperature bias (deg C)', fontsize=14)
a1.set_xlim(0,1)
a1.set_ylim(-3,3)
a1.grid (True,color='LightGray')
a1.legend(loc=0)
#crmse = '%.2f' %(rmse)
#crmse = 'RMSE = '+crmse
#a1.annotate(crmse,xy=(0.3,1.0))

axloc = [0.8,0.11,0.19,0.82]
a1 = fig.add_axes(axloc)
m = Basemap(llcrnrlon=rlonb,llcrnrlat=rlatb,urcrnrlon=rlone,\
    urcrnrlat=rlate,projection='mill',resolution='l')
m.drawcountries(linewidth=0.3,color='Gray')
m.drawstates(linewidth=0.2,color='Gray')
try:
    m.drawcoastlines(linewidth=0.3,color='Gray')
except:
    print ('not near coast')

plot_title = 'boxed_bias_lead='+clead+'.pdf'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title)
print ('Plot done')
    
    
    
    
    
            

