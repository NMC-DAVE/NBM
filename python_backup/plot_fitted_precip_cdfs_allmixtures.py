"""
plot_fitted_precip_cdfs_allmixtures.py cmonth clead jy ix

"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import LSQUnivariateSpline, splrep, splev

import _pickle as cPickle
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

# =====================================================================
    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clead = sys.argv[2] # 03, 06, 12, etc.
cjy = sys.argv[3] # units of 0, 5, etc
cix = sys.argv[4] # units of 0, 5, etc
jy = int(cjy)
ix = int(cix)
cdomain = 'conus'

# --- load from cPickle file

infile = cmonth+'_'+cdomain+'_jy'+cjy+'_ix'+cix+\
    '_singlepoint_apcp_gamma_parameters_h' + clead + '.cPick'

print (jy,ix)
print ('reading from ', infile)
inf = open(infile, 'rb')
fraction_zero = cPickle.load(inf)
weight_save_1m = cPickle.load(inf)
alpha_save_1m = cPickle.load(inf)
beta_save_1m = cPickle.load(inf)
weight_save_2m = cPickle.load(inf)
alpha_save_2m = cPickle.load(inf)
beta_save_2m = cPickle.load(inf)
weight_save_3m = cPickle.load(inf)
alpha_save_3m = cPickle.load(inf)
beta_save_3m = cPickle.load(inf)
Dnstat1 = cPickle.load(inf)
Dnstat2 = cPickle.load(inf)
Dnstat3 = cPickle.load(inf)
precip_nonzero = cPickle.load(inf)
lon = cPickle.load(inf)
lat = cPickle.load(inf)
inf.close()

clat = str(lat)
clon = str(lon)

# ---- build empirical CDF

x = np.arange(0.0,50.01,0.1)
nx = len(x)
len_nonzero = len(precip_nonzero)
fnzero = float(len_nonzero)
nz = int(fnzero)
pnzsort = np.sort(precip_nonzero)
cdf_empirical = np.zeros((nx),dtype=np.float64)
for i, xi in enumerate(x):
    nbelow = (pnzsort < xi).sum() 
    cdf_empirical[i] = fraction_zero + (1.-fraction_zero)*float(nbelow)/fnzero

# --- build CDFs for 1-parameter

y0 = x / beta_save_1m
cdf0 = stats.gamma.cdf(y0, alpha_save_1m)
cdf_fitted1 = fraction_zero + (1.-fraction_zero)*cdf0 

# --- build CDFs for 2-parameter

y0 = x / beta_save_2m[0]
y1 = x / beta_save_2m[1]

cdf0 = stats.gamma.cdf(y0, alpha_save_2m[0])
cdf1 = stats.gamma.cdf(y1, alpha_save_2m[1])
cdf_fitted2 = fraction_zero + (1.-fraction_zero)* \
    (weight_save_2m[0]*cdf0 + weight_save_2m[1]*cdf1) 

# --- build CDFs for 3-parameter

y0 = x / beta_save_3m[0]
y1 = x / beta_save_3m[1]
y2 = x / beta_save_3m[2]

cdf0 = stats.gamma.cdf(y0, alpha_save_3m[0])
cdf1 = stats.gamma.cdf(y1, alpha_save_3m[1])
cdf2 = stats.gamma.cdf(y2, alpha_save_3m[2])
cdf_fitted3 = fraction_zero + (1.-fraction_zero)* \
    (weight_save_3m[0]*cdf0 + weight_save_3m[1]*cdf1 + \
    weight_save_3m[2]*cdf2)   

# ---- spline fit with the focus on knots at higher quantiles

query_these_indices = [ nz//10, nz//4, nz//2, (3*nz)//5, \
    (4*nz)//5, (17*nz)//20, (9*nz)//10, (19*nz)//20, (39*nz)//40]
empirical_precipvals = pnzsort[query_these_indices]
empirical_cdf = 1.0/(2.0*nz) + np.arange(nz)/nz    

spltemp = splrep(pnzsort, empirical_cdf, xb=0., task=-1, \
    t = empirical_precipvals)
print (len(spltemp[0]))
print (len(spltemp[1]))
print (spltemp[2])
joe = np.zeros((1,1,2,17), dtype=np.float64)
print (np.shape(joe))
joe[0,0,0,:] = spltemp[0]
joe[0,0,1,:] = spltemp[1]
print (joe)
print ('np.shape(joe) = ', np.shape(joe))

spline_fit = fraction_zero + (1.-fraction_zero)*splev(pnzsort, spltemp) 



# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('(a) Spline fit, '+cmonth+', '+clead+'-h forecast, '+\
    clon+' W '+clat+' N',fontsize=13)
ax.plot(pnzsort, spline_fit,color='Blue',lw=2,label='Cubic spline')
ax.plot(x,cdf_empirical,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.6,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_spline.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 
              
                
              
        
# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('(a) 1-component mixture, '+cmonth+', '+clead+'-h forecast, '+\
    clon+' W '+clat+' N',fontsize=13)
ax.plot(x,cdf_fitted1,color='Blue',lw=2,label='Fitted 1-component Gamma mixture')
ax.plot(x,cdf_empirical,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.6,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_1component.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()


# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('(b) 2-component mixture, '+cmonth+', '+clead+'-h forecast, '+\
    clon+' W '+clat+' N',fontsize=13)
ax.plot(x,cdf_fitted2,color='Blue',lw=2,label='Fitted 2-component Gamma mixture')
ax.plot(x,cdf_empirical,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.6,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_2component.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()


# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('(c) 3-component mixture, '+cmonth+', '+clead+'-h forecast, '+\
    clon+' W '+clat+' N',fontsize=13)
ax.plot(x,cdf_fitted3,color='Blue',lw=2,label='Fitted 3-component Gamma mixture')
ax.plot(x,cdf_empirical,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.6,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_3component.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()














