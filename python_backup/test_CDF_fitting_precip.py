"""
test_CDF_fitting_precip.py cmonth clead clon clat

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import _pickle as cPickle
import scipy.stats as stats
base = importr('base')
mixtools = importr('mixtools')


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

# ---- inputs from command line

cmonth = sys.argv[1] # e.g, Jan, Feb
clead = sys.argv[2] # 06, etc.
clon = sys.argv[3] # e.g., -185.0
clat = sys.argv[4] # e.g, -10.0 for 10 S.

already = False
rlon = float(clon)
rlat = float(clat)
master_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'
 

if already == False:
    
    # ---- load this grid point's precipitation totals from netCDF file.

    years = range(2000,2020)
    nyears = len(years)
    precip_nonzero = []
    number_zeros = 0
    nmembers = 5
    ndays = 31
    for iyear, year in enumerate(years):
        print (iyear, year)
        cyear = str(year)
        ncfile = master_directory + cmonth + cyear + '_h' + clead + '.nc'
        print (ncfile)
        nc = Dataset(ncfile)
        if iyear == 0:
            lons_1d = nc.variables['lons_fcst'][:]
            lats_1d = nc.variables['lats_fcst'][:]
            #lats_1d = np.flipud(lats_1d)
            ilon = find_nearest(lons_1d, rlon)
            ilat = find_nearest(lats_1d, rlat)
            print ('ilat,ilon = ', ilat,ilon)
        #precip_ens = np.squeeze(nc.variables['apcp_fcst'][:,:,ilat,ilon])
        precip_ens = nc.variables['apcp_fcst'][:,:,ilat,ilon]
        #if year == 2000:
        #    precip_ens_4D = nc.variables['apcp_fcst'][:,:,:,:]
        #    print ('precip_ens[29,2,ilat,:] = ', precip_ens_4D[29,2,ilat,:])
        #    sys.exit()
        nc.close()
        for imem in range(nmembers):
            for iday in range(ndays):
                print ('year, imem, iday, precip = ',year, imem, iday, precip_ens[iday,imem])
                if precip_ens[iday,imem] <= 0.0:
                    number_zeros = number_zeros+1
                else:
                    precip_nonzero.append(precip_ens[iday,imem])

    fraction_zero = float(number_zeros) / float(ndays*nmembers*nyears)
    precip_nonzero = np.array(precip_nonzero)

    outfile = 'precipitation_data.cPick'
    print ('writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(fraction_zero, ouf)
    cPickle.dump(precip_nonzero, ouf)
    ouf.close()
    
else:

    infile = 'precipitation_data.cPick'
    print ('reading from ', infile)
    inf = open(infile,'rb')
    fraction_zero = cPickle.load(inf)
    precip_nonzero = cPickle.load(inf)
    inf.close()


# --- convert to R vector, call R Gamma mixture function

precip_nonzero_R = robjects.FloatVector(precip_nonzero)
#rmixture = robjects.r['mixtools.gammamixEM']
result_R = mixtools.gammamixEM(precip_nonzero_R, k=3)
result_np = np.asarray(result_R, dtype=object)
#print (result_np)
#print (np.shape(result_np))
#print ('result_R = ', result_R)
result_weights_np = np.asarray(result_R[1])
#print ('result_weights_np = ', result_weights_np)
result_alpha_beta_np = np.asarray(result_R[2])
#print ('result_alpha_beta_np = ', result_alpha_beta_np)
# rows give alpha, beta

print ('alpha 0, 1 = ', result_alpha_beta_np[0,0], result_alpha_beta_np[0,1])
print ('beta 0, 1 = ', result_alpha_beta_np[1,0], result_alpha_beta_np[1,1])


x = np.arange(0.0,50.01,0.1)
nx = len(x)
y0 = x / result_alpha_beta_np[1,0]
y1 = x / result_alpha_beta_np[1,1]
y2 = x / result_alpha_beta_np[1,2]
#y3 = x / result_alpha_beta_np[1,3]

cdf0 = stats.gamma.cdf(y0, result_alpha_beta_np[0,0])
cdf1 = stats.gamma.cdf(y1, result_alpha_beta_np[0,1])
cdf2 = stats.gamma.cdf(y2, result_alpha_beta_np[0,2])
#cdf3 = stats.gamma.cdf(y3, result_alpha_beta_np[0,3])
cdf_fitted = fraction_zero + (1.-fraction_zero)* \
    (result_weights_np[0]*cdf0 + result_weights_np[1]*cdf1 + \
    result_weights_np[2]*cdf2) # + result_weights_np[3]*cdf3)

cdf_empirical = np.zeros((nx),dtype=np.float64)
len_nonzero = len(precip_nonzero)
fnzero = float(len_nonzero)
pnzsort = np.sort(precip_nonzero)
#print (np.shape(pnzsort))
#print (pnzsort[0:-1:10])
for i, xi in enumerate(x):
    #nbelow = len(np.where(pnzsort < x[i])) 
    #nbelow = np.sum(np.where(pnzsort < xi)) 
    nbelow = (pnzsort < xi).sum() 
    cdf_empirical[i] = fraction_zero + (1.-fraction_zero)*float(nbelow)/fnzero
    #print ('i,xi, nbelow,fnzero, cdf_em = ',i,xi, nbelow, fnzero, cdf_empirical[i])

# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,4.5))#
ax = f.add_axes([.15,.14,.8,.77])
ax.set_title('Example: 3 to 6-h forecast precipitation, fitted and empirical CDFs',fontsize=11)
ax.plot(x,cdf_fitted,color='Blue',lw=2,label='Fitted 3-component Gamma mixture')
ax.plot(x,cdf_empirical,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.7,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,20)
ax.set_xlabel('3-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()

# ---- Q-Q plot of forecast data

f = plt.figure(figsize=(5.5,5.5))

ax = f.add_axes([.15,.14,.8,.77])
ax.set_title('Q-Q plot',fontsize=13)
ax.plot([0,1],[0,1],color='Gray',lw=0.5)
ax.plot(cdf_empirical,cdf_fitted,color='Red',lw=2)
ax.set_ylim(0.7,1)
ax.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0.7,1)
ax.set_xlabel('Empirical forecast quantile',fontsize=11)
ax.set_ylabel('Fitted forecast quantile',fontsize=11)

figname = 'QQplot_example.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()



