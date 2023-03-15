"""
plot_fitted_precip_cdfs_allmixtures.py cmonth clead jy ix

"""

import os, sys
from datetime import datetime
from dateutils import daterange, dateshift
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import LSQUnivariateSpline, splrep, splev
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clon = sys.argv[3] # units of 0, 5, etc
clat = sys.argv[4] # units of 0, 5, etc

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cmonths_before = ['Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov']
cmonths_after = ['Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan']

# ---- read in lat/lon for the older GEFSv10 reforecasts.

gefsv10_directory = '/Volumes/NBM/refcstv2/'
infile = gefsv10_directory + 'refcstv2_precip_ccpav3_000_to_006.nc'
nc = Dataset(infile)
lons_gefsv10 = nc.variables['lons_fcst'][:,:]
lats_gefsv10 = nc.variables['lats_fcst'][:,:]
lons_1d = lons_gefsv10[0,:]
lats_1d = lats_gefsv10[:,0]
print ('lons_1d = ',lons_1d)
print ('lats_1d = ',lats_1d)
ilon_gefsv10 = find_nearest(lons_1d,float(clon))
jlat_gefsv10 = find_nearest(lats_1d,float(clat))
print ('nearest lat/lon = ',jlat_gefsv10, ilon_gefsv10)
print ('lat, lon of nearest = ', lats_gefsv10[jlat_gefsv10, ilon_gefsv10],\
    lons_gefsv10[jlat_gefsv10, ilon_gefsv10])

# ---- read in time series of GEFSv10 precipitation at this lat, lon

for ilead, clead in enumerate(['000_to_006', \
    '006_to_012', '012_to_018', '018_to_024']):
    
    # ---- read in data for the older GEFSv10 reforecasts.

    infile = gefsv10_directory + 'refcstv2_precip_ccpav3_000_to_006.nc'
    print (infile)
    nc = Dataset(infile)
    if ilead == 0:
        yyyymmddhh_begin_gefsv10 = nc.variables['yyyymmddhh_begin'][:]
        apcp_fcst_ens_gefsv10 = nc.variables['apcp_fcst_ens'] \
            [:,:,jlat_gefsv10, ilon_gefsv10]
    else:
        apcp_fcst_ens_gefsv10 = apcp_fcst_ens_gefsv10 + \
            nc.variables['apcp_fcst_ens'][:,:,jlat_gefsv10, ilon_gefsv10]
    if ilead == 3:
        yyyymmddhh_end_gefsv10 = nc.variables['yyyymmddhh_end'][:]
    nc.close()
    
apcp_fcst_ens_gefsv10_mean = np.mean(apcp_fcst_ens_gefsv10,axis=1)
print ('sample apcp_fcst_ens_gefsv10_mean = ', apcp_fcst_ens_gefsv10_mean[0:20])
print ('sample yyyymmddhh_end_gefsv10 = ', yyyymmddhh_end_gefsv10[0:20])
    
    
# ---- read in lat/lon for the newer GEFSv12 reforecasts.

gefsv12_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
infile = gefsv12_directory + 'Oct_conus_reforecast_precip_h006.nc'
nc = Dataset(infile)
lons_gefsv12 = nc.variables['lons_fcst'][:]
lats_gefsv12 = nc.variables['lats_fcst'][:]
print ('lons_gefsv12 = ',lons_gefsv12)
print ('lats_gefsv12 = ',lats_gefsv12)
ilon_gefsv12 = find_nearest(lons_gefsv12,float(clon))
jlat_gefsv12 = find_nearest(lats_gefsv12,float(clat))
print ('nearest lat/lon = ',jlat_gefsv10, ilon_gefsv10)
print ('lat, lon of nearest = ', lats_gefsv12[jlat_gefsv12, ilon_gefsv12],\
    lons_gefsv12[jlat_gefsv12, ilon_gefsv12])

# ---- read in time series of GEFSv12 precipitation at this lat, lon

imonth = cmonths.index(cmonth)
for imonth, cmonth in enumerate(cmonths_before[imonth], \
    cmonths[imonth], cmonths_after[imonth]):
    for ilead, clead in enumerate(['006','012', '018', '024']):
    
        # ---- read in lat/lon for the newer GEFSv12 reforecasts.

        infile = gefsv12_directory + cmonths + \
            '_conus_reforecast_precip_h' + clead + '.nc'
        print ('reading ',infile)
        nc = Dataset(infile)
        if ilead == 0:
            apcp_fcst_ens_gefsv12 = nc.variables['apcp_fcst']\
                [:,jlat_gefsv12, ilon_gefsv12]
        else:
            apcp_fcst_ens_gefsv12 = apcp_fcst_ens_gefsv12 + \
                nc.variables['apcp_fcst'][:,jlat_gefsv12, ilon_gefsv12]
        if ilead == 3:
            yyyymmddhh_end_gefsv12 = nc.variables['yyyymmddhh_fcst'][:]
        nc.close()
        
        # ---- append to end of record so far
        
        if imonth == 0:
            psamples_gefsv12 = apcp_fcst_ens_gefsv12
            yyyymmddhh_end_gefsv12_3month = yyyymmddhh_end_gefsv12
        else:
            psamples_gefsv12 = np.append(psamples_gefsv12, apcp_fcst_ens_gefsv12)
            yyyymmddhh_end_gefsv12_3month = \
                np.append(yyyymmddhh_end_gefsv12_3month, yyyymmddhh_end_gefsv12)
                

# ---- 5 samples every time, get mean of these.  Also subset dates.

nsamps = len(psamples_gefsv12)
yyyyymmddhh_end_gefsv12 = []
psamples_gefsv12_mean = []
for isamp in range(0,nsamps,5):
    yyyyymmddhh_end_gefsv12 = \
        yyyyymmddhh_end_gefsv12.append(yyyymmddhh_end_gefsv12_3month[isamp])
    psamples_gefsv12_mean.append \
        (np.mean(yyyymmddhh_end_gefsv12_3month[isamp:isamp+5])

# ---- strip GEFSv12 data to only those with years 2002 - 2015 for overlap with
#      GEFSv10

print ('yyyyymmddhh_end_gefsv12[0:10] = ', yyyyymmddhh_end_gefsv12[0:10])
yyyy_end_gefsv12 = yyyymmddhh_end_gefsv12 // 1000000
print ('yyyy_end_gefsv12[0:10] = ',yyyy_end_gefsv12[0:10] )
a = np.where(np.logical_and(\
    yyyy_end_gefsv12 >= 2002,  yyyy_end_gefsv12_3month <= 2015))
    
ktr = 0
psamples_thinned_gefsv12 = []
psamples_thinned_gefsv10 = []
for idx in a:
    yyyymmddhh = yyyymmddhh_end_gefsv12[idx]
    psamples_thinned_gefsv12 = psamples_thinned_gefsv12.append(psamples_gefsv12[idx])
    index = yyyymmddhh_end_gefsv10.index(yyyymmddhh)
    psamples_thinned_gefsv10 = psamples_thinned_gefsv10.append(psamples_gefsv12[index])
    
    




# ---- build empirical CDFs of both

a = where(psamples_thinned_gefsv10 > 0.001)
precip_nonzero_gefsv10 = psamples_thinned_gefsv10[a]
a = where(psamples_thinned_gefsv12 > 0.001)
precip_nonzero_gefsv12 = psamples_thinned_gefsv12[a]
fraction_zero_gefsv10 = float(len(precip_nonzero_gefsv10)) / float(len(psamples_thinned_gefsv10))
fraction_zero_gefsv12 = float(len(precip_nonzero_gefsv12)) / float(len(psamples_thinned_gefsv12))
print ('fraction_zero_gefsv10 = ',fraction_zero_gefsv10 )
print ('fraction_zero_gefsv12 = ',fraction_zero_gefsv12 )



len_nonzero_v10 = len(precip_nonzero_gefsv10)
len_nonzero_v12 = len(precip_nonzero_gefsv12)
fnzero_v10 = float(len_nonzero_v10)
fnzero_v12 = float(len_nonzero_v12)
nz = int(fnzero)
pnzsort_v10 = np.sort(precip_nonzero_gefsv10)
pnzsort_v12 = np.sort(precip_nonzero_gefsv12)
cdf_empirical_v10 = np.zeros((nx),dtype=np.float64)
cdf_empirical_v12 = np.zeros((nx),dtype=np.float64)

x = np.arange(0.0,50.01,0.1)
nx = len(x)
for i, xi in enumerate(x):
    nbelow_v10 = (pnzsort_v10 < xi).sum()
    cdf_empirical_v10[i] = fraction_zero_gefsv10 + \
        (1.-fraction_zero_gefsv10)*float(nbelow_v10)/fnzero_v10
    nbelow_v12 = (pnzsort_v12 < xi).sum()
    cdf_empirical_v12[i] = fraction_zero_gefsv12 + \
        (1.-fraction_zero_gefsv12)*float(nbelow_v12)/fnzero_v12


# ---- scatterplot v10 vs. v12, and plot the CDFs

fig1 = plt.figure(figsize=(9.,4.9))
a1 = fig1.add_axes([0.07,0.14,0.41,0.73])
a1.set_title('(a) GEFSv10 vs. GEFSv12 Day +1\nprecipitation at Ithaca Game Farm',\
    fontsize=16)
a1.plot([0,50],[0,50], linewidth=1, color='Gray')
a1.plot(psamples_thinned_gefsv10[:], psamples_thinned_gefsv12[:],\
    's',color='Blue',markersize=1.5,  rasterized=True)
r = stats.pearsonr(psamples_thinned_gefsv10[:], psamples_thinned_gefsv12[:]))
plot_units = ' (mm)'
a1.set_xlabel(r'GEFSv12 mean forecast'+plot_units, fontsize=13)
a1.set_ylabel(r'GEFSv10 mean forecast'+plot_units, fontsize=13)
a1.set_xlim(-0.01,50.0)
a1.set_ylim(-0.01,50.0)
a1.grid(color='Gray', linestyle='-', linewidth=0.5)
print (r[0])
cr = 'r ='+'{0:7.2f}'.format(r[0])
a1.text(2,45,cr,fontsize=17,color='Blue')

a2 = f.add_axes([.57,.14,.41,.75])
a2.set_title('(b) Forecast CDFs',fontsize=13)
a2.plot(x,cdf_empirical_v10,color='Red',lw=2,label='GEFSv10')
a2.plot(x,cdf_empirical_v12,color='Blue',lw=2,label='GEFSv12')
a2.set_ylabel('Non-exceedance probability',fontsize=11)
a2.legend(loc=0)
a2.set_ylim(0.3,1)
a2.grid(True,lw=0.25,color='LightGray')
a2.set_xlim(0,50)
ax.set_xlabel('Daily total precipitation (mm)',fontsize=11)

figname = 'Day_precip_GEFSv10_GEFSv12.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()



