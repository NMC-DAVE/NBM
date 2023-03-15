"""
python plot_fitted_forecast_precip_cdfs_spline.py cmonth clead clat clon

where cmonth = 01 to 12
clead = 006 to 240 every 006
clat = latitude
clon = longitude, e.g., -100.0

this will plot both the forecast CDF and the cumulative hazard function.

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

def set_domain_boundaries(cdomain):
    
    """ used grib file of 2.5-km blend output grid to determine bounding 
        lat and lon, and from that, the domain bounding indices for the 
        0.25 GEFSv12 reforecast data that will encompass the domain.    
    """
    if cdomain == 'conus': 
        jmin = 93
        jmax = 246
        imin = 368
        imax = 686
    elif cdomain == 'pr':
        jmin = 243
        jmax = 256
        imin = 649
        imax = 667   
    elif cdomain == 'ak':
        jmin = 19
        jmax = 161
        imin = 201
        imax = 967
    else:
        print ('invalid domain.  Exiting.')     
        sys.exit()    
 
    return jmin, jmax, imin, imax

# =====================================================================

def find_nearest(vec, value):
    
    """ given a vector vec and a particular value, find the index in vec
    that is nearest to value"""
    
    idx = np.abs(vec-value).argmin()
    return idx

# =====================================================================

def get_surrounding_months(cmonth):

    if cmonth == 'Jan':
        cmonth_middle = 'Jan'
        cmonth_early = 'Dec'
        cmonth_late = 'Feb'
    elif cmonth == 'Feb':
        cmonth_middle = 'Feb'
        cmonth_early = 'Jan'
        cmonth_late = 'Mar'
    elif cmonth == 'Mar':
        cmonth_middle = 'Mar'
        cmonth_early = 'Feb'
        cmonth_late = 'Apr'
    elif cmonth == 'Apr':
        cmonth_middle = 'Apr'
        cmonth_early = 'Mar'
        cmonth_late = 'May'
    elif cmonth == 'May':
        cmonth_middle = 'May'
        cmonth_early = 'Apr'
        cmonth_late = 'Jun'
    elif cmonth == 'Jun':
        cmonth_middle = 'Jun'
        cmonth_early = 'May'
        cmonth_late = 'Jul'
    elif cmonth == 'Jul':
        cmonth_middle = 'Jul'
        cmonth_early = 'Jun'
        cmonth_late = 'Aug'
    elif cmonth == 'Aug':
        cmonth_middle = 'Aug'
        cmonth_early = 'Jul'
        cmonth_late = 'Sep'
    elif cmonth == 'Sep':
        cmonth_middle = 'Sep'
        cmonth_early = 'Aug'
        cmonth_late = 'Oct'
    elif cmonth == 'Oct':
        cmonth_middle = 'Oct'
        cmonth_early = 'Sep'
        cmonth_late = 'Nov'
    elif cmonth == 'Nov':
        cmonth_middle = 'Nov'
        cmonth_early = 'Oct'
        cmonth_late = 'Dec'
    elif cmonth == 'Dec':
        cmonth_middle = 'Dec'
        cmonth_early = 'Nov'
        cmonth_late = 'Jan'
    else:
        print ('invalid month')
        sys.exit()
    
    return cmonth_middle, cmonth_early, cmonth_late

# =====================================================================

def fraczero_possamps(nsamps, precip_ens):
    
    """
    
    from the vector input sample precip_ens, define the fraction of
    samples with zero precipitation.   For the positive samples, add
    a small random number to deal with the fact that the data was 
    discretized to ~0.1 mm, so that when later creating CDFs we don't 
    have empirical values with lots of tied amounts.  Also, sort the 
    nonzero amounts and return.
    
    """
    number_zeros = 0
    
    precip_ens_nonzero = np.delete(precip_ens, \
        np.where(precip_ens <= 0.0))  # censor at 0.0 mm
    precip_ens_nonzero = precip_ens_nonzero + \
        np.random.uniform(low=-0.1,high=0.1,size=len(precip_ens_nonzero))
    precip_ens_nonzero = np.delete(precip_ens_nonzero, \
        np.where(precip_ens_nonzero <= 0.0))  # censor at 0.0 mm
    nz = len(precip_ens_nonzero)    
        
    precip_ens_nonzero = np.sort(precip_ens_nonzero)  
    ntotal = len(precip_ens)
    nzero = ntotal - len(precip_ens_nonzero)
    fraction_zero = float(nzero) / float(ntotal)
    return fraction_zero, precip_ens_nonzero, nz

# =====================================================================

def read_forecast_spline_info(ccmonth, \
    master_directory_forecast_spline):

    """ load spline information for forecast netCDF file.
        For initial times other than 00 UTC, we need to find the
        appropriate spline file to read in based on 00 UTC
        forecast initialization times.   Assuming that
        the CDF characteristics are primarily a function of the
        diurnal cycle, make this adjustment. """

    ccycle = '00'
    icycle = int(ccycle)
    ilead = int(clead)
    ilead_use = icycle + ilead # diurnal cycle adjust
    if ilead_use > 240: ilead_use = ilead_use - 24
    if ilead_use < 10: # make 3-character variable of ilead_use
        clead_use = '00' + str(ilead_use)
    elif ilead_use >= 10 and ilead_use < 100:
        clead_use = '0' + str(ilead_use)
    else:
        clead_use = str(ilead_use)

    infile = master_directory_forecast_spline + \
        ccmonth + '_conus_GEFSv12_spline_info_h' + clead_use + '.nc'
    print ('reading forecast spline information from ', infile)
    nc = Dataset(infile)
    lons_spline_gefsv12_1d = nc.variables['lons'][:]
    lats_spline_gefsv12_1d = nc.variables['lats'][:]
    spline_info_gefsv12 = nc.variables['spline_info'][:,:,:,:]
    fraction_zero_gefsv12 = nc.variables['fzero'][:,:]
    usegamma_gefsv12 = nc.variables['usegamma'][:,:]
    quantile_99_gefsv12 = nc.variables['quantile_99'][:,:]
    number_knots_gefsv12 = nc.variables['number_knots'][:,:]
    ny_gefsv12, nx_gefsv12 = np.shape(usegamma_gefsv12)
    nc.close()
    lons_fcst_2d, lats_fcst_2d = \
        np.meshgrid(lons_spline_gefsv12_1d,lats_spline_gefsv12_1d)

    return lons_spline_gefsv12_1d, lats_spline_gefsv12_1d, \
        spline_info_gefsv12, fraction_zero_gefsv12, \
        usegamma_gefsv12, quantile_99_gefsv12, \
        number_knots_gefsv12, ny_gefsv12, nx_gefsv12, \
        lons_fcst_2d, lats_fcst_2d, clead_use


# =====================================================================
    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clead = sys.argv[2] # 006, 012, 018, 024, etc.
cleade = clead
clat = sys.argv[3] 
clon = sys.argv[4] 
ileadb = int(clead)-6
if ileadb < 10:
    cleadb = '00'+str(ileadb)
elif ileadb >= 10 and ileadb < 100:
    cleadb = '0'+str(ileadb)
else:
    cleadb = str(ileadb)

# ---- define constants.

rlon_input = float(clon)
rlat_input = float(clat)
cdomain = 'conus'
jmin, jmax, imin, imax = set_domain_boundaries(cdomain)
pflag = True
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
cmonths_early = ['12','01','02','03','04','05','06','07','08','09','10','11']
cmonths_late =  ['02','03','04','05','06','07','08','09','10','11','12','01']
cmonths_middle =  ['01','02','03','04','05','06','07','08','09','10','11','12']
cmonths = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
yearstart = 2002 # CCPA only available starting 2002
yearend = 2020 # companion reforecasts end at end of 2019
cmonth_middle, cmonth_early, cmonth_late = get_surrounding_months(cmonth)
imonth = cmonths.index(cmonth)

# --- load forecast spline info.

master_directory_forecast_spline = \
    '/Volumes/NBM/conus_gefsv12/CDF_spline/'
    # where forecast spline info is stored
master_directory = '/Volumes/NBM/'+\
    cdomain+'_gefsv12/precip/netcdf/'
    
# ---- read forecast spline parameters from netCDF file

print ('reading forecast spline parameters')
lons_1d_realtime, lats_1d_realtime, \
    spline_info_gefsv12, fraction_zero_gefsv12, \
    usegamma_gefsv12, quantile_99_gefsv12, \
    number_knots_gefsv12, ny_gefsv12, nx_gefsv12, \
    lons_fcst_2d, lats_fcst_2d, clead_use = \
    read_forecast_spline_info( \
    cmonth, master_directory_forecast_spline)

# --- make sure longitudes are negative for west, and
#     create 2D arrays of lat/lon

if np.max(lons_1d_realtime) > 180.: \
    lons_1d_realtime = lons_1d_realtime - 180.
lons_2d, lats_2d = \
    np.meshgrid(lons_1d_realtime, lats_1d_realtime)
    
# ---- find index of nearest grid point for spline data

print ('rlon_input, rlat_input = ', rlon_input, rlat_input)
distx = np.abs(lons_2d - rlon_input)
disty = np.abs(lats_2d - rlat_input)
dist = distx**2 + disty**2
i = np.argmin(dist)
print (i)
jspline, ispline = np.unravel_index\
    (i, shape=(ny_gefsv12, nx_gefsv12), order='C')    

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time + surrounding months.
#      All members, dates for this month have been
#      smushed into one leading index, dimension
#      nsamps, since the date of the forecast within the month and
#      the member number is irrelevant for the distribution fitting.

ncfile = master_directory + cmonth_middle + \
    '_apcp_sfc_h' + clead + '.nc'
print (ncfile)
nc = Dataset(ncfile)
lons_fcst = nc.variables['lons_fcst'][:]
lats_fcst = nc.variables['lats_fcst'][:]
lons_2d_fcst, lats_2d_fcst = np.meshgrid(lons_fcst, lats_fcst)
nxgefs = len(lons_fcst)
nygefs = len(lats_fcst)

# ---- find index of nearest grid point for spline data

distx = np.abs(lons_2d_fcst - rlon_input)
disty = np.abs(lats_2d_fcst - rlat_input)
dist = distx**2 + disty**2
i = np.argmin(dist)
jgefs, igefs = np.unravel_index(i, shape=(nygefs, nxgefs), order='C')    
precip_middle = nc.variables['apcp_fcst'][:,jgefs, igefs]
nsamps_middle = len(precip_middle)
nc.close()

ncfile = master_directory + cmonth_early + '_apcp_sfc_h' + clead + '.nc'
print (ncfile)
nc = Dataset(ncfile)
precip_early = nc.variables['apcp_fcst'][:,jgefs, igefs]
nsamps_early = len(precip_early)
nc.close()

ncfile = master_directory + cmonth_late + '_apcp_sfc_h' + clead + '.nc'
print (ncfile)
nc = Dataset(ncfile)
precip_late = nc.variables['apcp_fcst'][:,jgefs, igefs]
nsamps_late = len(precip_late)
nc.close()

# --- form a vector for this grid point of all 3 months of data

nsamps = nsamps_middle + nsamps_early + nsamps_late
precip_tseries = np.concatenate((precip_middle, \
    precip_early, precip_late))    
    
tp = np.min( [np.abs(np.min(precip_tseries)), 0.0] )
teeny_precip = tp*np.ones(nsamps)
precip_tseries = precip_tseries - teeny_precip[:]

# ---- determine the fraction zero, the number of nonzeros (nz)
#      and the sorted sample of nonzero values (precip_ens_nonzero)

fraction_zero, precip_ens_nonzero, nz = \
    fraczero_possamps(nsamps, precip_tseries) # return sorted

# ---- build empirical CDF for this grid point's nonzero data

x = np.arange(0.0,50.01,0.1)
nx = len(x)
len_nonzero = len(precip_ens_nonzero)
fnzero = float(len_nonzero)
nz = int(fnzero)
print ('nz = ', nz)
empirical_cdf = 1.0 / (2.0*nz) + \
    np.arange(nz)/nz  # "Hazen" plotting function
empirical_cdf_with_zeros = \
    fraction_zero + (1.-fraction_zero)*empirical_cdf
empirical_hazard_fn = -np.log(1.0 - empirical_cdf)

# ---- get the spline fit.

numknots = number_knots_gefsv12[jspline, ispline]  # now this number includes boundary
knots_in = spline_info_gefsv12[jspline, ispline,0,0:numknots]
bspline_coef_in = spline_info_gefsv12[jspline, ispline,1,0:numknots]
precip_regular = 0.0001+np.arange(0.,401.,1.)/10.

npr = len(precip_regular)
nthree = 3
ier = 0
cdf_spline = np.zeros((nz), dtype=np.float64)

splines_tuple = (spline_info_gefsv12[jspline, ispline,0,0:numknots], \
    spline_info_gefsv12[jspline, ispline,1,0:numknots], 3)
spline_hazard = splev(precip_regular, splines_tuple)
qpositive = 1.0 - np.exp(-spline_hazard) 
cdf_spline = qpositive   
hazard_fn = -np.log(1.0 - cdf_spline)

interior_knots_precip = knots_in[4:13]
cdf_at_knots = np.zeros((9), dtype=np.float32)
haz_at_knots = np.zeros((9), dtype=np.float32)
precip_at_knots = interior_knots_precip 
for i in range(9):
    idx = find_nearest(precip_regular,interior_knots_precip[i])
    haz_at_knots[i] = hazard_fn[idx]
    cdf_at_knots[i] = cdf_spline[idx]

cdf_spline_with_zeros = fraction_zero + (1.-fraction_zero)*qpositive


# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('Spline fit to hazard function,\n'+\
    cmonth+', '+clead+' h, 00 UTC IC, '+\
    clon+' W '+clat+' N',fontsize=13)
ax.plot(precip_ens_nonzero, empirical_hazard_fn,\
    color='Red',lw=2,label='Empirical')
ax.plot(precip_regular, spline_hazard,\
    color='Blue',lw=2,label='Cubic spline')

plt.ylabel('Cumulative hazard function',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.0,10.0)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,25)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'hazard_precip_forecast_example_spline.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 

# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('Spline fit to positive forecast precipitation values,\n'+\
    cmonth+', '+clead+' h, 00 UTC IC, '+\
    clon+' W '+clat+' N',fontsize=13)

ax.plot(precip_ens_nonzero, empirical_cdf_with_zeros,\
    color='Red',lw=2,label='Empirical')
ax.plot(precip_regular, cdf_spline_with_zeros,\
    color='Blue',lw=2,label='Cubic spline')

plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.0,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,25)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_forecast_example_spline.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 






# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(8,4.3))#
ax = f.add_axes([.07,.14,.41,.75])
ax.set_title('(a) Forecast cumulative hazard function',fontsize=12.5)
ax.plot(precip_ens_nonzero, empirical_hazard_fn,'--',\
    color='Red',lw=2,label='Empirical')
ax.plot(precip_regular, spline_hazard,\
    color='RoyalBlue',lw=2,label='Cubic spline')
ax.plot(precip_at_knots, haz_at_knots, 'k.')     

plt.ylabel('Cumulative hazard function',fontsize=9)
ax.legend(loc=0)
ax.set_ylim(0.0,10.0)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=9)

ax = f.add_axes([.57,.14,.41,.75])
ax.set_title('(b) Forecast CDF',fontsize=12.5)
ax.plot(precip_ens_nonzero, empirical_cdf,\
    color='Red',lw=1.6,label='Empirical')
ax.plot(precip_regular, cdf_spline,\
    color='Blue',lw=1.6,label='Cubic spline')
ax.plot(precip_at_knots, cdf_at_knots, 'k.')    

plt.ylabel('Non-exceedance probability',fontsize=9)
ax.legend(loc=0)
ax.set_ylim(0.0,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=9)
figname = 'CDF_precip_forecast_example_spline_ITH.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 





                
outfile = 'forecast_data_cdf.cPick'
ouf = open(outfile, 'wb')
cPickle.dump(precip_ens_nonzero, ouf)
cPickle.dump(empirical_hazard_fn, ouf)
cPickle.dump(spline_hazard, ouf)
cPickle.dump(empirical_cdf, ouf)
cPickle.dump(precip_regular, ouf)
cPickle.dump(cdf_spline, ouf)
cPickle.dump(precip_at_knots, ouf)
cPickle.dump(haz_at_knots, ouf)
cPickle.dump(cdf_at_knots, ouf)
ouf.close()



