"""
python plot_fitted_precip_cdfs_spline_v3.py cmonth clead clat clon

where cmonth = 01 to 12
clead = 006 to 240 by 006
clat = latitude
clon = longitude (west is negative)

plot analyzed CDF with spline fit, and cumulative hazard function

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
from control_splev import control_splev

import _pickle as cPickle
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'


# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx
# =====================================================================

def get_surrounding_months(cmonth):

    if cmonth == 'Jan':
        cmonth_middle = '01'
        cmonth_early = '12'
        cmonth_late = '02'
    elif cmonth == 'Feb':
        cmonth_middle = '02'
        cmonth_early = '01'
        cmonth_late = '03'
    elif cmonth == 'Mar':
        cmonth_middle = '03'
        cmonth_early = '02'
        cmonth_late = '04'
    elif cmonth == 'Apr':
        cmonth_middle = '04'
        cmonth_early = '03'
        cmonth_late = '05'
    elif cmonth == 'May':
        cmonth_middle = '05'
        cmonth_early = '04'
        cmonth_late = '06'
    elif cmonth == 'Jun':
        cmonth_middle = '06'
        cmonth_early = '05'
        cmonth_late = '07'
    elif cmonth == 'Jul':
        cmonth_middle = '07'
        cmonth_early = '06'
        cmonth_late = '08'
    elif cmonth == 'Aug':
        cmonth_middle = '08'
        cmonth_early = '07'
        cmonth_late = '09'
    elif cmonth == 'Sep':
        cmonth_middle = '09'
        cmonth_early = '08'
        cmonth_late = '10'
    elif cmonth == 'Oct':
        cmonth_middle = '10'
        cmonth_early = '09'
        cmonth_late = '11'
    elif cmonth == 'Nov':
        cmonth_middle = '11'
        cmonth_early = '10'
        cmonth_late = '12'
    elif cmonth == 'Dec':
        cmonth_middle = '12'
        cmonth_early = '11'
        cmonth_late = '01'
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

def read_analysis_spline_info (clead_use, master_directory_panal_spline,\
    cmonth):

    """ read netCDF file here for inverse spline parameters for
    the combined CCPA/MSWEP precip analysis CDFs on the 
    Lambert-conformal NDFD 2.5 km grid surrounding the CONUS.  
  
    **** Note that if we are applying to cycles other than 
    forecasts with 00UTC init, we'll need adjust the dates 
    of the spline files that we use to sync up. 
    """ 
   
    ndays = int(clead_use) // 24
    ilead = int(clead_use)-ndays*24
    if ilead == 0:
        cleada = '00'
    elif ilead == 6:
        cleada = '06'
    elif ilead == 12:
        cleada = '12'
    elif ilead == 18:
        cleada = '18'
        
    infile = master_directory_panal_spline + cmonth + \
        '_conus_CCPA_spline_info_h' + cleada + 'UTC.nc'
    print ('reading from ', infile)
    nc = Dataset(infile)
    spline_info_inv = nc.variables['spline_info_inv'][:,:,:,:]
    fraction_zero_ndfd = nc.variables['fzero'][:,:]
    number_knots = nc.variables['number_knots'][:,:]
    usegamma_ndfd = nc.variables['usegamma'][:,:]
    lons_ndfd = nc.variables['lons'][:,:]
    lons_ndfd = lons_ndfd
    lats_ndfd = nc.variables['lats'][:,:]
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
    nc.close()     
    
    return spline_info_inv, fraction_zero_ndfd, usegamma_ndfd, \
        number_knots, lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd


# =====================================================================
    
# ---- inputs from command line

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
cend_hour = sys.argv[2] # 06, 12, 18, 00
clat = sys.argv[3] 
clon = sys.argv[4] 

# ---- define constants.

rlon_input = float(clon)
rlat_input = float(clat)
cdomain = 'conus'
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
imonth = cmonths.index(cmonth)

# --- load analysis spline info.

master_directory_panal_spline = \
    '/Volumes/NBM/conus_panal/CDF_spline/'
spline_info_inv, fraction_zero_ndfd, \
    usegamma_ndfd, number_knots, \
    lons_ndfd, lats_ndfd, ny_ndfd, nx_ndfd = \
    read_analysis_spline_info (cend_hour, \
    master_directory_panal_spline, cmonth)
if np.max(lons_ndfd) > 180.: lons_ndfd = lons_ndfd - 180.
    
# ---- find index of nearest ndfd point

distx = np.abs(lons_ndfd - rlon_input)
disty = np.abs(lats_ndfd - rlat_input)
dist = distx**2 + disty**2
i = np.argmin(dist)
jndfd, indfd = np.unravel_index(i, shape=(ny_ndfd, nx_ndfd), order='C')   

# ---- determine the overall number of daily precipitation 
#      samples across all years for this month and the surrounding
#      two months

iearly = int(cmonths_early[imonth])-1
ilate = int(cmonths_late[imonth])-1

if imonth != 1:  # not Feb
    nsamps_mid = ndaysomo[imonth]*18
else:
    nsamps_mid = 4*ndaysomo_leap[imonth] + 14*ndaysomo[imonth]

if iearly != 1:  # not Feb    
    nsamps_early = ndaysomo[iearly]*20
else:
    nsamps_early = 4*ndaysomo_leap[iearly] + 14*ndaysomo[iearly]
if ilate != 1:  # not Feb    
    nsamps_late = ndaysomo[ilate]*20
else:
    nsamps_late = 4*ndaysomo_leap[ilate] + 14*ndaysomo[ilate]
nsamps = nsamps_mid + nsamps_early + nsamps_late

# ---- read in the previously generated analysis netCDF file with 
#      precipitation for this month and lead time as well as the 
#      surrounding two months.  All dates for this month have
#      been smushed into one leading index, dimension nsamps,
#      since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution 
#      fitting.

ktr = 0
master_directory = '/Volumes/NBM/conus_panal/'
for iyear in range(yearstart, yearend):

    # --- loop over the month in question and the surrounding 2 months,
    #     and read in precipitation analysis for this month.

    for cmo in [cmonths_middle[imonth], \
    cmonths_early[imonth], cmonths_late[imonth]]:
        imo = int(cmo)-1
        if iyear%4 == 0:
            ndays = ndaysomo_leap[imo]
        else:
            ndays = ndaysomo[imo]
        cyear = str(iyear)    
        infile = master_directory + cyear + cmo + \
            '_ccpa_on_ndfd_grid_6hourly.nc'
        nc = Dataset(infile)
        yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
        for iday in range(1,ndays+1):
            if iday < 10:
                cday = '0'+str(iday)
            else:
                cday = str(iday)
            iyyyymmddhh = int(str(iyear)+cmo+cday+cend_hour)
            idx = np.where(yyyymmddhh_end == iyyyymmddhh)[0]
            precip_in = np.squeeze(nc.variables['apcp_anal'][idx,jndfd,indfd])
            if iyear == 2002 and iday == 1 and cmo == cmonths_middle[imonth]:
                precip_tseries = np.zeros((nsamps), \
                    dtype=np.float64)
                missingv = -99.99
                lons = nc.variables['lons'][:,:]
                lats = nc.variables['lats'][:,:]
            precip_in = np.where(precip_in < 500., precip_in, missingv)
            precip_tseries[ktr] = precip_in
            ktr = ktr+1
        nc.close()

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
empirical_cdf = 1.0 / (2.0*nz) + np.arange(nz)/nz  # "Hazen" plotting function
hazard_fn_in = -np.log(1.0 - empirical_cdf) 
empirical_cdf_withzeros = fraction_zero_ndfd[jndfd, indfd] + \
    (1.0-fraction_zero_ndfd[jndfd, indfd])*empirical_cdf

# ---- for diagnostic purposes, also generate a spline fit to CDF
#      with many more percentiles.

nz5000 = 5000
empirical_cdf_hires = 1.0 / (2.0*nz5000) + np.arange(nz5000)/nz5000
hazard_fn_hires = -np.log(1.0 - empirical_cdf_hires) 

# ---- get the spline fit.

numknots = number_knots[jndfd, indfd]  # now this number includes boundary
knots_in = spline_info_inv[jndfd, indfd,0,0:numknots]
bspline_coef_in = spline_info_inv[jndfd, indfd,1,0:numknots]
interior_knots_chf=knots_in[4:13]
precip_at_knots = np.zeros((9), dtype=np.float32)
cdf_at_knots = np.zeros((9), dtype=np.float32)
haz_at_knots = np.zeros((9), dtype=np.float32)
for i in range(9):
    idx = find_nearest(hazard_fn_hires,interior_knots_chf[i])
    haz = hazard_fn_hires[idx]
    haz_at_knots[i] = haz
    cd = 1.-np.exp(-haz)
    cdf_at_knots[i] = cd
    idx = find_nearest(empirical_cdf, cd)
    precip_at_knots[i] = precip_ens_nonzero[idx]
    
nthree = 3
ier = 0
quantile_mapped_precip = np.zeros((nz), dtype=np.float64)
quantile_mapped_precip_hires = np.zeros((nz5000), dtype=np.float64)

ier, quantile_mapped_precip = control_splev(knots_in, \
    bspline_coef_in, hazard_fn_in, numknots, nz)
    
ier, quantile_mapped_precip_hires = control_splev(knots_in, \
    bspline_coef_in, hazard_fn_hires, numknots, nz5000)

# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('Spline fit to analysis positive precipitation values,\n'+\
    cmonth+', '+cend_hour+' UTC, '+\
    clon+' W '+clat+' N',fontsize=13)

ax.plot(precip_ens_nonzero, empirical_cdf_withzeros,\
    color='Red',lw=2,label='Empirical')
ax.plot(quantile_mapped_precip, empirical_cdf_withzeros,\
    color='Blue',lw=2,label='Cubic spline')

plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.0,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_analysis_example_spline.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 
              
              
              
# ---- plot the hazard function, empirical and best fitted.

f = plt.figure(figsize=(6.5,6.5))#
ax = f.add_axes([.13,.14,.84,.75])
ax.set_title('Spline fit to hazard function of positive precip values,\n'+\
    cmonth+', '+cend_hour+' UTC, '+\
    clon+' W '+clat+' N',fontsize=13)

ax.plot(precip_ens_nonzero, hazard_fn_in,\
    color='Red',lw=1.5,label='Empirical')
ax.plot(quantile_mapped_precip, hazard_fn_in,\
    color='Blue',lw=1.5,label='Cubic spline')

plt.ylabel('Cumulative hazard function',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.0,10.0)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'hazard_precip_analysis_example_spline.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()               
                
                
                
                
# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(8,4.3))#

ax = f.add_axes([.07,.14,.41,.75])
ax.set_title('(a) Analyzed cumulative hazard function',fontsize=12.5)
ax.plot(precip_ens_nonzero, hazard_fn_in,'--',\
    color='Red',lw=1.6,label='Empirical')
ax.plot(quantile_mapped_precip, hazard_fn_in,\
    color='RoyalBlue',lw=1.6,label='Cubic spline')
ax.plot(precip_at_knots, haz_at_knots, 'k.')     
plt.ylabel('Cumulative hazard function',fontsize=9)
ax.legend(loc=0)
ax.set_ylim(0.0,10.0)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=9)


ax = f.add_axes([.57,.14,.41,.75])
ax.set_title('(b) Analyzed CDF',fontsize=12.5)
ax.plot(precip_ens_nonzero, empirical_cdf,'--',\
    color='Red',lw=1.6,label='Empirical')
ax.plot(quantile_mapped_precip, empirical_cdf,\
    color='RoyalBlue',lw=1.6,label='Cubic spline')
ax.plot(precip_at_knots, cdf_at_knots, 'k.')     
    
plt.ylabel('Non-exceedance probability',fontsize=9)
ax.legend(loc=0)
ax.set_ylim(0.0,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=9)
figname = 'CDF_Ithaca_example_00UTC_Aug.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 
              
# --- write to cPickle file
                
outfile = 'analysis_data_cdf.cPick'
ouf = open(outfile, 'wb')
cPickle.dump(precip_ens_nonzero, ouf)
cPickle.dump(hazard_fn_in, ouf)
cPickle.dump(quantile_mapped_precip, ouf)
cPickle.dump(empirical_cdf, ouf)
cPickle.dump(precip_at_knots, ouf)
cPickle.dump(haz_at_knots, ouf)
cPickle.dump(cdf_at_knots, ouf)
cPickle.dump
ouf.close()








