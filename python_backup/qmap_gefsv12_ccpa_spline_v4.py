"""
qmap_gefsv12_ccpa_spline_v4.py cyyyymmddhh clead cmem

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
import pygrib
from quantile_mapping_gamma_mixture_v2_f90 import \
    quantile_mapping_gamma_mixture_v2_f90
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import LSQUnivariateSpline, splrep, splev

import _pickle as cPickle
import scipy.stats as stats


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

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

def fraczero_possamps(nsamps, precip):
    """
    from the vector input sample precip_ens, define the fraction of
    samples with zero precipitation.   For the positive samples, add
    a small random number to deal with the fact that the data was 
    discretized to 0.1 mm, so that when later creating CDFs we don't 
    have values with lots of tied amounts.   Sort the nonzero amounts 
    and return.
    """
    number_zeros = 0
    precip_nonzero = np.delete(precip, \
        np.where(precip <= 0.0))  # censor at 0.1 mm
    nz = len(precip_nonzero)
    # data discretized, so add random component of this magnitude
    precip_nonzero = precip_nonzero + \
        np.random.uniform(low=-0.005,high=0.005,size=nz) 
    precip_nonzero = np.sort(precip_nonzero)  
    #print (precip_ens_nonzero[0:10]) 
    ntotal = len(precip)
    nzero = ntotal - len(precip_nonzero)
    fraction_zero = float(nzero) / float(ntotal)
    return fraction_zero, precip_nonzero, nz


# =====================================================================

# --- read grib data on a single level

def read_gribdata(gribfilename, endStep):
    istat = -1
    fexist_grib = False
    fexist_grib = os.path.exists(gribfilename)
    print (gribfilename, endStep)
    if fexist_grib:
        try:
            fcstfile = pygrib.open(gribfilename)
            grb = fcstfile.select(shortName='tp',endStep=endStep)[0]
            precip_realtime = grb.values
            lats_full, lons_full = grb.latlons()
            istat = 0
            fcstfile.close()
        except IOError:
            print ('   IOError in read_gribdata reading ', \
                gribfilename, validityDate, validityTime)
            istat = -1
        except ValueError:
            print ('   ValueError in read_gribdata reading ', \
                gribfilename, validityDate, validityTime)
            istat = -1
        except RuntimeError:
            print ('   RuntimeError in read_gribdata reading ', \
                gribfilename, validityDate, validityTime)
            istat = -1
    return istat, precip_realtime, lats_full, lons_full
    
# =====================================================================

def get_qmapped_precip(gefs_quantile, precip_amount, jy, ix, \
    spline_info_ndfd_inv, indices_to_query_ndfd,\
    fraction_zero_ndfd):

    """ this gets the analyzed precipitation associated with
        this quantile in the distribution by using (for very 
        light precipitation climatologies) the inverse CDF
        (percent point function) of the Gamma distribution, and
        for heavier amounts using a spline fitted precipitation 
        based on the Hazard function from 
        https://doi.org/10.1175/MWR-D-20-0096.1
    """

    if precip_amount[jy,ix] == 0.0:

        # ---- arbitrarile assign the CDF to zero if precip is zero.
        
        qmp = 0.0
    elif gefs_quantile[jy,ix] < fraction_zero_ndfd[jy,ix]:
        qmp = 0.0
    else:
        qpositive = (gefs_quantile[jy,ix] - fraction_zero_ndfd[jy,ix])/ \
            (1.0 - fraction_zero_ndfd[jy,ix])
        if indices_to_query_ndfd[jy,ix,0] == -1: 

            # ---- flagged as a dry point that estimated CDF with a Gamma.

            alpha = spline_info_ndfd[jy,ix,0,0] # previously stored here
            beta = spline_info_ndfd[jy,ix,1,0] # previously stored here
            scale = 1.0 / beta
            qmp = stats.gamma.ppf(qpositive, alpha, \
                loc=0, scale=scale)
        else:
                
            # ---- flagged as a wet-enough point to estimate the CDF with 
            #      the spline fit to a hazard function. 
            
            splines_tuple_inv = (spline_info_ndfd_inv[jy,ix,0,:], \
                spline_info_ndfd_inv[jy,ix,1,:], 3)
            hazard_fn = -np.log(1.0 - qpositive)    
            qmp = splev(hazard_fn, splines_tuple_inv)

    return qmp

# =====================================================================

def get_quantile_gefsv12(precip_amount, jy, ix, \
    spline_info_gefsv12, fraction_zero_gefsv12):

    """ this gets the quantile associated with a given precipitation 
    amount for GEFSv12 data, this month and 6-hour period. """
    
    if precip_amount[jy,ix] == 0.0:
        
        # ---- arbitrarile assign the CDF to zero if precip is zero.
        
	    quantile = 0.0
    else:	
	
	    # ---- flagged as a wet-enough point to estimate the CDF with 
        #      the spline fit to a hazard function. 
            
        splines_tuple = (spline_info_gefsv12[jy,ix,0,:], \
            spline_info_gefsv12[jy,ix,1,:], 3)
        spline_hazard = splev(precip_amount[jy,ix], splines_tuple)
        spline_cdf = 1.0 - np.exp(-spline_hazard)
        quantile = fraction_zero_gefsv12[jy,ix] + \
            (1.0 - fraction_zero_gefsv12[jy,ix])*spline_cdf

    return quantile
    
# =====================================================================
    
def get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead, jmin, \
    jmax, imin, imax):

    # ---- get the desired 2021 GEFSv12 forecast as grib file downloaded
    #      from NOMADS server. Then subset to big domain encompassing
    #      all NDFD domains, and subset again to processing domain of
    #      interest.

    nib1 = 518 # ~lon 220E
    nie1 = 1440 # up to ~lon 310E
    nib2 = 0 # ~lon 220E
    nie2 = 45 # up to ~lon 310E
    njb = 38 # lat ~ 80.5
    nje = 483 # down to lat ~ -30.5
    nj = nje - njb 
    ni1 = nie1 - nib1
    ni2 = nie2 - nib2
    ni = ni1 + ni2 
    
    # ---- get the whole globe.
    
    input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
    if int(clead) > 100:
        infile = input_directory + cyyyymmddhh + \
            '_ge'+cmem+'.t00z.pgrb2s.0p25.f' + clead 
    else:
        infile = input_directory + cyyyymmddhh + \
            '_ge'+cmem+'.t00z.pgrb2s.0p25.f0' + clead
    print ('reading from ', infile)
    endStep = int(clead)
    istat, precip_realtime, lats_full, lons_full = \
        read_gribdata(infile, endStep)
    
    # ---- subset for the big grid encompassing all NBM domains.
    #      The eastern boundary crosses Greenwich meridian, so
    #      fiddle to make sure longitudes in ascending order.
    
    precip_realtime_biggrid = np.hstack((\
        precip_realtime[njb:nje,nib1:nie1], \
        precip_realtime[njb:nje,nib2:nie2]))
        
    print ('np.shape(lons_full) = ',np.shape(lons_full))
    lons_1D_biggrid = \
        np.hstack((lons_full[0,nib1:nie1], lons_full[0,nib2:nie2]))
    lons_1D_biggrid = np.where(lons_1D_biggrid > 90.0, \
        lons_1D_biggrid - 360., lons_1D_biggrid)
        
    lats_1D_biggrid = lats_full[njb:nje,0]
    
    # ---- further subset the data to the specific NBM domain

    precip_realtime = precip_realtime_biggrid[jmin:jmax, imin:imax]
    lons_1D_realtime = lons_1D_biggrid[imin:imax]
    lats_1D_realtime = lats_1D_biggrid[jmin:jmax]

    return precip_realtime, lons_1D_realtime, lats_1D_realtime
                
# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2]
cmem = sys.argv[3]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
jmin, jmax, imin, imax = set_domain_boundaries(cdomain)


# ---- read in the previously generated netCDF file with precipitation
#      to get lat/lon of NDFD grid

ccpa_directory = '/Volumes/Backup Plus/ccpa/'
infile = ccpa_directory + '200201_ccpa_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile)
lons_ndfd = nc.variables['lons'][:,:]
lons_ndfd = lons_ndfd
lats_ndfd = nc.variables['lats'][:,:]
ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
nc.close()
   
# --- load from cPickle file the spline parameters
#     for the GEFSv12 forecasts

master_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'
infile = master_directory + ccmonth+'_'+cdomain+\
    '_GEFSv12_spline_info_h' + clead + '.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
spline_info_gefsv12 = cPickle.load(inf)
fraction_zero_gefsv12 = cPickle.load(inf)
indices_to_query_gefsv12 = cPickle.load(inf)
empirical_precipvals_gefsv12 = cPickle.load(inf)
inf.close()
   
# --- load from cPickle file the spline/Gamma parameters
#     for the ccpa precipitation analyses
#	  np.shape(spline_info_ndfd) =  (1597, 2345, 2, 17) # CCPA on NDFD
#	  np.shape(indices_to_query_ndfd) =  (1597, 2345, 9)
   
infile = ccpa_directory + cmonth+'_'+cdomain+\
    '_CCPA_spline_info_h' + clead + '.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
spline_info_ndfd = cPickle.load(inf)
spline_info_inv = cPickle.load(inf)
fraction_zero_ndfd = cPickle.load(inf)
indices_to_query_ndfd = cPickle.load(inf)
inf.close()   

# ---- extract the real-time precipitation for this domain.
#      set GEFSv12 grid dimensions and 1/4-deg lon and lat

precip_realtime, lons_1d_realtime, lats_1d_realtime = \
    get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead, jmin, \
    jmax, imin, imax) 
nx_gefsv12 = len(lons_1d_realtime)
ny_gefsv12 = len(lats_1d_realtime)
lons_fcst_2d, lats_fcst_2d = \
    np.meshgrid(lons_1d_realtime,lats_1d_realtime)

# ---- flip upside down, as subsequent interpolation requires
#      S to N with increasing j index.  Then
#      interpolate precipitation to NDFD grid. 

precip_realtime_flipud = np.flipud(precip_realtime) 
lats_1d_realtime_flipud = np.flipud(lats_1d_realtime)
print ('interpolating GEFSv12 forecast of interest to NDFD grid')
precip_gefsv12_on_ndfd = interp(precip_realtime_flipud, \
    lons_1d_realtime, lats_1d_realtime_flipud,  \
    lons_ndfd, lats_ndfd, checkbounds=False, \
    masked=False, order=1)  
    
# ---- interpolate the GEFSv12 fraction zeros to NDFD grid
    
fraction_zero_gefsv12_flipped = np.flipud(fraction_zero_gefsv12)
fraction_zero_gefsv12_on_ndfd = interp(fraction_zero_gefsv12_flipped, \
    lons_1d_realtime, lats_1d_realtime_flipud, \
    lons_ndfd, lats_ndfd, checkbounds=False, \
    masked=False, order=1)
    
# ---- now loop over grid points and obtain the forecast quantile

gefsv12_quantiles = np.zeros((ny_gefsv12, nx_gefsv12), dtype=np.float64)
for jy in range(ny_gefsv12):
    for ix in range(nx_gefsv12):
        gefsv12_quantiles[jy,ix] = get_quantile_gefsv12(precip_realtime, jy, \
            ix, spline_info_gefsv12, fraction_zero_gefsv12)

# ---- interpolate the forecast quantile to the NDFD grid.

gefsv12_quantiles_flipud = np.flipud(gefsv12_quantiles) 
lats_1d_realtime_flipud = np.flipud(lats_1d_realtime)
print ('interpolating GEFSv12 forecast of interest to NDFD grid')
gefsv12_quantiles_on_ndfd = interp(gefsv12_quantiles_flipud, \
    lons_1d_realtime, lats_1d_realtime_flipud,  \
    lons_ndfd, lats_ndfd, checkbounds=False, \
    masked=False, order=1)  

# ---- apply quantile mapping procedure

qmapped_precip = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
now = datetime.now()
begin_time = now.strftime("%H:%M:%S")
for jy in range(ny_ndfd):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    if jy%10 == 0: print ('********** jy = ',jy,' of ',ny_ndfd,\
        '  elapsed time = ',current_time,' begin time = ',\
        begin_time) 
    for ix in range(nx_ndfd):

        # ---- get the quantile associated with today's GEFSv12 forecast
        #      through either use of Gamma fitting (very dry points)
        #      or spline fitting.

        qmapped_precip[jy,ix] = get_qmapped_precip(\
            gefsv12_quantiles_on_ndfd, \
            precip_gefsv12_on_ndfd, jy, ix, spline_info_inv, \
            indices_to_query_ndfd, fraction_zero_ndfd)
        

# ======================================================================


# ---- plot the quantile-mapped forecast.

m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ndfd, lats_ndfd)

clevs = [0.0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Quantile-mapped forecast, IC = '+cyyyymmddhh+\
    ' lead = '+clead+' h'
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, qmapped_precip, clevs,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'qmapped_precip_'+cyyyymmddhh+'_lead'+clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



# ---- first plot GEFS control forecast amount

m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ndfd, lats_ndfd)

# ---- make plots of fraction positive precip

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

clevs = [0.0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20]
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Interpolated control forecast precipitation amount (mm) for '+\
    clead+'-h IC = '+cyyyymmddhh
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, precip_gefsv12_on_ndfd, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Mean precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'forecast_precip_'+clead+'_h_IC'+cyyyymmddhh+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



# ---- first plot GEFS control forecast amount

m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ndfd, lats_ndfd)

# ---- make plots of fraction positive precip

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

clevs = [0.0,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.8,0.9,0.95]
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Interpolated GEFSv12 1-fraction zero for '+clead+\
    '-h IC = '+cyyyymmddhh
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, 1.-fraction_zero_gefsv12_on_ndfd, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.06,0.07,0.88,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Fraction zero',fontsize=9)

# ---- set plot title

plot_title = 'one_minus_fraczero_GEFSv12_'+clead+\
    '_h_IC'+cyyyymmddhh+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



