"""
qmap_gefsv12_mswep_spline_v2.py cyyyymmddhh clead

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
            #print ('opened ', gribfilename)
            #grb = fcstfile.select(shortName='tp',\
            #    validityDate=validityDate, \
            #    validityTime=validityTime)[0]
            grb = fcstfile.select(shortName='tp',endStep=endStep)[0]
            #print ('selected grb')
            precip_realtime = grb.values
            #print ('read values')
            lats_full, lons_full = grb.latlons()
            #print ('got lat/lon ',lats_full[0,0], lons_full[0,0])
            istat = 0
            fcstfile.close()
            #print ('closed file ')
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
    spline_info_mswep_inv, indices_to_query_mswep,\
    fraction_zero_mswep):

    """ this gets the analyzed precipitation associated with
        this quantile in the distribution by using (for very 
        light precipitation climatologies) the inverse CDF
        (percent point function) of the Gamma distribution, and
        for heavier amounts using a spline fitted precipitation 
        based on the Hazard function from 
        https://doi.org/10.1175/MWR-D-20-0096.1
    """

    if precip_amount == 0.0:

        # ---- arbitrarile assign the CDF to zero if precip is zero.
        
        qmp = 0.0
    elif gefs_quantile < fraction_zero_mswep[jy,ix]:
        qmp = 0.0
    else:
        qpositive = (gefs_quantile - fraction_zero_mswep[jy,ix])/ \
            (1.0 - fraction_zero_mswep[jy,ix])
        if indices_to_query_mswep[jy,ix,0] == -1: 

            # ---- flagged as a dry point that estimated CDF with a Gamma.

            alpha = spline_info_mswep[jy,ix,0,0] # previously stored here
            beta = spline_info_mswep[jy,ix,1,0] # previously stored here
            scale = 1.0 / beta
            qmp = stats.gamma.ppf(qpositive, alpha, \
                loc=0, scale=scale)
        else:
                
            # ---- flagged as a wet-enough point to estimate the CDF with 
            #      the spline fit to a hazard function. 
            
            splines_tuple_inv = (spline_info_mswep_inv[jy,ix,0,:], \
                spline_info_mswep_inv[jy,ix,1,:], 3)
            hazard_fn = -np.log(1.0 - qpositive)    
            qmp = splev(hazard_fn, splines_tuple_inv)

    return qmp

# =====================================================================

def get_quantile_gefsv12(precip_amount, jy, ix, \
    spline_info_gefsv12, fraction_zero_gefsv12):

    """ this gets the quantile associated with a given precipitation 
    amount for GEFSv12 data, this month and 6-hour period. """
    
    if precip_amount == 0.0:
        
        # ---- arbitrarile assign the CDF to zero if precip is zero.
        
	    quantile = 0.0
    else:	
	
	    # ---- flagged as a wet-enough point to estimate the CDF with 
        #      the spline fit to a hazard function. 
            
        splines_tuple = (spline_info_gefsv12[jy,ix,0,:], \
            spline_info_gefsv12[jy,ix,1,:], 3)
        spline_hazard = splev(precip_amount, splines_tuple)
        spline_cdf = 1.0 - np.exp(-spline_hazard)
        quantile = fraction_zero_gefsv12[jy,ix] + \
            (1.0 - fraction_zero_gefsv12[jy,ix])*spline_cdf

    return quantile

# =====================================================================

def get_bounding_indices(gefsv12_quantile, anal_quantiles):

    ihigh = np.searchsorted(anal_quantiles, gefsv12_quantile)
    ilow = ihigh-1
    return ilow, ihigh

# ======================================================================

# ---- inputs from command line

nstride = 1
cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
mswep_directory = '/Volumes/Backup Plus/mswep/'
#pamounts = [0.000001, 0.1, 0.5, 1.0, 3.0,  5.0, 10.0, 25.0, 50.0, 500.0]
conv_criterion = 0.001

# ---- read in the previously generated netCDF file with precipitation
#      to get lat/lon of MSWEP grid

infile = mswep_directory + '200001_on_ndfd_grid_6hourly.nc'
nc = Dataset(infile)
lons_mswep = nc.variables['lons'][:,:]
lons_mswep = lons_mswep - 360.0
lats_mswep = nc.variables['lats'][:,:]
ny_mswep, nx_mswep = np.shape(lons_mswep)
#print ('lats_mswep[0], [-1] = ', lats_mswep[0], lats_mswep[-1])
#print ('mswep lons dtype = ',lons_mswep.dtype)
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
#     for the mxswep forecasts
#	  np.shape(spline_info_mswep) =  (1597, 2345, 2, 17) # MSWEP
#	  np.shape(indices_to_query_mswep) =  (1597, 2345, 9)
   
mswep_directory = '/Volumes/Backup Plus/mswep/'
infile = mswep_directory + cmonth+'_'+cdomain+\
    '_MSWEP_spline_info_h' + clead + '.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
spline_info_mswep = cPickle.load(inf)
spline_info_inv = cPickle.load(inf)
fraction_zero_mswep = cPickle.load(inf)
indices_to_query_mswep = cPickle.load(inf)
inf.close()   

# ---- get the GEFSv12 domain subgrid lat/lon
#	   np.shape(spline_info) =  (153, 318, 2, 17) # GEFSv12
#	   np.shape(indices_to_query) =  (153, 318, 9)
   
jmin, jmax, imin, imax = set_domain_boundaries(cdomain)
ncfile = master_directory + ccmonth + '_apcp_h' + clead + '.nc'
nc = Dataset(ncfile)
print ('reading from ', infile)
lons_1d_gefsv12 = nc.variables['lons_fcst'][imin:imax]
lats_1d_gefsv12 = nc.variables['lats_fcst'][jmin:jmax]
lons_1d_gefsv12_full = nc.variables['lons_fcst'][:]
lats_1d_gefsv12_full = nc.variables['lats_fcst'][:]
nx_gefsv12 = len(lons_1d_gefsv12)
ny_gefsv12 = len(lats_1d_gefsv12)
nc.close()
lons_fcst_2d, lats_fcst_2d = \
    np.meshgrid(lons_1d_gefsv12,lats_1d_gefsv12)

# ---- get the desired 2021 GEFSv12 forecast as grib file downloaded
#      from NOMADS server. Flip upside down if arrays are not S to N.
#      This is needed for basemap.interp

input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
infile = input_directory + cyyyymmddhh + \
    '_gec00.t00z.pgrb2s.0p25.f0' + clead 
print ('reading from ', infile)
endStep = int(clead)
istat, precip_realtime, lats_full, lons_full = \
    read_gribdata(infile, endStep)
if lats_full[0,0] > lats_full[-1,0]: 
    lats_full = np.flipud(lats_full)
    precip_realtime = np.flipud(precip_realtime)  
lons_full = lons_full - 360.
    
# ---- interpolate precip to MSWEP on NDFD grid.   

print ('interpolating forecast of interest to MSWEP NDFD grid')
precip_gefsv12_on_mswep = interp(precip_realtime, \
    lons_full[0,:], lats_full[:,0], \
    lons_mswep, lats_mswep, checkbounds=False, \
    masked=False, order=1)  
    
fraction_zero_flipped = np.flipud(fraction_zero_gefsv12)
fraction_zero_gefsv12_on_mswep = interp(np.flipud(fraction_zero_gefsv12), \
    lons_1d_gefsv12, np.flipud(lats_1d_gefsv12), \
    lons_mswep, lats_mswep, checkbounds=False, \
    masked=False, order=1)
    
# ---- determine the GEFSv12 indices of closest grid point to each MSWEP point

readclosest = True
iclose = np.zeros((ny_mswep, nx_mswep), dtype=np.int)
jclose = np.zeros((ny_mswep, nx_mswep), dtype=np.int)
if readclosest == False:
    print ('Determining closest GEFSv12 grid points')
    for jy in range(ny_mswep):
        if jy%10 == 0: print('Processing jy ',jy,' of ',ny_mswep)
        for ix in range(nx_mswep):
            rlon = lons_mswep[jy,ix]
            frac = (rlon - lons_1d_gefsv12[0]) / \
                (lons_1d_gefsv12[-1] - lons_1d_gefsv12[0])
            idxx = int(nx_gefsv12*frac)
            if idxx < 0 or idxx > nx_gefsv12-1:
                if idxx < 0: idxx = 0
                if idxx > nx_gefsv12-1: idxx = nx_gefsv12-1
            iclose[jy,ix] = idxx

            rlat = lats_mswep[jy,ix]
            frac = (rlat-lats_1d_gefsv12[0]) / \
                (lats_1d_gefsv12[-1] - lats_1d_gefsv12[0])
            idxy = int(ny_gefsv12*frac)

            if idxy < 0 or idxy > ny_gefsv12-1:
                if idxy < 0: idxy = 0
                if idxy > ny_gefsv12-1: idxy = ny_gefsv12-1
            jclose[jy,ix] = idxy

    outfile = 'ijclosest.cPick'
    ouf = open(outfile, 'wb')
    cPickle.dump(iclose, ouf)
    cPickle.dump(jclose, ouf)
    ouf.close()
else:
    infile = 'ijclosest.cPick'
    inf = open(infile, 'rb')
    iclose = cPickle.load(inf)
    jclose = cPickle.load(inf)
    inf.close()

# ---- now loop over grid points and quantile map each

anal_quantiles = np.zeros((10), dtype=np.float64)
qmapped_precip = np.zeros((ny_mswep, nx_mswep), dtype=np.float64)
now = datetime.now()
begin_time = now.strftime("%H:%M:%S")
for jy in range(ny_mswep):
#for jy in range(291,292):
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    if jy%10 == 0: print ('********** jy = ',jy,' of ',ny_mswep,\
        '  elapsed time = ',current_time,' begin time = ',\
        begin_time) 
    #for ix in range(1183,1184):
    for ix in range(nx_mswep):

        # ---- get the quantile associated with today's GEFSv12 forecast
        #      at the nearest GEFSv12 grid point
        
        ixgefs = iclose[jy,ix]
        jygefs = jclose[jy,ix]
        pgefs = precip_gefsv12_on_mswep[jy,ix]
        gefsv12_quantile = get_quantile_gefsv12(pgefs, jygefs, \
            ixgefs, spline_info_gefsv12, fraction_zero_gefsv12)
        qmapped_precip[jy,ix] = get_qmapped_precip(gefsv12_quantile, \
            pgefs, jy, ix, spline_info_inv, indices_to_query_mswep, \
            fraction_zero_mswep)          
        
# ---- plot the quantile-mapped forecast.

m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_mswep, lats_mswep)

clevs = [0.0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Quantile mapped forecast, IC = '+cyyyymmddhh+' lead = '+clead+' h'
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
x, y = m(lons_mswep, lats_mswep)

# ---- make plots of fraction positive precip

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

clevs = [0.0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20]
fig = plt.figure(figsize=(9,6.5))
axloc = [0.02,0.1,0.96,0.84]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = 'Interpolated control forecast precipitation amount (mm) for '+clead+\
    '-h IC = '+cyyyymmddhh
ax1.set_title(title, fontsize=14,color='Black')
CS2 = m.contourf(x, y, precip_gefsv12_on_mswep, clevs,\
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
x, y = m(lons_mswep, lats_mswep)

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
CS2 = m.contourf(x, y, 1.-fraction_zero_gefsv12_on_mswep, clevs,\
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

plot_title = 'one_minus_fraczero_GEFSv12_'+clead+'_h_IC'+cyyyymmddhh+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



