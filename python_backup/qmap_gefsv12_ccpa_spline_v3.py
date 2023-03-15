"""
qmap_gefsv12_ccpa_spline_v2.py cyyyymmddhh clead cmem

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

    if precip_amount == 0.0:

        # ---- arbitrarile assign the CDF to zero if precip is zero.
        
        qmp = 0.0
    elif gefs_quantile < fraction_zero_ndfd[jy,ix]:
        qmp = 0.0
    else:
        qpositive = (gefs_quantile - fraction_zero_ndfd[jy,ix])/ \
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

# ======================================================================

def bilinear_interpolate_qmapped(jy,ix, jlow, ilow, jhigh, ihigh, \
    yfrac, xfrac, precip_realtime, spline_info_gefsv12, \
    fraction_zero_gefsv12, spline_info_inv, indices_to_query_ndfd, \
    fraction_zero_ndfd):


    """ perform a quantile mapping for each surrounding GEFSv12 point
    and bilinearly interpolate """
    
            
    pgefs_jminus_iminus = precip_realtime[jlow[jy,ix],ilow[jy,ix]]
    pgefs_jplus_iminus = precip_realtime[jhigh[jy,ix],ilow[jy,ix]]
    pgefs_jminus_iplus = precip_realtime[jlow[jy,ix],ihigh[jy,ix]]
    pgefs_jplus_iplus = precip_realtime[jhigh[jy,ix],ihigh[jy,ix]]
        
    gefsv12_quantile_jminus_iminus = get_quantile_gefsv12(pgefs_jminus_iminus, jlow[jy,ix], \
        ilow[jy,ix], spline_info_gefsv12, fraction_zero_gefsv12)
    gefsv12_quantile_jplus_iminus = get_quantile_gefsv12(pgefs_jplus_iminus, jhigh[jy,ix], \
        ilow[jy,ix], spline_info_gefsv12, fraction_zero_gefsv12)
    gefsv12_quantile_jminus_iplus = get_quantile_gefsv12(pgefs_jminus_iplus, jlow[jy,ix], \
        ihigh[jy,ix], spline_info_gefsv12, fraction_zero_gefsv12)
    gefsv12_quantile_jplus_iplus = get_quantile_gefsv12(pgefs_jplus_iplus, jhigh[jy,ix], \
        ihigh[jy,ix], spline_info_gefsv12, fraction_zero_gefsv12)
                
    qm_jminus_iminus = get_qmapped_precip(gefsv12_quantile_jminus_iminus, \
        pgefs_jminus_iminus, jy, ix, spline_info_inv, indices_to_query_ndfd, \
        fraction_zero_ndfd)
    qm_jplus_iminus = get_qmapped_precip(gefsv12_quantile_jplus_iminus, \
        pgefs_jplus_iminus, jy, ix, spline_info_inv, indices_to_query_ndfd, \
        fraction_zero_ndfd)
    qm_jminus_iplus = get_qmapped_precip(gefsv12_quantile_jminus_iplus, \
        pgefs_jminus_iplus, jy, ix, spline_info_inv, indices_to_query_ndfd, \
        fraction_zero_ndfd)
    qm_jplus_iplus = get_qmapped_precip(gefsv12_quantile_jplus_iplus, \
        pgefs_jplus_iplus, jy, ix, spline_info_inv, indices_to_query_ndfd, \
        fraction_zero_ndfd)                
                
    qmp = (1.0 -xfrac[jy,ix])*(1.0-yfrac[jy,ix])*qm_jminus_iminus + \
        xfrac[jy,ix]*(1.0-yfrac[jy,ix])*qm_jminus_iplus + \
        (1.0 -xfrac[jy,ix])*yfrac[jy,ix]*qm_jplus_iminus + \
        xfrac[jy,ix]*yfrac[jy,ix]*qm_jplus_iplus 
        
    return qmp

# =====================================================================

def find_nearest(vec, value):
    
    """ given a vector vec and a particular value, find the index in vec
        that is nearest to value"""
    
    idx = np.abs(vec-value).argmin()
    return idx
                
# =====================================================================
# =====================================================================


# ---- inputs from command line

nstride = 1
cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2]
cmem = sys.argv[3]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
jmin, jmax, imin, imax = set_domain_boundaries(cdomain)

#pamounts = [0.000001, 0.1, 0.5, 1.0, 3.0,  5.0, 10.0, 25.0, 50.0, 500.0]
conv_criterion = 0.001

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
fraction_zero_gefsv12_on_ccpa = interp(fraction_zero_gefsv12_flipped, \
    lons_1d_realtime, lats_1d_realtime_flipud, \
    lons_ndfd, lats_ndfd, checkbounds=False, \
    masked=False, order=1)
    
# ---- determine the GEFSv12 indices of closest grid point 
#      to each NDFD grid point

readclosest = True
iclose = np.zeros((ny_ndfd, nx_ndfd), dtype=np.int)
jclose = np.zeros((ny_ndfd, nx_ndfd), dtype=np.int)
ilow = np.zeros((ny_ndfd, nx_ndfd), dtype=np.int)
ihigh = np.zeros((ny_ndfd, nx_ndfd), dtype=np.int)
jlow = np.zeros((ny_ndfd, nx_ndfd), dtype=np.int)
jhigh = np.zeros((ny_ndfd, nx_ndfd), dtype=np.int)
xfrac = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
yfrac = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
if readclosest == False:
    print ('Determining closest GEFSv12 grid points')
    for jy in range(ny_ndfd):
    #for jy in range(ny_ndfd//2, ny_ndfd//2+1):
        if jy%10 == 0: print('Processing jy ',jy,' of ',ny_ndfd)
        for ix in range(nx_ndfd):
        #for ix in range(nx_ndfd//2, nx_ndfd//2+1):
            
            # ---- determine the indices of the bounding longitude box,
            #      and the fraction from lower to upper index in x dir
            rlon = lons_ndfd[jy,ix]
            inear = find_nearest(lons_1d_realtime, rlon)
            
            #print ('rlon, inear, nearest gefs = ', rlon, inear, lons_1d_realtime[inear])
            if np.abs(lons_1d_realtime[inear] > rlon):
                ilowg = max([0,inear-1])
                ihighg = min([inear,nx_gefsv12-1])
            else: # np.abs(lons_1d_realtime[inear] <= rlon):
                ilowg = max([0,inear])
                ihighg = min([inear+1,nx_gefsv12-1])
            #print (ilowg, ihighg)
            #print ('ilowg, ihighg, bounding lons = ', ilowg, ihighg, \
            #    lons_1d_realtime[ilowg], lons_1d_realtime[ihighg])          
            
            if ihighg != ilowg:
                xfrac[jy,ix] = (rlon - lons_1d_realtime[ilowg]) / \
                    (lons_1d_realtime[ihighg] - lons_1d_realtime[ilowg])
            else:
                xfrac[jy,ix] = 1.0 # at boundary
            iclose[jy,ix] = inear
            ilow[jy,ix] = ilowg
            ihigh[jy,ix] = ihighg
            #print ('xfrac = ', xfrac[jy,ix])
            
            
            # --- determine the indices of the bounding latitude box,
            #     and the fraction from lower to upper index in y dir.
            #     indexing is a bit tricky below as grid is oriented 
            #     N to S (90 to -90).
            
            rlat = lats_ndfd[jy,ix]
            jnear = find_nearest(lats_1d_realtime, rlat)
            
            #print ('rlat, jnear, nearest gefs = ', rlat, jnear, lats_1d_realtime[jnear])
            if np.abs(lats_1d_realtime[jnear] > rlat):
                jlowg = jnear
                jhighg = min([jnear+1, ny_gefsv12-1])
            else: # np.abs(lats_1d_realtime[jnear] <= rlat):
                jlowg = max([0,jnear-1])
                jhighg = min([jnear, ny_gefsv12-1])     
            #print ('jlowg, jhighg, bounding lats = ', jlowg, jhighg, \
            #    lats_1d_realtime[jlowg], lats_1d_realtime[jhighg])   
            
            if jhighg != jlowg:
                yfrac[jy,ix] = (rlat - lats_1d_realtime[jlowg]) / \
                    (lats_1d_realtime[jhighg] - lats_1d_realtime[jlowg])
            else:
                yfrac[jy,ix] = 1.0 # at boundary
        
            jclose[jy,ix] = jnear
            jlow[jy,ix] = jlowg
            jhigh[jy,ix] = jhighg
            #print ('yfrac = ', yfrac[jy,ix])
            
            #sys.exit()
            

    outfile = 'ijclosest.cPick'
    ouf = open(outfile, 'wb')
    cPickle.dump(iclose, ouf)
    cPickle.dump(jclose, ouf)
    cPickle.dump(ilow, ouf)
    cPickle.dump(jlow, ouf)
    cPickle.dump(ihigh, ouf)
    cPickle.dump(jhigh, ouf)
    cPickle.dump(xfrac, ouf)
    cPickle.dump(yfrac, ouf)    
    ouf.close()
else:
    infile = 'ijclosest.cPick'
    inf = open(infile, 'rb')
    iclose = cPickle.load(inf)
    jclose = cPickle.load(inf)
    ilow = cPickle.load(inf)
    jlow = cPickle.load(inf)
    ihigh = cPickle.load(inf)
    jhigh = cPickle.load(inf)
    xfrac = cPickle.load(inf)
    yfrac = cPickle.load(inf)
    inf.close()
    
# ---- now loop over grid points and quantile map each

anal_quantiles = np.zeros((10), dtype=np.float64)
qmapped_precip = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
now = datetime.now()
begin_time = now.strftime("%H:%M:%S")
interpolate = True
for jy in range(ny_ndfd):
#for jy in range(7*ny_ndfd//8,7*ny_ndfd//8+1):
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    if jy%10 == 0: print ('********** jy = ',jy,' of ',ny_ndfd,\
        '  elapsed time = ',current_time,' begin time = ',\
        begin_time) 
    #for ix in range(nx_ndfd//8, nx_ndfd//8+1):
    for ix in range(nx_ndfd):

        # ---- get the quantile associated with today's GEFSv12 forecast
        #      at the nearest GEFSv12 grid point
        
        
        if interpolate == False:
            ixgefs = iclose[jy,ix]
            jygefs = jclose[jy,ix]
            pgefs = precip_gefsv12_on_ndfd[jy,ix]
            gefsv12_quantile = get_quantile_gefsv12(pgefs, jygefs, \
                ixgefs, spline_info_gefsv12, fraction_zero_gefsv12)
            qmapped_precip[jy,ix] = get_qmapped_precip(gefsv12_quantile, \
                pgefs, jy, ix, spline_info_inv, indices_to_query_ndfd, \
                fraction_zero_ndfd)
        else:
            
            ixgefs = iclose[jy,ix]
            jygefs = jclose[jy,ix]
            #print ('NDFD jy, ix, lat, lon = ', jy, ix, lats_ndfd[jy,ix], lons_ndfd[jy,ix])
            #print ('GEFSv12 closest [lower] jy,ix index, lat, lon: ',\
            #    jygefs, ixgefs, lats_fcst_2d[jygefs,ixgefs], lons_fcst_2d[jygefs,ixgefs])
            #print ('surrounding latitudes: lats_fcst_2d[jygefs-1: jygefs+2,ixgefs] = ', \
            #    lats_fcst_2d[jygefs-1:jygefs+2,ixgefs])
            #print ('surrounding longitude: lons_fcst_2d[jygefs,ixgefs-1:ixgefs+2] = ', \
            #    lons_fcst_2d[jygefs,ixgefs-1:ixgefs+2])
            
            qmapped_precip[jy,ix] = bilinear_interpolate_qmapped(\
                jy,ix,jlow,ilow,jhigh,ihigh, yfrac, xfrac, precip_realtime,\
                spline_info_gefsv12, fraction_zero_gefsv12, \
                spline_info_inv, indices_to_query_ndfd, \
                fraction_zero_ndfd)
                
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
title = 'Quantile-mapped forecast, IC = '+cyyyymmddhh+' lead = '+clead+' h'
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
title = 'Interpolated control forecast precipitation amount (mm) for '+clead+\
    '-h IC = '+cyyyymmddhh
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
CS2 = m.contourf(x, y, 1.-fraction_zero_gefsv12_on_ccpa, clevs,\
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



