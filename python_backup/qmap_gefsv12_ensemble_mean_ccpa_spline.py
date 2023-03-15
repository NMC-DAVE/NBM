"""
qmap_gefsv12_ensemble_mean_ccpa_spline.py cyyyymmddhh clead 

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
from qmapping_spline import qmapping_spline
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
    #print (gribfilename, endStep)
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
                gribfilename, endStep)
            istat = -1
        except ValueError:
            print ('   ValueError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1
        except RuntimeError:
            print ('   RuntimeError in read_gribdata reading ', \
                gribfilename, endStep)
            istat = -1
    return istat, precip_realtime, lats_full, lons_full
    
# =====================================================================

def get_qmapped_precip(gefs_quantile, precip_amount, jy, ix, \
    spline_info_ndfd_inv, fraction_zero_ndfd):

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
    spline_info_gefsv12, fraction_zero_gefsv12, usegamma_gefsv12):

    """ this gets the quantile associated with a given precipitation 
    amount for GEFSv12 data, this month and 6-hour period. """
    
    if precip_amount[jy,ix] == 0.0:
        
        # ---- arbitrarile assign the CDF to zero if precip is zero.
        
	    quantile = 0.0
    else:	
	
        if usegamma_gefsv12[jy,ix] == 0:
            
	        # ---- flagged as a wet-enough point to estimate the CDF with 
            #      the spline fit to a hazard function. 
            
            splines_tuple = (spline_info_gefsv12[:,0,jy,ix], \
                spline_info_gefsv12[:,1,jy,ix], 3)
            spline_hazard = splev(precip_amount[jy,ix], splines_tuple)
            spline_cdf = 1.0 - np.exp(-spline_hazard)
            quantile = fraction_zero_gefsv12[jy,ix] + \
                (1.0 - fraction_zero_gefsv12[jy,ix])*spline_cdf
        else:
            if usegamma_gefsv12[jy,ix] = -1:
                # --- flagged as basically no training data.
                quantile = 0.0
            else:  # --- flagged as minimal training data - use Gamma
                alpha_hat = spline_info_fcst[0,0,jy,ix] 
                beta_hat = spline_info_fcst[0,1,jy,ix] 
                y0 = precip_amount[jy,ix] / beta_hat
                quantile = stats.gamma.cdf(y0, alpha_hat)

    return quantile
    
# =====================================================================
    
def get_domain_subset_of_gefs(cyyyymmddhh, group, clead, jmin, \
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
    
    
    isub = imax - imin 
    jsub = jmax - jmin
    precip_mean = np.zeros((jsub, isub), dtype=np.float64)
    print ('jsub, isub = ', jsub, isub)
    first = True
    for cmem in group:
    
        # ---- get the whole globe.
    
        input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
        if int(clead) > 100:
            infile = input_directory + cyyyymmddhh + \
                '_ge'+cmem+'.t00z.pgrb2s.0p25.f' + clead 
        else:
            infile = input_directory + cyyyymmddhh + \
                '_ge'+cmem+'.t00z.pgrb2s.0p25.f0' + clead
        print ('  reading from ', infile)
        endStep = int(clead)
        #print ('endStep = ', endStep)
        istat, precip_realtime, lats_full, lons_full = \
            read_gribdata(infile, endStep)
    
        # ---- subset for the big grid encompassing all NBM domains
        #      including Gaum, Hawaii, CONUS, AK, Puerto Rico.
        #      The eastern boundary crosses Greenwich meridian, so
        #      fiddle to make sure longitudes in ascending order.
    
        precip_realtime_biggrid = np.hstack((\
            precip_realtime[njb:nje,nib1:nie1], \
            precip_realtime[njb:nje,nib2:nie2]))
        
        if first == True:
            lons_1D_biggrid = \
                np.hstack((lons_full[0,nib1:nie1], lons_full[0,nib2:nie2]))
            lons_1D_biggrid = np.where(lons_1D_biggrid > 90.0, \
                lons_1D_biggrid - 360., lons_1D_biggrid)
            lats_1D_biggrid = lats_full[njb:nje,0]
            first = False
    
        # ---- further subset the data to the specific NBM domain

        precip_realtime = precip_realtime_biggrid[jmin:jmax, imin:imax]
        lons_1D_realtime = lons_1D_biggrid[imin:imax]
        lats_1D_realtime = lats_1D_biggrid[jmin:jmax]
        precip_mean = precip_mean + precip_realtime
        
    precip_mean = precip_mean / float(len(group))

    return precip_mean, lons_1D_realtime, lats_1D_realtime
                
# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2] # 3-digit number here
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
cdomain = 'conus'
jmin, jmax, imin, imax = set_domain_boundaries(cdomain)

groups = [['c00','p01', 'p02','p03','p04'],['p05','p06','p07','p08','p09'],\
    ['p10','p11', 'p12','p13','p14'],['p15','p16','p17','p18','p19'],\
    ['p20','p21', 'p22','p23','p24'],['p25','p26','p27','p28','p29','p30']]
groupnames = ['c00_to_p04', 'p05_to_p09','p10_to_p14',\
    'p15_to_p20','p20_to_p24','p24_to_p30']

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
   
# ---- load spline information for forecast netCDF file.  
#      Should be same lat/lon as real-time, but check

master_directory_forecast_spline = \
    '/Volumes/NBM/'+cdomain+'_gefsv12/CDF_spline/'
infile = master_directory_forecast_spline + \
    cmonth + '_' + cdomain + \
    '_GEFSv12_mean_spline_info_h' + clead + '.nc' 
nc = Dataset(infile)    
lons_spline_gefsv12_1d = nc.variables['lons'][:]
print ('lons_spline_gefsv12_1d = ', lons_spline_gefsv12_1d[0:-1:4])
lats_spline_gefsv12_1d = nc.variables['lats'][:]
print ('lats_spline_gefsv12_1d = ', lats_spline_gefsv12_1d[0:-1:4])
spline_info_gefsv12 = nc.variables['spline_info'][:,:,:,:]
fraction_zero_gefsv12 = nc.variables['fzero'][:,:]
usegamma_gefsv12 = nc.variables['usegamma'][:,:]
ny_gefsv12, nx_gefsv12 = np.shape(usegamma_gefsv12)
nc.close()
   
# --- load from netCDF file the spline/Gamma parameters
#     for the ccpa precipitation analyses
#	  np.shape(spline_info_ndfd) =  (1597, 2345, 2, 17) # CCPA on NDFD
#	  np.shape(indices_to_query_ndfd) =  (1597, 2345, 9)
   
if int(clead) > 18:
    ndays = int(clead) // 24
    ilead = int(clead)-ndays*24
    if ilead == 0:
        cleada = '00'
    elif ilead == 6:
        cleada = '06'
    elif ilead == 12:
        cleada = '12'
    elif ilead == 18:
        cleada = '18'
else:
    cleada = clead
    master_directory_panal_spline = '/Volumes/NBM/'+cdomain+'_panal/CDF_spline/'
    infile = master_directory_panal_spline + cmonth+'_'+cdomain+\
        '_CCPA_spline_info_h' + cleada + 'UTC.nc' 
    nc = Dataset(infile)
    spline_info_inv = nc.variables['spline_info_inv'][:,:,:,:]
    fraction_zero_ndfd = nc.variables['fzero'][:,:]
    usegamma_ndfd = nc.variables['usegamma'][:,:]
    lons_ndfd = nc.variables['lons'][:,:]
    lons_ndfd = lons_ndfd
    lats_ndfd = nc.variables['lats'][:,:]
    ny_ndfd, nx_ndfd = np.shape(lons_ndfd)
    nc.close()   

firstgroup = True
ngroups = len(groups)
print ('ngroups = ', ngroups)
qmapped_precip_mean  = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
precip_realtime_grand_mean = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
for group, groupname in zip(groups, groupnames):
    
    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print ('processing group = ',group,time)
    
    # ---- read in & extract the real-time precipitation for this domain.
    #      Set GEFSv12 grid dimensions and 1/4-deg lon and lat.

    precip_realtime_mean, lons_1d_realtime, lats_1d_realtime = \
        get_domain_subset_of_gefs(cyyyymmddhh, group, clead, jmin, \
        jmax, imin, imax) 
    nx_gefsv12 = len(lons_1d_realtime)
    ny_gefsv12 = len(lats_1d_realtime)
    lons_fcst_2d, lats_fcst_2d = \
        np.meshgrid(lons_1d_realtime,lats_1d_realtime)
        
    # ---- interpolate the GEFSv12 climatological fraction zeros to NDFD grid
    
    if firstgroup == True:
        fraction_zero_gefsv12_flipped = np.flipud(fraction_zero_gefsv12)
        lats_1d_realtime_flipud = np.flipud(lats_1d_realtime)
        fraction_zero_gefsv12_on_ndfd = interp(fraction_zero_gefsv12_flipped, \
            lons_1d_realtime, lats_1d_realtime_flipud, \
            lons_ndfd, lats_ndfd, checkbounds=False, \
            masked=False, order=1)
        firstgroup = False

    # ---- flip upside down, as subsequent interpolation requires
    #      S to N with increasing j index.  Then bilinear
    #      interpolate precipitation to NDFD grid. 

    precip_realtime_flipud = np.flipud(precip_realtime_mean) 
    precip_gefsv12_on_ndfd = interp(precip_realtime_flipud, \
        lons_1d_realtime, lats_1d_realtime_flipud,  \
        lons_ndfd, lats_ndfd, checkbounds=False, \
        masked=False, order=1)   
    
    # ---- now loop over grid points and obtain the forecast quantile

    gefsv12_quantiles = np.zeros((ny_gefsv12, nx_gefsv12), \
        dtype=np.float64)
    for jy in range(ny_gefsv12):
        for ix in range(nx_gefsv12):
            gefsv12_quantiles[jy,ix] = get_quantile_gefsv12(\
                precip_realtime_mean, jy, ix, spline_info_gefsv12, \
                fraction_zero_gefsv12, usegamma_gefsv12)

    # ---- interpolate the forecast quantile to the NDFD grid.

    gefsv12_quantiles_flipud = np.flipud(gefsv12_quantiles) 
    gefsv12_quantiles_on_ndfd = interp(gefsv12_quantiles_flipud, \
        lons_1d_realtime, lats_1d_realtime_flipud,  \
        lons_ndfd, lats_ndfd, checkbounds=False, \
        masked=False, order=1)  

    # ---- apply quantile mapping procedure to this group of members

    qmapped_precip = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
    qmapped_precip = qmapping_spline(gefsv12_quantiles_on_ndfd, \
    	precip_gefsv12_on_ndfd, spline_info_inv, \
    	fraction_zero_ndfd, usegamma_ndfd, ny_ndfd, nx_ndfd)
              
    # ---- write the raw mean and quantile mapped mean arrays to file
    
    # ---- form grand means of raw and qmapped
    
    qmapped_precip_mean  = qmapped_precip_mean + qmapped_precip
    precip_realtime_grand_mean = precip_realtime_grand_mean + precip_gefsv12_on_ndfd 
    
qmapped_precip_mean  = qmapped_precip_mean / float(ngroups)
precip_realtime_grand_mean = precip_realtime_grand_mean / float(ngroups)




# ---- store to file


    
# ======================================================================


# ---- plot the quantile-mapped forecast.

# 233 for whole US
m = Basemap(llcrnrlon=263.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ndfd, lats_ndfd)

clevs = [0.0, 0.25, 0.5, 0.75,1,1.5,2,3,4,5,6,8,10,15,20,25,30]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
fig = plt.figure(figsize=(5,8))
axloc = [0.02,0.1,0.96,0.8]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = '(b) Quantile-mapped mean forecast,\nIC = '+cyyyymmddhh+\
    ' lead = '+clead+' h'
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, qmapped_precip_mean, clevs,\
    cmap=None, colors=colorst, extend='both')
    
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.1,0.07,0.8,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Precipitation estimate (mm)',fontsize=9)

# ---- set plot title

plot_title = 'qmapped_precip_mean_'+cyyyymmddhh+'_lead'+clead+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



# ---- first plot GEFS control forecast amount

# 233 for whole US
m = Basemap(llcrnrlon=263.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_ndfd, lats_ndfd)

# ---- make plots of fraction positive precip

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

#clevs = [0.0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20]
fig = plt.figure(figsize=(5,8))
axloc = [0.02,0.1,0.96,0.8]
ax1 = fig.add_axes(axloc)
cleadb = str(int(clead)-6)
title = '(a) Ensemble-mean forecast precipitation\n amount (mm) for '+\
    clead+'-h IC = '+cyyyymmddhh
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, precip_realtime_grand_mean, clevs,\
    cmap=None, colors=colorst, extend='both')

m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
# ---- use axes_grid toolkit to make colorbar axes.

cax = fig.add_axes([0.1,0.07,0.8,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label('Mean precipitation (mm)',fontsize=9)

# ---- set plot title

plot_title = 'forecast_mean_precip_'+clead+'_h_IC'+cyyyymmddhh+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



