"""
qmap_gefsv12_ensemble_ccpa_spline_nc.py cyyyymmddhh clead 

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
        based on the Cumulative Hazard function from 
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
    spline_info_gefsv12, fraction_zero_gefsv12,\
    usegamma_gefsv12, quantile_98, use98):

    """ this gets the quantile associated with a given precipitation 
    amount for GEFSv12 data, this month and 6-hour period. """
    
    offset_out = 0.0
    if precip_amount[jy,ix] == 0.0:
        
        # ---- arbitrarily assign the CDF to zero if precip is zero.
        
	    quantile = 0.0
    else:	
	
        if usegamma_gefsv12[jy,ix] == 0:
            
	        # ---- flagged as a wet-enough point to estimate the CDF with 
            #      the spline fit to a hazard function. 
            
            splines_tuple = (spline_info_gefsv12[jy,ix,0,:], \
                spline_info_gefsv12[jy,ix,1,:], 3)
            spline_hazard = splev(precip_amount[jy,ix], splines_tuple)
            spline_cdf = 1.0 - np.exp(-spline_hazard)
            quantile = fraction_zero_gefsv12[jy,ix] + \
                (1.0 - fraction_zero_gefsv12[jy,ix])*spline_cdf
        else:
            
            if usegamma_gefsv12[jy,ix] == -1:
                # --- flagged as basically no training data.
                quantile = 0.0
            else:  # --- flagged as minimal training data - use Gamma
                alpha_hat = spline_info_fcst[jy,ix,0,0] 
                beta_hat = spline_info_fcst[jy,ix,1,0] 
                y0 = precip_amount[jy,ix] / beta_hat
                quantile = stats.gamma.cdf(y0, alpha_hat)   
                
        if use98 == True and quantile > 0.98:
            offset_out = precip_amount[jy,ix] - quantile_98[jy,ix]
            if quantile > 0.98: quantile = 0.98         

    return quantile, offset_out
    
# =====================================================================
    
def get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead):

    # ---- read in the global forecast.  Subset to CONUS.

    nib = 886 
    nie = nib+318 
    njb = 131 
    nje = njb+153 
    
    # ---- get the whole globe.
    
    input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
    infile = input_directory + cyyyymmddhh + \
        '_ge'+cmem+'.t00z.pgrb2s.0p25.f' + clead 
    print ('  reading from ', infile)
    endStep = int(clead)
    istat, precip_realtime, lats_full, lons_full = \
        read_gribdata(infile, endStep)
    
    # ---- subset for the big grid encompassing all NBM domains.
    #      The eastern boundary crosses Greenwich meridian, so
    #      fiddle to make sure longitudes in ascending order.
    
    precip_realtime = precip_realtime[njb:nje,nib:nie]
    lons_1D_realtime = lons_full[0,nib:nie]-360.
    lats_1D_realtime = lats_full[njb:nje,0]
    #print ('njb,nje, nib,nie =  ', njb,nje, nib,nie)
    #print ('min, max lons_1D_realtime = ',\
    #    np.min(lons_1D_realtime), np.max(lons_1D_realtime))
    #print ('min, max lats_1D_realtime = ',\
    #    np.min(lats_1D_realtime), np.max(lats_1D_realtime))
    return precip_realtime, lons_1D_realtime, lats_1D_realtime    
    
                
# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
use98 = True
cmembers = ['c00','p01', 'p02','p03','p04','p05','p06','p07','p08','p09','p10',\
    'p11', 'p12','p13','p14','p15','p16','p17','p18','p19','p20',\
    'p21', 'p22','p23','p24','p25','p26','p27','p28','p29','p30']
nmembers = len(cmembers)
    
master_directory_panal_spline = '/Volumes/NBM/conus_panal/CDF_spline/'
master_directory_forecast_spline = '/Volumes/NBM/conus_gefsv12/CDF_spline/'
   
# ---- load spline information for forecast netCDF file.  
#      Should be same lat/lon as real-time, but check

infile = master_directory_forecast_spline + \
    ccmonth + '_' + cdomain + \
    '_GEFSv12_spline_info_h' + clead + '.nc' 
nc = Dataset(infile)    
lons_spline_gefsv12_1d = nc.variables['lons'][:]
lats_spline_gefsv12_1d = nc.variables['lats'][:]
spline_info_gefsv12 = nc.variables['spline_info'][:,:,:,:]
fraction_zero_gefsv12 = nc.variables['fzero'][:,:]
usegamma_gefsv12 = nc.variables['usegamma'][:,:]
quantile_98 = nc.variables['quantile_98'][:,:]
ny_gefsv12, nx_gefsv12 = np.shape(usegamma_gefsv12)
nc.close()

# --- read netCDF file here for inverse spline parameters for
#     the combined CCPA/MSWEP precip analysis CDFs.   
#     **** Note that if we are applying to 
#     cycles other than 00UTC init, we'll need to substitute
#     clead for the appropriate other period time.   
   
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

# ---- set up output grids

precip_ens_raw_ndfd = np.zeros((nmembers, ny_ndfd, nx_ndfd), dtype=np.float32)
precip_ens_qmapped_ndfd = np.zeros((nmembers, ny_ndfd, nx_ndfd), dtype=np.float32)
prob_POP = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float32)
prob_1mm = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float32)
prob_5mm = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float32)
prob_10mm = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float32)
prob_25mm = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float32)
ny_ndfd, nx_ndfd = np.shape(lons_ndfd)

for imem, cmem in enumerate(cmembers): 
    
    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print ('processing member = ',cmem,time)
    
    # ---- read in & extract the real-time precipitation for this domain.
    #      Set GEFSv12 grid dimensions and 1/4-deg lon and lat.

    precip_realtime, lons_1d_realtime, lats_1d_realtime = \
        get_domain_subset_of_gefs(cyyyymmddhh, cmem, clead) 
    if cmem == 'c00':
        nx_gefsv12 = len(lons_1d_realtime)
        ny_gefsv12 = len(lats_1d_realtime)
    lons_fcst_2d, lats_fcst_2d = \
        np.meshgrid(lons_1d_realtime,lats_1d_realtime)
        
    # ---- interpolate the GEFSv12 climatological fraction zeros 
    #      and latitude array to NDFD grid. Set up output grids.
    
    if cmem == 'c00':
        fraction_zero_gefsv12_flipped = np.flipud(fraction_zero_gefsv12)
        lats_1d_realtime_flipud = np.flipud(lats_1d_realtime)
        fraction_zero_gefsv12_on_ndfd = interp(fraction_zero_gefsv12_flipped, \
            lons_1d_realtime, lats_1d_realtime_flipud, \
            lons_ndfd, lats_ndfd, checkbounds=False, \
            masked=False, order=1)

    # ---- flip upside down, as subsequent interpolation requires
    #      S to N with increasing j index.  Then
    #      interpolate precipitation to NDFD grid. 

    print ('   interpolating and flipping real-time forecast')
    precip_realtime_flipud = np.flipud(precip_realtime) 
    precip_gefsv12_on_ndfd = interp(precip_realtime_flipud, \
        lons_1d_realtime, lats_1d_realtime_flipud,  \
        lons_ndfd, lats_ndfd, checkbounds=False, \
        masked=False, order=1)  
    precip_ens_raw_ndfd[imem,:,:] = precip_gefsv12_on_ndfd[:,:]
        
    # ---- now loop over grid points and obtain the forecast quantile
    #      associated with this 0.25 degree forecast grid point

    print ('   getting quantiles of forecast')
    gefsv12_quantiles = np.zeros((ny_gefsv12, nx_gefsv12), \
        dtype=np.float64)
    offset = np.zeros((ny_gefsv12, nx_gefsv12), \
        dtype=np.float64)
    for jy in range(ny_gefsv12):
        for ix in range(nx_gefsv12):
            gefsv12_quantiles[jy,ix], offset[jy,ix] = get_quantile_gefsv12( \
                precip_realtime, jy, ix, spline_info_gefsv12, \
                fraction_zero_gefsv12, usegamma_gefsv12, quantile_98, use98)

    # ---- interpolate the forecast quantile to the NDFD grid.

    print ('   interpolating quantiles to NDFD')
    gefsv12_quantiles_flipud = np.flipud(gefsv12_quantiles) 
    gefsv12_quantiles_on_ndfd = interp(gefsv12_quantiles_flipud, \
        lons_1d_realtime, lats_1d_realtime_flipud,  \
        lons_ndfd, lats_ndfd, checkbounds=False, \
        masked=False, order=1)  
    offset_flipud = np.flipud(offset)    
    offset_on_ndfd = interp(offset_flipud, \
        lons_1d_realtime, lats_1d_realtime_flipud,  \
        lons_ndfd, lats_ndfd, checkbounds=False, \
        masked=False, order=1)

    # ---- apply quantile mapping procedure to this member

    print ('   applying quantile mapping')
    qmapped_precip = np.zeros((ny_ndfd, nx_ndfd), dtype=np.float64)
    qmapped_precip = qmapping_spline(gefsv12_quantiles_on_ndfd, \
    	precip_gefsv12_on_ndfd, spline_info_inv, \
    	fraction_zero_ndfd, usegamma_ndfd, use98, offset_on_ndfd, \
        ny_ndfd, nx_ndfd)
    precip_ens_qmapped_ndfd[imem,:,:] = qmapped_precip[:,:]
              
    # ---- write the quantile mapped array to file
    
    input_directory = '/Volumes/Backup Plus/gefsv12/2021/'
    if int(clead) > 100:
        outfile = input_directory + cyyyymmddhh + \
            '_qmap_'+cmem+'.t00z.pgrb2s.0p25.f' + clead +'.cPick'
    else:
        outfile = input_directory + cyyyymmddhh + \
            '_qmap_'+cmem+'.t00z.pgrb2s.0p25.f0' + clead+'.cPick'
    print ('   writing to ', outfile)
    ouf = open(outfile,'wb')
    cPickle.dump(precip_gefsv12_on_ndfd, ouf)
    cPickle.dump(qmapped_precip, ouf)
    ouf.close()
    
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



