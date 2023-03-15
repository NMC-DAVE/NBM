# python max_precip_byhuc.py

import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
from netCDF4 import Dataset
import scipy.stats as stats
import pygrib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

# --- read masks by huc, 1-18 for CONUS.  Generated 
#     from code create_mask_byhuc.py

infile = 'owp/HUC2_mask.cPick'
inf = open(infile, 'rb')
mask_HUC2 = cPickle.load(inf)
xgrid = cPickle.load(inf) #2D
ygrid = cPickle.load(inf) #2D
print ('max, min xgrid = ', np.max(xgrid), np.min(xgrid))
print ('ygrid[0], [-1] = ', ygrid[0,0], ygrid[-1,0])
inf.close()

# ---- read in global mask for 0.25 degree grid.  
#      subset to large NBM domain, then to CONUS.
#      Global mask obtained from https://
#      noaa-gefs-retrospective.s3.amazonaws.com/index.html

infile = 'landsfc.pgrb2.0p25'
lsmaskfile = pygrib.open(infile)
grb = lsmaskfile.select()[0]
lsmask_global = grb.values
lats_global, lons_global = grb.latlons()
print ('shape lats_global, lons_global = ', \
    np.shape(lats_global), np.shape(lons_global))
lsmaskfile.close()

# ---- following other processing Hamill does for 
#      PPGC-funded project to integrate GEFSv12 to NBM,
#      first carve out a domain that includes AK, HI,
#      Guam, PR, CONUS.

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
    
lsmask_nbm = lsmask_global[njb:nje,nib1:nie1]
lons_nbm = lons_global[njb:nje, nib1:nie1]
lats_nbm = lats_global[njb:nje, nib1:nie1]

# ---- next carve out the further subset surrounding the
#      CONUS NBM domain

jmin = 93
jmax = 246
imin = 368
imax = 686
lsmask_conus = lsmask_nbm[jmin:jmax, imin:imax]
lats_conus = lats_nbm[jmin:jmax, imin:imax]
lons_conus = lons_nbm[jmin:jmax, imin:imax]

# ---- process warm season separately from cool season.

master_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
cseason = sys.argv[1]
if cseason == 'warm':
    cmonths = ['Apr','May','Jun','Jul','Aug','Sep']
else:
    cmonths = ['Oct','Nov','Dec','Jan','Feb','Mar']
    
# ---- loop over all HUCs, where the HUC number is ihuc+1. 
#      Here, when ihuc+1 = 19 this denotes the full conus domain.

for ihuc in range(19):
#for ihuc in range(1,2):

    if ihuc < 18:
        mask = mask_HUC2[ihuc,:,:]
    else:
        mask = lsmask_conus[:,:]
        
    huc_precip_mean = []
    huc_precip_max20 = []
    huc_precip_max20_meanlon = []
    huc_precip_max20_meanlat = []
    yyyymmddhh_init_thisseason = []
    for cmonth in cmonths:
    
        # --- read in the time series of GEFSv12 gridded CONUS precipitation 
        #     tallied over the full first 240h of the forecast. Previously
        #     generated on Hamill's computer by reforecast_2netcdf_240h.py
        
        ncfile = master_directory + cmonth+ \
            '_conus_reforecast_ens_mean_0_to_240h.nc'
        print (ncfile)
        nc = Dataset(ncfile)
        apcp_fcst = nc.variables['apcp_fcst'][:,:,:]
        yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
        ndates = len(yyyymmddhh_init)
        if ihuc == 0:
            lons_1d = nc.variables['lons_fcst'][:]
            lats_1d = nc.variables['lats_fcst'][:]
            ny = len(lats_1d)
            nx = len(lons_1d)
        nc.close()

        # --- for each date, determine the mean precipitation over this HUC
        #     (a simple average over all grid points, no cos(lat) weighting)
        #     as well as the precipitation at the subset of 20 grid points
        #     with the highest precipitation, so as to account for smaller-
        #     scale, very heavy events.
        
        for idate in range(ndates):
            pflag = False
            #if yyyymmddhh_init[idate] == 2006083100: pflag = True
            pmean = np.sum(mask.astype(float)*apcp_fcst[idate,:,:]) / \
                np.sum(mask.astype(float))
            product = (mask.astype(float)*apcp_fcst[idate,:,:]).flatten()
            if pflag == True:
                print ('Product = ', product)
            a = np.where(product > 0.0)
            psorted = np.sort(product[a])
            if pflag == True:
                print ('psorted = ', psorted)
            pmax_20 = np.mean(psorted[-20:])
            if pflag == True:
                print ('pmax_20 = ', pmax_20)
            huc_precip_mean.append(pmean)
            huc_precip_max20.append(pmax_20)
            yyyymmddhh_init_thisseason.append(yyyymmddhh_init[idate])
            
            if pflag == True:
            
                # ---- plot this mask

                m = Basemap(llcrnrlon=xgrid[0,0], llcrnrlat=ygrid[-1,0],\
                    urcrnrlon = xgrid[-1,-1], urcrnrlat = ygrid[0,-1],\
                    projection='mill',resolution ='l')
                title = 'Mask for HUC '+str(ihuc+1)
                figtitle = 'HUC_mask_situational.png'
                fig = plt.figure(figsize=(9,6.))
                axloc = [0.02,0.1,0.96,0.81]
                ax1 = fig.add_axes(axloc)
                ax1.set_title(title, fontsize=14,color='Black')
                for jy in range(ny):
                    for ix in range(nx):
                        if mask[jy,ix] > 0:
                            xdot, ydot = m(xgrid[jy,ix], ygrid[jy,ix])
                            m.plot(xdot,ydot,marker='.',markersize=0.2,\
                                color='Black',markerfacecolor='Black')
                m.drawcoastlines(linewidth=0.8,color='Gray')
                m.drawcountries(linewidth=0.8,color='Gray')
                m.drawstates(linewidth=0.8,color='Gray')
    
                fig.savefig(figtitle, dpi=300)
                print ('saving plot to file = ',figtitle)
                print ('Done!')
                plt.close()
            
                # ---- plot the mean precipitation for this case

                colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
                    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
                    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
                clevs = [0.0,20,40,70,100,150,200,300,400,500,600,700,800,900,1200, 1500]
                x, y = m(xgrid, ygrid)
                title = 'Ensemble-mean 0-10 day forecast precipitation, IC = '+\
                    str(yyyymmddhh_init[idate])
                figtitle = 'precip_forecast_IC'+str(yyyymmddhh_init[idate])+'.png'
                fig = plt.figure(figsize=(9,6.))
                axloc = [0.02,0.1,0.96,0.81]
                ax1 = fig.add_axes(axloc)
                ax1.set_title(title, fontsize=14,color='Black')
                CS2 = m.contourf(x, y, apcp_fcst[idate,:,:], clevs,\
                    cmap=None, colors=colorst, extend='both')
                m.drawcoastlines(linewidth=0.8,color='Gray')
                m.drawcountries(linewidth=0.8,color='Gray')
                m.drawstates(linewidth=0.8,color='Gray')
    
                # ---- use axes_grid toolkit to make colorbar axes.

                cax = fig.add_axes([0.02,0.07,0.96,0.02])
                cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
                    drawedges=True,ticks=clevs,format='%g')
                cb.ax.tick_params(labelsize=7)
                cb.set_label('Ensemble mean 0-10 day total precipitation (mm)',fontsize=9)
    
                fig.savefig(figtitle, dpi=300)
                print ('saving plot to file = ',figtitle)
                print ('Done!')
                plt.close()
    
    
    # --- reorder data in terms of ascending date
    
    indices = np.argsort(np.array(yyyymmddhh_init_thisseason))
    yyyymmddhh_init_thisseason_sorted = np.array(yyyymmddhh_init_thisseason, dtype=np.int32)[indices]
    huc_precip_mean_sorted = np.array(huc_precip_mean)[indices]
    huc_precip_max20_sorted = np.array(huc_precip_max20)[indices]
    huc_precip_max20_meanlon_sorted = np.array(huc_precip_max20_meanlon)[indices]
    huc_precip_max20_meanlat_sorted = np.array(huc_precip_max20_meanlat)[indices]
    
    # --- write the resulting lists for this HUC-2 to ascii file 
    #     for later use.
    
    outfile = cseason + '_precip_stats_huc2number'+str(ihuc+1)+'.cPick'
    print ('writing to ', outfile)
    ouf = open(outfile, 'wb')
    cPickle.dump(yyyymmddhh_init_thisseason_sorted, ouf)
    cPickle.dump(huc_precip_mean_sorted, ouf)
    cPickle.dump(huc_precip_max20_sorted, ouf)
    ouf.close()
    #np.savetxt(outfile,arrout, fmt=["%d","%6.2f","%7.2f"])