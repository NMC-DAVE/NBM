# python max_precip_byhuc_savelatlon_roebber.py

import os, sys
from datetime import datetime
from dateutils import daterange
import numpy as np
import _pickle as cPickle
from netCDF4 import Dataset
import scipy.stats as stats
import pygrib
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ---- process warm season separately from cool season.

master_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
cseason = sys.argv[1]
if cseason == 'cool' or cseason == 'cold':
    cmonths = ['Apr','May','Jun','Jul','Aug','Sep']
else:
    cmonths = ['Oct','Nov','Dec','Jan','Feb','Mar']

# --- read masks by huc, 1-18 for CONUS.  Generated 
#     from code create_mask_byhuc.py

infile = 'owp/HUC2_mask_roebber.cPick'
inf = open(infile, 'rb')
mask_HUC2 = cPickle.load(inf)
xgrid = cPickle.load(inf) #2D
ygrid = cPickle.load(inf) #2D
print ('max, min xgrid = ', np.max(xgrid), np.min(xgrid))
print ('ygrid[0], [-1] = ', ygrid[0,0], ygrid[-1,0])
print ('xgrid[0], [-1] = ', xgrid[0,0], xgrid[0,-1])
inf.close()
nhucs, ny, nx = np.shape(mask_HUC2)
print ('nhucs, ny, nx = ', nhucs, ny, nx)
lsmask_conus = np.zeros((ny,nx), dtype=np.int32)
for ihuc in range(nhucs):
    lsmask_conus[:,:] = lsmask_conus[:,:] + mask_HUC2[ihuc,:,:]
print ('max lsmask_conus = ', np.max(lsmask_conus))


# ---- loop over all HUCs, where the HUC number is ihuc+1. 
#      Here, when ihuc+1 = 19 this denotes the full conus domain.

for ihuc in range(19):
#for ihuc in range(10,11):

    if ihuc < 18:
        mask = mask_HUC2[ihuc,:,:]
    else:
        mask = lsmask_conus[:,:]
        
    print ('ihuc, sum(mask) = ',ihuc,np.sum(mask))
    huc_precip_mean = []
    huc_precip_max20 = []
    huc_precip_max20_meanlon = []
    huc_precip_max20_meanlat = []
    yyyymmddhh_init_thisseason = []
    
    for iyear in range(2000,2020):
    
        # --- read in the time series of GEFSv12 gridded CONUS precipitation 
        #     tallied over the full first 240h of the forecast. Previously
        #     generated on Hamill's computer by reforecast_2netcdf_240h.py
            
        cyear = str(iyear)
        ncfile = '/Volumes/NBM/python/roebber/ReforecastThinning_UWM/'+\
            'DATA/AccumulatedPrecip10_F_A/apcp.'+cyear+'.nc'
        date_list_fullyear = daterange(cyear+'010100', cyear+'123100',24)
        ndates_fullyear = len(date_list_fullyear)
        iyyyymmddhh_list = [int(i) for i in date_list_fullyear]
        if cseason == 'cool':
            cyyyymmddhh_begin1 = cyear+'010100'
            cyyyymmddhh_end1 = cyear+'033100'
            cyyyymmddhh_begin2 = cyear+'100100'
            cyyyymmddhh_end2 = cyear+'123100'
            date_list1 = daterange(cyyyymmddhh_begin1,cyyyymmddhh_end1, 24)
            date_list2 = daterange(cyyyymmddhh_begin2,cyyyymmddhh_end2, 24)
            date_list = date_list1 + date_list2
        else:
            cyyyymmddhh_begin = cyear+'070100'
            cyyyymmddhh_end = cyear+'093000'
            date_list = daterange(cyyyymmddhh_begin,cyyyymmddhh_end, 24)
        ndates = len(date_list)

        print (ncfile)
        nc = Dataset(ncfile)
        apcp_fcst = 1000.*nc.variables['apcp'][:,:,:]+0.001
        apcp_fcst 
        if ihuc == 0:
            lons_2d = nc.variables['lon'][:,:]
            lats_2d = nc.variables['lat'][:,:]
            ny, nx = np.shape(lons_2d)
            lons_2d = np.where(lons_2d > 0., lons_2d-360., lons_2d)
        nc.close()
        zeros = np.zeros((ndates_fullyear, ny, nx), dtype=np.float64)
        apcp_fcst = np.where(apcp_fcst != apcp_fcst, zeros, apcp_fcst)
        

        # --- for each date, determine the mean precipitation over this HUC
        #     (a simple average over all grid points, no cos(lat) weighting)
        #     as well as the precipitation at the subset of 20 grid points
        #     with the highest precipitation, so as to account for smaller-
        #     scale, very heavy events.
        
        for idate in range(ndates):
            pflag = False
            if pflag == True: print ('apcp_fcst[idate, ny//2,0:nx:4] = ', apcp_fcst[idate, ny//2,0:nx:4])
            if pflag == True: print ('mask[ny//2,0:nx:4] = ', mask[ny//2,0:nx:4])
            if pflag == True: print ('apcp_fcst[idate, 0:ny:4, nx//2] = ', apcp_fcst[idate, 0:ny:4, nx//2])
            if pflag == True: print ('mask[0:ny:4,ny//2] = ', mask[0:ny:4,ny//2])
            idate = iyyyymmddhh_list.index(int(date_list[idate]))

            pmean = np.sum(mask.astype(float)*apcp_fcst[idate,:,:]) / \
                np.sum(mask.astype(float))
            product = (mask.astype(float)*apcp_fcst[idate,:,:]).flatten()
            psorted_args = np.argsort(product)
            psorted = product[psorted_args]
            psorted_lons = lons_2d.flatten()[psorted_args]
            psorted_lats = lats_2d.flatten()[psorted_args]
            
            if pflag == True:
                print ('psorted = ', psorted)
            pmax_20 = np.mean(psorted[-20:])
            pmax_20_meanlon = np.mean(psorted_lons[-20:])
            pmax_20_meanlat = np.mean(psorted_lats[-20:])
            if pflag == True:
                print ('pmax_20 = ', pmax_20)
            huc_precip_mean.append(pmean)
            huc_precip_max20.append(pmax_20)
            huc_precip_max20_meanlon.append(pmax_20_meanlon)
            huc_precip_max20_meanlat.append(pmax_20_meanlat)
            yyyymmddhh_init_thisseason.append(iyyyymmddhh_list[idate])
            
            if pflag == True: print ('psorted[-20:] = ', psorted[-20:])
            if pflag == True: print ('ihuc, date, pmean, pmax20, \
                pmax20_meanlon, pmax20_meanlat = ',\
                ihuc, iyyyymmddhh_list[idate], pmean, \
                pmax_20, pmax_20_meanlon, pmax_20_meanlat)
            
            if pflag == True:
            
                # ---- plot this mask

                #m = Basemap(llcrnrlon=xgrid[0,0], llcrnrlat=ygrid[-1,0],\
                #    urcrnrlon = xgrid[-1,-1], urcrnrlat = ygrid[0,-1],\
                #    projection='mill',resolution ='l')
                    
                m = Basemap(llcrnrlon=-130.0, llcrnrlat=20.0,\
                    urcrnrlon = -60.0, urcrnrlat = 52.0,\
                    projection='mill',resolution ='l')
                x, y = m(lons_2d, lats_2d)
                
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
                clevs = [0,5,10,20,40,60,80,100,130,160,200,250,300,400]
                x, y = m(xgrid, ygrid)
                title = 'Ensemble-mean 0-10 day forecast precipitation, IC = '+\
                    str(iyyyymmddhh_list[idate])
                figtitle = 'precip_forecast_IC'+str(iyyyymmddhh_list[idate])+'.png'
                fig = plt.figure(figsize=(9,6.))
                axloc = [0.02,0.1,0.96,0.81]
                ax1 = fig.add_axes(axloc)
                ax1.set_title(title, fontsize=14,color='Black')
                print ('max, min conus apcp_fcst', \
                    np.max(apcp_fcst[idate,:,:]*lsmask_conus[:,:]), \
                    np.min(apcp_fcst[idate,:,:]*lsmask_conus[:,:]))
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
    
    outfile = cseason + '_precip_stats_roebber_huc2number'+str(ihuc+1)+'.cPick'
    print ('writing to ', outfile)
    ouf = open(outfile, 'wb')
    cPickle.dump(yyyymmddhh_init_thisseason_sorted, ouf)
    cPickle.dump(huc_precip_mean_sorted, ouf)
    cPickle.dump(huc_precip_max20_sorted, ouf)
    cPickle.dump(huc_precip_max20_meanlon_sorted, ouf)
    cPickle.dump(huc_precip_max20_meanlat_sorted, ouf)
    ouf.close()
    
    
    
    #np.savetxt(outfile,arrout, fmt=["%d","%6.2f","%7.2f"])