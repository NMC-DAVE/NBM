
"""

plot_thinned_forecast_anal_1mm.py cyyyymm

"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.stats as stats
import pygrib
from netCDF4 import Dataset
from dateutils import hrs_since_day1CE_todate, \
    dateto_hrs_since_day1CE, hrstodate, datetohrs, dateshift
from mpl_toolkits.basemap import Basemap, interp

# ---- get the month and end time from the commmand line

cyyyymm = sys.argv[1] # 01 etc
clead = sys.argv[2] # 018 etc, 3 digit
cyyyy = cyyyymm[0:4]
cmm = cyyyymm[4:6]
imonth = int(cmm)-1
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
cmonths = ['Jan','Feb','Mar','Apr','May',\
    'Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
if int(cyyyy)%4 == 0:
    ndays = ndaysomo_leap[imonth]
else:
    ndays = ndaysomo[imonth]
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0]
nthresh = len(thresholds)
rthresh = 1.0
cthresh = str(rthresh)

# ---- loop over forecast samples

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
for iday in range(ndays):
    if iday+1 < 10:
        cday = '0'+str(iday+1)
    else:
        cday = str(iday+1)
        
    cyyyymmddhh_init = cyyyymm + cday + '00'
    iyyyymmddhh_init = int(cyyyymmddhh_init) 
    print (cyyyymmddhh_init, int(clead))  
    cyyyymmddhh_end = dateshift(cyyyymmddhh_init, int(clead))
    print ('IC, FC = ', cyyyymmddhh_init, cyyyymmddhh_end)
    iyyyymmddhh_end = int(cyyyymmddhh_end)
        
    # ---- read the stored netCDF thinned-field probabilities, raw
    #      and quantile mapped

    master_directory_thinned_output = '/Volumes/NBM/conus_gefsv12/thinned/'
    infile = master_directory_thinned_output+\
        ccmonth+cyyyy+'_lead'+clead+'_probabilities_thinned.nc'
    print ('reading from ', infile)
    nc = Dataset(infile,'r')
    if iday == 0:
        lons_thinned = nc.variables['lons'][:,:]
        lats_thinned = nc.variables['lats'][:,:]
        iyyyymmddhh_ic = nc.variables['yyyymmddhh_init'][:]
        ny_thinned, nx_thinned = np.shape(lons_thinned)
        print ('lons_thinned[ny_thinned//2,0:nx_thinned//2:5] = ', \
            lons_thinned[ny_thinned//2,0:nx_thinned//2:5])
        m = Basemap(llcrnrlon=233.7234-360.,llcrnrlat=19.229,
            urcrnrlon = 300.95782-360., urcrnrlat = 54.37279,\
            projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
            resolution ='l',area_thresh=1000.)
        x, y = m(lons_thinned, lats_thinned)
    
    idd = int(np.where(iyyyymmddhh_ic == int(iyyyymmddhh_init))[0])
    itt = 1 # int(np.where(thresholds == rthresh)[0])
    prob_raw_thinned = nc.variables['probability_raw'][idd,itt,:,:]
    prob_qmapped_thinned = nc.variables['probability_qmapped'][idd,itt,:,:]
    nc.close()

    # ---- read the time-coincident precipitation analyses 
    
    cyyyymm_anal = cyyyymmddhh_end[0:6]
    master_directory = '/Volumes/NBM/conus_panal/'
    infile = master_directory + cyyyymm_anal + \
        '_ccpa_on_ndfd_grid_6hourly_thinned.nc'
    nc = Dataset(infile)
    yyyymmddhh_end_in = nc.variables['yyyymmddhh_end'][:]
    idx = np.where(yyyymmddhh_end_in == int(cyyyymmddhh_end))[0]
    conusmask_in = nc.variables['conusmask'][:,:]
    precip_anal = np.squeeze(nc.variables['apcp_anal'][idx,:,:])
    ny_anal, nx_anal = np.shape(precip_anal)
    if iday == 0:
        ones = np.ones((ny_anal, nx_anal), dtype=np.float64)
        zeros = np.zeros((ny_anal, nx_anal), dtype=np.float64)
        lons_anal = nc.variables['lons'][:,:]
        lats_anal = nc.variables['lats'][:,:]
        print ('lons_anal[ny_anal//2,0:nx_anal//2:5] = ', \
            lons_thinned[ny_anal//2,0:nx_anal//2:5])

    weight = np.where(conusmask_in == 1.0, \
        ones*np.cos(lats_anal*3.1415926/180.), zeros)
    nc.close()    
        
        
    # ---- plot the raw forecast probability
    
    plot_raw = False
    if plot_raw == True:
        clevs = [0.0,0.02,0.04,0.07,0.1,0.2,0.3,0.4,\
            0.5,0.6,0.7,0.8,0.9,1.0]
        fig = plt.figure(figsize=(9,6.5))
        axloc = [0.06,0.15,0.92,0.8]
        ax1 = fig.add_axes(axloc)
        cleadb = str(int(clead)-6)
        title = 'Thinned raw GEFSv12 P(obs > '+\
            cthresh+' mm), IC = '+cyyyymmddhh_init+', lead = '+clead+' h'
        ax1.set_title(title, fontsize=14,color='Black')
        CS2 = m.contourf(x, y, prob_raw_thinned[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')

        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')

        cax = fig.add_axes([0.06,0.07,0.88,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Probability',fontsize=9)

        plot_title = 'thinned_raw_probability_'+cthresh+'mm_'+\
            cyyyymmddhh_init+'_lead'+clead+'.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)    
    
    # ---- plot the quantile-mapped forecast probability
    
    plot_qmapped = False
    if plot_qmapped == True:
        clevs = [0.0,0.02,0.04,0.07,0.1,0.2,0.3,0.4,\
            0.5,0.6,0.7,0.8,0.9,1.0]
        fig = plt.figure(figsize=(9,6.5))
        axloc = [0.06,0.15,0.92,0.8]
        ax1 = fig.add_axes(axloc)
        cleadb = str(int(clead)-6)
        title = 'Thinned quantile-mapped GEFSv12 P(obs > '+\
            cthresh+' mm), IC = '+cyyyymmddhh_init+', lead = '+clead+' h'
        ax1.set_title(title, fontsize=14,color='Black')
        CS2 = m.contourf(x, y, prob_qmapped_thinned[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')

        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')

        cax = fig.add_axes([0.06,0.07,0.88,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Probability',fontsize=9)

        plot_title = 'thinned_qmapped_probability_'+cthresh+'mm_'+\
            cyyyymmddhh_init+'_lead'+clead+'.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)
    
    # ---- plot time-coincident precipitation analysis
 

    clevs = [0.0,0.1,0.3,0.6,1,2,3,5,7,10,15,20,25]
    fig = plt.figure(figsize=(9,6.5))
    axloc = [0.02,0.1,0.96,0.84]
    ax1 = fig.add_axes(axloc)

    CS2 = m.contourf(x, y, precip_anal, clevs,\
        cmap=None, colors=colorst, extend='both')

    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    
    countabove = 0
    countbelow = 0
    for jy in range(0,ny_anal):
        for ix in range(0,nx_anal):
            xdot, ydot = m(lons_thinned[jy,ix],lats_thinned[jy,ix])
            prob_err = np.abs(prob_qmapped_thinned[jy,ix]-1./31.)
            if prob_err < 0.01 and conusmask_in[jy,ix] == 1:
                if precip_anal[jy,ix] > 1.0:
                    m.plot(xdot,ydot,marker='.',markersize=0.6,\
                        color='Orange',markerfacecolor='Orange')
                    countabove = countabove + 1
                else:
                    m.plot(xdot,ydot,marker='.',markersize=0.6,\
                        color='Green',markerfacecolor='Green')
                    countbelow = countbelow + 1
    
    relia = countabove / (countabove+countbelow)
    crelia = str(relia)
    title = 'Thinned CCPA precipitation, 6 h ending '+\
        cyyyymmddhh_end+', 1/31 relia = '+crelia
    ax1.set_title(title, fontsize=14,color='Black')
    
    # ---- use axes_grid toolkit to make colorbar axes.

    cax = fig.add_axes([0.06,0.07,0.88,0.02])
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=5)
    cb.set_label('Precipitation (mm)',fontsize=9)

    # ---- set plot title

    plot_title = 'CCPA_thinned_'+cyyyymmddhh_end+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    print ('Done!')  
    plt.close() 
  
