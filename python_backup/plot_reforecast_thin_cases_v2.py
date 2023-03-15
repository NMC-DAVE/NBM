
# python plot_reforecast_thin_cases_v2.py cseason ctotal_ncases cexclude_hucs

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

from matplotlib.backends.backend_pdf import PdfPages


cseason = sys.argv[1]   # cool, warm
ctotal_ncases = sys.argv[2] 
cexclude_hucs = sys.argv[3] # set to 1 if exclude 9, 13,16
if cexclude_hucs == '1':
    cexclude = '_no9_13_16'
else:
    cexclude = '_allhucs'

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
jmin = 93
jmax = 246
imin = 368
imax = 686
lsmask_conus = lsmask_nbm[jmin:jmax, imin:imax]
lats_conus = lats_nbm[jmin:jmax, imin:imax]
lons_conus = lons_nbm[jmin:jmax, imin:imax]
print ('max, min lons_conus = ', np.max(lons_conus), np.min(lons_conus))



# ---- read in the list of case dates.

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
infile = 'case_list_'+cseason+'season_ncases'+ctotal_ncases+cexclude+'.txt'
print (infile)
inf = open(infile,'r')
yyyymmddhh_save = []
hucnumber_save = []
lons_save = []
lats_save = []
for line in inf.readlines():
    yyyymmddhh_in, hucnumber_in, lon_in, lat_in = line.split()
    print (yyyymmddhh_in, hucnumber_in)
    yyyymmddhh_save.append(yyyymmddhh_in)
    hucnumber_save.append(hucnumber_in)
    lons_save.append(float(lon_in))
    lats_save.append(float(lat_in))


# ---- plot mean precip for each case date    
    
ncases = len(hucnumber_save)
print ('number of cases = ', ncases) 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
#clevs = [0.0,20,40,70,100,150,200,300,400,500,600,700,800,900,1200, 1500]
clevs = [0.0,5,10,15,20,30,40,60,80,100,120,140,170,200,300,400]
m = Basemap(llcrnrlon=lons_conus[0,0], llcrnrlat=lats_conus[-1,0],\
    urcrnrlon = lons_conus[-1,-1], urcrnrlat = lats_conus[0,-1],\
    projection='mill',resolution ='l')
x, y = m(lons_conus, lats_conus)   

print ('yyyymmddhh_save[0:30] = ', yyyymmddhh_save[0:30] )
outfile = 'multipage_reforecast_samples_'+cseason+cexclude+'.pdf'
with PdfPages(outfile) as pdf:    

    for icase in range(0,ncases-1,2):
    #for icase in range(0,30,2):
    #for icase in range(0,4,2):
    
        print ('processing ',icase,icase+1)
        fig1 = plt.figure(figsize=(6.,9))

        # ---- process first case to plot on page

        yyyymmddhh_in = yyyymmddhh_save[icase] 
        hucnumber_in = int(hucnumber_save[icase])
        lonmax = lons_save[icase]
        latmax = lats_save[icase]
        print ('icase, latmax = ',icase,latmax)
        chuc = str(hucnumber_in)
        if hucnumber_in == 19:
            ctitle = 'IC = '+str(yyyymmddhh_in)+' optimized for whole CONUS'
            mask = 0*lsmask_conus
        else:
            if latmax > 0.0:
                mean_or_max =  ' max20 '
            else:
                mean_or_max =  ' mean '
            ctitle = 'IC = '+str(yyyymmddhh_in)+' optimized for HUC '+chuc+mean_or_max
            mask = mask_HUC2[hucnumber_in-1,:,:]
        cmm = yyyymmddhh_in[4:6]
        imm = int(cmm)-1    
    
        # ---- open the appropriate file and read the forecast for this case day
    
        master_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
        cmonth = cmonths[imm]
        ncfile = master_directory + cmonth+ \
            '_conus_reforecast_ens_mean_0_to_240h.nc'
        nc = Dataset(ncfile)
        yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
        idx = np.where(yyyymmddhh_init == int(yyyymmddhh_in))[0]
        apcp_fcst = np.squeeze(nc.variables['apcp_fcst'][idx,:,:])
        lons_1d = nc.variables['lons_fcst'][:]
        lats_1d = nc.variables['lats_fcst'][:]
        ny = len(lats_1d)
        nx = len(lons_1d)
        nc.close()


        ax = fig1.add_axes([0.02,.54,0.96,.4])

        ax.set_title(ctitle,fontsize=12)
        CS2 = m.contourf(x, y, apcp_fcst[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')
        for jy in range(ny):
            for ix in range(nx):
                if mask[jy,ix] > 0:
                    #print (jy,ix)
                    xdot, ydot = m(lons_conus[jy,ix], lats_conus[jy,ix])
                    m.plot(xdot,ydot,marker='.',markersize=0.25,\
                        color='Black',markerfacecolor='Black')
                        
        # --- if this case was selected by max20, plot a big black dot
        #     at the centroid of max precip       
        if latmax > 0.0:
            print ('Max20:  ',lonmax, latmax)
            xdot, ydot = m(lonmax+360., latmax)
            m.plot(xdot,ydot,'ko',markersize=3)
                
        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
        

        # ---- move on to the second day to plot on this page. 
        #      open the appropriate file and read the forecast for this case day
        
        # ---- process second case to plot on page

        yyyymmddhh_in = yyyymmddhh_save[icase+1] 
        hucnumber_in = int(hucnumber_save[icase+1])
        lonmax = lons_save[icase+1]
        latmax = lats_save[icase+1]
        print ('icase+1, latmax = ',icase+1,latmax)
        chuc = str(hucnumber_in)
        if hucnumber_in == 19:
            ctitle = 'IC = '+str(yyyymmddhh_in)+' optimized for whole CONUS'
            mask = 0*lsmask_conus
        else:
            if latmax > 0.0:
                mean_or_max =  ' max20 '
            else:
                mean_or_max =  ' mean '
            ctitle = 'IC = '+str(yyyymmddhh_in)+' optimized for HUC '+chuc+mean_or_max
            mask = mask_HUC2[hucnumber_in-1,:,:]
        cmm = yyyymmddhh_in[4:6]
        imm = int(cmm)-1    
    
        # ---- open the appropriate file and read the forecast for this case day
    
        master_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
        cmonth = cmonths[imm]
        ncfile = master_directory + cmonth+ \
            '_conus_reforecast_ens_mean_0_to_240h.nc'
        nc = Dataset(ncfile)
        yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
        idx = np.where(yyyymmddhh_init == int(yyyymmddhh_in))[0]
        apcp_fcst = np.squeeze(nc.variables['apcp_fcst'][idx,:,:])
        lons_1d = nc.variables['lons_fcst'][:]
        lats_1d = nc.variables['lats_fcst'][:]
        ny = len(lats_1d)
        nx = len(lons_1d)
        nc.close()


        ax = fig1.add_axes([0.02,.13,0.96,.4])
        ax.set_title(ctitle,fontsize=12)
        CS2 = m.contourf(x, y, apcp_fcst[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')
        for jy in range(ny):
            for ix in range(nx):
                if mask[jy,ix] > 0:
                    #print (jy,ix)
                    xdot, ydot = m(lons_conus[jy,ix], lats_conus[jy,ix])
                    m.plot(xdot,ydot,marker='.',markersize=0.25,\
                        color='Black',markerfacecolor='Black')
                        
        # --- if this case was selected by max20, plot a big black dot
        #     at the centroid of max precip       
        if latmax > 0.0:
            print ('Max20:  ',lonmax, latmax)
            xdot, ydot = m(lonmax+360., latmax)
            m.plot(xdot,ydot,'ko',markersize=3)
            
        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
        
        
        # ---- use axes_grid toolkit to make colorbar axes.

        cax = fig1.add_axes([0.02,0.07,0.96,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Ensemble-mean 0-10 day total precipitation (mm)',fontsize=9)
            
        
    
        pdf.savefig()
        plt.close()
        
print ('saved to ', outfile)