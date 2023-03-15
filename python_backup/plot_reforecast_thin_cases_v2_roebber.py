
# python plot_reforecast_thin_cases_v2_roebber.py cseason clon clat

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

from matplotlib.backends.backend_pdf import PdfPages


cseason = sys.argv[1]   # cool, warm
 
data_directory = '/Volumes/NBM/python/roebber/ReforecastThinning_UWM/DATA/'
data_directory_fcst = '/Volumes/NBM/python/roebber/ReforecastThinning_UWM/DATA/AccumulatedPrecip10_F_A/'

# ---- read in the list of case dates.

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
infile = data_directory + cseason + 'season520.txt'
print (infile)
inf = open(infile,'r')
yyyymmddhh_save = []
hucnumber_save = []
for line in inf.readlines():
    yyyymmddhh_in, hucnumber_in = line.split()
    print (yyyymmddhh_in, hucnumber_in)
    yyyymmddhh_save.append(yyyymmddhh_in)
    hucnumber_save.append(int(hucnumber_in))
    
# --- read masks by huc, 1-18 for CONUS.  Generated
#     from code create_mask_byhuc.py

infile = 'owp/HUC2_mask.cPick'
inf = open(infile, 'rb')
mask_HUC2 = cPickle.load(inf)
xgrid = cPickle.load(inf) #2D
ygrid = cPickle.load(inf) #2D
nhuc,nyconus,nxconus = np.shape(mask_HUC2)
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

# ---- plot mean precip for each case date    
    
ncases = len(hucnumber_save)
print ('number of cases = ', ncases) 
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
clevs = [0,5,10,20,40,60,80,100,130,160,200,250,300,400]
 

print ('yyyymmddhh_save[0:30] = ', yyyymmddhh_save[0:30] )
outfile = 'multipage_roebber_kratsov_samples_'+cseason+'.pdf'
with PdfPages(outfile) as pdf:    

    for icase in range(0,ncases-1,2):
    #for icase in range(0,8,2):

        print ('processing ',icase,icase+1)
        fig1 = plt.figure(figsize=(6.,8.))

        # ---- process first case to plot on page
    
        yyyymmddhh_in = int(yyyymmddhh_save[icase])
        cyyyy = str(yyyymmddhh_save[icase])[0:4]
        cyyyymmddhh_start = cyyyy + '010100'
        cyyyymmddhh_end = cyyyy + '123100'
        cyyyymmddhh_list = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)
        iyyyymmddhh_list = [int(i) for i in cyyyymmddhh_list]
        
        infile = data_directory_fcst + 'apcp.'+cyyyy+'.nc'
        nc = Dataset(infile)    
        idx = iyyyymmddhh_list.index(yyyymmddhh_in)
        apcp_fcst1 = 1000.*np.squeeze(nc.variables['apcp'][idx,:,:])
        lons = nc.variables['lon'][:,:]
        lats = nc.variables['lat'][:,:]
        lons = np.where(lons > 0., lons-360., lons)
        ny, nx = np.shape(lons) 
        nc.close()
        
        if icase == 0:
            m = Basemap(llcrnrlon=-130.0, llcrnrlat=20.0,\
                urcrnrlon = -60.0, urcrnrlat = 52.0,\
                projection='mill',resolution ='l')
            x, y = m(lons, lats)

        ax = fig1.add_axes([0.02,.56,0.96,.4])
        hucnumber_in = int(hucnumber_save[icase])
        if hucnumber_in <= 18:
            chuc = 'HUC '+str(hucnumber_in)
        else:
            chuc = 'full CONUS'
        ctitle = 'Roebber/Kratsov case '+str(icase+1)+\
            ', '+str(yyyymmddhh_save[icase])+\
            ', optimized for '+chuc
        ax.set_title(ctitle,fontsize=12)
        CS2 = m.contourf(x, y, apcp_fcst1[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')
            
            
        if hucnumber_in == 19:
            mask = 0*lsmask_conus
        else:
            mask = mask_HUC2[hucnumber_in-1,:,:]    
        for jy in range(nyconus):
            for ix in range(nxconus):
                if mask[jy,ix] > 0:
                    xdot, ydot = m(lons_conus[jy,ix]-360., lats_conus[jy,ix])
                    m.plot(xdot,ydot,marker='.',markersize=0.25,\
                        color='Black',markerfacecolor='Black')
                
        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
        

        # ---- move on to the second day to plot on this page. 
        #      open the appropriate file and read the forecast for this case day
    
        yyyymmddhh_in = int(yyyymmddhh_save[icase+1])
        cyyyy = str(yyyymmddhh_save[icase+1])[0:4]
        cyyyymmddhh_start = cyyyy + '010100'
        cyyyymmddhh_end = cyyyy + '123100'
        cyyyymmddhh_list = daterange(cyyyymmddhh_start, cyyyymmddhh_end, 24)
        iyyyymmddhh_list = [int(i) for i in cyyyymmddhh_list]
        infile = data_directory_fcst + 'apcp.'+cyyyy+'.nc'
        nc = Dataset(infile)    
        idx = iyyyymmddhh_list.index(yyyymmddhh_in)
        apcp_fcst2 = 1000.*np.squeeze(nc.variables['apcp'][idx,:,:])
        nc.close()   
    
        ax = fig1.add_axes([0.02,.11,0.96,.4])
        hucnumber_in = int(hucnumber_save[icase+1])
        if hucnumber_in <= 18:
            chuc = 'HUC '+str(hucnumber_in)
        else:
            chuc = 'full CONUS'
        ctitle = 'Roebber/Kratsov case '+str(icase+2)+\
            ', '+str(yyyymmddhh_save[icase+1])+\
            ', optimized for '+chuc
        ax.set_title(ctitle,fontsize=12)
        CS2 = m.contourf(x, y, apcp_fcst2[:,:], clevs,\
            cmap=None, colors=colorst, extend='both')
        if hucnumber_in == 19:
            mask = 0*lsmask_conus
        else:
            mask = mask_HUC2[hucnumber_in-1,:,:]
        for jy in range(nyconus):
            for ix in range(nxconus):
                if mask[jy,ix] > 0:
                    xdot, ydot = m(lons_conus[jy,ix]-360., lats_conus[jy,ix])
                    m.plot(xdot,ydot,marker='.',markersize=0.25,\
                        color='Black',markerfacecolor='Black')
            
        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
        
        # ---- use axes_grid toolkit to make colorbar axes.

        cax = fig1.add_axes([0.02,0.07,0.96,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Roebber/Kratsov 10-day total mean precipitation (mm)',fontsize=9)
    
        pdf.savefig()
        plt.close()
        
print ('saved to ', outfile)

