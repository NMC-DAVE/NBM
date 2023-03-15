"""
plot_reforecast_precip_conus.py cyyyymmddhh clead cmem

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


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

# =====================================================================

def define_month(imonth):
    cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    return cmonths[imonth]

# =====================================================================
    
# ---- inputs from command line

cyyyymmddhh_init = sys.argv[1] # year-month-day-hour of forecast initial time
clead = sys.argv[2] # lead time of end of precip accumulation period; 006, 012, etc.
cmem = sys.argv[3] # 0-4; 0=control, 1-4 = perturbed members

# ---- various setup

iyyyymmddhh_init = int(cyyyymmddhh_init)
imonth = int(cyyyymmddhh_init[4:6])-1 # month 0-11
ilead = int(clead)
ileadb = ilead-6
if ileadb < 10:
    cleadb = '00'+str(ileadb)
elif ileadb > 10 and ileadb < 100:
    cleadb = '0'+str(ileadb)
else:
    cleadb = str(ileadb)
print ('cleadb, ileadb ', cleadb, ileadb)
cmonth = define_month(imonth)
master_directory = '/Volumes/NBM/conus_gefsv12/precip/netcdf/'
imem = int(cmem)


ncfile = master_directory + cmonth + '_conus_reforecast_precip_h'+clead+'.nc'
print (ncfile)
fexist = os.path.exists(ncfile)
if fexist == True:
    
    # ---- read in data.
    
    nc = Dataset(ncfile)
    lons_1d = nc.variables['lons_fcst'][:]
    lats_1d = nc.variables['lats_fcst'][:]
    yyyymmddhh_init = nc.variables['yyyymmddhh_init'][:]
    yyyymmddhh_fcst = nc.variables['yyyymmddhh_init'][:]
    idx = np.where(yyyymmddhh_init == iyyyymmddhh_init)[0]
    idxuse = idx[0] + imem
    print ('idxuse = ', idxuse)
    precip = nc.variables['apcp_fcst'][idxuse,:,:]
    nc.close()
    lons_2d, lats_2d = np.meshgrid(lons_1d,lats_1d)

    # ---- make plot of the precipitation

    m = Basemap(llcrnrlon=lons_2d[0,0],llcrnrlat=lats_2d[-1,-1],\
        urcrnrlon=lons_2d[-1,-1],urcrnrlat=lats_2d[0,0],\
        resolution='l', projection='mill')
    x, y = m(lons_2d, lats_2d)
    clevs = [0.0,0.25,0.5,1,2,3,5,7,10,15,20,25,35,50]
    colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
        '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
        '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    fig = plt.figure(figsize=(9,6.))
    axloc = [0.05,0.1,0.9,0.8]
    ax1 = fig.add_axes(axloc)
    title = 'GEFSv12 reforecast precipitation, IC = '+\
        cyyyymmddhh_init+', member '+cmem+', '+cleadb+' to '+clead+' h forecast'
    ax1.set_title(title, fontsize=13,color='Black')
    CS2 = m.contourf(x, y, precip, clevs,\
        cmap=None, colors=colorst, extend='both')
    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawcounties(linewidth=0.1,color='LightGray')
    m.drawstates(linewidth=0.3,color='Gray')
    
    # ---- use axes_grid toolkit to make colorbar axes.

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("bottom", size="3%", pad=0.28)
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('precipitation accumulation (mm)',fontsize=9)

    # ---- set plot title

    plot_title = 'precip_fcst_'+cyyyymmddhh_init+'_mem'+cmem+'_h'+clead+'.png'
    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)
    print ('Done!')


else:
    print ('Invalid date/lead/mem ', cyyyymmddhh_init, clead, cmem)
    print ('Please check and try again.')