"""
gefsv12_stamp_raw_qmapped_offset.py cyyyymmddhh clead
"""

import os, sys
from datetime import datetime
from dateutils import dateshift
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cf
import _pickle as cPickle
import scipy.stats as stats
from pyproj import Proj
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'


reader = shpreader.Reader('countyl010g_shp_nt00964/countyl010g.shp')
counties = list(reader.geometries())
COUNTIES = cf.ShapelyFeature(counties, ccrs.PlateCarree())

reader = shpreader.Reader('statesl010g/statesp010g.shp')
states = list(reader.geometries())
STATES = cf.ShapelyFeature(states, ccrs.PlateCarree())

shpfilename = shpreader.natural_earth(resolution='10m',
    category='cultural', name='admin_0_countries')
reader = shpreader.Reader(shpfilename)
mexico = [country for country in reader.records() if country.attributes["NAME_LONG"] == "Mexico"][0]

# =====================================================================
# =====================================================================

# ---- inputs from command line

cyyyymmddhh = sys.argv[1] # 
clead = sys.argv[2] # 018 not 18
ilead = int(clead)
cyear = cyyyymmddhh[0:4]
cmonth = cyyyymmddhh[4:6]
imonth = int(cmonth)-1
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
ccmonth = cmonths[imonth]
master_directory = '/Volumes/NBM/conus_gefsv12/qmapped/'
jnd = 491
ind = 690

#latb = 36.7 # Colorado domain
#late = 41.2 # Colorado domain
#lonb = -110. # Colorado domain
#lone = -100.9 # Colorado domain



# ---- read ensemble data

infile = master_directory + cyyyymmddhh+'_use99_lead='+clead+'.cPick'
print ('reading raw and qmapped ens from ', infile) 
inf = open(infile,'rb')
precip_ens_raw = cPickle.load(inf)
precip_ens_qmapped = cPickle.load(inf)
lons_ndfd = cPickle.load(inf)
lats_ndfd = cPickle.load(inf)
offset_ens = cPickle.load(inf)
quantile_99 = cPickle.load(inf)
fraction_zero_gefsv12_on_ndfd = cPickle.load(inf)
fraction_zero_ndfd = cPickle.load(inf)
usegamma_ndfd = cPickle.load(inf)
forecast_quantiles_ndfd = cPickle.load(inf)
analysis_quantile_ens = cPickle.load(inf)
analysis_quantile_nozero_ens = cPickle.load(inf)
ny_ndfd, nx_ndfd = np.shape(usegamma_ndfd)

print ('   --- after reading from file ')
print ('   offset_ens[2,jnd,ind] = ', offset_ens[2,jnd,ind])
print ('   forecast_quantiles_ndfd[2,jnd,ind]', forecast_quantiles_ndfd[2,jnd,ind])
print ('   analysis_quantile_ens[2,jnd,ind]', analysis_quantile_ens[2,jnd,ind])
print ('   analysis_quantile_nozero_ens[2,jnd,ind]', analysis_quantile_nozero_ens[2,jnd,ind])
inf.close()

rlat_point = lats_ndfd[jnd,ind]
rlon_point = lons_ndfd[jnd,ind]
print  ('rlon_point, rlat_point = ',rlon_point, rlat_point)

# ---- plot the stamp map

clevs = [0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, \
    10.0, 12.5, 15.0, 17.5, 20.0, 25.0, 30.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']


# ---- plot the 99th percentile over Arizona  

#latb = 33. # Rockies
#late = 47.
#lonb = -117.
#lone = -101.9

latb = 27.
late = 41.
lonb = -120.
lone = -104.9

xdim = 6.
ydim = 8.
drawcoasts = False

proj = ccrs.LambertConformal(\
    central_latitude = (latb+late)/2.,
    central_longitude = (lonb+lone)/2,
    standard_parallels = (latb, late))
    
xbegin = [0.025,0.185,0.345,0.505,0.665,0.825]
xlen = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15]
ybegin = [.76, .59, .42, .25, .08]
ylen = [0.17, 0.17,0.17, 0.17,0.17]

process_misc = False
if process_misc == True:

    print ('processing misc figures')
    fig = plt.figure(figsize=(xdim, ydim))
    axloc = [0.02,0.1,0.96,0.81]
    ax = plt.axes(axloc,projection = proj)
    ax.set_extent([lonb,lone,latb,late])
    ax.coastlines(resolution='50m',lw=0.5)
    ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
    ax.add_feature(cf.BORDERS,lw=0.5)
    ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)

    title = '99th percentile of '+cmonth+\
        '\n6-h GEFSv12 forecast rainfall ending '+clead+' h'
    ax.set_title(title, fontsize=16,color='Black')
    CS = ax.contourf(lons_ndfd, lats_ndfd, \
        quantile_99, clevs, cmap=None, colors=colorst, \
        extend='both', transform=ccrs.PlateCarree())
    
    ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
        linestyle='',color='Black',transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.02,0.08,0.96,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=9)
    cb.set_label('99th percentile precipitation amount (mm)',\
        fontsize=11)

    plot_title = 'q99_'+cmonth+'_Arizona_lead'+clead+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)


    # ---- plot the fraction zero over Arizona  

    fig = plt.figure(figsize=(xdim, ydim))
    axloc = [0.02,0.1,0.96,0.81]
    ax = plt.axes(axloc,projection = proj)
    ax.set_extent([lonb,lone,latb,late])
    ax.coastlines(resolution='50m',lw=0.5)
    ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
    ax.add_feature(cf.BORDERS,lw=0.5)
    ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)

    title = 'GEFSv12 1 minus fraction zero ending '+clead+' h'
    clevs = [-100.0, 0.0,0.03,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    ax.set_title(title, fontsize=16,color='Black')
    print ('min, max fraction_zero_gefsv12_on_ndfd = ', np.min(fraction_zero_gefsv12_on_ndfd),\
         np.max(fraction_zero_gefsv12_on_ndfd))

    CS = ax.contourf(lons_ndfd, lats_ndfd, \
        1.-fraction_zero_gefsv12_on_ndfd, clevs, cmap=None, colors=colorst, \
        extend='both', transform=ccrs.PlateCarree())
    
    ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
        linestyle='',color='Black',transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.02,0.08,0.96,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=9)
    cb.set_label('1 - fraction zero (fraction rainy)',\
        fontsize=11)

    plot_title = 'gefsv12_1mfraction_zero_'+cmonth+'_Arizona_lead'+clead+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)


    # ---- plot the fraction zero over Arizona  

    plot_fzero_anal = False
    if plot_fzero_anal == True:
        fraction_zero_ndfd = np.where(fraction_zero_ndfd < 0.0, np.zeros((ny_ndfd, nx_ndfd),\
             dtype=np.float64),fraction_zero_ndfd )

        fig = plt.figure(figsize=(xdim, ydim))
        axloc = [0.02,0.1,0.96,0.81]
        ax = plt.axes(axloc,projection = proj)
        ax.set_extent([lonb,lone,latb,late])
        ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(cf.BORDERS,lw=0.5)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)

        title = 'Analysis 1 minus fraction zero ending '+clead+' h'
        ax.set_title(title, fontsize=16,color='Black')
        print ('min, max fraction_zero_ndfd = ', np.min(fraction_zero_ndfd), np.max(fraction_zero_ndfd))
        CS = ax.contourf(lons_ndfd, lats_ndfd, \
            1.-fraction_zero_ndfd, clevs, cmap=None, colors=colorst, \
            extend='both', transform=ccrs.PlateCarree())
        print ('after contour')  
        ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
            linestyle='',color='Black',transform=ccrs.PlateCarree())

        cax = fig.add_axes([0.02,0.08,0.96,0.02])
        cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=9)
        cb.set_label('1 - fraction zero (fraction rainy)',\
            fontsize=11)

        plot_title = 'analysis_1mfraction_zero_'+cmonth+'_Arizona_lead'+clead+'h.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)

    # ---- plot the usegamma over Arizona  

    fig = plt.figure(figsize=(xdim, ydim))
    axloc = [0.02,0.1,0.96,0.81]
    ax = plt.axes(axloc,projection = proj)
    ax.set_extent([lonb,lone,latb,late])
    ax.coastlines(resolution='50m',lw=0.5)
    ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
    ax.add_feature(cf.BORDERS,lw=0.5)
    ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)

    title = 'Points using Gamma distributions, '+clead+' h'
    ax.set_title(title, fontsize=16,color='Black')
    clevs = [0.0,1.0]
    colorst = ['white', 'black']
    CS = ax.contourf(lons_ndfd, lats_ndfd, \
        usegamma_ndfd, clevs, cmap=None, colors=colorst, \
        extend='both', transform=ccrs.PlateCarree())
    
    ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
        linestyle='',color='Black',transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.02,0.08,0.96,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=9)
    cb.set_label('Fraction zero',\
        fontsize=11)

    plot_title = 'usegamma_'+cmonth+'_Arizona_lead'+clead+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)




latb = 29.7 # Arizona domain
late = 34.2 # Arizona domain
lonb = -118. # Arizona domain
lone = -108.9 # Arizona domain

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
Mexico = ShapelyFeature([mexico.geometry], ccrs.PlateCarree(), \
    facecolor="none", edgecolor='black', lw=0.5)
    
proj = ccrs.LambertConformal(\
    central_latitude = (latb+late)/2.,
    central_longitude = (lonb+lone)/2.,
    standard_parallels = (latb, late))

# ----- plot analysis quantile data.

process_aquantile = False
if process_aquantile == True:
    print ('processing full analysis quantile ')
    fig = plt.figure(figsize=(9.,6.5))

    ilead = int(clead)

    fig.suptitle('CCPA quantile in full distribution, IC = '+\
        cyyyymmddhh+', 6 h period ending '+clead+' h',\
        fontsize=12,color='Black')
    
    clevs = [0.0, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 0.992, 0.995, 0.997]
    for imem in range(30):
        print ('processing ',imem)
        irow = imem//6
        icol = imem - irow*6
        axloc = [xbegin[icol],ybegin[irow],xlen[icol],ylen[irow]]
        ax = plt.axes(axloc, projection = proj)
        if imem == 0:
            ax.set_title('Control', fontsize=8,color='Black')
        else:
            ax.set_title('Perturbed member '+str(imem), fontsize=8,color='Black')
        ax.set_extent([lonb,lone,latb,late])
        #ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        ax.add_feature(Mexico)
        
        ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
            linestyle='',color='Black',transform=ccrs.PlateCarree())
        aq = analysis_quantile_ens[imem,:,:]        
        CS = ax.contourf(lons_ndfd, lats_ndfd, \
            aq, clevs, cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.01,0.065,0.98,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Analysis quantile in full distribution, including zeros',fontsize=9)

    plot_title = 'aquantile_stamp_'+cyyyymmddhh+'_lead'+clead+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)





process_aquantile_nonzero = False
if process_aquantile_nonzero == True:
    
    print ('processing nonzero analysis quantile ')
    
    # ----- plot analysis quantile data without zeros.

    fig = plt.figure(figsize=(9.,6.5))
    ilead = int(clead)
    fig.suptitle('CCPA quantile of nonzeros, IC = '+\
        cyyyymmddhh+', 6 h period ending '+clead+' h',\
        fontsize=12,color='Black')
    clevs = [0.0, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 0.992, 0.995, 0.997]
    
    for imem in range(30):
        print ('processing ',imem)
        irow = imem//6
        icol = imem - irow*6
        axloc = [xbegin[icol],ybegin[irow],xlen[icol],ylen[irow]]
        ax = plt.axes(axloc, projection = proj)
        if imem == 0:
            ax.set_title('Control', fontsize=8,color='Black')
        else:
            ax.set_title('Perturbed member '+str(imem), fontsize=8,color='Black')
        ax.set_extent([lonb,lone,latb,late])
        #ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        ax.add_feature(Mexico)
        
        ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
            linestyle='',color='Black',transform=ccrs.PlateCarree())
        aq = analysis_quantile_nozero_ens[imem,:,:]        
        CS = ax.contourf(lons_ndfd, lats_ndfd, \
            aq, clevs, cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.01,0.065,0.98,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Analysis quantile in full distribution, NOT including zeros',fontsize=9)

    plot_title = 'aquantile_nonzero_stamp_'+cyyyymmddhh+'_lead'+clead+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)



process_fquantile = True
if process_fquantile == True:
    
    print ('processing forecast quantile ')
    clevs = [0.0, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 0.99, 0.992, 0.995, 0.997]
    fig = plt.figure(figsize=(9.,6.5))
    ilead = int(clead)    
    fig.suptitle('Forecast quantile including zeros, IC = '+\
        cyyyymmddhh+', 6 h period ending '+clead+' h',\
        fontsize=12,color='Black')
            
    for imem in range(30):
        print ('processing ',imem)
        irow = imem//6
        icol = imem - irow*6
        axloc = [xbegin[icol],ybegin[irow],xlen[icol],ylen[irow]]
        ax = plt.axes(axloc, projection = proj)
        if imem == 0:
            ax.set_title('Control', fontsize=8,color='Black')
        else:
            ax.set_title('Perturbed member '+str(imem), fontsize=8,color='Black')
        ax.set_extent([lonb,lone,latb,late])
        #ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        ax.add_feature(Mexico)
        
        ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
            linestyle='',color='Black',transform=ccrs.PlateCarree())
    
        fq = forecast_quantiles_ndfd[imem,:,:]
        print ('imem, fq[jnd,ind] = ', imem, fq[jnd,ind])
        print ('max, min forecast quantile = ', np.max(fq), np.min(fq))
        CS = ax.contourf(lons_ndfd, lats_ndfd, \
            fq, clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.01,0.065,0.98,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Forecast quantile in full distribution, including zeros',fontsize=9)

    plot_title = 'fquantile_stamp_'+cyyyymmddhh+'_lead'+clead+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)

# ----- plot raw data.



process_raw = False
if process_raw == True:
    
    print ('processing raw')
    fig = plt.figure(figsize=(9.,6.5))
    ilead = int(clead)
    fig.suptitle('Raw GEFSv12 total precipitation, IC = '+\
        cyyyymmddhh+', 6 h period ending '+clead+' h',\
        fontsize=12,color='Black')
    clevs = [0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, \
        10.0, 12.5, 15.0, 17.5, 20.0, 25.0, 30.0]
    for imem in range(30):
        print ('processing ',imem)
        irow = imem//6
        icol = imem - irow*6
        axloc = [xbegin[icol],ybegin[irow],xlen[icol],ylen[irow]]
        ax = plt.axes(axloc, projection = proj)
        if imem == 0:
            ax.set_title('Control', fontsize=8,color='Black')
        else:
            ax.set_title('Perturbed member '+str(imem), fontsize=8,color='Black')
        ax.set_extent([lonb,lone,latb,late])
        #ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        ax.add_feature(Mexico)
        
        ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
            linestyle='',color='Black',transform=ccrs.PlateCarree())
    
        CS = ax.contourf(lons_ndfd, lats_ndfd, precip_ens_raw[imem,:,:], clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.01,0.065,0.98,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Precipitation amount (mm)',fontsize=9)

    plot_title = 'raw_stamp_'+cyyyymmddhh+'_lead'+clead+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)

# ----- plot quantile-mapped data.

process_qmap = True
if process_qmap == True:
    print ('processing qmapped')
    clevs = [0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, \
        10.0, 12.5, 15.0, 17.5, 20.0, 25.0, 30.0]
    fig = plt.figure(figsize=(9.,6.5))
    fig.suptitle('Quantile-mapped GEFSv12 total precipitation, IC = '+\
        cyyyymmddhh+', 6-h period ending '+clead+' h',\
        fontsize=12,color='Black')
    for imem in range(30):
        print ('processing ',imem)
        irow = imem//6
        icol = imem - irow*6
        axloc = [xbegin[icol],ybegin[irow],xlen[icol],ylen[irow]]
        ax = plt.axes(axloc, projection = proj)
        if imem == 0:
            ax.set_title('Control', fontsize=8,color='Black')
        else:
            ax.set_title('Perturbed member '+str(imem), fontsize=8,color='Black')
        ax.set_extent([lonb,lone,latb,late])
        #ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        ax.add_feature(Mexico)
        ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
            linestyle='',color='Black',transform=ccrs.PlateCarree())
    
        CS = ax.contourf(lons_ndfd, lats_ndfd, precip_ens_qmapped[imem,:,:], clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.01,0.065,0.98,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('Precipitation amount (mm)',fontsize=9)

    plot_title = 'qmapped_stamp_'+\
        cyyyymmddhh+'_lead'+clead+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)

# ----- plot offset.

process_offset = True
if process_offset == True:
    print ('processing offset')
    fig = plt.figure(figsize=(9.,6.5))
    fig.suptitle('Offset of GEFSv12 precipitation, IC = '+\
        cyyyymmddhh+', 6-h period ending '+clead+' h',\
        fontsize=12,color='Black')
    for imem in range(30):
        print ('processing ',imem)
        irow = imem//6
        icol = imem - irow*6
        axloc = [xbegin[icol],ybegin[irow],xlen[icol],ylen[irow]]
        ax = plt.axes(axloc, projection = proj)
        if imem == 0:
            ax.set_title('Control', fontsize=8,color='Black')
        else:
            ax.set_title('Perturbed member '+str(imem), fontsize=8,color='Black')
        ax.set_extent([lonb,lone,latb,late])
        #ax.coastlines(resolution='50m',lw=0.5)
        ax.add_feature(COUNTIES, facecolor='none', edgecolor='gray',lw=0.13)
        ax.add_feature(STATES, facecolor='none', edgecolor='black',lw=0.5)
        ax.add_feature(Mexico)
        ax.plot(rlon_point, rlat_point, markersize=1,marker='o',\
            linestyle='',color='Black',transform=ccrs.PlateCarree())
        CS = ax.contourf(lons_ndfd, lats_ndfd, offset_ens[imem,:,:], clevs,\
            cmap=None, colors=colorst, extend='both', \
            transform=ccrs.PlateCarree())

    cax = fig.add_axes([0.01,0.065,0.98,0.02])
    cb = plt.colorbar(CS,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label('GEFSv12 forecast offset from 99th percentile (mm)',fontsize=9)

    plot_title = 'offset_99_stamp_'+cyyyymmddhh+'_lead'+clead+'h.png'
    fig.savefig(plot_title, dpi=300)
    plt.close()
    print ('saving plot to file = ',plot_title)
    


    



