""" create bias correlation matrix
"""
import pygrib
import os, sys
from datetime import datetime
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.signal as signal
import scipy.stats as stats
from calculate_terrain_gradients_f90 import calculate_terrain_gradients_f90

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx
    
# =====================================================================    

hlocal = 10. # measured in grid points
vlocal = 300. # measured in meters
#tlocal = 10.0
lsmask_local = 0.05


# --- read terrain and land-sea masks for ERA5 grid over CONUS

infile = 'lsmask.grib'
grbfile = pygrib.open(infile)
grb = grbfile.select()[0] # shortName = lsm
lats, lons = grb.latlons()
lons = lons-360
lsmask = grb.values
grbfile.close()
ny, nx = np.shape(lsmask)
lons = np.flipud(lons)
lats = np.flipud(lats)
lsmask = np.flipud(lsmask)

lons_1d = lons[0,:]
lats_1d = lats[:,0]


infile = 'terrain_height.grib'
grbfile = pygrib.open(infile)
grb = grbfile.select()[0] # shortName = orog
terrain_height = grb.values
terrain_height = np.flipud(terrain_height)
print ('max, min terrain_height = ', np.max(terrain_height), np.min(terrain_height))
grbfile.close()

# --- calculate terrain height gradients

#earth_radius_meters = 6378100.
#print (calculate_terrain_gradients_f90.__doc__)
#terrain_gradx, terrain_grady = calculate_terrain_gradients_f90(\
#    earth_radius_meters, terrain_height, lats, ny, nx)
#txmax = np.max(terrain_gradx)
#txmin = np.min(terrain_gradx)
#tymax = np.max(terrain_grady)
#tymin = np.min(terrain_grady)
#print ('max terrain_gradx, min terrain_gradx = ', txmax, txmin)
#print ('max terrain_grady, min terrain_grady = ', tymax, tymin)
#txfactor = np.max([txmax, -1.*txmin])
#tyfactor = np.max([tymax, -1.*tymin])
#terrain_gradx = terrain_gradx / txfactor
#terrain_grady = terrain_grady / tyfactor
    
# --- develop a model of correlations as a function of horizontal distance,
#     land-sea mask, terrain_height difference, and terrain gradient difference

already = False
if already == False:
    
    correlation_bias = np.zeros((ny, nx, ny, nx), dtype=np.float64)
    for ix1 in range(nx):
        print ('processing ',ix1,' of ', nx)
        for jy1 in range(ny):
            for ix2 in range(nx):
                for jy2 in range(ny):
                
                    hdist_factor = np.sqrt( (ix1-ix2)**2 + (jy1-jy2)**2) / hlocal
                    vdist_factor = np.abs(terrain_height[jy1,ix1] - terrain_height[jy2,ix2]) / vlocal
                    #grad_factor = np.sqrt((terrain_gradx[jy1,ix1] - terrain_gradx[jy2,ix2])**2 + \
                    #    (terrain_grady[jy1,ix1] - terrain_grady[jy2,ix2])**2) / tlocal
                    lsmask_factor = np.abs(lsmask[jy1,ix1]-lsmask[jy2,ix2])**2 / lsmask_local
                    dist = np.sqrt(hdist_factor**2 + vdist_factor**2 + lsmask_factor**2)
                    if dist != dist:
                        print ('improper distance at ix1,jy1,ix2,jy2 = ', ix1,jy1,ix2,jy2)
                        print ('hdist_factor, vdist_factor, grad_factor, lsmask_factor = ', \
                            hdist_factor, vdist_factor, grad_factor, lsmask_factor)
                        sys.exit()
                    c = np.exp(-dist**2)
                    if ix1 == ix2 and jy1 == jy2:
                        correlation_bias[jy1,ix1,jy2,ix2] = 1.0
                    else:
                        correlation_bias[jy1,ix1,jy2,ix2] = 0.03*np.exp(-dist**2)
                
    # --- save correlation of bias to file

    cfile = 'correlation_bias_ERA5grid.cPick'
    ouf = open(cfile, 'wb')
    cPickle.dump(correlation_bias, ouf)
    ouf.close()
    
else:
    
    # --- save correlation of bias to file

    cfile = 'correlation_bias_ERA5grid.cPick'
    inf = open(cfile, 'rb')
    correlation_bias = cPickle.load(inf)
    inf.close()
    

# ---- make some sample plots

plot_lons = [-105.0, -123.0, -105.5, -100.0, -90.0, -122.5, -126.0]
plot_lats = [40.0, 45.0, 40.0, 40.0, 35.0, 45.0, 45.0 ]
nplots = len(plot_lons)

for rlon, rlat in zip(plot_lons, plot_lats):
    
    ilon = find_nearest(lons_1d, rlon)
    ilat = find_nearest(lats_1d, rlat)
    clon = str(rlon)
    clat = str(rlat)
    bias_error_corr_map = correlation_bias[:,:,ilat,ilon]
    print ('max, min bias_error_corr_map ',\
        np.max(bias_error_corr_map), np.min(bias_error_corr_map))

    # --- now plot the single-point correlation model for the selected point 

    clevs = [-0.99,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9,0.99]
    colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
        'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000'] 
    
    fig = plt.figure(figsize=(6.,4.2))
    axloc = [0.02,0.09,0.96,0.82]
    ax1 = fig.add_axes(axloc)
    title = r'Estimated T$_{2m}$ bias error correlation map, '+ \
        ' lon = '+clon+', lat = '+clat
    ax1.set_title(title, fontsize=11,color='Black')
    m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[0,0],\
        urcrnrlon=lons[-1,-1],urcrnrlat=lats[-1,-1],\
        resolution='l', projection='mill')
    x, y = m(lons, lats)
    CS2 = m.contourf(x,y,bias_error_corr_map,clevs,\
        cmap=None,colors=colorst,extend='both')
    xdot, ydot = m(rlon,rlat)
    m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.4,color='Gray')

    # ---- use axes_grid toolkit to make colorbar axes.

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("bottom", size="3%", pad=0.1)
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.set_label('Temperature bias correlation')

    # ---- set plot title

    plot_title = 't2m_bias_correlation'+'_lon'+clon+'_lat'+clat+'.png'
    fig.savefig(plot_title, dpi=300,fontsize=9)
    print ('saving plot to file = ',plot_title)
    
print ('Done!')



