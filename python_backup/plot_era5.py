from netCDF4 import Dataset
import numpy as np
from dateutils import daterange, datetohrs, dayofyear
import sys
import pygrib
import os
import os.path
from os import path
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'

cyyyymmddhh = sys.argv[1]
cleade = sys.argv[2]
infile = 'ecmwf/2t_'+cyyyymmddhh+'_f'+cleade+'.grib2'
print (infile)
fexist = path.exists(infile)
if fexist == True:
        
    # ---- read in the skin temperatures from the operational
        
    grbfile = pygrib.open(infile)
    grb = grbfile.select(shortName = '2t')[0]
    t2m = grb.values - 273.15
    print ('min, max t2m = ',np.min(t2m), np.max(t2m))
    lats, lons = grb.latlons()
    lons = lons-360.
    if cleade == '0':
        t2m = np.flipud(t2m)
        lats = np.flipud(t2m)
    print ('lats = ', lats[:])
    print ('lons = ', lons[:])
    nlats, nlons = np.shape(lats)
    print ('nlats, nlons = ', nlats, nlons)
    grbfile.close()
    #sys.exit()

else:
    print ('unable to read ', infile)
    sys.exit()

#lons = np.flipud(lons)
#lats = np.flipud(lats)

# ---- code to use for plotting

fig = plt.figure(figsize=(9.,5.6))
axloc = [0.07,0.11,0.9,0.82]
ax = fig.add_axes(axloc)
ax.set_title('ECMWF 2-m temperature, initial date =  '+cyyyymmddhh+', lead time = '+cleade+' h')
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000']
clevs = [-30,-20,-10,-5,5,10,20,30]
colorstblack='Black'
m = Basemap(llcrnrlon=-125,llcrnrlat=20,\
    urcrnrlon=-60,urcrnrlat=50.,\
    projection='mill',resolution='l')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,t2m,clevs,cmap=None,colors=colorst,extend='both')
m.drawcoastlines(linewidth=0.5,color='Gray')
m.drawcountries(linewidth=0.3,color='Gray')
m.drawstates(linewidth=0.3,color='Gray')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax)
cax = divider.append_axes("bottom", size="3%", pad=0.35)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.set_label('Temperature (deg C)')

# ---- set plot title

plot_title = 't2m_'+cyyyymmddhh+'_f'+cleade+'.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')

