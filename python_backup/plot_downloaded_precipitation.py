"""
plot_downloaded_precipitation.py

"""
import pygrib
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os, sys

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

date = sys.argv[1] # yyyymmddhh format
cleadb = sys.argv[2] # end lead time in hours
cleade = sys.argv[3]
ileadb = int(cleadb)
ileade = int(cleade)
cyear = date[0:4]

url = 'https://noaa-gefs-retrospective.s3.amazonaws.com/GEFSv12/reforecast/'+\
    cyear + '/'+date+ '/c00/Days:1-10/apcp_sfc_'+date+'_c00.grib2'
file = 'apcp_sfc_'+date+'_c00.grib2'
exist = os.path.exists(file)
if exist == False:
    cmd = 'wget '+url
    try:
        istat = os.system(cmd)
    except:
        print ('could not find data for ', date)
    
url = 'https://noaa-gefs-retrospective.s3.amazonaws.com/GEFSv12/reforecast/'+\
    cyear + '/'+date+ '/c00/Days:1-10/acpcp_sfc_'+date+'_c00.grib2'
file = 'acpcp_sfc_'+date+'_c00.grib2'   
exist = os.path.exists(file)
if exist == False:
    cmd = 'wget '+url
    try:
        istat = os.system(cmd)
    except:
        print ('could not find data for ', date)


local_gribfile1 = 'apcp_sfc_'+date+'_c00.grib2'
local_gribfile2 = 'acpcp_sfc_'+date+'_c00.grib2'

grbfile = pygrib.open(local_gribfile1)
grb = grbfile.select(startStep=ileadb, endStep=ileade)[0]
total_precip = grb.values
lats, lons = grb.latlons()
lons = lons-360.
grbfile.close()

grbfile = pygrib.open(local_gribfile2)
grb = grbfile.select(startStep=ileadb, endStep=ileade)[0]
convective_precip = grb.values
grbfile.close()
        
# --- now plot the forecast precipitation 

clevs = [0.0,0.01,0.05,0.1,0.2,0.5,0.75,1.0,1.5,2.0,3.0,5.0,10.0,25.0,50.0]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

fig = plt.figure(figsize=(6.,8.))
axloc = [0.02,0.56,0.96,0.4]
ax1 = fig.add_axes(axloc)
title = '(a) '+cleadb+' to '+cleade+' hour total precipitation, initial date = '+date
ax1.set_title(title, fontsize=10,color='Black')
m = Basemap(llcrnrlon=-130.0,llcrnrlat=15.0,\
    urcrnrlon=-55.0,urcrnrlat=55.0,\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,total_precip, clevs, \
    cmap=None, colors=colorst, extend='both') #
m.drawcoastlines(linewidth=0.5,color='Black')
m.drawcountries(linewidth=0.5,color='Black')
m.drawstates(linewidth=0.5,color='Black')

axloc = [0.02,0.12,0.96,0.4]
ax1 = fig.add_axes(axloc)
title = '(b) '+cleadb+' to '+cleade+' hour convective precipitation, initial date = '+date
ax1.set_title(title, fontsize=10,color='Black')
m = Basemap(llcrnrlon=-130.0,llcrnrlat=15.0,\
    urcrnrlon=-55.0,urcrnrlat=55.0,\
    resolution='l', projection='mill')
x, y = m(lons, lats)
CS2 = m.contourf(x,y,convective_precip, clevs, \
    cmap=None, colors=colorst, extend='both') #
m.drawcoastlines(linewidth=0.5,color='Black')
m.drawcountries(linewidth=0.5,color='Black')
m.drawstates(linewidth=0.5,color='Black')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
#cax = divider.append_axes("bottom", size="3%", pad=0.07)
cax = fig.add_axes([0.02,0.06,0.96,0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs, format='%g') # ticks=clevs,
cb.ax.tick_params(labelsize=7)
cb.set_label(r'precipitation amount (mm)',fontsize=9)

# ---- set plot title

plot_title = 'precipitation_'+date+'_'+cleadb+'_to_'+cleade+'.png'
fig.savefig(plot_title, dpi=300,fontsize=9)

print ('saving plot to file = ',plot_title)
print ('Done!')

