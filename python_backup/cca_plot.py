"""
cca_plot.py: user inputs the date and the forecast lead time. program
then loads in the HUC precipitation and the ERA-5 analysis data,
and it plots out the prediction from previously calculated CCA analysis.

"""

from netCDF4 import Dataset
import numpy as np
from dateutils import daterange, datetohrs, dateshift
from datetime import datetime
import sys
import os
from os import path
import numpy.ma as ma
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap
import scipy.stats as stats
from sklearn.cross_decomposition import CCA
import shapefile
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon

rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'
rcParams['legend.fontsize']='xx-small'
rcParams['legend.fancybox']=True

def load_era5(cyear, cmonth, cvar):
    
    infile = 'cca/'+cyear+cmonth+'_'+cvar+'_era5.cPick'
    #print (infile)
    inf = open(infile,'rb')
    input_data = cPickle.load(inf)
    ntimes, nlevels, ny, nx = np.shape(input_data)
    lon = cPickle.load(inf)
    lat = cPickle.load(inf)
    yyyymmddhh = cPickle.load(inf)
    inf.close()
    return yyyymmddhh, input_data, ntimes, nlevels, ny, nx, lon, lat
    
def convert_to_yyyymmddhh(precipDates):
    
    # --- convert Matt's HUC date array to yyyymmddhh format, making the 
    #     assumption that the 7am local time is close enough to 12Z.
    
    npdates, nymd = np.shape(precipDates)
    yyyymmddhh = []
    #print ('npdates, nymd = ', npdates, nymd)
    for i in range(npdates):
        
        yyyy = str(precipDates[i,0])
        imm = precipDates[i,1]
        if imm < 10:
            cmm = '0'+str(imm)
        else:
            cmm = str(imm)
        idd = precipDates[i,2]
        if idd < 10:
            cdd = '0'+str(idd)
        else:
            cdd = str(idd)
        yyyymmddhh.append(int(yyyy+cmm+cdd+'12'))
        #print (precipDates[i,0], precipDates[i,1], precipDates[i,2], int(yyyy+cmm+cdd+'12'))
        #if i == 1000: sys.exit()
    return yyyymmddhh
        
def find_color_index(levels, Y):
    foundit = False
    n = len(levels)
    #print ('number of levels = ', n)
    for i in range(0,n-1):
        if Y >= levels[i] and Y < levels[i+1] and foundit == False:
            idxcolor = i
            foundit = True
    return idxcolor

# --- get inputs from command line

cdate = sys.argv[1]
clead = sys.argv[2]
cmonth = cdate[4:6]
cyear = cdate[0:4]
#cmonth = sys.argv[1] # 01 to 12
#clead = sys.argv[2] # forecast lead time in days, e.g., 2.5   Use half days so the 00 UTC
cyears = ['1981', '1982', '1983', '1984','1985', '1986', '1987', '1988', '1989', '1990', \
    '1991', '1992', '1993', '1994', '1995', '1996', '1997', '1998', '1999', '2000', \
    '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011']
nyears = len(cyears)
    
ifhour = int(float(clead)*24)
print ('forecast lead in hours: ', ifhour)
imonth = int(cmonth)

# --- read in the ERA5 data for this month.


yyyymmddhh_era5, input_data, ntimes, nlevels, ny, nx, lon, lat = \
    load_era5(cyear, cmonth, 'air')
nlons = len(lon)
nlats = len(lat)
lon2d, lat2d = np.meshgrid(lon,lat)
zeros = np.zeros((nlats, nlons), dtype=np.float32)    

# --- read in the HUC data provided by Matt Switanek

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('Reading HUC data. Current time is ', current_time)
f1 = open("cca/PRISM_pr_hucs_19810101-20180930.pickle",'rb')
precipHUCs = cPickle.load(f1) #the daily precip accumulations in mm (days,hucs)
ndayhucs, nhucs = np.shape(precipHUCs)
precipDates = cPickle.load(f1) #the dates
hucLats = cPickle.load(f1) #centroid lats of hucs
hucLons = cPickle.load(f1) #centroid lons of hucs
hucShapes = cPickle.load(f1) #embedded lists of huc boundaries
hucDivision4 = cPickle.load(f1) #the division 4 numeric codes of the hucs
f1.close()

# ---- read in the previously calculated CCA output

infile = 'cca/PLSR_regression_data_month='+cmonth+'_lead='+clead+'days.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
Ypred_full = cPickle.load(inf)
yyyymmddhh_ICdates = cPickle.load(inf)
inf.close()
ys = np.shape(Ypred_full)
yzeros = np.zeros(ys,dtype=np.float32)
Ypred_full = np.where (Ypred_full < 0.0, yzeros, Ypred_full)

# ---- convert the HUC precipitation dates into yyyymmddhh format

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print ('Converting HUC precipitation dates.  Current time is ', current_time)
yyyymmddhh_hucs = convert_to_yyyymmddhh(precipDates)

# ---- for the chosen month, load the ERA5 analysis data and extract the analysis at the initial time

yyyymmddhh, temp, ntimes, nlevels, ny, nx, lon, lat = load_era5(cyear, cmonth, 'air')
yyyymmddhh, shum, ntimes, nlevels, ny, nx, lon, lat = load_era5(cyear, cmonth, 'shum')    
yyyymmddhh, uwnd, ntimes, nlevels, ny, nx, lon, lat = load_era5(cyear, cmonth, 'uwnd')    
yyyymmddhh, vwnd, ntimes, nlevels, ny, nx, lon, lat = load_era5(cyear, cmonth, 'vwnd')
yyyymmddhh_list = yyyymmddhh.tolist()

timeindex = yyyymmddhh_list.index(int(cdate))    
temp_day = temp[timeindex,:,:]  
shum_day = shum[timeindex,:,:]  
uwnd_day = uwnd[timeindex,:,:]  
vwnd_day = vwnd[timeindex,:,:]  
    
# --- for an n.5 -day forecast, pluck the HUC data offset by +n.5 days.
#     got to get this data individually for the given month of each year.
#     yyyymmddhh_store has the vector of initial condition dates.
#     yyyymmddhh_hucs has the vector of HUC verification period (end)
    
fcst_date = int(dateshift(cdate, ifhour))
timeindex_PLSR = np.where(yyyymmddhh_ICdates == int(cdate)) [0]
Ypred2 = np.squeeze(Ypred_full[timeindex_PLSR,:])
timeindex_HUC = yyyymmddhh_hucs.index(fcst_date)
Y = precipHUCs[timeindex_HUC,:]

# --- get the patches and colors for the analyzed and forecast precip

patches_analyzed = []
patches_forecast = []
colors_analyzed = []
colors_forecast = []
levels = [0.0,0.1,0.3,0.5,1.0, 1.5,2.0,3.0,5.0,10.0,20.0,30.0,40.0,50.0,200.]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray'] 
        
# --- now plot the analyzed

fig1 = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.12,0.96,0.77]
ax1 = fig1.add_axes(axloc)
ax1.set_title('HUC4 analyzed precipitation for 24 hours ending '+ \
    str(fcst_date), fontsize=16,color='Black')
map = Basemap(llcrnrlon=-125,llcrnrlat=25,urcrnrlon=-65,urcrnrlat=51.,
    resolution='l', projection='mill')
patches_analyzed = []
colors_analyzed = []
for v in range(0,202): 
    #print (v, Y[v], levels)
    idxcolor = find_color_index(levels, Y[v])
    for v2 in range(0,len(hucShapes[v])):   
        a = hucShapes[v][v2]
        xlon = a[:,0]
        ylat = a[:,1]
        xs, ys = map(xlon, ylat) 
        patches_analyzed.append(Polygon(np.column_stack([xs, ys]), True) ) 
        colors_analyzed.append(colorst[idxcolor])
ax1.add_collection(PatchCollection(patches_analyzed, facecolor=colors_analyzed, \
    edgecolor='Gray', linewidths=0.3, zorder=2)) 
xc, yc = map(lon2d, lat2d)
cs1 = map.contourf(xc, yc, zeros, levels, colors=colorst,extend='neither')   
map.drawcoastlines()
map.drawcountries()
map.drawstates()

cax = fig1.add_axes([0.02,0.07,0.96,0.02])
cbar = plt.colorbar(cs1, orientation='horizontal',\
    cax=cax, extend='both', ticks=levels[0:-1], format='%g') 
cbar.ax.tick_params(labelsize=9)
cbar.set_label('24-h accumulated precipitation amount (mm)')

# ---- set plot title

plot_title = 'analyzed_precipitation_HUCs_'+str(fcst_date)+'.png'
print ('saving plot to ', plot_title)
fig1.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


# --- now plot the analyzed

fig1 = plt.figure(figsize=(10,6.2))
axloc = [0.02,0.12,0.96,0.77]
ax1 = fig1.add_axes(axloc)
ax1.set_title('HUC4 PLSR forecast precipitation for 24 hours ending '+ str(fcst_date)+\
    ' lead = '+clead+' days', fontsize=16,color='Black')
map = Basemap(llcrnrlon=-125,llcrnrlat=25,urcrnrlon=-65,urcrnrlat=51.,
    resolution='l', projection='mill')
patches_forecast = []
colors_forecast = []
print ('max, min Ypred2 = ', np.max(Ypred2), np.min(Ypred2))
for v in range(0,202): 
    #print (v, Ypred2[v], levels)
    idxcolor = find_color_index(levels, Ypred2[v])
    for v2 in range(0,len(hucShapes[v])):   
        a = hucShapes[v][v2]
        xlon = a[:,0]
        ylat = a[:,1]
        xs, ys = map(xlon, ylat) 
        patches_forecast.append(Polygon(np.column_stack([xs, ys]), True) ) 
        colors_forecast.append(colorst[idxcolor])

ax1.add_collection(PatchCollection(patches_forecast, facecolor=colors_forecast, \
    edgecolor='Gray', linewidths=0.3, zorder=2)) 
xc, yc = map(lon2d, lat2d)
cs1 = map.contourf(xc, yc, zeros, levels, colors=colorst,extend='neither')   
map.drawcoastlines()
map.drawcountries()
map.drawstates()

cax = fig1.add_axes([0.02,0.07,0.96,0.02])
cbar = plt.colorbar(cs1, orientation='horizontal',\
    cax=cax, extend='both', ticks=levels[0:-1], format='%g') 
cbar.ax.tick_params(labelsize=9)
cbar.set_label('24-h accumulated precipitation amount (mm)')

# ---- set plot title

plot_title = 'forecast_precipitation_IC='+cdate+'_lead='+clead+'days_HUCs.png'
print ('saving plot to ', plot_title)
fig1.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')



      