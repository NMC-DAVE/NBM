"""
plot_conserved_quantities.py 
plot the dry pressure, the evaporation, the precipitation rate
""" 


from netCDF4 import Dataset
import numpy as np
from dateutils import daterange, datetohrs, dayofyear
import sys
import os
import os.path
from os import path
import numpy.ma as ma
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

date_list = daterange('2000010100','2019123100',24)
ndates = len(date_list)

# ---- read in from Sherrie's file the date and dry pressure.

a = np.loadtxt('reanalysis/surface_all.txt')
yyyymmddhh = a[:,0]
drypressure = a[:,2]
drypressure_masked = ma.masked_where(drypressure <= 0.0, drypressure)

for i in range(ndates):
    x = int(date_list[i])
    #if x != yyyymmddhh[i]:
    print (i, x, date_list[i], drypressure[i])
sys.exit()

# ---- read the evaporation and the precipitation rate.


infile = 'evap_prate.cPick'
print ('reading from ', infile)
inf = open(infile,'rb')
evap_timeseries = cPickle.load(inf)
prate_timeseries = cPickle.load(inf)
print (evap_timeseries[0:-1:100])
print (prate_timeseries[0:-1:100])

inf.close()
evap_timeseries_masked = ma.masked_where(evap_timeseries < -99., evap_timeseries)
prate_timeseries_masked = ma.masked_where(prate_timeseries <= -99., prate_timeseries)

print (evap_timeseries_masked[0:-1:100])
print (prate_timeseries_masked[0:-1:100])

# ---- convert the dates to decimal years.

decimalyear = np.zeros((ndates), dtype=np.float32)
hourssince2000 = np.zeros((ndates), dtype=np.float32)
hourssince1CE_2000 = datetohrs('2000010100')
datetime_list = []
for i in range(ndates):
    
    hourssince2000[i] = datetohrs(date_list[i]) - hourssince1CE_2000
    decimalyear[i] = 2000. + hourssince2000[i]/(24.*365.25)
#
# ---- plot (a) dry pressure, (b) precipitation and evaporation, and (c) p-e
#
f = plt.figure(figsize=(6.5,9.))

ax = f.add_axes([.12,.72,.84,.24])
plt.title('(a) Global dry pressure (hPa)',fontsize=16)
ax.plot(decimalyear, drypressure_masked, marker='o', color='Red',\
    lw=0.0, markersize=0.3, markerfacecolor='Red')
plt.ylabel('Dry pressure (hPa)',fontsize=13)
plt.grid(True, lw=0.15)
ax.set_xlim(2000, 2020)
ax.set_ylim(983., 984.)
ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
ax.set_xlabel('Date', fontsize=13)

ax = f.add_axes([.12,.39,.84,.24])
plt.title(r'(b) Precipitation and evaporation rate ($\times\ 10^5$)',fontsize=16)
ax.plot(decimalyear, evap_timeseries_masked*100000., marker='o', color='Red',\
    lw=0.00, markersize=0.3, markerfacecolor='Red',label='Evaporation')
ax.plot(decimalyear, prate_timeseries_masked*100000., marker='o', color='RoyalBlue',\
    lw=0.00, markersize=0.3, markerfacecolor='RoyalBlue',label='Precipitation') # factor of two for bug in data.
plt.ylabel(r'Rate ($kg m^{-2} s^{-1}$)',fontsize=13)
plt.grid(True, lw=0.15)
ax.set_xlim(2000, 2020)
ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
ax.set_xlabel('Date', fontsize=13)
ax.set_ylim(3.0, 5.0)
ax.legend(loc='upper center')

ax = f.add_axes([.12,.06,.84,.24])
plt.title(r'(c) Precipitation minus evaporation rate ($\times\ 10^5$)',fontsize=16)
ax.plot(decimalyear, evap_timeseries*100000. - prate_timeseries*100000., marker='o', color='Red',\
    lw=0.0, markersize=0.3, markerfacecolor='Red')
plt.ylabel(r'Rate ($kg m^{-2} s^{-1}$)',fontsize=13)
plt.grid(True, lw=0.15)
ax.set_xlim(2000, 2020)
ax.plot([2000,2020], [0,0], 'k-', lw=0.9)
ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
ax.set_xlabel('Date', fontsize=13)
ax.set_ylim(-0.6, 0.6)

imagefile = 'conserved_quantities_timeseries.png'
plt.savefig(imagefile, dpi=600)
print ('Plot done', imagefile)
plt.close()

