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


infile = 'precip_budget_and_pressure_reanalysis_timeseries.cPick'
inf = open(infile, 'rb')
decimalyear = cPickle.load(inf)
#cPickle.dump(total_precip, ouf)
surface_pressure = cPickle.load(inf)
precipitation_rate = cPickle.load(inf)
evaporation_rate = cPickle.load(inf)
inf.close()

print ('mean precipitation rate, evaporation rate = ', \
    ma.mean(precipitation_rate), ma.mean(evaporation_rate) )
print ('max, min surface_pressure = ', np.max(surface_pressure), np.min(surface_pressure))

# ---- plot (a) dry pressure, (b) precipitation and evaporation, and (c) p-e
#
f = plt.figure(figsize=(6.5,9.))

ax = f.add_axes([.12,.72,.84,.24])
plt.title('(a) Global dry pressure (hPa)',fontsize=16)
ax.plot(decimalyear, surface_pressure, marker='o', color='Red',\
    lw=0.0, markersize=0.5, markerfacecolor='Red')
plt.ylabel('Dry pressure (hPa)',fontsize=13)
plt.grid(True, lw=0.15)
ax.set_xlim(2000, 2020)
ax.set_ylim(983., 983.7)
ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
ax.set_xlabel('Date', fontsize=13)

ax = f.add_axes([.12,.39,.84,.24])
plt.title(r'(b) Precipitation and evaporation rate ($\times\ 10^5$)',fontsize=16)
ax.plot(decimalyear, evaporation_rate*100000., marker='o', color='Red',\
    lw=0.00, markersize=0.5, markerfacecolor='Red',label='Evaporation')
ax.plot(decimalyear, precipitation_rate*100000., marker='o', color='RoyalBlue',\
    lw=0.00, markersize=0.5, markerfacecolor='RoyalBlue',label='Precipitation') # factor of two for bug in data.
plt.ylabel(r'Rate ($kg m^{-2} s^{-1}$)',fontsize=13)
plt.grid(True, lw=0.15)
ax.set_xlim(2000, 2020)
ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
ax.set_xlabel('Date', fontsize=13)
ax.set_ylim(2.5, 4.5)
ax.legend(loc='upper center')

ax = f.add_axes([.12,.06,.84,.24])
plt.title(r'(c) Precipitation minus evaporation rate ($\times\ 10^5$)',fontsize=16)
ax.plot(decimalyear, evaporation_rate*100000. - precipitation_rate*100000., marker='o', \
    lw=0.0, markersize=0.5, color='Red')
plt.ylabel(r'Rate ($kg m^{-2} s^{-1}$)',fontsize=13)
plt.grid(True, lw=0.15)
ax.set_xlim(2000, 2020)
ax.plot([2000,2020], [0,0], 'k-', lw=0.9)
ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
ax.set_xlabel('Date', fontsize=13)
ax.set_ylim(-1.0, 1.0)

imagefile = 'conserved_quantities_timeseries.png'
plt.savefig(imagefile, dpi=600)
print ('Plot done', imagefile)
plt.close()

