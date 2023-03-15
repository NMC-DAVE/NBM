"""
plot_bias_timeseries_v2.py

"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from read_reanalysis_timeseries import read_reanalysis_timeseries
from read_forecast_timeseries_GEFSv12 import read_forecast_timeseries_GEFSv12

# --------------------------------------------------------------   

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx

# --------------------------------------------------------------  

def read_bias_corrections(infile):
    print ('   reading beta estimates from ', infile)
    inf = open(infile,'rb')
    beta_3d = cPickle.load(inf)  
    lats = cPickle.load(inf) 
    lons = cPickle.load(inf) 
    inf.close()
    return beta_3d

# --------------------------------------------------------------  

def initialize_date_lists(cyear, clead):

    start_date = cyear+'010100'
    end_date = cyear+'123100'
    date_list_anal = daterange(start_date, end_date, 24)
    ndates = len(date_list_anal)
    date_list_forecast = []
    for i in range(ndates):
        date_list_forecast.append(dateshift(date_list_anal[i], \
            int(clead)))
    return date_list_anal, date_list_forecast    

# --------------------------------------------------------------  


    

# ---- various initialization

clead = sys.argv[1]
clon = sys.argv[2]
clat = sys.argv[3]
rlon = float(clon)
rlat = float(clat)
iskip = int(clead)//24
cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_forecast = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_beta = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_random = '/Volumes/Backup Plus/gefsv12/t2m/'
cpath_gain = '/Volumes/Backup Plus/gefsv12/t2m/'

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print (current_time,'   &&&&&&&&&& LEAD TIME = ',clead,'  &&&&&&&&&&')
    
cyear = '2019'
date_list_anal, date_list_forecast = initialize_date_lists(cyear, clead)
ndates = len(date_list_forecast)
print (date_list_forecast)




# - get lat/lon indices

infile = cpath_era5 +'2001/t2m_era5_halfdegree_2001010100.cPick'
inf = open(infile, 'rb')
analysis = cPickle.load(inf) - 273.16
analysis = np.flipud(analysis)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
nlats, nlons = np.shape(lats)
lats = np.flipud(lats)
lons_1d = lons[0,:]
lats_1d = lats[:,0]
imin = find_nearest(lons_1d, rlon)
jmin = find_nearest(lats_1d, rlat)
inf.close()

print ('lons_1d = ',lons_1d)
print ('lats_1d = ',lons_1d)
imin = find_nearest(lons_1d, rlon)
jmin = find_nearest(lats_1d, rlat)
print ('jmin, imin = ', jmin, imin)


# ---- read climatology

if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
    infile = cpath_era5 + 'ERA5_temperature_climatology_12UTC.cPick'
else:
    infile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'

print ('reading from ', infile)
inf = open(infile,'rb')
climo_temps_estimated = cPickle.load(inf)
inf.close()
analysis_climatology_1d = climo_temps_estimated[:,jmin, imin]
print (analysis_climatology_1d[180:211])

#print (date_list_forecast)
analyses_3d, lats, lons = read_reanalysis_timeseries(cpath_era5, \
    date_list_forecast)
analyses_1d = analyses_3d[:,jmin, imin]
ndates, nlats, nlons = np.shape(analyses_3d)
print (analyses_1d[180:211])
    
# ---- read the forecast time series on the dates specified.   Pass in 
#      the lead time and the initial time of the analysis.

#print (date_list_anal)
forecast_3d, latsf, lonsf = read_forecast_timeseries_GEFSv12 ( \
    cpath_forecast, date_list_anal, clead)
forecast_1d = forecast_3d[:,jmin, imin]
print (forecast_1d[180:211])


infile = cpath_forecast + '2019_qmap_lead'+clead+'.cPick'
beta_qmap = read_bias_corrections(infile)
beta_qmap_1d = beta_qmap[:,jmin,imin]
r = np.corrcoef(beta_qmap_1d[0:-1],beta_qmap_1d[1:] )
cr_qmap = r"$r_1$ = {:.2f}".format(r[0,1])
print (cr_qmap)

infile = cpath_forecast + '2019_decayavg_lead'+clead+'.cPick'
beta_decayavg = read_bias_corrections(infile)
beta_decayavg_1d = beta_decayavg[:,jmin,imin]
r = np.corrcoef(beta_decayavg_1d[0:-1],beta_decayavg_1d[1:] )
cr_decay = r"$r_1$ = {:.2f}".format(r[0,1])

infile = cpath_forecast + '2019_MOS_lead'+clead+'.cPick'
beta_uMOS = read_bias_corrections(infile) 
beta_uMOS_1d = beta_uMOS[:,jmin,imin]
r = np.corrcoef(beta_uMOS_1d[0:-1],beta_uMOS_1d[1:] )
cr_uMOS = r"$r_1$ = {:.2f}".format(r[0,1])

infile = cpath_forecast + '2019_MOS_mvr_lead'+clead+'.cPick'
beta_mvMOS = read_bias_corrections(infile) 
beta_mvMOS_1d = beta_mvMOS[:,jmin,imin]
r = np.corrcoef(beta_mvMOS_1d[0:-1],beta_mvMOS_1d[1:] )
cr_mvMOS = r"$r_1$ = {:.2f}".format(r[0,1])
    
# ---- loop thru months

cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

f = plt.figure(figsize=(9.,6.5))

ax = f.add_axes([.08,.57,.9,.36])
ax.set_title('(a) 00 UTC forecasts and analyses',fontsize=15)
ax.plot(range(ndates),analyses_1d,'o-',color='Red',\
    lw=1, markersize=1.5, label='ERA5 analysis', markerfacecolor='Red')
ax.plot(range(ndates),analysis_climatology_1d,'-',color='Peru',\
    lw=2, markersize=1.5, label='ERA5 climatology')
ax.plot(range(ndates),forecast_1d,'o-',color='RoyalBlue',\
    lw=1, markersize=1.5, label='Forecast', markerfacecolor='RoyalBlue')
ax.set_ylim(-20,40)
ax.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,365)
ax.set_xticks([0,30,58,89,119,150,180, 211,242,272,303,333])
ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
ax.set_xlabel('Julian day number',fontsize=11)
ax.set_ylabel('Temperature (deg C)',fontsize=11)
ax.legend(loc=0)

ax = f.add_axes([.08,.07,.9,.36])
ax.set_title('(b) bias estimates',fontsize=15)
ax.plot([0,365],[0,0],lw=1,color='Black')
ax.plot(range(ndates),beta_decayavg_1d,'o-',color='Peru',\
    lw=0.5, markersize=1.5, label='DAV, '+cr_decay, markerfacecolor='Peru')
ax.plot(range(ndates),beta_qmap_1d,'o-',color='Blue',\
    lw=0.5, markersize=1.5, label='QM, '+cr_qmap, markerfacecolor='Blue')    
ax.plot(range(ndates),beta_uMOS_1d,'o-',color='Green',\
    lw=0.5, markersize=1.5, label='uMOS, '+cr_uMOS, markerfacecolor='Green')
ax.plot(range(ndates),beta_mvMOS_1d,'o-',color='Gray',\
    lw=0.5, markersize=1.5, label='mvMOS, '+cr_mvMOS, markerfacecolor='Gray') 
ax.set_xticks([0,30,58, 89,119,150, 180,211,242, 272,303,333])
ax.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
ax.legend(loc=0)  
    
ax.set_xlim(0,365)
ax.grid(True,lw=0.25,color='LightGray')
ax.set_ylim(-3,2)
ax.set_xlabel('Julian day number',fontsize=11)
ax.set_ylabel('Bias estimate (deg C)',fontsize=11)


figname = 'bias_tseries_lon='+clon+'_lat='+clat+'_lead='+clead+'.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()
    

    
    
    