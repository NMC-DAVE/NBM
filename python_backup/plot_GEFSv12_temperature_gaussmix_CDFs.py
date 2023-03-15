"""
plot_GEFSv12_temperature_CDFs.py

"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
import numpy as np
import numpy.ma as ma
import scipy
import scipy.stats as stats
import matplotlib.pyplot as plt

# --------------------------------------------------------------   


def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx
    
# =====================================================================
    
# ---- various initialization

clead = sys.argv[1]
cmonth = sys.argv[2] 
clon = sys.argv[3]
clat = sys.argv[4]
rlon = float(clon)
rlat = float(clat)


if cmonth == 'Jan':
    imonth = 0
elif cmonth == 'Feb':
    imonth = 1
elif cmonth == 'Mar':
    imonth = 2
elif cmonth == 'Apr':
    imonth = 3
elif cmonth == 'May':
    imonth = 4
elif cmonth == 'Jun':
    imonth = 5
elif cmonth == 'Jul':
    imonth = 6
elif cmonth == 'Aug':
    imonth = 7
elif cmonth == 'Sep':
    imonth = 8
elif cmonth == 'Oct':
    imonth = 9
elif cmonth == 'Nov':
    imonth = 10
elif cmonth == 'Dec':
    imonth = 11

iskip = int(clead)//24
cvariable = '2t'
cpath_era5 = '/Volumes/Backup Plus/ecmwf/'
cpath_gefsv12 = '/Volumes/Backup Plus/gefsv12/t2m/'

date_end = dateshift('2018123100',-2*int(clead))
date_list_anal = daterange('2000010100',date_end,24)

ndates = len(date_list_anal)
date_list_forecast = []
for i in range(ndates):
    date_list_forecast.append(dateshift(date_list_anal[i], int(clead)))
    
date_list_anal_ver = []
for i in range(ndates):
    date_list_anal_ver.append(dateshift(date_list_forecast[i], int(clead)))

# ---- loop through dates and process day by day

ndatektr = 0
ndatektr_yearly = np.zeros((18), dtype=np.int32)
for idate, date in enumerate(date_list_anal):
    
    # ---- read reanalysis appropriate at the time.
    
    rem = idate%30
    if rem == 0: print ('processing date = ', date)
    cyear = date[0:4]
    cmm = date[4:6]
    cmmdd = date[4:8]
    imm = int(cmm)
    iyear = int(cyear)-2000
    datef = date_list_forecast[idate]
    cyearf = datef[0:4]
    
    infile = cpath_era5 +cyearf+'/t2m_era5_halfdegree_'+datef+'.cPick'
    fexist1 = os.path.exists(infile)
    if fexist1 == True:
        inf = open(infile, 'rb')
        analysis = cPickle.load(inf) - 273.16
        analysis = np.flipud(analysis)
        if idate == 0: 
            lats = cPickle.load(inf)
            lons = cPickle.load(inf)
            #print ('min, max lons = ', np.min(lons), np.max(lons))
            nlats, nlons = np.shape(lats)
            lats = np.flipud(lats)
            forecast_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            analysis_3d = np.zeros((ndates,nlats,nlons), dtype=np.float64)
            lons_1d = lons[0,:]
            lats_1d = lats[:,0]
            #print ('lats_1d = ', lats_1d)
            imin = find_nearest(lons_1d, rlon)
            jmin = find_nearest(lats_1d, rlat)
    else:
        print ('1. did not find file ', infile)
            
    # ---- read the forecast information .
            
    cyear = date[0:4]    
    cpath_forecast = cpath_gefsv12+cyear+'/'
    infile = cpath_forecast + date + '_lead'+clead+'_conus_0.5deg_hour'+clead+'.cPick'
    fexist2 = os.path.exists(infile)
    if fexist2 == True:
        inf = open(infile,'rb')
        forecast = cPickle.load(inf)
        inf.close()
        if cmmdd == '0101':
            infile = '/Volumes/Backup Plus/gefsv12/t2m/gefsv12_latlon_subset.cPick'
            inf = open(infile,'rb')
            latsf = cPickle.load(inf)
            lonsf = cPickle.load(inf)
            #print ('latsf[:,0] = ',latsf[:,0])
            #print ('min, max lonsf = ', np.min(lonsf), np.max(lonsf))
            inf.close()
    else:
        print ('2. did not find file ', infile)
        
    # ---- enter into 3D array 
    
    forecast_3d[idate,:,:] = forecast[:,:] 
    analysis_3d[idate,:,:] = analysis[:,:]
        
# ---- thin down the data to just the month of interest and surrounding months

datevalid_indices = np.zeros(ndates,dtype=np.int32)
for idate, date in enumerate(date_list_anal):
    thismonth = int(str(date)[4:6]) - 1
    if imonth == 0:
        if thismonth == imonth or thismonth == 11 or \
            thismonth == imonth+1: datevalid_indices[idate] = 1
    elif imonth == 11:
        if thismonth == imonth or thismonth == 0 or \
            thismonth == imonth-1: datevalid_indices[idate] = 1
    else:
        if thismonth == imonth or thismonth == imonth+1 or \
            thismonth == imonth-1: datevalid_indices[idate] = 1
                
ndates_valid = np.sum(datevalid_indices)
print ('ndates_valid, nlats, nlons', ndates_valid, nlats, nlons)
forecast_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float32)
analysis_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float32)

ktr = 0
for idate, date in enumerate(date_list_anal):
    if datevalid_indices[idate] == 1:
        forecast_validdates[ktr,:,:] = forecast_3d[idate,:,:]
        analysis_validdates[ktr,:,:] = analysis_3d[idate,:,:]
        ktr = ktr+1

print ('min, max forecast_validdates = ', np.min(forecast_validdates), np.max(forecast_validdates))
print ('min, max analysis_validdates = ', np.min(analysis_validdates), np.max(analysis_validdates))

# ---- load the information for generating CDFs, the climatological distributions 
    
infile = cpath_gefsv12+'GEFSv12_forecast_Gaussmix2_parameters_f'+clead+'.cPick'
inf = open(infile, 'rb')
weights_forecast = cPickle.load(inf)
means_forecast = cPickle.load(inf)
stddevs_forecast = cPickle.load(inf)
inf.close()

infile = cpath_era5+'ERA5_analyzed_Gaussmix2_parameters_f'+clead+'.cPick'
inf = open(infile, 'rb')
weights_analysis = cPickle.load(inf)
means_analysis = cPickle.load(inf)
stddevs_analysis = cPickle.load(inf)
inf.close()

# ---- determine the empirical CDF for this grid point  

print ('determining forecast empirical quantiles')
sample = forecast_validdates[:,jmin, imin]
forecast_sample_sorted = np.sort(sample)
nsamps = len(forecast_sample_sorted)
empirical_CDF_forecast = 0.5/float(nsamps) + np.arange(nsamps)/float(nsamps)

print ('determining analyzed empirical quantiles')
sample = analysis_validdates[:,jmin, imin]
analysis_sample_sorted = np.sort(sample)
nsamps = len(analysis_sample_sorted)
empirical_CDF_analysis = 0.5/float(nsamps) + np.arange(nsamps)/float(nsamps)

fmin = forecast_sample_sorted[0]
fmax = forecast_sample_sorted[-1]
amin = analysis_sample_sorted[0]
amax = analysis_sample_sorted[-1]
print ('fmin, fmax, amin, amax = ', fmin, fmax, amin, amax)

# --- using the pre-determined power transformation for this grid point and the sample mean,
#     develop a 

fweights = weights_forecast[imonth,:,jmin,imin]
fmeans = means_forecast[imonth,:,jmin,imin]
fstds = stddevs_forecast[imonth,:,jmin,imin]
print ('fweights, fmeans, fstds = ', fweights, fmeans, fstds )
#Zval_forecast1 = (temps_forecast - fmeans[0]) / fstds[0]
#Zval_forecast2 = (temps_forecast - fmeans[1]) / fstds[1]
#Zval_forecast3 = (temps_forecast - fmeans[2]) / fstds[2]
Zval_forecast1 = (forecast_sample_sorted - fmeans[0]) / fstds[0]
Zval_forecast2 = (forecast_sample_sorted - fmeans[1]) / fstds[1]
Zval_forecast3 = (forecast_sample_sorted - fmeans[2]) / fstds[2]
CDF_forecast1 = scipy.stats.norm.cdf(Zval_forecast1, loc=0., scale=1.)
CDF_forecast2 = scipy.stats.norm.cdf(Zval_forecast2, loc=0., scale=1.)
CDF_forecast3 = scipy.stats.norm.cdf(Zval_forecast3, loc=0., scale=1.)
CDF_forecast = fweights[0]*CDF_forecast1 + \
    fweights[1]*CDF_forecast2 + fweights[2]*CDF_forecast3



aweights = weights_analysis[imonth,:,jmin,imin]
ameans = means_analysis[imonth,:,jmin,imin]
astds = stddevs_analysis[imonth,:,jmin,imin]
print ('aweights, ameans, astds = ', aweights, ameans, astds )
#Zval_analysis1 = (temps_analysis - ameans[0]) / astds[0]
#Zval_analysis2 = (temps_analysis - ameans[1]) / astds[1]
#Zval_analysis3 = (temps_analysis - ameans[2]) / astds[2]
Zval_analysis1 = (analysis_sample_sorted - ameans[0]) / astds[0]
Zval_analysis2 = (analysis_sample_sorted - ameans[1]) / astds[1]
Zval_analysis3 = (analysis_sample_sorted - ameans[2]) / astds[2]
CDF_analysis1 = scipy.stats.norm.cdf(Zval_analysis1, loc=0., scale=1.)
CDF_analysis2 = scipy.stats.norm.cdf(Zval_analysis2, loc=0., scale=1.)
CDF_analysis3 = scipy.stats.norm.cdf(Zval_analysis3, loc=0., scale=1.)
CDF_analysis = aweights[0]*CDF_analysis1 + \
    aweights[1]*CDF_analysis2 + aweights[2]*CDF_analysis3

# ---- plot the CDFs, forecast and analyzed, empirical and best fitted. 

f = plt.figure(figsize=(4.5,9.))

ax = f.add_axes([.15,.57,.8,.38])
#ax.set_title('(a) Forecast CDFs for '+cmonth+', lon = '+clon+' lat = '+clat+' lead = '+clead+'h',fontsize=13)
ax.set_title('(a) Forecast CDFs',fontsize=16)
ax.plot(forecast_sample_sorted,empirical_CDF_forecast,':',color='Red',lw=2,label='Empirical')
ax.plot(forecast_sample_sorted,CDF_forecast,color='Blue',lw=2,label='Fitted')
plt.ylabel('Non-exceedance probability',fontsize=13)
ax.legend(loc=0)
ax.set_ylim(0,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(int(fmin)-1, int(fmax)+1)
ax.set_xlabel('Temperature (degrees C)',fontsize=13)

ax = f.add_axes([.15,.07,.8,.38])
#ax.set_title('(b) Analysis CDFs for '+cmonth+', lon = '+clon+' lat = '+clat+', lead = '+clead+'h',fontsize=13)
ax.set_title('(b) Analysis CDFs',fontsize=16)
ax.plot(analysis_sample_sorted,empirical_CDF_analysis,':',color='Red',lw=2,label='Empirical')
ax.plot(analysis_sample_sorted,CDF_analysis,color='Blue',lw=2,label='Fitted')
plt.ylabel('Non-exceedance probability',fontsize=13)
ax.legend(loc=0)
ax.set_ylim(0,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(int(fmin)-1, int(fmax)+1)
ax.set_xlabel('Temperature (degrees C)',fontsize=13)        

figname = 'CDF_'+cmonth+'_lon='+clon+'_lat='+clat+'_lead='+clead+'.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()


# ---- plot the CDFs, forecast and analyzed, empirical and best fitted. 

f = plt.figure(figsize=(6.5,4.5))

ax = f.add_axes([.11,.14,.88,.74])
ax.set_title('Forecast and analyzed CDFs for '+cmonth+',\nlon = '+clon+' lat = '+clat+', lead = '+clead+'h',fontsize=12)
ax.plot(forecast_sample_sorted,CDF_forecast,color='Blue',lw=2,label='Forecast')
ax.plot(analysis_sample_sorted,CDF_analysis,color='Red',lw=2,label='Analyzed')
plt.ylabel('Non-exceedance probability',fontsize=13)
ax.legend(loc=0)
ax.set_ylim(0,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(int(fmin)-1, int(fmax)+1)
ax.set_xlabel('Temperature (C)',fontsize=13)

figname = 'CDF_together_'+cmonth+'_lon='+clon+'_lat='+clat+'_lead='+clead+'.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()



# ---- P-P plot of forecast data

f = plt.figure(figsize=(4.5,9.))

ax = f.add_axes([.15,.56,.81,.38])
#ax.set_title('(a) Forecast Q-Q plot for '+cmonth+',\n lon = '+clon+' lat = '+clat+' lead = '+clead+'h',fontsize=13)
ax.set_title('(a) Forecast P-P plot',fontsize=13)
ax.plot([0,1],[0,1],color='Gray',lw=0.5)
ax.plot(empirical_CDF_forecast,CDF_forecast,color='Red',lw=2)
ax.set_ylim(0,1)
ax.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,1)
ax.set_xlabel('Empirical forecast quantile',fontsize=11)
ax.set_ylabel('Fitted forecast quantile',fontsize=11)

ax = f.add_axes([.15,.08,.81,.38])
#ax.set_title('(b) Analysis Q-Q plot for '+cmonth+',\n lon = '+clon+' lat = '+clat+', lead = '+clead+'h',fontsize=13)
ax.set_title('(b) Analysis P-P plot',fontsize=13)
ax.plot([0,1],[0,1],color='Gray',lw=0.5)
ax.plot(empirical_CDF_analysis,CDF_analysis,color='Red',lw=2,label='Empirical')
ax.set_ylim(0,1)
ax.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,1)
ax.set_xlabel('Empirical analyzed quantile',fontsize=11)   
ax.set_ylabel('Fitted analyzed quantile',fontsize=11)     

figname = 'PPplot_'+cmonth+'_lon='+clon+'_lat='+clat+'_lead='+clead+'.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()
    