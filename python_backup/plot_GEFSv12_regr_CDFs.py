"""
plot_GEFSv12_regr_CDFs.py

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
import statsmodels.api as sm

# --------------------------------------------------------------   


def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx
    
# =====================================================================
    
# ---- various initialization

clead = '24' # sys.argv[1]
cmonth = 'Feb' # 'Jul' # sys.argv[2] 
clon = -105.0 # sys.argv[3]
clat = 40.0 # sys.argv[4]
rlon = float(clon)
rlat = float(clat)

#infile = 'boulder_july_data_lead='+clead+'h.cPick'
infile = 'boulder_feb_data_lead='+clead+'h.cPick'
inf = open(infile,'rb')
analysis_bou = cPickle.load(inf)
forecast_bou = cPickle.load(inf)
ktr = len(analysis_bou)
inf.close() 
print ('min, max analysis_bou = ', np.min(analysis_bou), np.max(analysis_bou))
print ('min, max forecast_bou = ', np.min(forecast_bou), np.max(forecast_bou))

#infile = 'boulder_july_data_2019_lead='+clead+'h.cPick'
infile = 'boulder_feb_data_2019_lead='+clead+'h.cPick'
inf = open(infile,'rb')
analysis_bou_2019 = cPickle.load(inf)
forecast_bou_2019 = cPickle.load(inf)
inf.close() 
print ('min, max analysis_bou_2019 = ', np.min(analysis_bou_2019), np.max(analysis_bou_2019))
print ('min, max forecast_bou_2019 = ', np.min(forecast_bou_2019), np.max(forecast_bou_2019))

X = np.zeros((ktr,2), dtype=np.float64)
X[:,0] = 1.0
X[:,1] = forecast_bou[:]
y = np.zeros((ktr), dtype=np.float64)
y[:] = analysis_bou[:]
model = sm.OLS(y, X)
results = model.fit()
regr_coefs = results.params 
#x = np.arange(10.,40.5,0.5)
x = np.arange(-10.,20.5,0.5)
predicted = regr_coefs[0] + regr_coefs[1]*x

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

# ---- read climatology

if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
    infile = cpath_era5 + 'ERA5_temperature_climatology_12UTC.cPick'
else:
    infile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'

print ('reading from ', infile)
inf = open(infile,'rb')
climo_temps_estimated = cPickle.load(inf)
inf.close()

# ---- loop through dates and process day by day

ndatektr = 0
ndatektr_yearly = np.zeros((18), dtype=np.int32)
for idate, date in enumerate(date_list_anal):
    
    # ---- read reanalysis appropriate at the time.
    
    rem = idate%30
    #if rem == 0: print ('processing date = ', date)
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
            print ('jmin, imin = ',jmin, imin)
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
#print ('ndates_valid, nlats, nlons', ndates_valid, nlats, nlons)
forecast_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float32)
analysis_validdates = np.zeros((ndates_valid, nlats, nlons), dtype=np.float32)

ktr = 0
for idate, date in enumerate(date_list_anal):
    if datevalid_indices[idate] == 1:
        forecast_validdates[ktr,:,:] = forecast_3d[idate,:,:]
        analysis_validdates[ktr,:,:] = analysis_3d[idate,:,:]
        ktr = ktr+1

#print ('min, max forecast_validdates = ', np.min(forecast_validdates), np.max(forecast_validdates))
#print ('min, max analysis_validdates = ', np.min(analysis_validdates), np.max(analysis_validdates))

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
#print ('fmin, fmax, amin, amax = ', fmin, fmax, amin, amax)

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
PDF_forecast1 = scipy.stats.norm.pdf(Zval_forecast1, loc=0., scale=1.)
PDF_forecast2 = scipy.stats.norm.pdf(Zval_forecast2, loc=0., scale=1.)
PDF_forecast3 = scipy.stats.norm.pdf(Zval_forecast3, loc=0., scale=1.)
PDF_forecast = fweights[0]*PDF_forecast1 + \
    fweights[1]*PDF_forecast2 + fweights[2]*PDF_forecast3



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
PDF_analysis1 = scipy.stats.norm.pdf(Zval_analysis1, loc=0., scale=1.)
PDF_analysis2 = scipy.stats.norm.pdf(Zval_analysis2, loc=0., scale=1.)
PDF_analysis3 = scipy.stats.norm.pdf(Zval_analysis3, loc=0., scale=1.)
PDF_analysis = aweights[0]*PDF_analysis1 + \
    aweights[1]*PDF_analysis2 + aweights[2]*PDF_analysis3


# ---- plot the Boulder CDFs and regression data. 

f = plt.figure(figsize=(4.5,9.))

ax = f.add_axes([.16,.57,.8,.37])
ax.set_title('(a) CDFs',fontsize=15)
ax.plot(forecast_sample_sorted,CDF_forecast,color='Blue',lw=2,label='Forecast')
ax.plot(analysis_sample_sorted,CDF_analysis,'--',color='Red',lw=2,label='Analyzed')
ax.set_ylabel('Non-exceedance probability',fontsize=13)
ax.legend(loc=0)
ax.set_ylim(0,1)
plt.grid(True,lw=0.25,color='LightGray')
#ax.set_xlim(10,40)
ax.set_xlim(-10,20)
ax.set_xlabel('Temperature (C)',fontsize=13)

ax = f.add_axes([.16,.07,.8,.37])
ax.set_title('(b) Forecasts and analyses',fontsize=15)
print (forecast_bou[0:-1:10])
print (analysis_bou[0:-1:10])
ax.scatter(forecast_bou,analysis_bou,color='Gray',marker='o',s=0.1, label='2000-2018')
ax.scatter(forecast_bou_2019,analysis_bou_2019,color='Red',marker='o',s=2.5,label='2019')
ax.plot(x,predicted,'r-',lw=2,label='uMOS')
#ax.plot([10,40],[10,40],lw=1,color='Black')
ax.plot([-10,20],[-10,20],lw=1,color='Black')
#ax.plot(10.+PDF_analysis*15.0,analysis_sample_sorted,lw=1,color='Gray')
#ax.plot(forecast_sample_sorted,10.+PDF_forecast*15.0,lw=1,color='Gray')
ax.plot(-10.+PDF_analysis*15.0,analysis_sample_sorted,lw=1,color='Gray')
ax.plot(forecast_sample_sorted,-10.+PDF_forecast*15.0,lw=1,color='Gray')
print ('forecast_sample_sorted[0:-1:10] = ', forecast_sample_sorted[0:-1:10])
print ('analysis_sample_sorted[0:-1:10] = ', analysis_sample_sorted[0:-1:10])
print ('PDF_forecast[0:-1:10] = ', PDF_forecast[0:-1:10])
print ('PDF_analysis[0:-1:10] = ', PDF_analysis[0:-1:10])
ax.set_ylabel('Analyzed temperature (deg C)',fontsize=13)
#ax.set_xlim(10,40)
#ax.set_ylim(10,40)
ax.set_xlim(-10,20)
ax.set_ylim(-10,20)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlabel('Forecast temperature (deg C)',fontsize=13)
ax.legend(loc=0)

figname = 'boulder_24h_'+cmonth+'.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()



    