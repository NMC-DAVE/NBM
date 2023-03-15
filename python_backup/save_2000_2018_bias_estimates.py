"""
save_2019_bias_estimates.py

"""

import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
from datetime import datetime
import _pickle as cPickle
import numpy as np
import numpy.ma as ma

# --------------------------------------------------------------

def set_alpha(clead):
    if clead == '24':
        #alpha = 0.16  # will yield approx 0.08 alpha coefficient if R=Bx=1.0
        alpha = 0.1  # will yield approx 0.08 alpha coefficient if R=Bx=1.0
    elif clead == '48':
        alpha = 0.12  
    elif clead == '72':
        alpha = 0.08  
    elif clead == '96':
        alpha = 0.07  
    elif clead == '120':
        alpha = 0.06  

    return alpha

# --------------------------------------------------------------   

# ---- various initialization

clead = sys.argv[1]
alpha = set_alpha(clead)
iskip = int(clead)//24
cvariable = '2t'
rmse_raw = np.float64(0.0)
mae_raw = np.float64(0.0)
bia_raw = np.float64(0.0)
rmse_bcorr = np.float64(0.0)
mae_bcorr = np.float64(0.0)
bia_bcorr = np.float64(0.0)

rmse_raw_yearly = np.zeros((20), dtype=np.float64)
mae_raw_yearly = np.zeros((20), dtype=np.float64)
bia_raw_yearly = np.zeros((20), dtype=np.float64)
rmse_bcorr_yearly = np.zeros((20), dtype=np.float64)
mae_bcorr_yearly = np.zeros((20), dtype=np.float64)
bia_bcorr_yearly = np.zeros((20), dtype=np.float64)


cpath_era5 = '/Volumes/Backup Plus/ecmwf/'

date_end = dateshift('2017123100',-2*int(clead))

date_list_anal = daterange('2000010100',date_end,24)
#date_list_anal = daterange('2000010100','2000013100',24)

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
    
    # ---- read reanalysis appropriate at the time of this forecast for bias corr.
    
    cyear = date[0:4]
    iyear = int(cyear)-2000
    datef = date_list_forecast[idate]
    cyearf = datef[0:4]
    infile = cpath_era5 +cyearf+'/t2m_era5_halfdegree_'+datef+'.cPick'
    fexist1 = os.path.exists(infile)
    #print (infile, fexist1)
    if fexist1 == True:
        inf = open(infile, 'rb')
        analysis = cPickle.load(inf) - 273.16
        if idate == 0: 
            lats = cPickle.load(inf)
            lons = cPickle.load(inf)
            nlats, nlons = np.shape(lats)
            bias_3d = np.zeros((ndates,nlats,nlons), dtype=np.float32)
            random_3d = np.zeros((ndates,nlats,nlons), dtype=np.float32)
            #print ('analysis lats[:,0] = ', lats[:,0])
            #print ('analysis lons[0,:] = ', lons[0,:])
            #print ('analysis[0,0] = ', analysis[0,0])
    else:
        print ('1. did not find file ', infile)
            
    # ---- read the forecast information for bias corr.
            
    cyear = date[0:4]    
    cpath_forecast = '/Volumes/Backup Plus/gefsv12/t2m/'+cyear+'/'
    infile = cpath_forecast + date + '_lead'+clead+'_conus_0.5deg_hour'+clead+'.cPick'
    fexist2 = os.path.exists(infile)
    #print (infile, fexist2)
    if fexist2 == True:
        inf = open(infile,'rb')
        forecast = cPickle.load(inf)
        inf.close()
        if idate == 0:
            infile = '/Volumes/Backup Plus/gefsv12/t2m/gefsv12_latlon_subset.cPick'
            inf = open(infile,'rb')
            latsf = cPickle.load(inf)
            lonsf = cPickle.load(inf)
            inf.close()
            #print ('forecast latsf[:,0] = ', latsf[:,0])
            #print ('forecast lonsf[0,:] = ', lonsf[0,:])
            #print ('forecast[0,0] = ', forecast[0,0])
    else:
        print ('2. did not find file ', infile)
    
    # ---- perform decaying average bias correction
    
    if fexist1 == True and fexist2 == True:
        obsinc = forecast - analysis # - forecast
        if idate == 0:
            bias_3d[idate,:,:] = alpha*obsinc[:,:]
            random_3d[idate,:,:] = 0.0
        else:
            bias_3d[idate,:,:] = (1.-alpha)*bias_3d[idate-1,:,:] + \
                alpha*obsinc[:,:]
            random_3d[idate,:,:] = (forecast[:,:] - bias_3d[idate-1,:,:]) - analysis[:,:]
    else:
        if idate > 0: bias_3d[idate,:,:] = bias_3d[idate-1,:,:] 
    
    # ---- read reanalysis appropriate at the time of the verification
    
    datef = date_list_anal_ver[idate]
    cyearf = datef[0:4]
    infile = cpath_era5+cyearf+'/t2m_era5_halfdegree_'+datef+'.cPick'
    fexist2 = os.path.exists(infile)
    #print (infile, fexist2)
    if fexist2 == True:
        inf = open(infile, 'rb')
        analysis = cPickle.load(inf) - 273.16
        inf.close()
    else:
        print ('3. did not find file ', infile)
        
    # ---- read the forecast information at the time of the verification.
    
    cyearf = datef[0:4]        
    cpath_forecast = '/Volumes/Backup Plus/gefsv12/t2m/'+cyearf+'/'
    infile = cpath_forecast + datef + '_lead'+clead+\
        '_conus_0.5deg_hour'+clead+'.cPick'
    fexist1 = os.path.exists(infile)
    #print (infile, fexist1)
    if fexist1 == True:
        inf = open(infile,'rb')
        forecast = cPickle.load(inf)
        inf.close()
    else:
        print ('4. did not find file ', infile)
        
        
    meanf = np.mean(forecast)
    meana = np.mean(analysis)
    if meanf > 50. or meana > 50:
        print ('date, meanf, meana = ', date, meanf, meana)  

    if fexist1 == True and fexist2 == True:
        
        # ---- verify raw
    
        difference = analysis - forecast 
        d0 = np.mean(difference)
        if d0 > 10.0: print ('abs difference > 10 for ',idate,date)
        rmse_raw = rmse_raw + np.sum(difference**2)
        mae_raw = mae_raw + np.sum(np.abs(difference))
        bia_raw = bia_raw + np.sum(difference)
        rmse_raw_yearly[iyear] = rmse_raw_yearly[iyear]  + np.sum(difference**2)
        mae_raw_yearly[iyear]  = mae_raw_yearly[iyear]  + np.sum(np.abs(difference))
        bia_raw_yearly[iyear]  = bia_raw_yearly[iyear]  + np.sum(difference)
        
        if np.abs(np.sum(difference/(nlats*nlons))) > 4.0:
            print ('problematic data at ',datef)
            print (infile)
            print ('mae = ', mae)
            sys.exit()
        
        # ---- verify bias_corrected
    
        b0 = np.mean(bias_3d[idate,:,:])
        if b0 == 0.0:   print ('bias = 0 for ',idate,date)
        if np.abs(b0) > 10.0: print ('abs bias > 10 for ',idate,date)
        difference = analysis - (forecast - bias_3d[idate,:,:])
        rmse_bcorr = rmse_bcorr + np.sum(difference**2)
        mae_bcorr = mae_bcorr + np.sum(np.abs(difference))
        bia_bcorr = bia_bcorr + np.sum(difference)
        rmse_bcorr_yearly[iyear] = rmse_bcorr_yearly[iyear]  + np.sum(difference**2)
        mae_bcorr_yearly[iyear]  = mae_bcorr_yearly[iyear]  + np.sum(np.abs(difference))
        bia_bcorr_yearly[iyear]  = bia_bcorr_yearly[iyear]  + np.sum(difference)
        
        ndatektr = ndatektr + 1
        ndatektr_yearly[iyear] = ndatektr_yearly[iyear] + 1
        
    else:
        print ('missing date', date)
        

# ---- compute the final stats and print 

nsamps = ndatektr*nlats*nlons
bia_raw = bia_raw / float(nsamps)
mae_raw = mae_raw / float(nsamps)
rmse_raw = np.sqrt( rmse_raw / float(nsamps))

bia_bcorr = bia_bcorr / float(nsamps)
mae_bcorr = mae_bcorr / float(nsamps)
rmse_bcorr = np.sqrt( rmse_bcorr / float(nsamps))

print (ndatektr, ndatektr/365.)
print ('raw rmse, bias, mae = ', rmse_raw, bia_raw, mae_raw)
print ('bias-corrected rmse, bias, mae = ', rmse_bcorr, bia_bcorr, mae_bcorr)

for iyear in range(18):
    nsamps = ndatektr_yearly[iyear]*nlats*nlons
    bia_raw_yearly[iyear] = bia_raw_yearly[iyear] / float(nsamps)
    mae_raw_yearly[iyear] = mae_raw_yearly[iyear] / float(nsamps)
    rmse_raw_yearly[iyear] = np.sqrt( rmse_raw_yearly[iyear] / float(nsamps))

    bia_bcorr_yearly[iyear] = bia_bcorr_yearly[iyear] / float(nsamps)
    mae_bcorr_yearly[iyear] = mae_bcorr_yearly[iyear] / float(nsamps)
    rmse_bcorr_yearly[iyear] = np.sqrt( rmse_bcorr_yearly[iyear] / float(nsamps))
    
    print (iyear+2000, bia_raw_yearly[iyear], mae_raw_yearly[iyear], \
        rmse_raw_yearly[iyear],  bia_bcorr_yearly[iyear],\
        mae_bcorr_yearly[iyear], rmse_bcorr_yearly[iyear])

# ---- save the bias correction

bias_file = 'bias_correction_decayavg_lead'+clead+'.cPick'
ouf = open(bias_file, 'wb')
cPickle.dump(bias_3d, ouf)
cPickle.dump(date_list_anal, ouf)
ouf.close()

# ---- save the random data

bias_file = 'random_error_lead'+clead+'.cPick'
ouf = open(bias_file, 'wb')
cPickle.dump(random_3d, ouf)
cPickle.dump(date_list_anal, ouf)
ouf.close()
    






    
    

    
    
    