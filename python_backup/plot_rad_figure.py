#!/usr/local/bin/python2.7
import numpy as np
import numpy.ma as ma
import dateutils
from netCDF4 import Dataset
from datetime import datetime, date, time, timedelta
import matplotlib as mpl
mpl.use('Agg') #for web 
from dateutils import datetohrs
import matplotlib.pyplot as plt
import sys
import scipy.signal as signal
from astropy.convolution import convolve
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.interpolate import UnivariateSpline
from sklearn.linear_model import LinearRegression


# ---- output directory for images
imagedir='/Volumes/Backup Plus/python/'

# input directory with RAD_*.nc files
datapath='/Volumes/Backup Plus/gefsv12/obs/'

if len(sys.argv) < 2:
    raise SystemExit('python plot_rad_figure.py <begindate> <enddate>, eg. python plot_rad_figure.py 2000010100 2019123118')
begindate=sys.argv[1]
enddate=sys.argv[2]

# always use GEFS & CFSR, global for amsua_n15

modelstreams = ['GEFS', 'CFSR']
region='GLOBL'
instrmnt='amsua'
satlite='n15'
channels = [ 1, 2, 3, 4, 5, 6, 7, 8,  9, 10, 11, 12, 13, 15]
#channels = [7]
for channel in channels:
    
    imagefile = imagedir+begindate+'_'+enddate+'_'+instrmnt+'_'+\
        satlite+'_ch'+str(channel)+'.png'
    channel=int(channel)
    enddate=str(enddate)
    begindate=str(begindate)
    nummodel=len(modelstreams)
    beginyr=str(begindate)[0:4]
    endyr=str(enddate)[0:4]
    years = range(int(beginyr), int(endyr)+1)
    strdates = dateutils.daterange(begindate,enddate,300)
    decimalyear = np.zeros(len(strdates), dtype=np.float32)
    dates = np.zeros(len(strdates))
    hourssince2000 = np.zeros(len(strdates), dtype=np.float32)
    hourssince1CE_2000 = datetohrs('2000010100')
    datetime_list = []
    for i in range(len(dates)):
        dates[i] = int(strdates[i])
        hourssince2000[i] = datetohrs(strdates[i]) - hourssince1CE_2000
        datetime_list.append(datetime.strptime(strdates[i], "%Y%m%d%H"))
        decimalyear[i] = 2000. + hourssince2000[i]/(24.*365.25)
    timediff = datetime_list[len(dates)-1] - datetime_list[0]

    channum_str = str(channel)
    nobs_all  = ma.zeros((len(dates), nummodel));      nobs_all.mask  = True
    nobs_qcd  = ma.zeros((len(dates), nummodel));      nobs_qcd.mask  = True
    nobs_used = ma.zeros((len(dates), nummodel));      nobs_used.mask = True
    mean_obs_all  = ma.zeros((len(dates), nummodel));  mean_obs_all.mask  = True
    mean_obs_used = ma.zeros((len(dates), nummodel));  mean_obs_used.mask = True
    mean_obs_qcd  = ma.zeros((len(dates), nummodel));  mean_obs_qcd.mask  = True
    mean_omf_ctrl = ma.zeros((len(dates), nummodel));  mean_omf_ctrl.mask = True
    mean_oma_ctrl = ma.zeros((len(dates), nummodel));  mean_oma_ctrl.mask = True
    std_omf_ctrl  = ma.zeros((len(dates), nummodel));  std_omf_ctrl.mask  = True
    std_oma_ctrl  = ma.zeros((len(dates), nummodel));  std_oma_ctrl.mask  = True
    mean_omf_ens  = ma.zeros((len(dates), nummodel));  mean_omf_ens.mask  = True
    spread_f      = ma.zeros((len(dates), nummodel));  spread_f.mask      = True
    spread_obserr_f = ma.zeros((len(dates), nummodel));  spread_obserr_f.mask = True
    std_omf_ens   = ma.zeros((len(dates), nummodel));  std_omf_ens.mask   = True
    mean_biascor  = ma.zeros((len(dates), nummodel));  mean_biascor.mask  = True
    std_biascor   = ma.zeros((len(dates), nummodel));  std_biascor.mask   = True
    X = np.ones((len(dates), 4), dtype=np.float32)
    Y_gefsv12 = np.ones((len(dates)), dtype=np.float32)
    Y_cfsr = np.ones((len(dates)), dtype=np.float32)
    yfit_gefsv12 = ma.ones((len(dates)), dtype=np.float32)
    yfit_cfsr = ma.ones((len(dates)), dtype=np.float32)
    Y_gefsv12_rmse = ma.ones((len(dates)), dtype=np.float32)
    Y_cfsr_rmse = ma.ones((len(dates)), dtype=np.float32)
    Y_gefsv12_bias = ma.ones((len(dates)), dtype=np.float32)
    Y_cfsr_bias = ma.ones((len(dates)), dtype=np.float32)
    
    for modct in range(nummodel):
        for year in years:
    
            modelfile = datapath + '/RAD_'+modelstreams[modct] + '_' +\
                str(year) + '_' + instrmnt + '_' + satlite + '_' + region + '.nc'
            anndata  = Dataset(modelfile, 'r')
            chans = anndata['Channels'][:].tolist()
            chanindx =chans.index(channel)
            thisdates = anndata['All_Dates'][:]
            indx_in = np.where(np.in1d(thisdates, dates)) [0]
            indx_out = np.where(np.in1d(dates, thisdates)) [0]
            nobs_all[indx_out,modct]  = anndata['nobs_all'][indx_in,chanindx]
            nobs_qcd[indx_out,modct]  = anndata['nobs_qcd'][indx_in,chanindx]
            nobs_used[indx_out,modct] = anndata['nobs_used'][indx_in,chanindx]
            mean_obs_all[indx_out,modct]  = anndata['mean_obs_all'][indx_in,chanindx]
            mean_obs_used[indx_out,modct] = anndata['mean_obs_used'][indx_in,chanindx]
            mean_obs_qcd[indx_out,modct]  = anndata['mean_obs_qcd'][indx_in,chanindx]
            mean_omf_ctrl[indx_out,modct] = anndata['mean_omf_ctrl'][indx_in,chanindx]
            mean_oma_ctrl[indx_out,modct] = anndata['mean_oma_ctrl'][indx_in,chanindx]
            std_omf_ctrl[indx_out,modct]  = anndata['std_omf_ctrl'][indx_in,chanindx]
            std_oma_ctrl[indx_out,modct]  = anndata['std_oma_ctrl'][indx_in,chanindx]
            mean_omf_ens[indx_out,modct]  = anndata['mean_omf_ens'][indx_in,chanindx]
            spread_f[indx_out,modct]      = anndata['spread_f'][indx_in,chanindx]
            spread_obserr_f[indx_out,modct] = anndata['spread_obserr_f'][indx_in,chanindx]
            std_omf_ens[indx_out,modct]   = anndata['std_omf_ens'][indx_in,chanindx]
            mean_biascor[indx_out,modct]  = anndata['mean_biascor'][indx_in,chanindx]
            std_biascor[indx_out,modct]   = anndata['std_biascor'][indx_in,chanindx]

            anndata.close()
            weight_gefsv12 = np.ones((len(dates)), dtype=np.float32)
            weight_cfsr = np.ones((len(dates)), dtype=np.float32)


    #print ('ndates = ', len(dates))
    for i in range(len(dates)):
        
        #print ('i, date, O-F(GEFS), O-F(CFSR) = ',i, dates[i],\
        #    std_omf_ctrl[i,0],std_omf_ctrl[i,1])
            
        if std_omf_ctrl.mask[i,0] == True or std_omf_ctrl[i,0] != std_omf_ctrl[i,0]: 
            weight_gefsv12[i] = 0.0001 
            Y_gefsv12[i] = 0.2
        else:
            Y_gefsv12[i] = std_omf_ctrl[i,0]
        if std_omf_ctrl.mask[i,1] == True or std_omf_ctrl[i,1] != std_omf_ctrl[i,1]: 
            weight_cfsr[i] = 0.0001
            Y_cfsr[i] = 0.2
        else:    
            Y_cfsr[i] = std_omf_ctrl[i,1]
        X[i,0] = decimalyear[i]
        f = decimalyear[i] - int(decimalyear[i])
        X[i,1] = np.cos(2.*3.1415926*f)
        X[i,2] = np.sin(2.*3.1415926*f)
        X[i,3] = decimalyear[i]**2
        #print ('YG, YC, X0, X1, X2, X3, wg, wc = ', Y_gefsv12[i], Y_cfsr[i], \
        #    X[i,0], X[i,1], X[i,2], X[i,3], weight_gefsv12[i], weight_cfsr[i])

    rscale = 10.0
    for i in range(len(dates)):
    #for i in [100]:
        sumwt1 = 0.0
        sumprod1 = 0.0
        sumwt2 = 0.0
        sumprod2 = 0.0
        for j in range(len(dates)):
            idiff = np.float(np.abs(i-j))
            wt = np.exp(-idiff**2/rscale**2)
            #print (i,j,idiff,wt)
            if mean_omf_ctrl.mask[i,0] == False and mean_omf_ctrl[j,0] == mean_omf_ctrl[j,0]: 
                sumwt1 = sumwt1 + wt
                sumprod1 = sumprod1 + wt*mean_omf_ctrl[j,0]
            if mean_omf_ctrl.mask[i,1] == False and mean_omf_ctrl[j,1] == mean_omf_ctrl[j,1]: 
                sumwt2 = sumwt2 + wt
                sumprod2 = sumprod2 + wt*mean_omf_ctrl[j,1]
        if sumwt1 > 5:
            Y_gefsv12_bias[i] = sumprod1 / sumwt1
        else:
            Y_gefsv12_bias[i] = 0.0
            Y_gefsv12_bias[i] = ma.masked
            
        if sumwt2 > 5:
            Y_cfsr_bias[i] = sumprod2 / sumwt2
        else:
            Y_cfsr_bias[i] = 0.0
            Y_cfsr_bias[i] = ma.masked


    for i in range(len(dates)):
        sumwt1 = 0.0
        sumprod1 = 0.0
        sumwt2 = 0.0
        sumprod2 = 0.0
        for j in range(len(dates)):
            idiff = np.float(np.abs(i-j))
            wt = np.exp(-idiff**2/rscale**2)
            #print (i,j,idiff,wt)
            if std_omf_ctrl.mask[i,0] == False and std_omf_ctrl[j,0] == std_omf_ctrl[j,0]: 
                sumwt1 = sumwt1 + wt
                sumprod1 = sumprod1 + wt*std_omf_ctrl[j,0]
            if std_omf_ctrl.mask[i,1] == False and std_omf_ctrl[j,1] == std_omf_ctrl[j,1]: 
                sumwt2 = sumwt2 + wt
                sumprod2 = sumprod2 + wt*std_omf_ctrl[j,1]
        if sumwt1 > 5:
            Y_gefsv12_rmse[i] = sumprod1 / sumwt1
        else:
            Y_gefsv12_rmse[i] = 0.0
            Y_gefsv12_rmse[i] = ma.masked
            
        if sumwt2 > 5:
            Y_cfsr_rmse[i] = sumprod2 / sumwt2
        else:
            Y_cfsr_rmse[i] = 0.0
            Y_cfsr_rmse[i] = ma.masked


    # ---- linear regression
   
    reg = LinearRegression().fit(X, Y_gefsv12)
    yfit_gefsv12 = reg.predict(X)
    reg = LinearRegression().fit(X, Y_cfsr)
    yfit_cfsr = reg.predict(X)   
   
    # ---- plot data
    
    f = plt.figure(figsize=(10.,8.))

    ax = f.add_axes([.09,.57,.88,.38])
    plt.title('(a) NOAA 15 AMSU-A channel '+str(channel)+' background RMSE',fontsize=16)
    ax.plot(decimalyear, std_omf_ctrl[:,0], marker='o', color='Red',\
        lw=0, markersize=0.6, label='GEFSv12 reanalysis', markerfacecolor='Red')
    ax.plot(decimalyear, std_omf_ctrl[:,1], marker='o', color='Blue', \
        markersize=0.6, lw=0, label='CFSR', markerfacecolor='Blue')
    ax.plot(decimalyear[0:-1], Y_gefsv12_rmse[0:-1], color='Red', lw=1)
    ax.plot(decimalyear[0:-60], Y_cfsr_rmse[0:-60] , color='Blue', lw=1)
    #ax.plot(decimalyear[0:-1], yfit_gefsv12[0:-1], color='Red', lw=1)
    #ax.plot(decimalyear[0:-60], yfit_cfsr[0:-60] , color='Blue', lw=1)
    plt.ylabel('Observed - forecast RMSE (deg C)',fontsize=13)
    ax.legend(loc='upper center')
    plt.grid(True, lw=0.25)
    ax.set_xlim(2000, 2020)
    ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
    ax.set_xlabel('Date', fontsize=13)
    
    ax = f.add_axes([.09,.07,.88,.38])
    plt.title('(b) NOAA 15 AMSU-A channel '+str(channel)+' background bias',fontsize=16)
    ax.plot(decimalyear, -mean_omf_ctrl[:,0], marker='o', color='Red',\
        lw=0, markersize=0.6, label='GEFSv12 reanalysis', markerfacecolor='Red')
    ax.plot(decimalyear, -mean_omf_ctrl[:,1], marker='o', color='Blue', \
        markersize=0.6, lw=0, label='CFSR', markerfacecolor='Blue')
    ax.plot(decimalyear[0:-1], -Y_gefsv12_bias[0:-1], color='Red', lw=1)
    ax.plot(decimalyear[0:-60], -Y_cfsr_bias[0:-60] , color='Blue', lw=1)
    #ax.plot(decimalyear[0:-1], yfit_gefsv12[0:-1], color='Red', lw=1)
    #ax.plot(decimalyear[0:-60], yfit_cfsr[0:-60] , color='Blue', lw=1)
    plt.ylabel('Forecast - observed bias (deg C)',fontsize=13)
    ax.legend(loc='upper center')
    plt.grid(True, lw=0.25)
    ax.set_xlim(2000, 2020)
    ax.plot([2000,2020],[0,0],lw=1.5,color='LightGray')
    ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
    ax.set_xlabel('Date', fontsize=13)
    
    
    plt.savefig(imagefile, dpi=400)
    print ('Plot done', imagefile)
    plt.close()