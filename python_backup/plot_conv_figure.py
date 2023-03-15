#!/usr/local/bin/python2.7
import numpy as np
import numpy.ma as ma
import dateutils
from netCDF4 import Dataset
from datetime import datetime, date, time
from dateutils import datetohrs
import matplotlib as mpl
mpl.use('Agg') #for web 
import matplotlib.pyplot as plt
import sys
import scipy.signal as signal
from astropy.convolution import convolve
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.interpolate import UnivariateSpline
from sklearn.linear_model import LinearRegression

# output path to images
imagedir='/Volumes/Backup Plus/python/'

# input path with CONV_*.nc data
datapath='/Volumes/Backup Plus/gefsv12/obs/'

if len(sys.argv) < 2:
    raise SystemExit('python plot_conv_figure.py <begindate> <enddate> eg. python plot_conv_figure.py 2000010100 2019123118')
begindate=sys.argv[1]
enddate=sys.argv[2]

# always use GEFS + CFSR; global, temperature
modelstreams = ['GEFS', 'CFSR']
region='GLOBL'
varnames=['t']
plevels=range(0,1000,100)

print ('processing from ', begindate, ' to ', enddate)

for var in varnames:
    for plevel in plevels:
        
        figname = imagedir+begindate+'_'+enddate+'_'+var+'_'+str(plevel)+'.png'

        nummodel=len(modelstreams)
        enddate=str(enddate)
        begindate=str(begindate)


        # ---- first and last year
         
        beginyr = str(begindate)[0:4]
        endyr = str(enddate)[0:4]
        years = range(int(beginyr), int(endyr)+1)

        # ---- all dates (last parameter -- step between datapoints in hours)

        strdates = dateutils.daterange(begindate,enddate,300)
        dates = np.zeros(len(strdates))
        hourssince2000 = np.zeros(len(strdates), dtype=np.float32)
        hourssince1CE_2000 = datetohrs('2000010100')
        decimalyear = np.zeros(len(strdates), dtype=np.float32)
        datetime_list = []
        for i in range(len(dates)):
            dates[i] = int(strdates[i])
            datetime_list.append(datetime.strptime(strdates[i], "%Y%m%d%H"))
            hourssince2000[i] = datetohrs(strdates[i]) - hourssince1CE_2000
            decimalyear[i] = 2000. + hourssince2000[i]/(24.*365.25)
            
        timediff = datetime_list[len(dates)-1] - datetime_list[0]

        # ---- declare arrays
         
        nobs_all  = ma.zeros((len(dates), nummodel));      nobs_all.mask  = True
        nobs_used = ma.zeros((len(dates), nummodel));      nobs_used.mask = True
        mean_obs_all  = ma.zeros((len(dates), nummodel));  mean_obs_all.mask  = True
        mean_obs_used = ma.zeros((len(dates), nummodel));  mean_obs_used.mask = True
        mean_omf_ctrl = ma.zeros((len(dates), nummodel));  mean_omf_ctrl.mask = True
        mean_oma_ctrl = ma.zeros((len(dates), nummodel));  mean_oma_ctrl.mask = True
        std_omf_ctrl  = ma.zeros((len(dates), nummodel));  std_omf_ctrl.mask  = True
        std_oma_ctrl  = ma.zeros((len(dates), nummodel));  std_oma_ctrl.mask  = True
        mean_omf_ens  = ma.zeros((len(dates), nummodel));  mean_omf_ens.mask  = True
        spread_f      = ma.zeros((len(dates), nummodel));  spread_f.mask      = True
        spread_obserr_f = ma.zeros((len(dates), nummodel));  spread_obserr_f.mask = True
        std_omf_ens   = ma.zeros((len(dates), nummodel));  std_omf_ens.mask   = True
        std_omf_ctrl_mean  = ma.zeros((len(dates), nummodel));  std_omf_ctrl.mask  = True
        weight_gefsv12 = np.ones((len(dates)), dtype=np.float32)
        weight_cfsr = np.ones((len(dates)), dtype=np.float32)
        X = np.ones((len(dates), 3), dtype=np.float32)
        Y_gefsv12 = ma.ones((len(dates)), dtype=np.float32)
        Y_cfsr = ma.ones((len(dates)), dtype=np.float32)
        Y_gefsv12_mean = ma.ones((len(dates)), dtype=np.float32)
        Y_cfsr_mean = ma.ones((len(dates)), dtype=np.float32)
        Y_gefsv12_bias = ma.ones((len(dates)), dtype=np.float32)
        Y_cfsr_bias = ma.ones((len(dates)), dtype=np.float32)

        # ---- read all data
         
        for modct in range(nummodel):
            for year in years:
                modelfile = datapath + 'CONV_' + modelstreams[modct] + '_' + \
                    str(year) + '_' + var + '_' + region + '.nc'
                anndata  = Dataset(modelfile, 'r')
                plevs = anndata['Plevels'][:].tolist()
                plevindex = plevs.index(int(plevel))
                thisdates = anndata['All_Dates'][:]
                indx_in = np.where(np.in1d(thisdates, dates)) [0]
                indx_out = np.where(np.in1d(dates, thisdates)) [0]

                nobs_all[indx_out,modct]  = anndata['nobs_all'][indx_in,plevindex]
                nobs_used[indx_out,modct] = anndata['nobs_used'][indx_in,plevindex]
                mean_obs_all[indx_out,modct]  = anndata['mean_obs_all'][indx_in,plevindex]
                mean_obs_used[indx_out,modct] = anndata['mean_obs_used'][indx_in,plevindex]
                mean_omf_ctrl[indx_out,modct] = anndata['mean_omf_ctrl'][indx_in,plevindex]
                mean_omf_ens[indx_out,modct]  = anndata['mean_omf_ens'][indx_in,plevindex]
                mean_oma_ctrl[indx_out,modct] = anndata['mean_oma_ctrl'][indx_in,plevindex]
                spread_f[indx_out,modct]      = anndata['spread_f'][indx_in,plevindex]
                std_omf_ens[indx_out,modct]   = anndata['std_omf_ens'][indx_in,plevindex]
                std_omf_ctrl[indx_out,modct]  = anndata['std_omf_ctrl'][indx_in,plevindex]
                spread_obserr_f[indx_out,modct] = anndata['spread_obserr_f'][indx_in,plevindex]
                std_oma_ctrl[indx_out,modct]  = anndata['std_oma_ctrl'][indx_in,plevindex]
                anndata.close()
        
        
        # ---- set up data for linear regression of RMSE

        for i in range(len(dates)):
            if std_omf_ctrl.mask[i,0] > 1.4: 
                weight_gefsv12[i] = 0.8 
            if std_omf_ctrl.mask[i,1] > 1.4: 
                weight_cfsr[i] = 0.8
                            
            if std_omf_ctrl.mask[i,0] == True: 
                weight_gefsv12[i] = 0.0001 
                Y_gefsv12[i] = 1.4
            else:
                Y_gefsv12[i] = std_omf_ctrl[i,0]
            if std_omf_ctrl.mask[i,1] == True: 
                weight_cfsr[i] = 0.0001
                Y_cfsr[i] = 1.4
            else:
                Y_cfsr[i] = std_omf_ctrl[i,1]
                
            X[i,0] = decimalyear[i]
            f = decimalyear[i] - int(decimalyear[i])
            X[i,1] = np.cos(2.*3.1415926*f)
            X[i,2] = np.sin(2.*3.1415926*f)
            
                
        # ---- set up data for kernel fitting of RMSE bias

        rscale = 10.0
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
                Y_gefsv12_mean[i] = sumprod1 / sumwt1
            else:
                Y_gefsv12_mean[i] = 0.0
                Y_gefsv12_mean[i] = ma.masked
            
            if sumwt2 > 5:
                Y_cfsr_mean[i] = sumprod2 / sumwt2
            else:
                Y_cfsr_mean[i] = 0.0
                Y_cfsr_mean[i] = ma.masked

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

        reg = LinearRegression().fit(X, Y_gefsv12, sample_weight=weight_gefsv12)
        yfit_gefsv12 = reg.predict(X)
        reg = LinearRegression().fit(X, Y_cfsr, sample_weight=weight_cfsr)
        yfit_cfsr = reg.predict(X)
        print (datetime_list[0], datetime_list[1], datetime_list[-1])
        
        f = plt.figure(figsize=(10.,8.))
        
        ax = f.add_axes([.08,.57,.88,.38])
        ax.plot(decimalyear,std_omf_ctrl[:,0],marker='o',color='Red',\
            lw=0, markersize=0.8, label='GEFSv12 reanalysis', markerfacecolor='Red')
        ax.plot(decimalyear,std_omf_ctrl[:,1],marker='o',color='Blue', \
            markersize=0.8, lw=0, label='CFSR', markerfacecolor='Blue')
        #ax.plot(decimalyear,Y_gefsv12_mean,color='Red',lw=1)
        #ax.plot(decimalyear,Y_cfsr_mean,color='Blue',lw=1)
        ax.plot(decimalyear[0:-1],yfit_gefsv12[0:-1],color='Red',lw=1)
        ax.plot(decimalyear[0:-60],yfit_cfsr[0:-60] ,color='Blue',lw=1)
        plt.ylabel('Forecast - observed RMSE (deg C)',fontsize=13)
        ax.legend(loc='upper center')
        plt.grid(True,lw=0.25)
        plt.title('(a) Conventional temperature background forecast errors, '+\
            str(plevel)+' hPa to '+str(int(plevel)+100)+' hPa',fontsize=16)
        ax.set_xlim(2000,2020)
        #ax.set_xlim(2000-01-01, 2021-01-01)
        ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
        ax.set_xlabel('Date',fontsize=13)
    
        ax = f.add_axes([.08,.07,.88,.38])
        ax.plot(decimalyear,-mean_omf_ctrl[:,0],marker='o',color='Red',\
            lw=0, markersize=0.8, label='GEFSv12 reanalysis', markerfacecolor='Red')
        ax.plot(decimalyear,-mean_omf_ctrl[:,1],marker='o',color='Blue', \
            markersize=0.8, lw=0, label='CFSR', markerfacecolor='Blue')
        ax.plot(decimalyear,-Y_gefsv12_bias,color='Red',lw=1)
        ax.plot(decimalyear,-Y_cfsr_bias,color='Blue',lw=1)
        plt.ylabel('Forecast - observed bias (deg C)',fontsize=13)
        ax.legend(loc='upper center')
        plt.grid(True,lw=0.25)
        plt.title('(b) Conventional temperature background forecast biases, '+\
            str(plevel)+' hPa to '+str(int(plevel)+100)+' hPa',fontsize=16)
        ax.set_xlim(2000,2020)
        ax.plot([2000,2020],[0,0],lw=1.5,color='LightGray')
        #ax.set_xlim(2000-01-01, 2021-01-01)
        ax.set_xticks([2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020])
        ax.set_xlabel('Date',fontsize=13)

        plt.savefig(figname,dpi=400)
        print ('Plot done', figname)
        plt.close()        
        
        
        
        