
""" estimate_cov_params_analyses_10pct.py : estimate the parameters for a analysis-error covariance model. """

import numpy as np
import numpy.ma as ma
import math
import sys
import scipy.stats.mstats as stats
import _pickle as cPickle
import time
import matplotlib.pyplot as plt
from netCDF4 import Dataset, chartostring
import math
from scipy.optimize import curve_fit
from numba import jit
from cressman_update_background import cressman_update_background # python not fortran
from dateutils import dayofyear
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
    
# --------------------------------------------------------------------------    
    
@jit(nopython=True)
def gamma_exponential(data_array, rho_horiz, vdconst, cpconst, gamma):
    term = np.sqrt(data_array[0,:]**2/rho_horiz**2 + \
        data_array[1,:]**2/vdconst**2 + data_array[2,:]**2/cpconst**2)
    expterm_gamma = np.exp(-term**(gamma))      
    return expterm_gamma   
    
#@jit(nopython=True)
def gamma_exponential_nugget(data_array, rho_horiz, vdconst, cpconst, gamma, nugget):
    term = np.sqrt(data_array[0,:]**2/rho_horiz**2 + \
        data_array[1,:]**2/vdconst**2 + data_array[2,:]**2/cpconst**2)
    expterm_gamma = (1.-nugget)*np.exp(-term**(gamma))      
    return expterm_gamma       
    
    
def gamma_exponential_nugget_horiz(data_array, rho_horiz, gamma, nugget):
    term = np.sqrt(data_array[0,:]**2/rho_horiz**2)
    expterm_gamma = (1.-nugget)*np.exp(-term**(gamma))      
    return expterm_gamma  
    
# ---- initialize blank list to hold information on correlation between stations and 
#      the associated horizontal distance, vertical distance, and coastal proximity
#      difference.

makeplots = True

cmonthname = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#cmonthname = ['Mar']

#chours = ['00','01','02','03','04','05','06','07','08','09','10',\
#    '11','12','13','14','15','16','17','18','19','20','21','22','23']
chours = ['00']
    
dates_analysis_error_variance = ['2017123123','2018013023', '2018030123', \
    '2018033123', '2018043023', '2018053023', '2018062923', '2018072923', \
    '2018082823', '2018092723', '2018102723', '2018112623', '2018122623']
    
juldates_avariance = []
for idate in range(13):
    yyyy = int(dates_analysis_error_variance[idate][0:4])
    mm = int(dates_analysis_error_variance[idate][4:6])
    dd = int(dates_analysis_error_variance[idate][6:8])
    juldates_avariance.append(dayofyear(yyyy,mm,dd))  
    
center_julday = np.zeros((12), dtype=np.int16)
for imonth in range(12):
    center_julday[imonth] = dayofyear(2018,imonth+1,15)

# ---- read observation locations in prism grid coordinates 
#      that are both inside prism and hrrr grids.

infilename = 'obslocns_closest_prism_lons_indices.npy'
print (infilename)
inearest = np.load(infilename)

infilename = 'obslocns_closest_prism_lats_indices.npy'
print (infilename)
jnearest = np.load(infilename)

# ---- read in prism lat, lon, mask

infilename = \
    '/data/thamill/Rf2_tests/sfcdata/prism/diurnal_cycle_temp_CONUS_julian_day_0.nc'
nc = Dataset(infilename)
climo_temperature_prism = nc.variables['temperature'][0,:,:]
prism_lons = nc.variables['lons_anal'][:,:]
prism_lats = nc.variables['lats_anal'][:,:]
nyprism, nxprism = np.shape(prism_lons)
zeros = np.zeros((nyprism, nxprism), dtype=np.int16)
ones = np.ones((nyprism, nxprism), dtype=np.int16)
prism_mask = np.ones((nyprism, nxprism), dtype=np.int16)
prism_mask = np.where(climo_temperature_prism < -999, zeros, ones)
nc.close()

# ---- loop thru hours and days to get correlations 

ktr = 0
for chour in chours:
    ihour = int(chour)
    for imonth, cmonth in zip(range(12),cmonthname):
        
        print (' ------------- processing hour, month ', \
            chour, cmonth, time.asctime()) 
            
        # ---- read in the station-derived correlations and station variances, calculated
        #      previously in station_pair_correlations_statfirstg.py
                    
        #infile = '/data/thamill/sfcdata/station_pair_correlations_00UTC_'+cmonth+'.cpick'
        infile = '/data/thamill/sfcdata/station_pair_correlations_00UTCanal_10pct_'+cmonth+'.cpick'
        
        print ('reading from ', infile)
        inf = open(infile,'rb')
        corr = np.squeeze(cPickle.load(inf))
        hdist = np.squeeze(cPickle.load(inf))
        vdist = np.squeeze(cPickle.load(inf))
        cpdiff = np.squeeze(cPickle.load(inf))
        varis = np.squeeze(cPickle.load(inf))
        ones = np.ones(len(corr), dtype=np.float32)
        print ('mean of station variances = ', np.mean(varis))
        sample_var = 1.0 - corr**2
        sample_var = np.where(corr < 0.25, 1.0 - 0.0625*ones, sample_var)
                
        lons_stations = np.squeeze(cPickle.load(inf))
        lats_stations = np.squeeze(cPickle.load(inf))
        nstations = len(lons_stations)
        inf.close()
        nsamps = len(hdist)
        
        # ---- determine which analysis-error variance file to use as a first guess
        #      in making an estimate of background-error variance.   Then read in that 
        #      estimate from file 

        diff = np.abs(center_julday[imonth] - np.array(juldates_avariance))
        indices = np.argmin(diff)
        if hasattr(indices, "__len__"):
            idx = indices[0]
        else:
            idx = indices

        infile = '/data/thamill/sfcdata/analysis_error_variance_degC2_'+\
            dates_analysis_error_variance[idx]+'.npy'
        print ('reading ', infile)
        background = np.load(infile, allow_pickle=True) *3.24  # convert to deg F
        print ('mean gridded 1st guess for background error variance, deg F = ', \
            np.sum(background*prism_mask)/np.sum(prism_mask)) 

        # ---- perform 1-pass cressman analysis to update the analysis error variance
        #      to the station 1st-guess error variances.
        
        print ('min, max varis = ',np.min(varis), np.max(varis))
        niter = 2
        roi = np.array([500., 250.])
        background_error_variance = cressman_update_background (background, \
            prism_lons, prism_lats, prism_mask, nyprism, nxprism, \
            varis, jnearest, inearest, lats_stations, \
            lons_stations, nstations, niter, roi)
        print ('after updating with station variance data, background_error_variance[ny//2, 0:nx:10] = ', \
            background_error_variance[nyprism//2, 0:nxprism:10] )

        # ---- convert lists to numpy vectors
               
        data_array = np.zeros((3,len(hdist)),dtype=np.float32)
        data_array[0,:] = hdist[:]
        data_array[1,:] = vdist[:]
        data_array[2,:] = cpdiff[:]
        print ('max, min hdist = ', np.max(hdist), np.min(hdist))
        print ('max, min vdist = ', np.max(vdist), np.min(vdist))
        print ('max, min cpdiff = ', np.max(cpdiff), np.min(cpdiff))
        a = np.where(hdist != hdist)
        print ('nans for hdist = ', a)
        a = np.where(vdist != vdist)
        print ('nans for vdist = ', a)
        a = np.where(cpdiff != cpdiff)
        print ('nans for cpdiff = ', a)
        a = np.where(corr != corr)
        print ('nans for corr = ', a)

        # ------ estimate covariance parameters
        
        rho_horiz = 100.  # set to conservative values determined from Jul data.
        vdconst = 500.
        cpconst = 400. 
        gamma = 1.8  
        nugget = 1.0-np.max(corr)
            
        popt, pcov = curve_fit(gamma_exponential_nugget_horiz, \
            data_array, corr, p0=(rho_horiz, gamma, nugget), \
            check_finite=True,  method='trf', diff_step=0.000001, \
            sigma = sample_var, absolute_sigma=False, \
            bounds = ([0.,1.0,0.0], [2500.,2.0,1.0]))
        #popt, pcov = curve_fit(gamma_exponential_nugget, \
        #    data_array, corr, p0=(rho_horiz, vdconst, cpconst, gamma, nugget), \
        #    check_finite=True,  method='trf', diff_step=0.000001, \
        #    bounds = ([0.,0.,0.,1.0,0.0], \
        #    [2500., 5000., 10000.,2.0,1.0]))
            
        vari0 = np.mean(np.array(varis))
        rho_horiz = popt[0]
        #vdconst = popt[1]
        #cpconst = popt[2]
        gamma = popt[1]
        nugget = popt[2]
        
        rho_horiz_sprd = np.sqrt(pcov[0,0])
        #vdconst_sprd = np.sqrt(pcov[1,1])
        #cpconst_sprd = np.sqrt(pcov[2,2])
        gamma_sprd = np.sqrt(pcov[1,1])
        nugget_sprd = np.sqrt(pcov[2,2])
        
        #print ('rho_horiz, vdconst, cpconst, gamma, nugget = ', \
        #    rho_horiz, vdconst, cpconst, gamma, nugget)
        print ('rho_horiz, gamma, nugget = ', \
            rho_horiz,  gamma, nugget)

        # ----   make plot of horizontal covariance functions, 
        #        now for gamma-exponential with the correct data
        
        if makeplots == True:
            title = cmonth + ' '+chour+' UTC horizontal correlation from station pairs'
            plot_title = 'gamma_exponential_horizontal_'+cmonth+'_00UTC_00UTCanal.pdf'
            fig = plt.figure(figsize=(6., 4.5))
            a1 = fig.add_axes([.14,.11,.82,.8])
            a1.set_title(title,fontsize=14)
            dist = np.arange(0,501,1)
            term = dist/popt[0]
            corrfn = (1.-nugget)*np.exp(-term**gamma) # gamma-exponential
            a1.scatter(hdist, corr, marker='o', s=0.1, color='Gray',zorder=8)
            a1.plot(range(501),corrfn,'-',color='Red',linewidth=3,zorder=11)
            a1.plot([0,500],[0,0],'-',color='Black',linewidth=1,zorder=10)
            a1.set_xlim(0,500)
            a1.set_ylim(-0.4,1.02)
            a1.set_xlabel('Distance (km)',fontsize=13)
            a1.set_ylabel('Correlation',fontsize=13)
            a1.grid(color='LightGray',linewidth=0.3)
            print ('saving plot to file = ',plot_title)
            plt.savefig(plot_title)
            plt.close()
            print ('Plot done')
        
        # ---- save data assimilation parameter estimates to file
        
        outfile = 'data/00UTC_anal_gamma_exponential_nugget_horiz_corr_month='+\
            cmonth+'_hour='+chour+'.pydat'
        ouf = open(outfile,'wb')
        cPickle.dump(rho_horiz, ouf)
        #cPickle.dump(vdconst, ouf)
        #cPickle.dump(cpconst, ouf)
        cPickle.dump(gamma, ouf)
        cPickle.dump(nugget, ouf)
        cPickle.dump(background_error_variance, ouf)
        cPickle.dump(prism_lons, ouf)
        cPickle.dump(prism_lats, ouf)
        cPickle.dump(prism_mask, ouf)
        #cPickle.dump(background, ouf)
        cPickle.dump(rho_horiz_sprd, ouf) 
        #cPickle.dump(vdconst_sprd, ouf)
        #cPickle.dump(cpconst_sprd, ouf)
        cPickle.dump(gamma_sprd, ouf)
        cPickle.dump(nugget_sprd, ouf)         
        ouf.close()
        
        # ---- generate a 4-panel plot    

        fig1 = plt.figure(figsize=(9.,5.6))

        # ---- add colorbar
        
        m = Basemap(llcrnrlon=np.min(prism_lons),\
            llcrnrlat=np.min(prism_lats),urcrnrlon=np.max(prism_lons),\
            urcrnrlat=np.max(prism_lats),projection='mill',resolution='l')
        colorstblack='Black'

        x, y = m(prism_lons, prism_lats)  
        colorst = ['White','#ECFFFF','#D9F7FF','#C4E8FF','#E8FBE8','#C7F4C7','#92F592','Yellow',\
            'Orange','#FFB2B2','#EC5B71','Magenta','DarkOrchid','Black']
        levels = [0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,4.0]
        ax = fig1.add_axes([0.02,0.11,0.96,0.82])
        title = cmonth+ r' T$_{2m}$ first guess uncertainty $diag({\bf{P}}^a )^{1/2}$'
        ax.set_title(title,fontsize=10.)
        cf = ax.contourf(x,y,(np.sqrt(background_error_variance)/1.8)*prism_mask,cmap=None,colors=colorst,\
            vmin=0., vmax=4,levels=levels)    
        m.drawcoastlines(linewidth=.5,color='Gray') 
        m.drawstates(linewidth=.5,color='Gray')
        m.drawcountries(linewidth=.5,color='Gray')

        ax = fig1.add_axes([0.02,0.07,0.96,0.82])
        ax.set_frame_on(False)
        ax.set_visible(False)
        cb = m.colorbar(cf,"bottom", size="3%", pad="1%",\
            ticks=levels) # im1
        cb.set_label(r'Background uncertainty estimate $diag({\bf{P}}^a )$ (deg C)',fontsize=9)
        cb.ax.tick_params(labelsize=7)

        plot_title = 'anal_uncertainty_'+cmonth+'.png'
        print ('saving plot to file = ',plot_title)
        fig1.savefig(plot_title,dpi=300)
        plt.close()