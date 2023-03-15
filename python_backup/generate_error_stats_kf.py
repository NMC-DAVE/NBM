"""
generate_error_stats_kf.py 

"""
import pygrib
from dateutils import daterange, dateshift, dayofyear, splitdate
import os, sys
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'
rcParams["contour.negative_linestyle"]='solid'


# =====================================================================

clead = sys.argv[1]  # lead time, e.g., 12, 72, 120 (in hours)
alpha = 0.02
calpha = str(alpha)
makeplots = True

ilead = int(clead)
datadir = '/Users/Tom/python/ecmwf/'
cvariable = '2t'
datestart = dateshift('2019010100',ilead)
#date_list_anal = daterange(datestart,'2019123100',24)
date_list_anal = daterange(datestart,'2019013100',24)
ndates = len(date_list_anal)
rmse_kf = np.zeros((ndates), dtype=np.float32)
rmse_decayavg = np.zeros((ndates), dtype=np.float32)
rmse_raw = np.zeros((ndates), dtype=np.float32)
date_list_fcst = []
date_list_bcorr = []
frac2019 = np.zeros((ndates), dtype=np.float32)
for idate in range(ndates):
    date_list_fcst.append(dateshift(date_list_anal[idate],-ilead)) # initial times of fcst
    date_list_bcorr.append(dateshift(date_list_anal[idate],-2*ilead))

for idate, datea in enumerate(date_list_anal):
    
    datef = date_list_fcst[idate]
    dateb = date_list_bcorr[idate]
    if datea == '2019010100': dstart = idate
    #print ('------ processing analysis, forecast dates = ', datea, datef)
    yyyy,mm,dd,hh = splitdate(datef) 
    doy = dayofyear(yyyy,mm,dd)
    fracyear = doy/365.

    # ---- read the ECMWF ERA5 reanalysis at this analysis date.
    
    infile = datadir + 't2m_era5_halfdegree_'+datea+'.cPick'
    inf = open(infile, 'rb')
    analysis = cPickle.load(inf)
    if idate == 0:
        lats = cPickle.load(inf)
        lons = cPickle.load(inf)
        nlats, nlons = np.shape(lats)
        npts = nlats*nlons 
    inf.close()
    
    # ---- read the ECMWF control forecast at this lead time and initial date
 
    infile = datadir + cvariable+'_'+datef+'_f'+clead+'.grib2'  
    grbfile = pygrib.open(infile) 
    grb = grbfile.select()[0] 
    forecast = grb.values
    grbfile.close()
    
    # ---- read the decaying average and Kalman filter bias corrections estimates
    #      for this date.
    
    infilename = datadir + 'bias_est_'+dateb+'_f'+clead+'.cPick'
    inf = open(infilename, 'rb')
    bias_decayavg = cPickle.load(inf)
    bias_estimate = cPickle.load(inf)
    inf.close()
    
    frac2019[idate] = fracyear = doy/365.
    
    rmse_raw[idate] = np.sqrt(np.sum((forecast-analysis)**2)/(npts-1.))
    rmse_kf[idate] = np.sqrt(np.sum(((forecast-bias_estimate)-analysis)**2)/(npts-1.))
    rmse_decayavg[idate] = np.sqrt(np.sum(((forecast-bias_decayavg)-analysis)**2)/(npts-1.))
    

    if makeplots == True:
        
        # ---- plot standard bias correction approach

        fig = plt.figure(figsize=(6.5,9.0))
        fig.suptitle('Lead time = '+clead+' h, IC = '+datef, fontsize=18)
        colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#ccccff','#e6e6ff', \
            'White', '#ffe6e6', '#ffcccc', '#ffb2b2', '#ff7373', '#ff0000']
        clevs = [-4.0,-3.0,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,3.0,4.0]

        axloc = [0.02,0.53,0.96,0.4]
        a1 = fig.add_axes(axloc)
        a1.set_title('(a) Standard decaying-average bias correction',fontsize=14)
        m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],\
            urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
            resolution='l', projection='mill')
        x, y = m(lons, lats)
        CS2 = m.contourf(x,y,bias_decayavg,clevs,cmap=None,colors=colorst,extend='both')
        CS1 = m.contour(x,y,bias_decayavg,clevs,linestyle='-',colors='white',linewidths=0.5)
        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
    
        axloc = [0.02,0.1,0.96,0.4]
        a2 = fig.add_axes(axloc)
        a2.set_title('(b) Simultaneous CONUS-wide Kalman-filter bias correction',fontsize=14)
        m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],\
            urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],
            resolution='l', projection='mill')
        CS2 = m.contourf(x,y,bias_estimate,clevs,cmap=None,colors=colorst,extend='both')
        CS1 = m.contour(x,y,bias_estimate,clevs,linestyle='-',colors='white',linewidths=0.5)
        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')

        # ---- use axes_grid toolkit to make colorbar axes.

        divider = make_axes_locatable(a2)
        cax = cax = fig.add_axes([0.02,0.06,0.96,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.set_label('Bias estimate (deg C)',fontsize=12)

        plot_title = 'figures/CONUS_bias_estimates_lead'+clead+'_IC'+datef+'.png'
        print ('saving plot to file = ',plot_title)
        plt.savefig(plot_title,dpi=300)
        print ('Plot done')
        #sys.exit()


rmse_raw_mean = np.mean(rmse_raw)
rmse_kf_mean = np.mean(rmse_kf)
rmse_decayavg = np.mean(rmse_decayavg)
print ('rmse raw kf decay = ',rmse_raw_mean, rmse_kf_mean, rmse_decayavg )
    
    
    
    
            

