"""
analyzed_climatological_probability.py cmonth cend_hour ctype

"""

import os, sys
from datetime import datetime
import numpy as np
import numpy.ma as ma
import _pickle as cPickle
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.stats as stats
rcParams['legend.fontsize']='medium'

# =====================================================================
    
    
# ---- inputs from command line

cmonth = sys.argv[1] # '01', '02' etc.
cend_hour = sys.argv[2] # 06, 12, 18, 00 -- end hour of 6-h period
ctype = sys.argv[3] # fullfield, thinned, upscaled
imonth = int(cmonth) - 1
nstride = 1
plotit = True

# ---- set parameters

pflag = False # for print statements
master_directory = '/Volumes/NBM/conus_panal/'
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cmonths_early = ['12','01','02','03','04','05','06','07','08','09','10','11']
cmonths_late =  ['02','03','04','05','06','07','08','09','10','11','12','01']
thresholds = [0.254, 1.0, 5.0, 10.0, 25.0] # mm
#thresholds = [0.254]
nthresh = len(thresholds)

# ---- determine the overall number of daily precipitation 
#      samples across all years for this month

iearly = int(cmonths_early[imonth])-1
ilate = int(cmonths_late[imonth])-1

if imonth != 1:  # not Feb
    nsamps_mid = ndaysomo[imonth]*18
else:
    nsamps_mid = 4*ndaysomo_leap[imonth] + 14*ndaysomo[imonth]
    
if iearly != 1:  # not Feb    
    nsamps_early = ndaysomo[iearly]*20
else:
    nsamps_early = 4*ndaysomo_leap[iearly] + 14*ndaysomo[iearly]
if ilate != 1:  # not Feb    
    nsamps_late = ndaysomo[ilate]*20
else:
    nsamps_late = 4*ndaysomo_leap[ilate] + 14*ndaysomo[ilate]
nsamps = nsamps_mid + nsamps_early + nsamps_late

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time as well as the surrounding
#      two months.  All dates for this month have
#      been smushed into one leading index, dimension nsamps,
#      since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution 
#      fitting.

ktr = 0
for iyear in range(2002,2020):
    #print (iyear)
    for cmo in [cmonth, cmonths_early[imonth], cmonths_late[imonth]]:
        imo = int(cmo)-1
        if iyear%4 == 0:
            ndays = ndaysomo_leap[imo]
        else:
            ndays = ndaysomo[imo]
        cyear = str(iyear)  
        if ctype != 'fullfield':  
            infile = master_directory + cyear + cmo + \
                '_ccpa_on_ndfd_grid_6hourly_'+ctype+'.nc' 
        else:
            infile = master_directory + cyear + cmo + \
                '_ccpa_on_ndfd_grid_6hourly.nc'       
        print (infile)
        nc = Dataset(infile)
        
        for iday in range(1,ndays+1):
            precip_in = np.squeeze(nc.variables['apcp_anal'][iday-1,:,:])
            if iyear == 2002 and iday == 1 and cmo == cmonth:
                nyin, nxin = np.shape(precip_in)
                print ('nyin, nxin = ', nyin, nxin)
                precip_tseries = np.zeros((nsamps,nyin,nxin), \
                    dtype=np.float64)
                missingv = -99.99*np.ones((nyin, nxin), dtype=np.float64)
                lons = nc.variables['lons'][:,:]
                lats = nc.variables['lats'][:,:]
                ones = np.ones((nyin, nxin), dtype=np.float64)
                zeros = np.zeros((nyin, nxin), dtype=np.float64)
            precip_in = np.where(precip_in < 500., precip_in, missingv)
            precip_tseries[ktr,:,:] = precip_in[:,:]
            ktr = ktr+1
        nc.close()

# --- set up for making plots

colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,\
    urcrnrlon = 300.95782, urcrnrlat = 54.37,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons, lats)
print ('Done computing m(lons, lats)')

# --- loop over thresholds, calculate probability, and plot

bine_numerator = np.zeros((nsamps,nyin,nxin), dtype=np.int32)
bine_denominator = np.zeros((nsamps,nyin,nxin), dtype=np.int32)
climatological_probability = np.zeros((nthresh, nyin, nxin), dtype=np.float64)

#print (np.shape(precip_tseries))
for ithresh, thresh in enumerate(thresholds):
    
    # thresholds = [0.254, 1.0, 5.0, 10.0, 25.0] # mm
    if thresh == 0.254 or thresh == 1.0:
        clevs = [0,0.01, 0.03, 0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8 ]
    elif thresh == 5.0 or thresh == 10.0:
        clevs = [0,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.12,0.14,0.16]
    else:
        clevs = [0,0.002,0.004,0.006,0.008,0.01,0.013,\
            0.016,0.02,0.025,0.03,0.04,0.06,0.08]
    
    # --- calculate probability
    
    print ('thresh = ', thresh)
    bine_numerator[:,:] = 0
    bine_denominator[:,:] = 0
    a = np.where(precip_tseries >= thresh)
    apos = np.where(precip_tseries >= 0.0)
    bine_numerator[a] = 1
    bine_denominator[apos] = 1
    bine_numerator_count = np.sum(bine_numerator, axis=0)
    bine_denominator_count = np.sum(bine_denominator, axis=0)
    climatological_probability[ithresh,:,:] = \
        np.where(bine_denominator_count > 0,
        bine_numerator_count.astype(float) / bine_denominator_count.astype(float),\
        -99.99*ones)
    
    
    if plotit == True:
        
        # --- make plots of probability.

        fig = plt.figure(figsize=(8.,6.5))
        axloc = [0.02,0.1,0.96,0.81]
        ax1 = fig.add_axes(axloc)
        title = cmonths[imonth]+' CCPA/MSWEP probability of exceeding '+\
            str(thresh)+' mm, 6 hour period ending '+cend_hour+' UTC'
        ax1.set_title(title, fontsize=13,color='Black')
        CS2 = m.contourf(x, y, climatological_probability[ithresh,:,:], clevs,\
            cmap=None, colors=colorst, extend='both')
    
        m.drawcoastlines(linewidth=0.8,color='Gray')
        m.drawcountries(linewidth=0.8,color='Gray')
        m.drawstates(linewidth=0.8,color='Gray')
    
        # ---- use axes_grid toolkit to make colorbar axes.

        cax = fig.add_axes([0.06,0.07,0.88,0.02])
        cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
            drawedges=True,ticks=clevs,format='%g')
        cb.ax.tick_params(labelsize=7)
        cb.set_label('Probability',fontsize=9)

        # ---- set plot title

        plot_title = 'ccpa_probability_'+ctype+'_'+str(thresh)+\
            'mm_month'+cmonth+'_'+cend_hour+'UTC.png'
        fig.savefig(plot_title, dpi=300)
        print ('saving plot to file = ',plot_title)
    
# ---- save array to cPickle file

outfile = master_directory + cmonth + '_'+cend_hour+\
    'UTC_ccpa_on_ndfd_grid_6hourly_climo_probability_'+ctype+'.cPick'
print ('writing to ', outfile) 
ouf = open(outfile, 'wb')
cPickle.dump(climatological_probability, ouf) 
cPickle.dump(lons, ouf)
cPickle.dump(lats, ouf)
cPickle.dump(thresholds, ouf)
ouf.close()

print ('Done!')