# plot_mswep_spline_cdf.py mm lead lat lon

import numpy as np
import _pickle as cPickle
import sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.stats as stats
from scipy.interpolate import LSQUnivariateSpline, splrep, splev

rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'

# ===========================================================

def find_nearest(vec, value):
    idx = np.abs(vec-value).argmin()
    return idx
    
# ===========================================================

def fraczero_possamps(nsamps, precip):
    """
    from the vector input sample precip_ens, define the fraction of
    samples with zero precipitation.   For the positive samples, add
    a small random number to deal with the fact that the data was 
    discretized to 0.1 mm, so that when later creating CDFs we don't 
    have values with lots of tied amounts.   Sort the nonzero amounts 
    and return.
    """
    number_zeros = 0
    precip_nonzero = np.delete(precip, \
        np.where(precip <= 0.0))  # censor at 0.1 mm
    nz = len(precip_nonzero)
    # data discretized, so add random component of this magnitude
    precip_nonzero = precip_nonzero + \
        np.random.uniform(low=-0.005,high=0.005,size=nz) 
    precip_nonzero = np.sort(precip_nonzero)  
    #print (precip_ens_nonzero[0:10]) 
    ntotal = len(precip)
    nzero = ntotal - len(precip_nonzero)
    fraction_zero = float(nzero) / float(ntotal)
    return fraction_zero, precip_nonzero, nz
    
# ===========================================================


cmonth = sys.argv[1] # 01, etc.
clead = sys.argv[2]
clat = sys.argv[3]
clon = sys.argv[4]
rlon = float(clon)
rlat = float(clat)

cclon = str(rlon)
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
imonth = int(cmonth)-1
ccmonth = cmonths[imonth]
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
cmonths_early = ['12','01','02','03','04','05','06','07','08','09','10','11']
cmonths_late =  ['02','03','04','05','06','07','08','09','10','11','12','01']
master_directory = '/Volumes/Backup Plus/mswep/'
cdomain = 'conus'

# ---- determine the overall number of daily precipitation 
#      samples across all years for this month

iearly = int(cmonths_early[imonth])-1
ilate = int(cmonths_late[imonth])-1

if imonth != 1:  # not Feb
    nsamps_mid = ndaysomo[imonth]*20
else:
    nsamps_mid = 6*ndaysomo_leap[imonth] + 14*ndaysomo[imonth]
if iearly != 1:  # not Feb    
    nsamps_early = ndaysomo[iearly]*20
else:
    nsamps_early = 6*ndaysomo_leap[iearly] + 14*ndaysomo[iearly]
if ilate != 1:  # not Feb    
    nsamps_late = ndaysomo[ilate]*20
else:
    nsamps_late = 6*ndaysomo_leap[ilate] + 14*ndaysomo[ilate]
nsamps = nsamps_mid + nsamps_early + nsamps_late
print ('nsamps = ', nsamps)

# --- load spline information from file.   These also store Gamma
#     parameters for situations where there are inadequate samples.

infile = master_directory + cmonth+'_'+cdomain+\
    '_MSWEP_spline_info_h' + clead + '.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
spline_info = cPickle.load(inf)
print ('np.shape(spline_info) = ', np.shape(spline_info))
indices_to_query = cPickle.load(inf)
print ('np.shape(indices_to_query) = ', np.shape(indices_to_query))
inf.close()

infile = master_directory + cmonth+'_'+cdomain+\
    '_MSWEP_Dnstat_h' + clead + '.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
Dnstat = cPickle.load(inf)
lons_mswep = cPickle.load(inf)
nlats, nlons = np.shape(lons_mswep)
print ('min, max lons_mswep = ', np.min(lons_mswep), \
    np.max(lons_mswep))
lats_mswep = cPickle.load(inf)
inf.close()

# ---- find the index of the minimum

difflon = np.abs(lons_mswep - rlon)
difflat = np.abs(lats_mswep - rlat)
diffboth = np.sqrt(difflon**2 + difflat**2)
index_array = np.argmin(diffboth)
print ('index_array = ', index_array)
unraveled = np.unravel_index(index_array, (nlats,nlons))
jy = unraveled[0]
ix = unraveled[1]
print ('jmin, imin = ', jy, ix)
print ('lon, lat at minimum = ', lons_mswep[jy, ix], \
    lats_mswep[jy,ix])

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time as well as the surrounding
#      two months.  All dates for this month have
#      been smushed into one leading index, dimension nsamps,
#      since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution 
#      fitting.
   
ktr = 0
for iyear in range(2000,2020):
    for cmo in [cmonth, cmonths_early[imonth], cmonths_late[imonth]]:
        imo = int(cmo)-1
        if iyear%4 == 0:
            ndays = ndaysomo_leap[imo]
        else:
            ndays = ndaysomo[imo]
        cyear = str(iyear)    
        infile = master_directory + cyear + cmo + \
            '_on_ndfd_grid_6hourly.nc'
        print (iyear, infile, ndays)
        nc = Dataset(infile)
        yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
        for iday in range(1,ndays+1):
            if iday < 10:
                cday = '0'+str(iday)
            else:
                cday = str(iday)
            iyyyymmddhh = int(str(iyear)+cmo+cday+clead)
            idx = np.where(yyyymmddhh_end == iyyyymmddhh)[0]
            precip_in = nc.variables['apcp_anal'][idx,jy,ix]
            if iyear == 2000 and iday == 1 and cmo == cmonth:
                precip_tseries = np.zeros((nsamps), \
                    dtype=np.float64)
            precip_tseries[ktr] = precip_in
            ktr = ktr+1
        nc.close()

# ---- determine the fraction zero, the positive samples, and sort

fraction_zero, precip_ens_nonzero, nz = \
    fraczero_possamps(nsamps, precip_tseries) # return sorted
    
empirical_cdf = 1.0/(2.0*nz) + np.arange(nz)/nz   
        
# ---- indices_to_query, if set to -1, indicates that there aren't
#      enough samples to use spline.   In this case, the spline
#      info arrays contain alpha and beta parameters.   Fit CDF 
#      accordingly

if indices_to_query[jy,ix,0] == -1:
    
    # ---- fit Gamma CDF
    
    alpha = spline_info[jy,ix,0]
    beta = spline_info[jy,ix,1]
    y = precip_ens_nonzero / beta
    fitted_cdf = stats.gamma.cdf(y, alpha)
    spg = 'Gamma'
else:
    
    # ---- fit with splines to hazard function. 
    
    splines_tuple = (spline_info[jy,ix,0,:], spline_info[jy,ix,1,:], 3)
    spline_hazard = splev(precip_ens_nonzero, splines_tuple)
    fitted_cdf = 1.0 - np.exp(-spline_hazard)
    spg = 'spline'
        
# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,4.5))#
ax = f.add_axes([.13,.12,.84,.75])
ax.set_title(r'Fitted and empirical MSWEP CDFs, '+\
    ccmonth+',\n6-hour period ending '+clead+' UTC, '+\
    cclon+'$^{\circ}$ W '+clat+'$^{\circ}$ N',fontsize=14)
ax.plot(precip_ens_nonzero,fitted_cdf,color='Blue',lw=2,\
    label='fitted '+spg)
ax.plot(precip_ens_nonzero,empirical_cdf,color='Red',\
    lw=2,label='empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=4)
ax.set_ylim(0.0,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,40)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_MSWEP_'+ccmonth+'_'+\
    cclon+'W_'+clat+'N_'+clead+'UTC.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()


# ---- make plots of Dnstats.   

#Dnstat_plot = (1-fraction_zero)*Dn

title = r'MSWEP D$_n$ for '+\
    ccmonth+', 6-h period ending '+clead+' UTC'    
plot_title = 'Dnstat_mswep_'+ccmonth+'_'+clead+'UTC.png' 

clevs = [0.0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
m = Basemap(llcrnrlon=233.7234-360.,llcrnrlat=19.229,
    urcrnrlon = 300.95782-360., urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_mswep, lats_mswep)    

fig = plt.figure(figsize=(8,6.5))
axloc = [0.08,0.1,0.87,0.88]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=16,color='Black')
CS2 = m.contourf(x, y, Dnstat, clevs,\
    cmap=None, colors=colorst, extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
xdot, ydot = m(rlon,rlat)
print ('xdot, ydot = ',xdot, ydot)
m.plot(xdot,ydot,marker='.',markersize=3,color='Black')
   
# draw parallels and meridians.
# label on left and bottom of map.
parallels = np.arange(20.,60.01,5.)
m.drawparallels(parallels,labels=[1,0,0,0],color='Black',fontsize=8,linewidth=0.3)
meridians = np.arange(220.,360.,5.)
m.drawmeridians(meridians,labels=[0,0,0,1],color='Black',fontsize=8,linewidth=0.3)
   
# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
#cax = divider.append_axes("bottom", size="3%", pad=0.25)
cax = fig.add_axes([0.08, 0.09, 0.88, 0.02])
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label(r'D$_n$ statistic',fontsize=11)

# ---- set plot title

fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')






