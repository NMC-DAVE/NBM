"""
plot_fitted_anal_precip_cdfs.py cmonth clead jy ix

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

import _pickle as cPickle
import scipy.stats as stats


rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

# =====================================================================

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


# =====================================================================
    
# ---- inputs from command line

nstride = 1
cmonth = sys.argv[1] # '01', '12', etc.
imonth = int(cmonth)-1
clead = sys.argv[2] # 06, 12, 18, 00
cjy = sys.argv[3] 
cix = sys.argv[4] 
jy = int(cjy)
ix = int(cix)
cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ccmonth = cmonths[imonth]
cdomain = 'conus'
mswep_directory = '/Volumes/Backup Plus/mswep/'

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time.  All members, dates for this 
#      month have been smushed into one leading index, dimension
#      nsamps, since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution fitting.

if imonth == 1:  # Feb
    nsamps = ndaysomo[imonth]*20
else:
    nsamps = 6*ndaysomo_leap[imonth] + 14*ndaysomo[imonth]
print ('nsamps = ', nsamps)


ktr = 0
for iyear in range(2000,2020):
    if iyear%4 == 0:
        ndays = ndaysomo_leap[imonth]
    else:
        ndays = ndaysomo[imonth]
    cyear = str(iyear)
    infile = mswep_directory + cyear + cmonth + '_on_ndfd_grid_6hourly.nc'
    print (iyear, infile)
    nc = Dataset(infile)
    yyyymmddhh_end = nc.variables['yyyymmddhh_end'][:]
    for iday in range(1,ndays+1):
        if iday < 10:
            cday = '0'+str(iday)
        else:
            cday = str(iday)
        iyyyymmddhh = int(str(iyear)+cmonth+cday+clead)
        idx = np.where(yyyymmddhh_end == iyyyymmddhh)[0]
        precip_in = np.squeeze(nc.variables['apcp_anal'][idx,::nstride,::nstride])
        #precip_in = np.flipud(precip_in)
        if iyear == 2000 and iday == 1:
            nyin, nxin = np.shape(precip_in)
            precip_tseries = np.zeros((nsamps,nyin,nxin), dtype=np.float32)
            lons = nc.variables['lons'][::nstride,::nstride]
            lats = nc.variables['lats'][::nstride,::nstride]
            #lats = np.flipud(lats)
        precip_tseries[ktr,:,:] = precip_in[:,:]
        
        ktr = ktr+1
    nc.close()
    
clat = '%.2f' %(lats[jy,ix]) 
clon = '%.2f' %(lons[jy,ix]) 

# --- determine the empirical fraction_zero, the sorted precip_nonzero

precip = precip_tseries[:,jy,ix]
fraction_zero, precip_nonzero, nz = fraczero_possamps(nsamps, precip)

# --- load the fitted Gamma precipitation parameters from cPickle file

cdomain = 'conus'
infile = mswep_directory + cmonth+'_'+cdomain+\
    '_MSWEP_apcp_gamma_parameters_h' + clead + '.cPick'
print ('reading from ', infile)
inf = open(infile, 'rb')
weights = cPickle.load(inf)
alpha = cPickle.load(inf)
beta = cPickle.load(inf)
fzero = cPickle.load(inf)
Dnstat1 = cPickle.load(inf)
Dnstat2 = cPickle.load(inf)
Dnstat3 = cPickle.load(inf)
nmixture = cPickle.load(inf)

weights = weights[:,::nstride,::nstride]
alpha = alpha[:,::nstride,::nstride]
beta = beta[:,::nstride,::nstride]
fzero = fzero[::nstride,::nstride]
Dnstat1 = Dnstat1[::nstride,::nstride]
Dnstat2 = Dnstat2[::nstride,::nstride]
Dnstat3 = Dnstat3[::nstride,::nstride]
nmixture = nmixture[::nstride,::nstride]
inf.close()

# --- build CDF

x = np.arange(0.0,50.01,0.1)
nx = len(x)
y0 = x / beta[0,jy,ix]
y1 = x / beta[1,jy,ix]
y2 = x / beta[2,jy,ix]

cdf0 = stats.gamma.cdf(y0, alpha[0,jy,ix])
cdf1 = stats.gamma.cdf(y1, alpha[1,jy,ix])
cdf2 = stats.gamma.cdf(y2, alpha[2,jy,ix])
cdf_fitted = fzero[jy,ix] + (1.-fzero[jy,ix])* \
    (weights[0,jy,ix]*cdf0 + weights[1,jy,ix]*cdf1 + \
    weights[2,jy,ix]*cdf2) 

cdf_empirical = np.zeros((nx),dtype=np.float64)
len_nonzero = len(precip_nonzero)
fnzero = float(len_nonzero)
nz = int(fnzero)
pnzsort = np.sort(precip_nonzero)
diffmax = 0.
iofmax = 0
for i, xi in enumerate(x):
    nbelow = (pnzsort < xi).sum() 
    cdf_empirical[i] = fraction_zero + \
        (1.-fraction_zero)*float(nbelow)/fnzero
    absdiff = np.abs(cdf_empirical[i] - cdf_fitted[i])
    if absdiff > diffmax:
        diffmax = absdiff
        iofmax = i
        
print ('maximum difference at ', x[iofmax],\
    '  absdiff, cdf_empir, cdf_fitted = ',absdiff,\
    cdf_empirical[iofmax], cdf_fitted[iofmax])
print ('Dnstat1, 2, 3 here: ', Dnstat1[jy,ix], \
    Dnstat2[jy,ix], Dnstat3[jy,ix])

        
# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,4.5))#
ax = f.add_axes([.13,.12,.84,.75])
ax.set_title(r'Fitted and empirical MSWEP CDFs, '+\
    ccmonth+',\n6-hour period ending '+clead+' UTC, '+\
    clon+'$^{\circ}$ W '+clat+'$^{\circ}$ N',fontsize=14)
ax.plot(x,cdf_fitted,color='Blue',lw=2,\
    label='Fitted 3-component Gamma mixture')
ax.plot(x,cdf_empirical,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=4)
ax.set_ylim(0.7,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,20)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_y'+cjy+'_x'+cix+'_'+clead+'UTC.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()



m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons, lats)

#m = Basemap(llcrnrlon=lons[0,0],llcrnrlat=lats[-1,-1],\
#    urcrnrlon=lons[-1,-1],urcrnrlat=lats[0,0],\
#    resolution='l', projection='mill')
#x, y = m(lons, lats)
nj, ny, nx = np.shape(alpha)

# ---- make plots of Dnstats.   ymask out negatives

d1 = ma.minimum((1-fzero)*Dnstat1, (1-fzero)*Dnstat2)
Dnstat_plot = ma.minimum(d1,(1-fzero)*Dnstat3)
print ('Dnstat_plot, 1.-fzero = ', Dnstat_plot[jy,ix], 1.-fzero[jy,ix])

title = r'D$_n$ statistic overall, MSWEP, '+\
    ccmonth+' 6-h period ending '+clead+' UTC'    
plot_title = 'Dnstat_mswep_y'+cjy+'_x'+cix+'_'+clead+'UTC.png' 

clevs = [0.0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

fig = plt.figure(figsize=(8,6.5))
axloc = [0.06,0.1,0.91,0.82]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=13,color='Black')
CS2 = m.contourf(x, y, Dnstat_plot, clevs,\
    cmap=None, colors=colorst, extend='both')
m.drawcoastlines(linewidth=0.8,color='Gray')
m.drawcountries(linewidth=0.8,color='Gray')
m.drawstates(linewidth=0.8,color='Gray')
    
xdot, ydot = m(lons[jy,ix],lats[jy,ix])
m.plot(xdot,ydot,marker='.',markersize=2,color='Black')
    
# draw parallels and meridians.
# label on left and bottom of map.
#parallels = np.arange(-20.,80,10.)
#m.drawparallels(parallels,labels=[1,0,0,0],color='LightGray')
#meridians = np.arange(0.,360.,20.)
#m.drawmeridians(meridians,labels=[0,0,0,1],color='LightGray')
       
ax.set_xlim(0,nx)
ax.set_ylim(0,ny)
ax.set_xticks(range(0,nx,25))
ax.set_yticks(range(0,ny,25))
ax.set_xlabel('Thinned x grid point')
ax.set_ylabel('Thinned y grid point')

# ---- use axes_grid toolkit to make colorbar axes.

divider = make_axes_locatable(ax1)
cax = divider.append_axes("bottom", size="3%", pad=0.25)
cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
    drawedges=True,ticks=clevs,format='%g')
cb.ax.tick_params(labelsize=7)
cb.set_label(r'D$_n$ statistic',fontsize=11)

# ---- set plot title

fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)

print ('Done!')









