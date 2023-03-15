"""
plot_fitted_precip_cdfs.py cmonth clead jy ix

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
rcParams['legend.fontsize']='medium'

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

cmonth = sys.argv[1] # 'Jan', 'Feb', etc.
clead = sys.argv[2] # 03, 06, 12, etc.
cjy = sys.argv[3] # units of 0, 5, etc
cix = sys.argv[4] # units of 0, 5, etc
jy = int(cjy)
ix = int(cix)
#ixd5 = ix//5
#jyd5 = jy//5


master_directory = '/Volumes/Backup Plus/gefsv12/precip/netcdf/'

# ---- read in the previously generated netCDF file with precipitation
#      for this month and lead time.  All members, dates for this 
#      month have been smushed into one leading index, dimension
#      nsamps, since the date of the forecast within the month and 
#      the member number is irrelevant for the distribution fitting.

jmin = 93
jmax = 246
imin = 368
imax = 686
        
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
ncfile = master_directory + cmonth + '_apcp' '_h' + clead + '.nc'
print (ncfile)
nc = Dataset(ncfile)
lons_1d = nc.variables['lons_fcst'][imin:imax]
lats_1d = nc.variables['lats_fcst'][jmin:jmax]
nlats = len(lats_1d)
nlons = len(lons_1d)
print ('nlats, lats[0], lats[-1] = ', nlats, lats_1d[0], lats_1d[-1])
print ('nlons, lons[0], lons[-1] = ', nlons, lons_1d[0], lons_1d[-1])
print ('sample lat, lons = ', lats_1d[jy], lons_1d[ix])
print ('nlats, nlons = ', nlats,nlons)
precip_full = nc.variables['apcp_fcst'][:,jmin:jmax,imin:imax]
precip = precip_full[:,jy,ix]
nsamps = len(precip)
teeny_precip = 0.06*np.ones(nsamps) 
precip = precip - teeny_precip
nc.close()
clat = str(lats_1d[jy])
clon = str(lons_1d[ix])

# --- determine the empirical fraction_zero, the sorted precip_nonzero

fraction_zero, precip_nonzero, nz = fraczero_possamps(nsamps, precip)

# --- load from cPickle file

#infile = master_directory + cmonth+ '_apcp_gamma_parameters_v3_h' + clead + '.cPick'
#infile = master_directory + cmonth+ '_apcp_gamma_parameters_v3_eps003_maxit60_h' + clead + '.cPick'

cdomain = 'conus'
infile = master_directory + cmonth+'_'+cdomain+'_apcp_gamma_parameters_h' + clead + '.cPick'


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
inf.close()

nstride = 1
lats_thin = lats_1d[0::nstride]
lons_thin = lons_1d[0::nstride]
weights_thin = weights[0::nstride, 0::nstride]
alpha_thin = alpha[0::nstride, 0::nstride]
beta_thin = beta[0::nstride, 0::nstride]
fzero_thin = fzero[0::nstride, 0::nstride]
Dnstat1_thin = Dnstat1[0::nstride, 0::nstride]
Dnstat2_thin = Dnstat2[0::nstride, 0::nstride]
Dnstat3_thin = Dnstat3[0::nstride, 0::nstride]
nmixture_thin = nmixture[0::nstride, 0::nstride]
fzero_thin = fzero[0::nstride, 0::nstride]
lons_thin_2D, lats_thin_2D = np.meshgrid(lons_thin, lats_thin)

# --- build CDF

x = np.arange(0.0,50.01,0.1)
nx = len(x)
y0 = x / beta[0,jy,ix]
y1 = x / beta[1,jy,ix]
y2 = x / beta[2,jy,ix]
#y3 = x / beta[3,jy,ix]

print ('weights[:,jy,ix] = ', weights[:,jy,ix])
print ('alpha[:,jy,ix] = ', alpha[:,jy,ix])
print ('beta[:,jy,ix] = ', beta[:,jy,ix])

cdf0 = stats.gamma.cdf(y0, alpha[0,jy,ix])
cdf1 = stats.gamma.cdf(y1, alpha[1,jy,ix])
cdf2 = stats.gamma.cdf(y2, alpha[2,jy,ix])
#cdf3 = stats.gamma.cdf(y3, alpha[3,jy,ix])
cdf_fitted = fzero[jy,ix] + (1.-fzero[jy,ix])* \
    (weights[0,jy,ix]*cdf0 + weights[1,jy,ix]*cdf1 + \
    weights[2,jy,ix]*cdf2) 

cdf_empirical = np.zeros((nx),dtype=np.float64)
len_nonzero = len(precip_nonzero)
fnzero = float(len_nonzero)
nz = int(fnzero)
pnzsort = np.sort(precip_nonzero)
for i, xi in enumerate(x):
    nbelow = (pnzsort < xi).sum() 
    cdf_empirical[i] = fraction_zero + (1.-fraction_zero)*float(nbelow)/fnzero
    
    
    
query_these_indices = [ nz//20, nz//10, (3*nz)//20, nz//5, nz//4, (3*nz)//10, \
    (7*nz)//20, (2*nz)//5, (9*nz)//20, nz//2, (11*nz)//20, (3*nz)//5, (13*nz)//20, \
    (7*nz)//10, (3*nz)//4, (4*nz)//5, (17*nz)//20, (35*nz)//40,(9*nz)//10, \
    (37*nz)//40, (19*nz)//20, (39*nz)//40, (79*nz)//80]
        
print ('   nz, query_these_indices = ', nz, query_these_indices)
empirical_CDF = fzero[jy,ix] + (1.-fzero[jy,ix])*np.array(query_these_indices, dtype=np.float) / float(nz) 
print ('   empirical_CDF = ',empirical_CDF)
empirical_precipvals = pnzsort[query_these_indices]     
y0 = empirical_precipvals / beta[0,jy,ix]
y1 = empirical_precipvals / beta[1,jy,ix]
y2 = empirical_precipvals / beta[2,jy,ix]
fitted_CDF0 = stats.gamma.cdf(y0, alpha[0,jy,ix])
fitted_CDF1 = stats.gamma.cdf(y1, alpha[1,jy,ix])
fitted_CDF2 = stats.gamma.cdf(y2, alpha[2,jy,ix])
fitted_CDF = fzero[jy,ix] + (1.-fzero[jy,ix])*(weights[0,jy,ix]*fitted_CDF0 + \
    weights[1,jy,ix]*fitted_CDF1 + weights[2,jy,ix]*fitted_CDF2)
print ('fitted_CDF every 0.05 = ', fitted_CDF)            
        
# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,4.5))#
ax = f.add_axes([.13,.14,.84,.77])
ax.set_title(r'Fitted and empirical GEFSv12 CDFs, '+cmonth+', '+clead+'-h forecast, '+\
    clon+'$^{\circ}$ W '+clat+'$^{\circ}$ N',fontsize=11)
ax.plot(x,cdf_fitted,color='Blue',lw=2,label='Fitted 3-component Gamma mixture')
#ax.plot(x,fzero[jy,ix]+(1.-fzero[jy,ix])*cdf0,color='Black',\
#    lw=1,label='Gamma mixture 0, weight = '+str(weights[0,jy,ix]))
#ax.plot(x,fzero[jy,ix]+(1.-fzero[jy,ix])*cdf1,color='Gray',\
#    lw=1,label='Gamma mixture 1, weight = '+str(weights[1,jy,ix]))
#ax.plot(x,fzero[jy,ix]+(1.-fzero[jy,ix])*cdf2,color='LightGray',\
#    lw=1,label='Gamma mixture 2, weight = '+str(weights[2,jy,ix]))
ax.plot(x,cdf_empirical,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.7,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,20)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()


# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,4.5))#
ax = f.add_axes([.15,.14,.8,.77])
ax.set_title(r'Fitted and empirical CDFs every 0.05, '+clon+'$^{\circ}$ W '+clat+'$^{\circ}$ N',fontsize=11)
ax.plot(empirical_precipvals,fitted_CDF,color='Blue',lw=2,label='Fitted 3-component Gamma mixture')

ax.plot(empirical_precipvals,empirical_CDF,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=0)
ax.set_ylim(0.7,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,20)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_v2.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()

# ---- Q-Q plot of forecast data

f = plt.figure(figsize=(5.5,5.5))
ax = f.add_axes([.15,.14,.8,.77])
ax.set_title('Q-Q plot, '+clon+' '+clat, fontsize=13)
ax.plot([0,1],[0,1],color='Gray',lw=0.5)
ax.plot(cdf_empirical,cdf_fitted,color='Red',lw=2)
ax.set_ylim(0.7,1)
ax.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0.7,1)
ax.set_xlabel('Empirical forecast quantile',fontsize=11)
ax.set_ylabel('Fitted forecast quantile',fontsize=11)

figname = 'QQplot_example.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()



m = Basemap(llcrnrlon=lons_thin_2D[0,0],llcrnrlat=lats_thin_2D[-1,-1],\
    urcrnrlon=lons_thin_2D[-1,-1],urcrnrlat=lats_thin_2D[0,0],\
    resolution='l', projection='mill')
x, y = m(lons_thin_2D, lats_thin_2D)
nj, ny, nx = np.shape(alpha)

# ---- make plots of Dnstats.   ymask out negatives

for idn in range(4):
    if idn == 0:
        Dnstat1_thin = ma.masked_where(Dnstat1_thin <=0, Dnstat1_thin)
        Dnstat_plot = (1-fzero_thin)*Dnstat1_thin
        title = r'D$_n$ for single Gamma'
        plot_title = 'Dnstat1_eps003_maxit40.png'
    elif idn == 1:
        Dnstat2_thin = ma.masked_where(Dnstat2_thin <=0, Dnstat2_thin)
        Dnstat_plot = (1-fzero_thin)*Dnstat2_thin
        title = r'D$_n$ for mixture of two Gammas'
        plot_title = 'Dnstat2_eps003_maxit40.png'
    elif idn == 2:
        Dnstat3_thin = ma.masked_where(Dnstat3_thin <=0, Dnstat3_thin)
        Dnstat_plot = (1-fzero_thin)*Dnstat3_thin
        title = r'D$_n$ for mixture of three Gammas'
        plot_title = 'Dnstat3_eps003_maxit40.png'
    elif idn == 3:
        d1 = ma.minimum((1-fzero_thin)*Dnstat1_thin, (1-fzero_thin)*Dnstat2_thin)
        Dnstat_plot = ma.minimum(d1,(1-fzero_thin)*Dnstat3_thin)
        print ('jy, ix, ny, nx', jy, ix, np.shape(Dnstat2_thin))
        print ('Dnstat1 2 3 o = ', (1.-fzero_thin[jy,ix])*Dnstat1_thin[jy,ix], \
            (1.-fzero_thin[jy,ix])*Dnstat2_thin[jy,ix],\
            (1.-fzero_thin[jy,ix])*Dnstat3_thin[jy,ix],Dnstat_plot[jy,ix])
        
        title = r'D$_n$ overall'    
        plot_title = 'Dnstat__eps003_maxit40_overall.png' 

    clevs = [0.0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05]
    colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
        '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
        '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']

    fig = plt.figure(figsize=(9,6.))
    axloc = [0.06,0.11,0.91,0.85]
    ax1 = fig.add_axes(axloc)
    ax1.set_title(title, fontsize=11,color='Black')
    CS2 = m.contourf(x, y, Dnstat_plot, clevs,\
        cmap=None, colors=colorst, extend='both')
    m.drawcoastlines(linewidth=0.8,color='Gray')
    m.drawcountries(linewidth=0.8,color='Gray')
    m.drawstates(linewidth=0.8,color='Gray')
    
    xdot, ydot = m(lons_thin_2D[jy,ix],lats_thin_2D[jy,ix])
    m.plot(xdot,ydot,marker='.',markersize=5,color='Black')
    
    # draw parallels and meridians.
    # label on left and bottom of map.
    parallels = np.arange(-20.,80,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],color='LightGray')
    meridians = np.arange(0.,360.,20.)
    m.drawmeridians(meridians,labels=[0,0,0,1],color='LightGray')
    
        
    ax.set_xlim(0,nx)
    ax.set_ylim(0,ny)
    ax.set_xticks(range(0,nx,25))
    ax.set_yticks(range(0,ny,25))
    ax.set_xlabel('x grid point')
    ax.set_ylabel('y grid point')

    # ---- use axes_grid toolkit to make colorbar axes.

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("bottom", size="3%", pad=0.25)
    cb = plt.colorbar(CS2,orientation='horizontal',cax=cax,\
        drawedges=True,ticks=clevs,format='%g')
    cb.ax.tick_params(labelsize=7)
    cb.set_label(r'D$_n$ statistic',fontsize=9)

    # ---- set plot title

    fig.savefig(plot_title, dpi=300)
    print ('saving plot to file = ',plot_title)

print ('Done!')









