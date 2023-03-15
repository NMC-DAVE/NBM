import numpy as np
import _pickle as cPickle
import sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.stats as stats

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


cmonth = sys.argv[1]
clead = sys.argv[2]
clon = sys.argv[3]
clat = sys.argv[4]
rlon = float(clon)
rlat = float(clat)
cclon = str(rlon-360)
cmonths = ['Jan','Feb','Mar','Apr','May','Jun',\
    'Jul','Aug','Sep','Oct','Nov','Dec']
imonth = int(cmonth)-1
ccmonth = cmonths[imonth]
ndaysomo = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]
ndaysomo_leap = [31, 28, 31,  30, 31, 30,  31, 31, 30,  31, 30, 31]

# ---- read Gamma parameter file

data_directory = '/Volumes/Backup Plus/mswep/'
infile = data_directory + cmonth+'_conus'+\
    '_MSWEP_apcp_gamma_parameters_h'+clead+'.cPick'
    
#infile = data_directory + cmonth+'_conus'+\
#    '_MSWEP_apcp_gamma_parameters_h'+clead+\
#    '_108_to_109.cPick'       
    
print ('reading from ', infile)
inf = open(infile, 'rb')
weights = cPickle.load(inf)
alpha = cPickle.load(inf)
beta = cPickle.load(inf)
fzero = cPickle.load(inf)
Dn = cPickle.load(inf)
nmixture = cPickle.load(inf)
ny, nx = np.shape(nmixture)
print (ny,nx)
inf.close()

# ---- read in the mswep lat/lon on ndfd grid

ncfile = data_directory + '200001_on_ndfd_grid_6hourly.nc'
print (ncfile)
nc = Dataset(ncfile)
lons_mswep = nc.variables['lons'][:,:]
lats_mswep = nc.variables['lats'][:,:]
print ('min, max lon = ', np.min(lons_mswep), np.max(lons_mswep))
ny, nx = np.shape(lats_mswep)
print (ny,nx)
nc.close()

# ---- make plots of Dnstats.   ymask out negatives

Dnstat_plot = (1-fzero)*Dn

title = r'MSWEP D$_n$ for '+\
    ccmonth+', 6-h period ending '+clead+' UTC'    
plot_title = 'Dnstat_mswep_'+ccmonth+'_'+clead+'UTC.png' 

clevs = [0.0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05]
colorst = ['White','#E4FFFF','#C4E8FF','#8FB3FF','#D8F9D8',\
    '#A6ECA6','#42F742','Yellow','Gold','Orange','#FCD5D9','#F6A3AE',\
    '#FA5257','Orchid','#AD8ADB','#A449FF','LightGray']
    
m = Basemap(llcrnrlon=233.7234,llcrnrlat=19.229,
    urcrnrlon = 300.95782, urcrnrlat = 54.37279,\
    projection='lcc',lat_1=25.,lat_2=25.,lon_0=265.,\
    resolution ='l',area_thresh=1000.)
x, y = m(lons_mswep, lats_mswep)
print ('x[0,0],x[0,-1],x[-1,0],x[-1,-1] = ', x[0,0],x[0,-1],x[-1,0],x[-1,-1])
print ('y[0,0],y[0,-1],y[-1,0],y[-1,-1] = ', y[0,0],y[0,-1],y[-1,0],y[-1,-1])    

fig = plt.figure(figsize=(8,6.5))
axloc = [0.08,0.1,0.87,0.88]
ax1 = fig.add_axes(axloc)
ax1.set_title(title, fontsize=16,color='Black')
CS2 = m.contourf(x, y, Dnstat_plot, clevs,\
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
m.drawparallels(parallels,labels=[1,0,0,0],color='LightGray',fontsize=8)
meridians = np.arange(220.,360.,5.)
m.drawmeridians(meridians,labels=[0,0,0,1],color='LightGray',fontsize=8)
   
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


# =====================================================================
# PART 2: plot CDFs, empirical and fitted, at selected lon,lat
# =====================================================================

# --- so, approximate to start the coordinates in the m() projection
#     x[0,0],x[0,-1],x[-1,0],x[-1,-1] =  0.4061070242896676 5952874.848155265 
#     -0.524982794187963 5952875.206909264
#     y[0,0],y[0,-1],y[-1,0],y[-1,-1] =  -0.08487691776826978 0.20302832452580333 
#     4053237.109193937 4053237.2167218113
# 

dx = 5952874.848155265/nx
dy = 4053237.109193937/ny
xcoords = np.arange(0, 5952874.848155265, dx)
ycoords = np.arange(0, 4053237.109193937, dy)
ixcoord = find_nearest(xcoords, xdot)
jycoord = find_nearest(ycoords, ydot)
cix = str(ixcoord)
cjy = str(jycoord)
print ('location of sample point jy, ix = ', jycoord,ixcoord)

# ---- read in the precip at this point

ktr = 0
for iyear in range(2000,2020):
    if iyear%4 == 0:
        ndays = ndaysomo_leap[imonth]
    else:
        ndays = ndaysomo[imonth]
    cyear = str(iyear)
    infile = data_directory + cyear + cmonth + '_on_ndfd_grid_6hourly.nc'
    print (iyear, infile)
    nc = Dataset(infile)
    if iyear == 2000:
        precip = np.squeeze(nc.variables['apcp_anal'][:,jycoord,ixcoord])
    else:
        precip = np.append(precip, np.squeeze(nc.variables['apcp_anal'][:,jycoord,ixcoord]) )  
    nc.close()
    
# ---- determine the number of samples, fraction zero, 
#      precip_nonzero vector, and number of zeros

nsamps = len(precip)
fraction_zero, precip_nonzero, nz = \
    fraczero_possamps(nsamps, precip)

# --- build CDFs, fitted and empirical

x = np.arange(0.0,50.01,0.1)
nx = len(x)
y0 = x / beta[0,jycoord,ixcoord]
y1 = x / beta[1,jycoord,ixcoord]
y2 = x / beta[2,jycoord,ixcoord]

cdf0 = stats.gamma.cdf(y0, alpha[0,jycoord,ixcoord])
cdf1 = stats.gamma.cdf(y1, alpha[1,jycoord,ixcoord])
cdf2 = stats.gamma.cdf(y2, alpha[2,jycoord,ixcoord])
cdf_fitted = fzero[jycoord,ixcoord] + (1.-fzero[jycoord,ixcoord])* \
    (weights[0,jycoord,ixcoord]*cdf0 + weights[1,jycoord,ixcoord]*cdf1 + \
    weights[2,jycoord,ixcoord]*cdf2) 

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
        
# ---- plot the CDFs, forecast and analyzed, empirical and best fitted.

f = plt.figure(figsize=(6.5,4.5))#
ax = f.add_axes([.13,.12,.84,.75])
ax.set_title(r'Fitted and empirical MSWEP CDFs, '+\
    ccmonth+',\n6-hour period ending '+clead+' UTC, '+\
    cclon+'$^{\circ}$ W '+clat+'$^{\circ}$ N',fontsize=14)
ax.plot(x,cdf_fitted,color='Blue',lw=2,\
    label='Fitted 3-component Gamma mixture')
ax.plot(x,cdf_empirical,color='Red',lw=2,label='Empirical')
plt.ylabel('Non-exceedance probability',fontsize=11)
ax.legend(loc=4)
#ax.set_ylim(0.7,1)
ax.set_ylim(0.0,1)
plt.grid(True,lw=0.25,color='LightGray')
ax.set_xlim(0,20)
#ax.set_xlim(0,5)
ax.set_xlabel('6-hourly total precipitation (mm)',fontsize=11)
figname = 'CDF_precip_example_MSWEP_'+ccmonth+'_'+\
    cclon+'W_'+clat+'N_'+clead+'UTC.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()






