"""
variogram_bias_random.py

"""
import os, sys
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.optimize import curve_fit

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'

#@jit(nopython=True)
def exponential(dist, dist_lengthscale, asymptote):
    term = dist / dist_lengthscale
    #print ('np.shape(term) = ',np.shape(term))
    #print (term[0:-1:20])
    expterm = asymptote*(1.0 - np.exp(-term))
    return expterm

clead = sys.argv[1] # hours
                
# ---- read 2018 data

infile = 'covstats_bias_random_ecmwf2018_lead='+clead+'.cPick'
inf = open(infile,'rb')
print ('reading from ',infile)
cov_bias_2018 = cPickle.load(inf)
print (np.shape(cov_bias_2018))
cov_random_2018 = cPickle.load(inf)
print (np.shape(cov_random_2018))
var_bias_2018 = cPickle.load(inf)
var_random_2018 = cPickle.load(inf)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
nlats, nlons = np.shape(lats)
print ('nlats, nlons = ', nlats, nlons)
inf.close()              

# ---- read 2019 data

infile = 'covstats_bias_random_ecmwf2019_lead='+clead+'.cPick'
inf = open(infile,'rb')
print ('reading from ',infile)
cov_bias_2019 = cPickle.load(inf)
print (np.shape(cov_bias_2019))
cov_random_2019 = cPickle.load(inf)
var_bias_2019 = cPickle.load(inf)
var_random_2019 = cPickle.load(inf)
lats = cPickle.load(inf)
lons = cPickle.load(inf)
inf.close()  


figname = 'variogram_bias_'+clead+'.pdf'
f = plt.figure(figsize=(6.5,4.5))

ax = f.add_axes([.13,.11,.82,.78])
ax.set_title('Bias spatial correlation absolute difference, 2019 vs. 2018',fontsize=14)

# ---- loop thru representative sample of point

mean_diff = np.zeros(40,dtype=np.float64)
ktr_diff = np.zeros(40,dtype=np.float64)
mean_diff_list = []
dist_list = []
for ix1 in range(10,nlons-10,5):
    print ('ix1 = ',ix1,' of ',nlons)
    for jy1 in range(10,nlats-10,5):
        iran = np.random.randint(low=0,high=nlons-1, size=100)
        jran = np.random.randint(low=0,high=nlats-1, size=100)
        for k in range(100):
            corr19 = cov_bias_2019[jy1,ix1,jran[k],iran[k]] / \
                (np.sqrt(var_bias_2019[jy1,ix1]) * np.sqrt(var_bias_2019[jran[k],iran[k]]))
            corr18 = cov_bias_2018[jy1,ix1,jran[k],iran[k]] / \
                (np.sqrt(var_bias_2018[jy1,ix1]) * np.sqrt(var_bias_2018[jran[k],iran[k]]))
            diff_bias = np.abs(corr19-corr18)
            dist = np.sqrt(float(ix1-iran[k])**2 + float(jy1-jran[k])**2)
            distidx = np.int(dist)
            if distidx <40:
                mean_diff[distidx] = mean_diff[distidx] + diff_bias
                ktr_diff[distidx] = ktr_diff[distidx] + 1
                mean_diff_list.append(diff_bias)
                dist_list.append(dist)
                ax.plot(dist,diff_bias,marker='o',color='Red',\
                    markersize=0.2, markerfacecolor='Red')

        iran = np.random.randint(low=-9,high=9, size=100)
        jran = np.random.randint(low=-9,high=9, size=100)        
        for k in range(100):
            corr19 = cov_bias_2019[jy1,ix1,jy1-jran[k],ix1-iran[k]] / \
                (np.sqrt(var_bias_2019[jy1,ix1]) * np.sqrt(var_bias_2019[jy1-jran[k],ix1-iran[k]]))
            corr18 = cov_bias_2018[jy1,ix1,jy1-jran[k],ix1-iran[k]] / \
                (np.sqrt(var_bias_2018[jy1,ix1]) * np.sqrt(var_bias_2018[jy1-jran[k],ix1-iran[k]]))
            diff_bias = np.abs(corr19-corr18)
            dist = np.sqrt((ix1-(ix1-iran[k]))**2 + (jy1-(jy1-jran[k]))**2)
            distidx = np.int(dist)
            if distidx <40:
                mean_diff[distidx] = mean_diff[distidx] + diff_bias
                ktr_diff[distidx] = ktr_diff[distidx] + 1
                mean_diff_list.append(diff_bias)
                dist_list.append(dist)
                ax.plot(dist,diff_bias,marker='o',color='Red',\
                    markersize=0.2, markerfacecolor='Red')

nsamps = len(mean_diff_list)
mean_diff_arr = np.array(mean_diff_list)
dist_arr = np.array(dist_list)
bigdist = np.argwhere(dist_arr > 25.0)
fives = 0.05*np.ones((nsamps), dtype=np.float32)
print (np.shape(fives), np.shape(dist_arr))
sample_var = 1.0 - np.exp(-dist_arr[:] / 5.)
sample_var = np.where(sample_var < 0.05, fives, sample_var)
sample_var = np.ones((nsamps), dtype=np.float32)
print ('min, sample_var = ', np.min(sample_var), np.max(sample_var))
dist_lengthscale = 5.  # set to conservative values determined from Jul data.
asymptote = 0.3
popt, pcov = curve_fit(exponential, dist_arr, mean_diff_arr,  p0=(dist_lengthscale, asymptote), \
    check_finite=True,  method='trf', diff_step=0.000001, \
    sigma = sample_var, absolute_sigma=False, \
    bounds = ([0.0,0.0], [20.0,1.0]))
rho_bias = popt[0]
asymptote = popt[1]
print ('bias correlation length scale = ', rho_bias)

#mean_diff = mean_diff / ktr_diff  
#ax.plot(range(40), mean_diff,'k-',lw=3)   

ax.plot(range(40),asymptote*(1.0 - np.exp(-np.arange(40)/ rho_bias)),'k-',lw=3)             
ax.set_ylabel('Absolute difference of correlation',fontsize=13)
#ax.legend(loc=0)
ax.grid(True,lw=0.25)
ax.set_ylim(0,1)
ax.set_xlim(0,40)
ax.set_xticks([0,5,10,15,20,25,30,35,40])
ax.set_xlabel('Distance (grid points)',fontsize=13)

plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close()          

# ---- random error plot


figname = 'variogram_randerr'+clead+'.pdf'
f = plt.figure(figsize=(6.5,4.5))

ax = f.add_axes([.13,.11,.82,.77])
ax.set_title('Random error spatial correlation absolute difference,\n2019 vs. 2018',fontsize=13)

# ---- loop thru representative sample of point

mean_diff = np.zeros(40,dtype=np.float64)
ktr_diff = np.zeros(40,dtype=np.float64)
mean_diff_list = []
dist_list = []
for ix1 in range(10,nlons-10,5):
    print ('ix1 = ',ix1,' of ',nlons)
    for jy1 in range(10,nlats-10,5):
        iran = np.random.randint(low=0, high=nlons-1, size=100)
        jran = np.random.randint(low=0, high=nlats-1, size=100)
        for k in range(100):
            
            corr19 = cov_random_2019[jy1,ix1,jran[k],iran[k]] / \
                (np.sqrt(var_random_2019[jy1,ix1]) * np.sqrt(var_random_2019[jran[k],iran[k]]))
            corr18 = cov_random_2018[jy1,ix1,jran[k],iran[k]] / \
                (np.sqrt(var_random_2018[jy1,ix1]) * np.sqrt(var_random_2018[jran[k],iran[k]]))
            diff_random = np.abs(corr19-corr18)
            dist = np.sqrt(float(ix1-iran[k])**2 + float(jy1-jran[k])**2)
            distidx = np.int(dist)
            if distidx <40:
                mean_diff[distidx] = mean_diff[distidx] + diff_random
                ktr_diff[distidx] = ktr_diff[distidx] + 1
                ax.plot(dist,diff_random,marker='o',color='Red',\
                    markersize=0.2, markerfacecolor='Red')
                mean_diff_list.append(diff_random)
                dist_list.append(dist)
                
        iran = np.random.randint(low=-9,high=9, size=100)
        jran = np.random.randint(low=-9,high=9, size=100)
        for k in range(100):
            corr19 = cov_random_2019[jy1,ix1,jy1-jran[k],ix1-iran[k]] / \
                (np.sqrt(var_random_2019[jy1,ix1]) * np.sqrt(var_random_2019[jy1-jran[k],ix1-iran[k]]))
            corr18 = cov_random_2018[jy1,ix1,jy1-jran[k],ix1-iran[k]] / \
                (np.sqrt(var_random_2018[jy1,ix1]) * np.sqrt(var_random_2018[jy1-jran[k],ix1-iran[k]]))
            diff_random = np.abs(corr19-corr18)
            dist = np.sqrt((ix1-(ix1-iran[k]))**2 + (jy1-(jy1-jran[k]))**2)
            distidx = np.int(dist)
            if distidx <40:
                mean_diff[distidx] = mean_diff[distidx] + diff_random
                ktr_diff[distidx] = ktr_diff[distidx] + 1
                ax.plot(dist,diff_random,marker='o',color='Red',\
                    markersize=0.2, markerfacecolor='Red')
                mean_diff_list.append(diff_random)
                dist_list.append(dist)
                
nsamps = len(mean_diff_list)
mean_diff_arr = np.array(mean_diff_list)
dist_arr = np.array(dist_list)
fives = 0.05*np.ones((nsamps), dtype=np.float32)
sample_var = 1.0 - np.exp(-dist_arr[:] / 5.)
sample_var = np.where(sample_var < 0.05, fives, sample_var)
sample_var = np.ones((nsamps), dtype=np.float32)
dist_lengthscale = 3.  # set to conservative values 
asymptote = 0.3
popt, pcov = curve_fit(exponential, dist_arr,  mean_diff_arr, p0=(dist_lengthscale, asymptote), \
    check_finite=True,  method='trf', diff_step=0.000001, \
    sigma = sample_var, absolute_sigma=False, \
    bounds = ([0.0,0.0], [20.0,1.0]))
rho_random = popt[0]
asymptote = popt[1]
print ('random correlation length scale = ', rho_random)
            
#mean_diff = mean_diff / ktr_diff  
#ax.plot(range(40), mean_diff,'k-',lw=3)  
ax.plot(range(40),asymptote*(1.0 - np.exp(-np.arange(40)/ rho_random)),'k-',lw=3)                
ax.set_ylabel('Absolute difference of correlation',fontsize=13)
#ax.legend(loc=0)
ax.grid(True,lw=0.25)
ax.set_ylim(0,0.3)
ax.set_xlim(0,40)
ax.set_xticks([0,5,10,15,20,25,30,35,40])
ax.set_xlabel('Distance (grid points)',fontsize=13)

plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 

# --- save correlation length scales to file

outfile = 'KF_localization_lengthscales_lead='+clead+'.cPick'
ouf = open(outfile, 'wb')
print (outfile)
cPickle.dump(rho_bias, ouf)
cPickle.dump(rho_random, ouf)
ouf.close()
