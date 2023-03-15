import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['legend.fontsize']='x-small'

idate = range(365) # Julian day number of simulation
bias = 2.0 # systematic error in the ECMWF system
phi_truth = 0.5 # autocorrelation of truth
phi_ecmwf = 0.7 # autocorrelation of ECMWF error
R = 0.15 # observation error variance, roughly 1/10th climatological 
random_ecmwf_mag = 0.5  # scaling factor for ECMWF random error
random_truth_mag = 0.7  # scaling factor for truth

# random error component of ECMWF synthetic analyses, N(0,1)
random_ecmwf = np.random.randn(365) 
# random error component of truth, N(0,1)
random_truth = np.random.randn(365) 
print ('random_truth[0:10] = ', random_truth[0:10])

truth = np.zeros((365), dtype=np.float32) # initialize vector to zero
# a climatological background, set to constant zero 
background = np.zeros((365), dtype=np.float32)  
# initialize a vector of autocorrelated ECMWF analysis error
autocorr_ecmwf = np.zeros((365), dtype=np.float32)

for i in range(1,365):
    # true state is AR1 time series
    truth[i] = phi_truth*truth[i-1] + random_truth[i]
    # autocorrelated ECMWF error, bias + random ecmwf error 
    #     + fraction of truth random component
    autocorr_ecmwf[i] = phi_ecmwf*autocorr_ecmwf[i-1] + \
        random_ecmwf_mag*random_ecmwf[i] + \
        random_truth_mag*random_truth[i] 
        
autocorr_ecmwf = autocorr_ecmwf + bias
# set the background error variance of the truth 
# relative to the constant (zero) background
B = np.var(truth)
# used to separate out the ECMWF error component that is time varying vs. bias
autocorr_mean_ecmwf = np.mean(autocorr_ecmwf) 
# simulate two observation time series, independent of each other
obs1 = np.sqrt(R)*np.random.randn(365) + truth
obs2 = np.sqrt(R)*np.random.randn(365) + truth

print ('np.var(truth) = ', B)
# assimilation of one observation
state_estimate_assim1 = np.zeros((365), dtype=np.float32)
state_estimate_assim2 = np.zeros((365), dtype=np.float32)

for i in range(365):
    # Kalman gain for first observation
    K1 = B / (B+R)
    # analysis error variance after assimilation of first observation
    Pa1 = B - K1*B
    # state estimate after assimilation of first observation
    state_estimate_assim1[i] = background[i] + K1*(obs1[i]-background[i])
    # Kalman gain for assimilation of second observation
    K2 = Pa1 / (Pa1+R)
    # state estimate after assimilation of second observation
    state_estimate_assim2[i] = state_estimate_assim1[i] + \
        K2*(obs2[i]-state_estimate_assim1[i])
    Pa2 = Pa1 - K2*Pa1
    
# ---- make 4-panel plot

fig = plt.figure(figsize=(6.5,9.0))

axloc = [0.11,0.83,0.85,0.14]
ax = fig.add_axes(axloc)
title = '(a) Truth, ECMWF 4D-Var, OI analyses'
ax.set_title(title, fontsize=11,color='Black')
ax.plot(idate,truth,'r-',lw=0.6,label='Truth')
ax.plot(idate,autocorr_ecmwf,'b-',lw=0.6,label='ECMWF 4D-Var')
ax.plot(idate,state_estimate_assim1,'g-',lw=0.6,label='OI, 1 obs')
ax.plot(idate,state_estimate_assim2,'-',color='Purple',lw=0.6,label='OI, 2 obs')
ax.set_ylabel('Temperature')
ax.legend(loc=0)
ax.set_xlim(0,365)
ax.set_ylim(-3,7)
ax.grid(True,lw=0.25,color='LightGray')

axloc = [0.11,0.63,0.85,0.14]
ax = fig.add_axes(axloc)
title = r'(b) $\delta_t^{\prime}$ relative to OI with 1 observation'
ax.set_title(title, fontsize=11,color='Black')
ax.plot(idate,truth,'r-',lw=0.6,label='Truth')
ax.plot(idate,autocorr_ecmwf-autocorr_mean_ecmwf,'b-',\
    lw=0.6,label='ECMWF 4D-Var random component')
ax.plot(idate,state_estimate_assim1,'g-',lw=0.6,label='OI, 1 obs')
ax.set_ylabel('Temperature')
for i in range(365):
    ax.plot([i,i],[autocorr_ecmwf[i]-autocorr_mean_ecmwf, \
        state_estimate_assim1[i]],color='Black',lw=1.3)
ax.legend(loc=0)
ax.set_xlim(0,365)
ax.set_ylim(-5,5)
ax.grid(True,lw=0.25,color='LightGray')

axloc = [0.11,0.43,0.85,0.14]
ax = fig.add_axes(axloc)
title = r'(c) $\delta_t^{\prime}$ relative to OI with 2 observations'
ax.set_title(title, fontsize=11,color='Black')
ax.plot(idate,truth,'r-',lw=0.6,label='Truth')
ax.plot(idate,autocorr_ecmwf-autocorr_mean_ecmwf,'b-',\
    lw=0.6,label=r'ECMWF 4D-Var random component')
ax.plot(idate,state_estimate_assim2,'-',color='Purple',\
    lw=0.6,label='OI, 2 obs')
ax.set_ylabel('Temperature')
for i in range(365):
    ax.plot([i,i],[autocorr_ecmwf[i]-autocorr_mean_ecmwf, \
        state_estimate_assim2[i]],color='Black',lw=1.3)
ax.legend(loc=0)
ax.set_xlim(0,365)
ax.set_ylim(-5,5)
ax.grid(True,lw=0.25,color='LightGray')

axloc = [0.11,0.23,0.85,0.14]
ax = fig.add_axes(axloc)
title = r'(d) OI error with 1 observation'
ax.set_title(title, fontsize=11,color='Black')
ax.plot(idate,truth,'r-',lw=0.6,label='Truth')
ax.plot(idate,state_estimate_assim1,'g-',lw=0.6,\
    label='OI, assimilate 1 obs')
ax.set_ylabel('Temperature')
for i in range(365):
    ax.plot([i,i],[truth[i], state_estimate_assim1[i]],\
        color='Black',lw=1.3)
ax.legend(loc=0)
ax.set_xlim(0,365)
ax.set_ylim(-5,5)
ax.grid(True,lw=0.25,color='LightGray')

axloc = [0.11,0.03,0.85,0.14]
ax = fig.add_axes(axloc)
title = r'(d) OI error with 2 observations'
ax.set_title(title, fontsize=11,color='Black')
ax.plot(idate,truth,'r-',lw=0.6,label='Truth')
ax.plot(idate,state_estimate_assim2,'-',color='Purple',\
    lw=0.6,label='OI, assimilate 2 obs')
for i in range(365):
    ax.plot([i,i],[truth[i], state_estimate_assim2[i]],\
        color='Black',lw=1.3)
ax.legend(loc=0)
ax.set_ylabel('Temperature')
ax.set_xlim(0,365)
ax.set_xlabel('Day number')
ax.set_ylim(-5,5)
ax.grid(True,lw=0.25,color='LightGray')

# ---- set plot title

plot_title = 'demo_sfcdata.png'
fig.savefig(plot_title,dpi=400)
print ('saving plot to file = ',plot_title)
print ('Done!')


print ('Var(truth - state_estimate_assim1) = ', \
    np.var(truth - state_estimate_assim1), Pa1)
print ('Var(truth - state_estimate_assim2) = ', \
    np.var(truth - state_estimate_assim2), Pa2)
print ('Var(autocorr_ecmwf -  state_estimate_assim1) = ', \
    np.var((autocorr_ecmwf-autocorr_mean_ecmwf) -  state_estimate_assim1))
print ('Var(autocorr_ecmwf -  state_estimate_assim2) = ', \
    np.var((autocorr_ecmwf-autocorr_mean_ecmwf) -  state_estimate_assim2))
        
    
    