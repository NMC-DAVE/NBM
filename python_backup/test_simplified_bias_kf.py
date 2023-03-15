""" apply Kalman filter bias correction to forecasts.  Note
        the mix of some arrays shapes; Kalman gain is shaped
        (nlats*nlons, nlats*nlons)
"""
import numpy as np
import sys
import os, sys
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='large'
    
R = np.zeros((100,100),dtype=np.float32)
for k in range(100):
    R[k,k] = 1.0
print ('R[0:2,0:2] = ',R[0:2,0:2] )
Bx1 = np.zeros((100,100), dtype=np.float32)
Bx4 = np.zeros((100,100), dtype=np.float32)
Bx16 = np.zeros((100,100), dtype=np.float32)
Bbeta = np.zeros((100,100), dtype=np.float32)
Beta_variance = np.zeros((100), dtype=np.float32)
for i in range(100):
    dist = np.amin(np.array([np.abs(30-i), np.abs(130-i), np.abs(-70-i)]),axis=0)
    Beta_variance[i] = 0.4 # 0.4*np.exp(-dist**2/2500.)
for i in range(100):
    for j in range(100):
        dist = np.amin(np.array([np.abs(j-i), np.abs(j+100-i), np.abs(j-100-i)]),axis=0)
        Bx1[j,i] = np.exp(-dist**2/100.)
        Bx4[j,i] = 4.*np.exp(-dist**2/100.)
        Bx16[j,i] = 16*np.exp(-dist**2/100.)
        Bbeta[j,i] = np.sqrt(Beta_variance[i])*np.sqrt(Beta_variance[j])*np.exp(-dist**2/2500.)

Bbeta_plus_Bx1_plus_R = R + Bx1 + Bbeta
Bbeta_plus_Bx4_plus_R = R + Bx4 + Bbeta
Bbeta_plus_Bx16_plus_R = R + Bx16 + Bbeta

Bbeta_plus_Bx1_plus_R_inv = \
    np.linalg.inv(Bbeta_plus_Bx1_plus_R)
Bbeta_plus_Bx4_plus_R_inv = \
    np.linalg.inv(Bbeta_plus_Bx4_plus_R)
Bbeta_plus_Bx16_plus_R_inv = \
    np.linalg.inv(Bbeta_plus_Bx16_plus_R)
Kalman_gain_beta1 = np.matmul(Bbeta, Bbeta_plus_Bx1_plus_R_inv)
Kalman_gain_beta4 = np.matmul(Bbeta, Bbeta_plus_Bx4_plus_R_inv)
Kalman_gain_beta16 = np.matmul(Bbeta, Bbeta_plus_Bx16_plus_R_inv)
       
beta_3d1 = np.zeros((100), dtype=np.float32)
beta_3d4 = np.zeros((100), dtype=np.float32)
beta_3d16 = np.zeros((100), dtype=np.float32)
obsinc = np.ones((100), dtype=np.float32)
    
for i in range(100):
    beta_3d1[i] = np.sum(Kalman_gain_beta1[i,:]*obsinc[:])
    beta_3d4[i] = np.sum(Kalman_gain_beta4[i,:]*obsinc[:])
    beta_3d16[i] = np.sum(Kalman_gain_beta16[i,:]*obsinc[:])
            
# ---- plot beta_3d[1,:]
    
fig = plt.figure(figsize=(6.5,8.5))

axloc = [0.12,0.72,0.8,0.23]
ax1 = fig.add_axes(axloc)
title = r'(a) Forecast-error covariance for i=50, $\sigma^b_x$ = 1 '
ax1.set_title(title, fontsize=16,color='Black')
ax1.plot(range(100), Bx1[50,:], 'r-', lw=2)
ax1.set_xlim(0,99)
ax1.set_xlabel('grid point number (i)')
ax1.set_ylim(0,1.1)
ax1.set_ylabel('Forecast-error covariance')
ax1.grid(True,lw=0.25,color='LightGray')

axloc = [0.12,0.39,0.8,0.23]
ax1 = fig.add_axes(axloc)
title = '(b) Bias estimate error covariance for i=50'
ax1.set_title(title, fontsize=16,color='Black')
ax1.plot(range(100), Bbeta[50,:], 'r-', lw=2)
ax1.set_xlim(0,99)
ax1.set_xlabel('grid point number (i)')
ax1.set_ylim(0,0.3)
ax1.set_ylabel('Bias error covariance')
ax1.grid(True,lw=0.25,color='LightGray')

axloc = [0.12,0.06,0.8,0.23]
ax1 = fig.add_axes(axloc)
title = '(c) Bias Kalman gain at i=50 for R=1'
ax1.set_title(title, fontsize=16,color='Black')
ax1.plot(range(100), Kalman_gain_beta1[:,50], 'r-', lw=2,label=r'$\sigma^b_x$ = 1')
ax1.plot(range(100), Kalman_gain_beta4[:,50], 'b-', lw=2,label=r'$\sigma^b_x$ = 4')
ax1.plot(range(100), Kalman_gain_beta16[:,50], 'g-', lw=2,label=r'$\sigma^b_x$ = 16')
ax1.set_xlim(0,99)
ax1.set_xlabel('grid point number (i)')
ax1.set_ylim(0,0.01)
ax1.set_ylabel('Bias Kalman gain')
ax1.grid(True,lw=0.25,color='LightGray')
ax1.legend(loc=0)

# ---- set plot title

plot_title = 'kalman_gain_bias_example.pdf'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')


 