#!/usr/bin/python
#
# This script reads in the closest-member histograms for the reforecast ensemble,
# plots the data.  It then fits a beta distribution to the histogram, using
# this to estimate the closest-member histogram for a larger, 51-member ensemble,
# plotting out these estimates.

import numpy as np
import sys
import os
import scipy.stats as st
import scipy.special as sp
import time as time2
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'


cp = '3.5' # cp = sys.argv[1]
cq = '1.0' # sys.argv[2]
cknots = sys.argv[1] # number of interior knots
rp = float(cp)
rq = float(cq)
nknots = int(cknots)

gamma_p_plus_q = -sp.gamma(rp+rq)
gamma_p = -sp.gamma(rp)
gamma_q = -sp.gamma(rq)
#print ('gamma_p_plus_q, gamma_p, gamma_q = ',\
#    gamma_p_plus_q, gamma_p, gamma_q)

xbeta = 0.0005 + np.arange(1000)/1000.
nxbeta = len(xbeta)
beta_distn = np.zeros((nxbeta),dtype=np.float32)

for i in range(nxbeta):
    beta_distn[i] = (gamma_p_plus_q / (gamma_p + gamma_q)) * \
        xbeta[i]**(rp-1.0) * (1.0 - xbeta[i])**(rq-1.0)

#print (beta_distn[0:-1:50])

# ---- make a plot of the quantile-mapped closest-member histograms and then overlay
#      the fitted beta distribution.

fig = plt.figure(figsize=(9.,3.8))
ctitle = r'Fitted beta distribution and knot locations, $\alpha$ =  '+cp+r' , $\beta$ = '+cq

ax = [0.1,0.15,0.85,0.76]


# ---- plot the closest-member histograms

a2 = fig.add_axes(ax,zorder=2)
a2.set_title(ctitle,fontsize=15)
#a2.plot(xp_mid, chist,color='Red',\
#    linewidth=1, alpha=1,linestyle='-')
#a2.vlines(xp_mid - (xp_mid[1]-xp_mid[0])/10., 0.0*chist, \
#    chist, colors='r', linestyles='solid',\
#    label='Closest-member histogram',lw=2)
#a2.vlines(xp_mid + (xp_mid[1]-xp_mid[0])/10., 0.0*chist_beta, \
#    chist_beta, colors='b', linestyles='solid',\
#    label='Beta distribution integral',lw=2)

# ---- overplot the fitted beta distribution
a2.plot(xbeta, beta_distn,color='Red',\
    linewidth=3, linestyle='-')
a2.set_xlim(0,1)
a2.set_ylim(0,3)
a2.set_xticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
a2.set_ylabel('Probability density',fontsize=12)
a2.set_xlabel('z',fontsize=12)
a2.grid(color='LightGray',linewidth=0.8)

print ()
for iknot in range (1,nknots+1):
    rknot = float(iknot)/(nknots+1)
    #rcdf = st.beta.cdf(rknot, rp, rq)
    xloc = st.beta.ppf(rknot, rp, rq, loc=0, scale=1)
    print ('rknot, xloc = ',rknot, xloc )
    a2.vlines(xloc, 0.0, 3.0, colors='Red', linestyles='dashed',lw=1)


#a2.legend(loc=9)

plot_title = 'beta.png'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title, dpi=400)
print ('Plot done')

