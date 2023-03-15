"""
python plot_dressing_example.py

"""

import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import norm

import _pickle as cPickle
import scipy.stats as stats

rcParams['xtick.labelsize']='medium'
rcParams['ytick.labelsize']='medium'
rcParams['legend.fontsize']='medium'

precip = np.array([0.2, 0.5, 0.8, 1.0, 1.25, 1.7, 2.0, 2.6, 3.3, 3.9])
sigma_p = 0.25 + 0.15*precip
x = np.arange(1001)/100.

f = plt.figure(figsize=(6.5,3.))
ax = f.add_axes([.11,.17,.87,.73])
ax.set_title('Dressing example',fontsize=13)
for i in range(10):
    ax.plot(x,norm.pdf(x,loc=precip[i],scale=sigma_p[i]),color='LightGray',lw=0.5)
    ax.plot([precip[i],precip[i]],[0,0],'r.')
ax.set_xlim(0,10)
ax.set_ylabel('Probability density',fontsize=11)
ax.set_xlabel('Precipitation (mm)',fontsize=11)
figname = 'dressing.png'
plt.savefig(figname,dpi=400)
print ('Plot done', figname)
plt.close() 
