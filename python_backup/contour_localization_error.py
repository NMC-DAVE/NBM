import csv
import sys
import numpy as np
from datetime import datetime
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams 
rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'

clead = sys.argv[1]
warm_or_cold = sys.argv[2]

cpath_errorstats = '/Volumes/Backup Plus/python/fcst_stats/'

rmse_array = np.zeros((4,7), dtype=np.float32)
bias_array = np.zeros((4,7), dtype=np.float32)

# ---- declare arrays

random_localizations = np.array([200.0, 400.0, 600.0, 800.0])
bias_localizations = np.array([400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0])
for ktrr, efold_random in enumerate(random_localizations):
    strr = str(efold_random)
    for ktrb, efold_bias in enumerate(bias_localizations):
        strb = str(efold_bias)
        if warm_or_cold == 'cold':
            statsfile = cpath_errorstats + '2018_KF_forecast_errorstats_'+\
                clead+'h_flocal'+strr+'_blocal'+strb+'.txt'
        else:
            statsfile = cpath_errorstats + '2018_KF_forecast_errorstats_'+\
                clead+'h_flocal'+strr+'_blocal'+strb+'_'+warm_or_cold+'.txt'
        values = np.loadtxt(statsfile)
        rmse_array[ktrr,ktrb] = values[2]
        bias_array[ktrr,ktrb] = values[0]

# ----- make plots

print ('min, max rmse = ', np.min(rmse_array), np.max(rmse_array))
print ('min, max bias = ', np.min(bias_array), np.max(bias_array))

fig = plt.figure(figsize=(6.5,3.4))

#clevels = [-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000'] 
fig.suptitle(clead+'-h forecasts, '+warm_or_cold+' season',fontsize=14)    
    
a1 = fig.add_axes([.11,.17,.37,.67])
a1.set_title('(a) RMSE',fontsize=13)
cf = a1.contour(bias_localizations,random_localizations,rmse_array,\
    colors = 'Black', linewidths=0.8)
a1.clabel(cf, inline=1, fontsize=9)
a1.set_ylabel('Random error localization\nlength scale (km)', fontsize=11)
a1.set_xlabel('Bias estimate localization\nlength scale (km)', fontsize=11)

a1 = fig.add_axes([.61,.17,.37,.67])
a1.set_title('(b) Bias',fontsize=13)
cf = a1.contour(bias_localizations,random_localizations,bias_array,\
    colors = 'Black', linewidths=0.8)
a1.clabel(cf, inline=1, fontsize=9)
a1.set_ylabel('Random error localization\nlength scale (km)', fontsize=11)
a1.set_xlabel('Bias estimate localization\nlength scale (km)', fontsize=11)

#cax = fig.add_axes([0.1,0.05,0.8,0.013])
#cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
#    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%",
#   pad="1%",
#cb.set_label('Bias (deg C)',fontsize=6)
#cb.ax.tick_params(labelsize=6)

plot_title = 'localization_rmse_bias_2018'+warm_or_cold+'_'+clead+'h.pdf'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title)
print ('Plot done')



        
        
