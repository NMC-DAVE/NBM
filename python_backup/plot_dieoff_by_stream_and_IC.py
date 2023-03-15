import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap, interp
from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib import rcParams
rcParams['xtick.labelsize']='small'
rcParams['ytick.labelsize']='small'
rcParams['legend.fontsize']='small'

infiles = ['sherrie/avg_cfsr_P500.txt','sherrie/avg_reanal_P250.txt',\
    'sherrie/avg_reanal_P850.txt','sherrie/avg_cfsr_P250.txt',\
    'sherrie/avg_cfsr_P850.txt','sherrie/avg_reanal_P500.txt']

AC_CFSR = np.zeros((6,3,4,5)) # lead time, level, region (GL, NH, TR, SH), stream
AC_GEFSv12 = np.zeros((6,3,4,5)) # lead time, level, region (GL, NH, TR, SH), stream
level_dict = {250:0, 500:1, 850:2}
lead_dict = {6:0, 24:1, 48:2, 72:3, 96:4, 120:5}
cregions = ['N. Hem. ','Tropics ','S. Hem. ']
clevels = ['250 hPa ','500 hPa ','850 hPa ']
cletters = ['(a) ','(b) ','(c) ','(d) ','(e) ','(f) ','(g) ','(h) ','(i) ']

# ----- read in anomaly correlation from text files, assign to 
#       array elements

for file in infiles:
    indata = np.loadtxt(file,delimiter=None)
    level = indata[:,0].astype(int)
    lead = indata[:,1].astype(int)
    streamno = indata[:,2].astype(int)
    region = indata[:,3].astype(int)
    correlation = indata[:,4]
    nsamps = len(level)
    ICtype = file[12:13]
    if ICtype == 'c':
        for isamp in range(nsamps):
            AC_CFSR[lead_dict[lead[isamp]],level_dict[level[isamp]],\
                region[isamp],streamno[isamp]] = correlation[isamp]
    else:
        for isamp in range(nsamps):
            AC_GEFSv12[lead_dict[lead[isamp]],level_dict[level[isamp]],\
                region[isamp], streamno[isamp]] = correlation[isamp]
                
# ---- make 9-panel plot

xstart = [0.09, 0.41, 0.73]
xlen = [0.24, 0.24, 0.24]
ystart = [0.72, 0.4, 0.08]
ylen = [0.22, 0.22, 0.22]


fig = plt.figure(figsize=(6.5,6.5))

ktr = 0
for jlevel in range(3):
    for ix in range(3):
        region = ix+1 # NH, TR, SH
        CFSR_data = AC_CFSR[:,jlevel,region,:]
        GEFSv12_data = AC_GEFSv12[:,jlevel,region,:]
        
        axloc = [xstart[ix],ystart[jlevel],xlen[ix],ylen[jlevel]]
        ax1 = fig.add_axes(axloc)
        title = cletters[ktr]+clevels[jlevel]+cregions[ix]
        ax1.set_title(title, fontsize=11,color='Black')
        for istream in range(5):
            if istream == 0:
                ax1.plot([0.25,1.0,2.0,3.0,4.0,5.0], CFSR_data[:,istream],\
                    color='Blue',linestyle='--',lw=0.6,label='CFSR')
                ax1.plot([0.25,1.0,2.0,3.0,4.0,5.0], GEFSv12_data[:,istream],\
                    color='Red',lw=0.6,label='GEFSv12')
            else:
                ax1.plot([0.25,1.0,2.0,3.0,4.0,5.0], CFSR_data[:,istream],\
                    color='Blue',linestyle='--',lw=0.6)
                ax1.plot([0.25,1.0,2.0,3.0,4.0,5.0], GEFSv12_data[:,istream],\
                    color='Red',lw=0.6)
        ax1.set_xlim(-0.1,5.1)
        ax1.set_ylim(.48,1.02)
        ax1.set_xticks([0,1,2,3,4,5])
        ax1.set_yticks([0.5,0.6,0.7,0.8,0.9,1.0])
        if ix == 0: ax1.set_ylabel('Anomaly correlation',fontsize=9)
        if jlevel == 2: ax1.set_xlabel('Lead time (days)',fontsize=9)
        if ktr == 1: ax1.legend(loc=0)
        ax1.grid(True,lw=0.25,color='LightGray')
        ktr = ktr+1

# ---- set plot title

plot_title = 'AC_bystream_temperature.png'
fig.savefig(plot_title, dpi=300)
print ('saving plot to file = ',plot_title)
print ('Done!')

    

