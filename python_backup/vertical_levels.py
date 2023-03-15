import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams 

filename = 'vertical_levels.txt'
x = np.loadtxt(filename, usecols=(0, 1))
preslevels = x[:,0] + x[:,1]*100000.
print (x)
print (len(preslevels))
plt.gca().invert_yaxis()
fig = plt.figure(figsize=(4.,6.))

fig.suptitle('Vertical levels',fontsize=15)
a1 = fig.add_axes([.17,.05,.28,.84])
a1.set_title('(a) logarithmic')
for i in range(len(preslevels)):
    print (x[i,0], x[i,1], preslevels[i]/100.)
    a1.semilogy([0,1],[preslevels[i]/100., preslevels[i]/100.],\
        '-',color='Gray',linewidth=0.7)
a1.set_ylabel('Pressure (hPa)')
a1.set_xlim(0,1)
a1.set_xticks([0,1])
a1.set_xticklabels([' ', ' '])
a1.set_ylim(1000,0.1)

a1 = fig.add_axes([.67,.05,.28,.84])
a1.set_title('(b) linear')
for i in range(len(preslevels)):
    print (x[i,0], x[i,1], preslevels[i]/100.)
    a1.plot([0,1],[preslevels[i]/100., preslevels[i]/100.],\
        '-',color='Gray',linewidth=0.7)
a1.set_ylabel('Pressure (hPa)')
a1.set_xlim(0,1)
a1.set_xticks([0,1])
a1.set_xticklabels([' ', ' '])
a1.set_ylim(1000,0.1)

plot_title = 'vertlevels.png'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title,dpi=600)
print ('Plot done')
