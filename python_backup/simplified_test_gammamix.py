"""
simplified_test_gammamix.py

testing calling the R function mixtools.gammamixEM

"""

import os, sys
import numpy as np
import _pickle as cPickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import _pickle as cPickle
base = importr('base')
mixtools = importr('mixtools')

infile = 'precipitation_data.cPick'
print ('reading from ', infile)
inf = open(infile,'rb')
fraction_zero = cPickle.load(inf)
precip_nonzero = cPickle.load(inf)
print ('fraction zero = ', fraction_zero)
print ('precip_nonzero[0:-1:10] = ', precip_nonzero[0:-1:10])
inf.close()

# --- convert to R vector, call R Gamma mixture function

precip_nonzero_R = robjects.FloatVector(precip_nonzero)

#rmixture = robjects.r['mixtools.gammamixEM']
#result_R = rmixture(precip_nonzero_R, alpha=NULL, beta=NULL,k=2)

result_R = mixtools.gammamixEM(precip_nonzero_R, k=2)

result_np = np.asarray(result_R)
print (result_np)



