# --- plot gamma distributions that illustrate sharpening of pdfs

import scipy.special as special
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt

def gamma_pdf(alpha, beta, x):
    
    gamma = special.gamma(alpha) 
    denom = beta * gamma
    pdf = (x/beta)**(alpha-1.0)*np.exp(-x/beta) / denom
    return pdf
    

mean_climatology = 6.1
var_climatology = 8.0
alpha_climatology = mean_climatology**2/var_climatology
beta_climatology = var_climatology/mean_climatology
print (alpha_climatology, beta_climatology)

mean_trend = 6.4
var_trend = 7.0
alpha_trend = mean_trend**2/var_trend
beta_trend = var_trend/mean_trend
print (alpha_trend, beta_trend)

mean_enso = 6.7
var_enso = 6.0
alpha_enso = mean_enso**2/var_enso
beta_enso = var_enso/mean_enso
print (alpha_enso, beta_enso)

mean_enso_subs = 6.9
var_enso_subs = 5.0
alpha_enso_subs = mean_enso_subs**2/var_enso_subs
beta_enso_subs = var_enso_subs/mean_enso_subs
print (alpha_enso_subs, beta_enso_subs)

mean_week2 = 7.5
var_week2 = 3
alpha_week2 = mean_week2**2/var_week2
beta_week2 = var_week2/mean_week2
print (alpha_week2, beta_week2)

mean_3_5 = 9.
var_3_5 = 1.4
alpha_3_5 = mean_3_5**2/var_3_5
beta_3_5 = var_3_5/mean_3_5
print (alpha_3_5, beta_3_5)

mean_1_2 = 11.0
var_1_2 = 0.05
alpha_1_2 = mean_1_2**2/var_1_2
beta_1_2 = var_1_2/mean_1_2
print (alpha_1_2, beta_1_2)

x = np.arange(0.0,25.01,0.1)
pdf_climatology = gamma_pdf(alpha_climatology, beta_climatology, x)
pdf_trend = gamma_pdf(alpha_trend, beta_trend, x)
pdf_enso = gamma_pdf(alpha_enso, beta_enso, x)
pdf_enso_subs = gamma_pdf(alpha_enso_subs, beta_enso_subs, x)
pdf_week2 = gamma_pdf(alpha_week2, beta_week2, x)
pdf_3_5 = gamma_pdf(alpha_3_5, beta_3_5, x)
#pdf_1_2 = gamma_pdf(alpha_1_2, beta_1_2, x)
pdf_1_2 = norm.pdf(x, loc=mean_1_2, scale=np.sqrt(var_1_2))
print (pdf_climatology[0:-1:10])
print (pdf_enso[0:-1:10])
print (pdf_enso_subs[0:-1:10])
print (pdf_3_5[0:-1:10])
print (pdf_1_2[0:-1:10])


fig = plt.figure(figsize=(6.5,5))
    
axloc = [0.13,0.13,0.8,0.8]
ax = fig.add_axes(axloc)
title = 'Refinement of predictive pdfs'
ax.set_title(title, fontsize=13,color='Black')
ax.plot(x,pdf_climatology,'k-',lw=2,label='Unconditional climatology')
ax.plot(x,pdf_trend,'-',color='Blue',lw=2,label='Climatology with trend')
ax.plot(x,pdf_enso,'-',color='Green',lw=2,label='Climatology with trend and SST anomaly')
ax.plot(x,pdf_enso_subs,'-',color='DarkRed',lw=2,label='Weeks 3-4')
ax.plot(x,pdf_week2,'-',color='DarkOrange',lw=2,label='Week 2')
ax.plot(x,pdf_3_5,'-',color='LightSkyBlue',lw=2,label='Day +5 forecast')
ax.plot(x,pdf_1_2,'-',color='DarkOrchid',lw=2,label='Day +1 forecast')
ax.legend(loc=0)
ax.set_xlim(0,13)
ax.set_ylim(0,1.0)
ax.set_xlabel('Temperature (deg C)')
ax.set_ylabel('Probability density')

# ---- set plot title

plot_title = 'pdf_sharpening.png'
fig.savefig(plot_title, dpi=400)
print ('saving plot to file = ',plot_title)
print ('Done!')


