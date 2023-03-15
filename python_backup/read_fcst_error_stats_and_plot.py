import csv
import sys
import _pickle as cPickle
import numpy as np
from datetime import datetime
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams 
rcParams['xtick.labelsize']='x-small'
rcParams['ytick.labelsize']='x-small'

# ---- declare arrays

vsum_T_gefs = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
vsum_Z_gefs = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
vsum_U_gefs = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
vsum_V_gefs = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
vsum_T_cfsr = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
vsum_Z_cfsr = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
vsum_U_cfsr = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
vsum_V_cfsr = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh

bsum_T_gefs = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
bsum_Z_gefs = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
bsum_U_gefs = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
bsum_V_gefs = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
bsum_T_cfsr = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
bsum_Z_cfsr = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
bsum_U_cfsr = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh
bsum_V_cfsr = np.zeros((14,6,3), dtype=np.float64) # nlevels, nleads, nh/tr/sh

ktr_T_gefs = np.zeros((14,6,3), dtype=np.int32) # nlevels, nleads, nh/tr/sh
ktr_Z_gefs = np.zeros((14,6,3), dtype=np.int32) # nlevels, nleads, nh/tr/sh
ktr_U_gefs = np.zeros((14,6,3), dtype=np.int32) # nlevels, nleads, nh/tr/sh
ktr_V_gefs = np.zeros((14,6,3), dtype=np.int32) # nlevels, nleads, nh/tr/sh
ktr_T_cfsr = np.zeros((14,6,3), dtype=np.int32) # nlevels, nleads, nh/tr/sh
ktr_Z_cfsr = np.zeros((14,6,3), dtype=np.int32) # nlevels, nleads, nh/tr/sh
ktr_U_cfsr = np.zeros((14,6,3), dtype=np.int32) # nlevels, nleads, nh/tr/sh
ktr_V_cfsr = np.zeros((14,6,3), dtype=np.int32) # nlevels, nleads, nh/tr/sh

# ---- loop thru levels, streams, mode, variable, load data and sum.

clevels = ['925','850','800','750','700', '600','500','400','300','250', '200','150','100','10']
for ilevel, input_level in enumerate(clevels):
    #for input_stream in ['1999','2003','2007','2011']:
    for input_stream in ['2007','2011','2015']:
    #for input_stream in ['1999','2003','2007','2011','2015']:
        for input_model in ['cfsr','reanalysis']:
            #for input_variable in ['T','U','V','Z']:
            for input_variable in ['T']:
                for ihr, input_hour in enumerate(['006','024','048','072','096','120']):

                    input_file = 'reanalysis/'+input_model+'_stream='+\
                        input_stream+'_lead='+input_hour+'_var='+\
                        input_variable+'_level='+input_level+'.cPick'
                    inf = open(input_file, 'rb')
                    print ('input_file = ',input_file)
                    ICdate_gl_arr = cPickle.load(inf)
                    bias_gl_arr = cPickle.load(inf)
                    rmse_gl_arr = cPickle.load(inf)
                    #print ('icdate = ', ICdate_gl_arr)
                    #print ('bias_gl_arr = ', bias_gl_arr)
                    #print ('rmse_gl_arr = ', rmse_gl_arr)
                    #sys.exit()
                    ICdate_nh_arr = cPickle.load(inf)
                    bias_nh_arr = cPickle.load(inf)
                    rmse_nh_arr = cPickle.load(inf) 
                    ICdate_sh_arr = cPickle.load(inf)
                    bias_sh_arr = cPickle.load(inf) 
                    rmse_sh_arr = cPickle.load(inf)
                    ICdate_tr_arr = cPickle.load(inf)
                    bias_tr_arr = cPickle.load(inf) 
                    rmse_tr_arr = cPickle.load(inf)
                    inf.close()
                        
                    bsum_nh = np.sum(bias_nh_arr)
                    vsum_nh = np.sum(rmse_nh_arr**2)
                    ktr_nh = len(bias_nh_arr)
                    bsum_tr = np.sum(bias_tr_arr)
                    vsum_tr = np.sum(rmse_tr_arr**2)
                    ktr_tr = len(bias_tr_arr)
                    bsum_sh = np.sum(bias_sh_arr)
                    vsum_sh = np.sum(rmse_sh_arr**2)
                    ktr_sh = len(bias_sh_arr)
                    #print ('bsum_nh, vsum_nh, ktr_nh = ', bsum_nh, vsum_nh, ktr_nh )
                    #print ('bsum_tr, vsum_tr, ktr_tr = ', bsum_tr, vsum_tr, ktr_tr )
                    #print ('bsum_sh, vsum_sh, ktr_sh = ', bsum_sh, vsum_sh, ktr_sh )
                        
                    if input_variable == 'T':
                        if input_model == 'reanalysis':
                            bsum_T_gefs[ilevel,ihr,0] = bsum_T_gefs[ilevel,ihr,0] + bsum_nh
                            vsum_T_gefs[ilevel,ihr,0] = vsum_T_gefs[ilevel,ihr,0] + vsum_nh
                            ktr_T_gefs[ilevel,ihr,0] = ktr_T_gefs[ilevel,ihr,0] + ktr_nh
                            bsum_T_gefs[ilevel,ihr,1] = bsum_T_gefs[ilevel,ihr,1] + bsum_tr
                            vsum_T_gefs[ilevel,ihr,1] = vsum_T_gefs[ilevel,ihr,1] + vsum_tr
                            ktr_T_gefs[ilevel,ihr,1] = ktr_T_gefs[ilevel,ihr,1] + ktr_tr
                            bsum_T_gefs[ilevel,ihr,2] = bsum_T_gefs[ilevel,ihr,2] + bsum_sh
                            vsum_T_gefs[ilevel,ihr,2] = vsum_T_gefs[ilevel,ihr,2] + vsum_sh
                            ktr_T_gefs[ilevel,ihr,2] = ktr_T_gefs[ilevel,ihr,2] + ktr_sh
                        else:
                            bsum_T_cfsr[ilevel,ihr,0] = bsum_T_cfsr[ilevel,ihr,0] + bsum_nh
                            vsum_T_cfsr[ilevel,ihr,0] = vsum_T_cfsr[ilevel,ihr,0] + vsum_nh
                            ktr_T_cfsr[ilevel,ihr,0] = ktr_T_cfsr[ilevel,ihr,0] + ktr_nh
                            bsum_T_cfsr[ilevel,ihr,1] = bsum_T_cfsr[ilevel,ihr,1] + bsum_tr
                            vsum_T_cfsr[ilevel,ihr,1] = vsum_T_cfsr[ilevel,ihr,1] + vsum_tr
                            ktr_T_cfsr[ilevel,ihr,1] = ktr_T_cfsr[ilevel,ihr,1] + ktr_tr
                            bsum_T_cfsr[ilevel,ihr,2] = bsum_T_cfsr[ilevel,ihr,2] + bsum_sh
                            vsum_T_cfsr[ilevel,ihr,2] = vsum_T_cfsr[ilevel,ihr,2] + vsum_sh
                            ktr_T_cfsr[ilevel,ihr,2] = ktr_T_cfsr[ilevel,ihr,2] + ktr_sh
                    elif input_variable == 'Z':
                        if input_model == 'reanalysis':
                            bsum_Z_gefs[ilevel,ihr,0] = bsum_Z_gefs[ilevel,ihr,0] + bsum_nh
                            vsum_Z_gefs[ilevel,ihr,0] = vsum_Z_gefs[ilevel,ihr,0] + vsum_nh
                            ktr_Z_gefs[ilevel,ihr,0] = ktr_Z_gefs[ilevel,ihr,0] + ktr_nh
                            bsum_Z_gefs[ilevel,ihr,1] = bsum_Z_gefs[ilevel,ihr,1] + bsum_tr
                            vsum_Z_gefs[ilevel,ihr,1] = vsum_Z_gefs[ilevel,ihr,1] + vsum_tr
                            ktr_Z_gefs[ilevel,ihr,1] = ktr_Z_gefs[ilevel,ihr,1] + ktr_tr
                            bsum_Z_gefs[ilevel,ihr,2] = bsum_Z_gefs[ilevel,ihr,2] + bsum_sh
                            vsum_Z_gefs[ilevel,ihr,2] = vsum_Z_gefs[ilevel,ihr,2] + vsum_sh
                            ktr_Z_gefs[ilevel,ihr,2] = ktr_Z_gefs[ilevel,ihr,2] + ktr_sh
                        else:
                            bsum_Z_cfsr[ilevel,ihr,0] = bsum_Z_cfsr[ilevel,ihr,0] + bsum_nh
                            vsum_Z_cfsr[ilevel,ihr,0] = vsum_Z_cfsr[ilevel,ihr,0] + vsum_nh
                            ktr_Z_cfsr[ilevel,ihr,0] = ktr_Z_cfsr[ilevel,ihr,0] + ktr_nh
                            bsum_Z_cfsr[ilevel,ihr,1] = bsum_Z_cfsr[ilevel,ihr,1] + bsum_tr
                            vsum_Z_cfsr[ilevel,ihr,1] = vsum_Z_cfsr[ilevel,ihr,1] + vsum_tr
                            ktr_Z_cfsr[ilevel,ihr,1] = ktr_Z_cfsr[ilevel,ihr,1] + ktr_tr
                            bsum_Z_cfsr[ilevel,ihr,2] = bsum_Z_cfsr[ilevel,ihr,2] + bsum_sh
                            vsum_Z_cfsr[ilevel,ihr,2] = vsum_Z_cfsr[ilevel,ihr,2] + vsum_sh
                            ktr_Z_cfsr[ilevel,ihr,2] = ktr_Z_cfsr[ilevel,ihr,2] + ktr_sh
                    elif input_variable == 'U':
                        if input_model == 'reanalysis':
                            bsum_U_gefs[ilevel,ihr,0] = bsum_U_gefs[ilevel,ihr,0] + bsum_nh
                            vsum_U_gefs[ilevel,ihr,0] = vsum_U_gefs[ilevel,ihr,0] + vsum_nh
                            ktr_U_gefs[ilevel,ihr,0] = ktr_U_gefs[ilevel,ihr,0] + ktr_nh
                            bsum_U_gefs[ilevel,ihr,1] = bsum_U_gefs[ilevel,ihr,1] + bsum_tr
                            vsum_U_gefs[ilevel,ihr,1] = vsum_U_gefs[ilevel,ihr,1] + vsum_tr
                            ktr_U_gefs[ilevel,ihr,1] = ktr_U_gefs[ilevel,ihr,1] + ktr_tr
                            bsum_U_gefs[ilevel,ihr,2] = bsum_U_gefs[ilevel,ihr,2] + bsum_sh
                            vsum_U_gefs[ilevel,ihr,2] = vsum_U_gefs[ilevel,ihr,2] + vsum_sh
                            ktr_U_gefs[ilevel,ihr,2] = ktr_U_gefs[ilevel,ihr,2] + ktr_sh
                        else:
                            bsum_U_cfsr[ilevel,ihr,0] = bsum_U_cfsr[ilevel,ihr,0] + bsum_nh
                            vsum_U_cfsr[ilevel,ihr,0] = vsum_U_cfsr[ilevel,ihr,0] + vsum_nh
                            ktr_U_cfsr[ilevel,ihr,0] = ktr_U_cfsr[ilevel,ihr,0] + ktr_nh
                            bsum_U_cfsr[ilevel,ihr,1] = bsum_U_cfsr[ilevel,ihr,1] + bsum_tr
                            vsum_U_cfsr[ilevel,ihr,1] = vsum_U_cfsr[ilevel,ihr,1] + vsum_tr
                            ktr_U_cfsr[ilevel,ihr,1] = ktr_U_cfsr[ilevel,ihr,1] + ktr_tr
                            bsum_U_cfsr[ilevel,ihr,2] = bsum_U_cfsr[ilevel,ihr,2] + bsum_sh
                            vsum_U_cfsr[ilevel,ihr,2] = vsum_U_cfsr[ilevel,ihr,2] + vsum_sh
                            ktr_U_cfsr[ilevel,ihr,2] = ktr_U_cfsr[ilevel,ihr,2] + ktr_sh      
                    elif input_variable == 'V':
                        if input_model == 'reanalysis':
                            bsum_V_gefs[ilevel,ihr,0] = bsum_V_gefs[ilevel,ihr,0] + bsum_nh
                            vsum_V_gefs[ilevel,ihr,0] = vsum_V_gefs[ilevel,ihr,0] + vsum_nh
                            ktr_V_gefs[ilevel,ihr,0] = ktr_V_gefs[ilevel,ihr,0] + ktr_nh
                            bsum_V_gefs[ilevel,ihr,1] = bsum_V_gefs[ilevel,ihr,1] + bsum_tr
                            vsum_V_gefs[ilevel,ihr,1] = vsum_V_gefs[ilevel,ihr,1] + vsum_tr
                            ktr_V_gefs[ilevel,ihr,1] = ktr_V_gefs[ilevel,ihr,1] + ktr_tr
                            bsum_V_gefs[ilevel,ihr,2] = bsum_V_gefs[ilevel,ihr,2] + bsum_sh
                            vsum_V_gefs[ilevel,ihr,2] = vsum_V_gefs[ilevel,ihr,2] + vsum_sh
                            ktr_V_gefs[ilevel,ihr,2] = ktr_V_gefs[ilevel,ihr,2] + ktr_sh
                        else:
                            bsum_V_cfsr[ilevel,ihr,0] = bsum_V_cfsr[ilevel,ihr,0] + bsum_nh
                            vsum_V_cfsr[ilevel,ihr,0] = vsum_V_cfsr[ilevel,ihr,0] + vsum_nh
                            ktr_V_cfsr[ilevel,ihr,0] = ktr_V_cfsr[ilevel,ihr,0] + ktr_nh
                            bsum_V_cfsr[ilevel,ihr,1] = bsum_V_cfsr[ilevel,ihr,1] + bsum_tr
                            vsum_V_cfsr[ilevel,ihr,1] = vsum_V_cfsr[ilevel,ihr,1] + vsum_tr
                            ktr_V_cfsr[ilevel,ihr,1] = ktr_V_cfsr[ilevel,ihr,1] + ktr_tr
                            bsum_V_cfsr[ilevel,ihr,2] = bsum_V_cfsr[ilevel,ihr,2] + bsum_sh
                            vsum_V_cfsr[ilevel,ihr,2] = vsum_V_cfsr[ilevel,ihr,2] + vsum_sh
                            ktr_V_cfsr[ilevel,ihr,2] = ktr_V_cfsr[ilevel,ihr,2] + ktr_sh                        
                        
                       
rmse_T_gefs = np.sqrt(vsum_T_gefs / ktr_T_gefs)
rmse_Z_gefs = np.sqrt(vsum_Z_gefs / ktr_Z_gefs)
rmse_U_gefs = np.sqrt(vsum_U_gefs / ktr_U_gefs)
rmse_V_gefs = np.sqrt(vsum_V_gefs / ktr_V_gefs)
rmse_T_cfsr = np.sqrt(vsum_T_cfsr / ktr_T_cfsr) 
rmse_Z_cfsr = np.sqrt(vsum_Z_cfsr / ktr_T_cfsr)
rmse_U_cfsr = np.sqrt(vsum_U_cfsr / ktr_T_cfsr)
rmse_V_cfsr = np.sqrt(vsum_V_cfsr / ktr_T_cfsr)

bia_T_gefs = bsum_T_gefs / ktr_T_gefs
bia_Z_gefs = bsum_Z_gefs / ktr_Z_gefs
bia_U_gefs = bsum_U_gefs / ktr_U_gefs
bia_V_gefs = bsum_V_gefs / ktr_V_gefs
bia_T_cfsr = bsum_T_cfsr / ktr_T_cfsr
bia_Z_cfsr = bsum_Z_cfsr / ktr_T_cfsr
bia_U_cfsr = bsum_U_cfsr / ktr_T_cfsr
bia_V_cfsr = bsum_V_cfsr / ktr_T_cfsr

# ----- make plots

fig = plt.figure(figsize=(6.5,7.))

fig.suptitle('2007-2019 Temperature forecast RMSE difference and biases',fontsize=13)


hours = [6,24,48,72,96,120]
plevels = [925,850,800,750,700,600,500,400,300,250,200,150,100,10]
#clevels = [-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.01,1.5,2.0,2.5]
clevels = [-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]
colorst = ['#0000ff', '#6666ff', '#b2b2ff', '#e6e6ff', \
    'White', '#ffe6e6', '#ffb2b2', '#ff7373', '#ff0000'] 
    
a1 = fig.add_axes([.10,.72,.24,.20])
a1.set_title('(a) NH RMSE diff (GEFSv12-CFSR)',fontsize=7)
cf = a1.contourf(hours,plevels,rmse_T_gefs[:,:,0]-rmse_T_cfsr[:,:,0],\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,rmse_T_gefs[:,:,0]-rmse_T_cfsr[:,:,0],\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=7)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

a1 = fig.add_axes([.42,.72,.24,.20])
a1.set_title('(b) TR RMSE diff(GEFSv12-CFSR)',fontsize=7)
a1.contourf(hours,plevels,rmse_T_gefs[:,:,1]-rmse_T_cfsr[:,:,1],\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,rmse_T_gefs[:,:,1]-rmse_T_cfsr[:,:,1],\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

a1 = fig.add_axes([.74,.72,.24,.20])
a1.set_title('(c) SH RMSE diff (GEFSv12-CFSR)',fontsize=7)
a1.contourf(hours,plevels,rmse_T_gefs[:,:,2]-rmse_T_cfsr[:,:,2],\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,rmse_T_gefs[:,:,2]-rmse_T_cfsr[:,:,2],\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

cax = fig.add_axes([0.1,0.67,0.8,0.013])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('RMSE, GEFSv12-CFSR (deg C)',fontsize=6)
cb.ax.tick_params(labelsize=6)

# ----- GEFS bias

clevels = [-2,-1.5,-1,-0.5,-0.2,0.2,0.5,1,1.5,2]

a1 = fig.add_axes([.10,.38,.24,.2])
a1.set_title('(d) NH bias, GEFSv12 init',fontsize=7.5)
cf = a1.contourf(hours,plevels,bia_T_gefs[:,:,0],\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,bia_T_gefs[:,:,0],\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=7)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

a1 = fig.add_axes([.42,.38,.24,.20])
a1.set_title('(e) TR bias, GEFSv12 initial',fontsize=7.5)
a1.contourf(hours,plevels,bia_T_gefs[:,:,1],\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,bia_T_gefs[:,:,1],\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])

a1 = fig.add_axes([.74,.38,.24,.20])
a1.set_title('(f) SH bias, GEFSv12 initial',fontsize=7.5)
a1.contourf(hours,plevels,bia_T_gefs[:,:,2],\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,bia_T_gefs[:,:,2],\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])



# ----- Bias CFSR

a1 = fig.add_axes([.10,.12,.24,.20])
a1.set_title('(g) NH bias, CFSR initial',fontsize=7.5)
cf = a1.contourf(hours,plevels,bia_T_cfsr[:,:,0],\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,bia_T_cfsr[:,:,0],\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
a1.set_ylabel('Pressure (hPa)', fontsize=7)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])
a1.set_xlabel('Forecast lead time (h)', fontsize=7)

a1 = fig.add_axes([.42,.12,.24,.2])
a1.set_title('(h) TR bias, CFSR initial',fontsize=7.5)
a1.contourf(hours,plevels,bia_T_cfsr[:,:,1],\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,bia_T_cfsr[:,:,1],\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])
a1.set_xlabel('Forecast lead time (h)', fontsize=7)

a1 = fig.add_axes([.74,.12,.24,.2])
a1.set_title('(i) SH bias, CFSR initial',fontsize=7.5)
a1.contourf(hours,plevels,bia_T_cfsr[:,:,2],\
    levels = clevels, colors=colorst,extend='both')
a1.contour(hours,plevels,bia_T_cfsr[:,:,2],\
    levels = clevels, colors='Black',linewidths=0.3)
a1.set_ylim(1000,10)
#a1.set_ylabel('Pressure (hPa)', fontsize=8)
a1.set_xlim(0,120)
a1.set_xticks([0,24,48,72,96,120])
a1.set_xlabel('Forecast lead time (h)', fontsize=7)

cax = fig.add_axes([0.1,0.05,0.8,0.013])
cb = fig.colorbar(cf, ticks=clevels,cax=cax, \
    orientation='horizontal',drawedges=True,format='%g') # im1 "bottom", size="3%", pad="1%",
cb.set_label('Bias (deg C)',fontsize=6)
cb.ax.tick_params(labelsize=6)


plot_title = 'temperature_rmse_difference_2007-2019.png'
print ('saving plot to file = ',plot_title)
plt.savefig(plot_title,dpi=600)
print ('Plot done')



        
        
