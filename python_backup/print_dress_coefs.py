# python print_dress_coefs.py


from netCDF4 import Dataset
import numpy as np
import sys, os
import matplotlib.pyplot as plt

cleads = ['012','024','036','048','060','072','084',\
    '096','108','120','132','144','156','168','180',\
    '192','204','216','228','240']
nleads = len(cleads)
nmonths = 12
cmonths = ['01','02','03','04','05','06',\
    '07','08','09','10','11','12']
cmonthnames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',\
    'Oct','Nov','Dec']

b0_mean_lowrank_all = np.zeros((nleads,nmonths), dtype=np.float64)
b1_mean_lowrank_all = np.zeros((nleads,nmonths), dtype=np.float64) 
b0_std_lowrank_all = np.zeros((nleads,nmonths), dtype=np.float64) 
b1_std_lowrank_all = np.zeros((nleads,nmonths), dtype=np.float64) 

b0_mean_midrank_all = np.zeros((nleads,nmonths), dtype=np.float64) 
b1_mean_midrank_all = np.zeros((nleads,nmonths), dtype=np.float64) 
b0_std_midrank_all = np.zeros((nleads,nmonths), dtype=np.float64)  
b1_std_midrank_all = np.zeros((nleads,nmonths), dtype=np.float64) 

b0_mean_highrank_all = np.zeros((nleads,nmonths), dtype=np.float64)  
b1_mean_highrank_all = np.zeros((nleads,nmonths), dtype=np.float64)  
b0_std_highrank_all = np.zeros((nleads,nmonths), dtype=np.float64)  
b1_std_highrank_all = np.zeros((nleads,nmonths), dtype=np.float64)  

# ---- read the stored netCDF thinned probabilities

master_directory_histogram = '/Volumes/NBM/chist/'

for imonth, cmonth in enumerate(cmonths):
    for ilead,clead in enumerate(cleads):
        infile = master_directory_histogram +\
            'closest_member_histogram_thinned'+\
            '_month='+cmonth+'_lead='+clead+'h.nc'
        print ('reading from ', infile)

        nc = Dataset(infile,'r')

        b0_mean_lowrank_all[ilead,imonth] = nc.variables['b0_mean_lowrank'][0]
        b1_mean_lowrank_all[ilead,imonth] =nc.variables['b1_mean_lowrank'][0]
        b0_std_lowrank_all[ilead,imonth] = nc.variables['b0_std_lowrank'][0]
        b1_std_lowrank_all[ilead,imonth] =nc.variables['b1_std_lowrank'][0]

        b0_mean_midrank_all[ilead,imonth] = nc.variables['b0_mean_midrank'][0]
        b1_mean_midrank_all[ilead,imonth] =nc.variables['b1_mean_midrank'][0]
        b0_std_midrank_all[ilead,imonth] = nc.variables['b0_std_midrank'][0]
        b1_std_midrank_all[ilead,imonth] = nc.variables['b1_std_midrank'][0]

        b0_mean_highrank_all[ilead,imonth] = nc.variables['b0_mean_highrank'][0]
        b1_mean_highrank_all[ilead,imonth] = nc.variables['b1_mean_highrank'][0]
        b0_std_highrank_all[ilead,imonth] = nc.variables['b0_std_highrank'][0]
        b1_std_highrank_all[ilead,imonth] = nc.variables['b1_std_highrank'][0]
        nc.close()
  
  
#print ('b0_std_highrank_all = ', b0_std_highrank_all)  
#print ('b1_std_highrank_all = ', b1_std_highrank_all)  
#sys.exit()
    
# ---- make plots

for imonth, cmonthname in enumerate(cmonthnames):
    for icat in range(3):
        if icat == 0:
            b0_mean = b0_mean_lowrank_all[:,imonth]
            b1_mean = b1_mean_lowrank_all[:,imonth]
            b0_std = b0_std_lowrank_all[:,imonth]
            b1_std = b1_std_lowrank_all[:,imonth]
            ct = 'Lowest member'
            ct2 = 'lowrank'
        elif icat == 1:
            b0_mean = b0_mean_midrank_all[:,imonth]
            b1_mean = b1_mean_midrank_all[:,imonth]
            b0_std = b0_std_midrank_all[:,imonth]
            b1_std = b1_std_midrank_all[:,imonth]
            ct = 'Middle members'
            ct2 = 'midrank'
        elif icat == 2:
            b0_mean = b0_mean_highrank_all[:,imonth]
            b1_mean = b1_mean_highrank_all[:,imonth]
            b0_std = b0_std_highrank_all[:,imonth]
            b1_std = b1_std_highrank_all[:,imonth]
            ct = 'Highest member'
            ct2 = 'highrank'

        
        fig = plt.figure(figsize=(9., 6.5))
        
        axlocn = [0.08,0.57,0.4,0.36]
        a1 = fig.add_axes(axlocn)
        a1.set_title('b0 mean for month = '+cmonthname+', '+ct,fontsize=11)
        a1.set_xlabel('Lead time (days)',fontsize=9)
        a1.set_ylabel('Coefficient value',fontsize=9)
        a1.set_xlim(0,241)
        #a1.set_ylim(0.001, xtop)
        a1.grid(color='Gray',lw=0.2,linestyle='--',axis='x')
        a1.plot(range(12,241,12),b0_mean,'-',color='Red',linewidth=1.5)
        a1.set_xticks([24,48,72,96,120,144,168,192,216,240])

        axlocn = [0.58,0.57,0.4,0.36]
        a1 = fig.add_axes(axlocn)
        a1.set_title('b1 mean for month = '+cmonthname+', '+ct,fontsize=11)
        a1.set_xlabel('Lead time (days)',fontsize=9)
        a1.set_ylabel('Coefficient value',fontsize=9)
        a1.set_xlim(0,241)
        #a1.set_ylim(0.001, xtop)
        a1.grid(color='Gray',lw=0.2,linestyle='--',axis='x')
        a1.plot(range(12,241,12),b1_mean,'-',color='Red',linewidth=1.5)
        a1.set_xticks([24,48,72,96,120,144,168,192,216,240])

        axlocn = [0.08,0.07,0.4,0.36]
        a1 = fig.add_axes(axlocn)
        a1.set_title('b0 std for month = '+cmonthname+', '+ct,fontsize=11)
        a1.set_xlabel('Lead time (days)',fontsize=9)
        a1.set_ylabel('Coefficient value',fontsize=9)
        a1.set_xlim(0,241)
        #a1.set_ylim(0.001, xtop)
        a1.grid(color='Gray',lw=0.2,linestyle='--',axis='x')
        a1.plot(range(12,241,12),b0_std,'-',color='Red',linewidth=1.5)
        a1.set_xticks([24,48,72,96,120,144,168,192,216,240])

        axlocn = [0.58,0.07,0.4,0.36]
        a1 = fig.add_axes(axlocn)
        a1.set_title('b1 std for month = '+cmonthname+', '+ct,fontsize=11)
        a1.set_xlabel('Lead time (days)',fontsize=9)
        a1.set_ylabel('Coefficient value',fontsize=9)
        a1.set_xlim(0,10)
        #a1.set_ylim(0.001, xtop)
        a1.grid(color='Gray',lw=0.2,linestyle='--',axis='x')
        a1.plot(range(12,241,12),b1_std,'-',color='Red',linewidth=1.5)
        a1.set_xticks([24,48,72,96,120,144,168,192,216,240])

        plot_title = 'dress_coeffs_'+cmonthname+'_'+ct2+'.png'
        fig.savefig(plot_title, dpi=300)
        plt.close()
        print ('saving plot to file = ',plot_title)




