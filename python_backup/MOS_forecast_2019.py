
def MOS_forecast_2019(npts, nlats, nlons, \
    forecast_3d, analyses_3d, beta_3d, lsmask, \
    date_list_forecast, clead, cpath_forecast, \
    cpath_era5):
    
    """ apply MOS forecast regression procedure in 2019. 
    """
    import numpy as np
    from datetime import datetime
    import sys
    import _pickle as cPickle
    from dateutils import daterange, dateshift, dayofyear, splitdate
    
    cmonths = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
         
    # ---- read the climatology file.
    
    if clead == '12' or clead == '36' or clead == '60' or clead == '84' or clead == '108':
        infile = cpath_era5 + 'ERA5_temperature_climatology_12UTC.cPick'
    else:
        infile = cpath_era5 + 'ERA5_temperature_climatology_00UTC.cPick'
    inf = open(infile,'rb')
    climo_temps_estimated = cPickle.load(inf)
    inf.close()
         
    # ---- sequentially loop through dates during the sample, updating
    #      the previous day's bias correction to the new days fcst vs. obs
    #      discrepancy.
    
    fcollect = []
    acollect = []
    fcollect_feb = []
    acollect_feb = []
    first = True
    ndates = int(len(date_list_forecast))
    obsinc_2d = np.zeros((nlats, nlons), dtype=np.float64)
    #for idate, date in enumerate(date_list_forecast[1:]):
    for idate, date in enumerate(date_list_forecast):

        cdd = date[6:8]
        cmm = date[4:6]
        cyearf = date[0:4]
        if first == True or cdd == '01':
            
            cmonth = cmonths[int(cmm)-1]
            first = False
            
            # ---- read the MOS regression coefficient data
            
            infile = cpath_forecast+'MOS_slope_intercept_'+cmonth+\
                '_lead='+clead+'.cPick'
            inf = open(infile,'rb')
            slope = cPickle.load(inf)
            intercept = cPickle.load(inf)
            inf.close()
            
        # ---- determine the julian day

        imm = int(cmm)
        idd = int(cdd)
        iyear_full = int(cyearf)
        julday = dayofyear(iyear_full, imm, idd) - 1
        if julday > 364: julday = 364 
         
        #ftoday = forecast_3d[idate,:,:] - climo_temps_estimated[julday,:,:]
        #atoday = analyses_3d[idate,:,:] - climo_temps_estimated[julday,:,:]
        ftoday = forecast_3d[idate,:,:] 
        atoday = analyses_3d[idate,:,:] 
        regressed = slope[:,:]*ftoday[:,:] + intercept[:,:]
        beta_3d[idate,:,:] = lsmask[:,:]*(ftoday[:,:] - regressed[:,:])
        
        if cmonth == 'Jul' and clead == '24':
            fcollect.append(ftoday[20,40])
            acollect.append(atoday[20,40])
            #print ('f,o collection for idate, date = ',ftoday[20,40], atoday[20,40], idate,date)
        if cmonth == 'Feb' and clead == '24':
            fcollect_feb.append(ftoday[20,40])
            acollect_feb.append(atoday[20,40])
            print ('f,o collection for idate, date = ',ftoday[20,40], atoday[20,40], idate,date)
        
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    
    if clead == '24':
        outfile = 'boulder_july_data_2019_lead='+clead+'h.cPick'
        ouf = open(outfile,'wb')
        cPickle.dump(np.array(acollect), ouf)
        cPickle.dump(np.array(fcollect), ouf)
        ouf.close()
        
        outfile = 'boulder_feb_data_2019_lead='+clead+'h.cPick'
        ouf = open(outfile,'wb')
        cPickle.dump(np.array(acollect_feb), ouf)
        cPickle.dump(np.array(fcollect_feb), ouf)
        ouf.close()
    
    return beta_3d