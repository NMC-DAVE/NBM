def fcst_error_stats(input_stream, input_model, input_variable, input_level, input_hour):
    import csv
    import sys
    import _pickle as cPickle
    import numpy as np
    from datetime import datetime

    input_file1 = input_stream+'_'+input_model+'_'+input_hour+'.txt'
    print ('input_file1 = ', input_file1)

    output_file = 'reanalysis/'+input_model+'_stream='+\
        input_stream+'_lead='+input_hour+'_var='+\
        input_variable+'_level='+input_level+'.cPick'


    ICdate_nh = []
    bias_nh = []
    rmse_nh = []

    ICdate_sh = []
    bias_sh = []
    rmse_sh = []

    ICdate_tr = []
    bias_tr = []
    rmse_tr = []

    ICdate_gl = []
    bias_gl = []
    rmse_gl = []

    #The format of the text files is as follows.
    #  Column 1  - date
    #  Column 2  - date also
    #  Column 3  - hour of data
    #  Column 4  - pressure level
    #  Column 5  - field
    #  Column 6  - mean value of the field
    #  Column 7  - mean vlaue of the ERA-interim data used as the analysis
    #  Column 8  - bias between model and ERA-interim
    #  Column 9  - RMSE between model and ERA-interim
    #  Column10  - anomaly correlation
    #  Column11  - standard deviation of the model
    #  Column12  - standard deviation of ERA-interim data
    #  Column13  - the number of grid points
    #  Column14  - the region
    #              GL  - Global data
    #              NH  - northern hemisphere data
    #              SH  - southern hemisphere data
    #              TR  - tropics data

    input_file1_full = 'reanalysis/'+input_file1
    print (input_file1_full)
    #sys.exit()
    with open(input_file1_full, 'r') as infile1:
        csv_reader1 = csv.reader(infile1, delimiter=' ', skipinitialspace=True)
        for line in csv_reader1:
            #print (line)
            cplevel = line[3]
            #print (cplevel)
            iplevel = int(cplevel[1:4])
            #print (iplevel, input_level)
            #print (line[13], line[7], line[8])
            if iplevel == int(input_level) and input_variable == line[4]:
                if line[13] == 'GL' and line[7] != 'NA' and line[8] != 'NA':
                    #print ('global append')
                    ICdate_gl.append(line[0])
                    bias_gl.append(float(line[7]))
                    rmse_gl.append(float(line[8]))
                elif line[13] == 'NH' and line[7] != 'NA' and line[8] != 'NA':
                    ICdate_nh.append(line[0])
                    bias_nh.append(float(line[7]))
                    rmse_nh.append(float(line[8]))
                elif line[13] == 'SH' and line[7] != 'NA' and line[8] != 'NA':
                    ICdate_sh.append(line[0])
                    bias_sh.append(float(line[7]))
                    rmse_sh.append(float(line[8]))
                elif line[13] == 'TR' and line[7] != 'NA' and line[8] != 'NA':
                    ICdate_tr.append(line[0])
                    bias_tr.append(float(line[7]))
                    rmse_tr.append(float(line[8]))
    infile1.close()   
        
    #print (bias_gl)
    #print (rmse_gl)
    ICdate_gl_arr = np.squeeze(np.array(ICdate_gl))
    bias_gl_arr = np.squeeze(np.array(bias_gl))
    rmse_gl_arr = np.squeeze(np.array(rmse_gl))
    ICdate_nh_arr = np.squeeze(np.array(ICdate_nh))
    bias_nh_arr = np.squeeze(np.array(bias_nh))
    rmse_nh_arr  = np.squeeze(np.array(rmse_nh))
    ICdate_sh_arr = np.squeeze(np.array(ICdate_sh))
    bias_sh_arr = np.squeeze(np.array(bias_sh))
    rmse_sh_arr = np.squeeze(np.array(rmse_sh))
    ICdate_tr_arr = np.squeeze(np.array(ICdate_tr))
    bias_tr_arr = np.squeeze(np.array(bias_tr))
    rmse_tr_arr = np.squeeze(np.array(rmse_tr))
    
    print ('min, max bias_gl_arr = ', np.min(bias_gl_arr), np.max(bias_gl_arr))
    print ('min, max rmse_gl_arr = ', np.min(rmse_gl_arr), np.max(rmse_gl_arr))
    print ('mean bias gl nh sh tr = ', np.mean(bias_gl_arr), \
        np.mean(bias_nh_arr), np.mean(bias_tr_arr), np.mean(bias_sh_arr))
    print ('mean rmse gl nh sh tr = ', np.mean(rmse_gl_arr), \
        np.mean(rmse_nh_arr), np.mean(rmse_tr_arr), np.mean(rmse_sh_arr))
    
    ouf = open(output_file, 'wb')
    print ('output_file = ',output_file)
    cPickle.dump(ICdate_gl_arr, ouf)
    cPickle.dump(bias_gl_arr, ouf)
    cPickle.dump(rmse_gl_arr, ouf) 
    cPickle.dump(ICdate_nh_arr, ouf)
    cPickle.dump(bias_nh_arr, ouf) 
    cPickle.dump(rmse_nh_arr, ouf) 
    cPickle.dump(ICdate_sh_arr, ouf) 
    cPickle.dump(bias_sh_arr, ouf) 
    cPickle.dump(rmse_sh_arr, ouf) 
    cPickle.dump(ICdate_tr_arr, ouf)
    cPickle.dump(bias_tr_arr, ouf) 
    cPickle.dump(rmse_tr_arr, ouf)
    ouf.close() 
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    #sys.exit()
    
    
    
    istat = 0
    return istat

        
        
