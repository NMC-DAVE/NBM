import sys
from fcst_error_stats import fcst_error_stats

clevels = ['925','850','800','750','700','600','500','400','300','250','200','150','100','10']
for input_hour in ['006','024','048','072','096','120']:
    #for input_stream in ['1999','2003','2007','2011']:
    for input_stream in ['2015']:
        for input_model in ['cfsr','reanalysis']:
            for input_variable in ['T','U','V','Z']:
                for input_level in clevels:
                    print (input_stream,input_model,input_variable,input_hour,input_level)
                    istat = fcst_error_stats(input_stream, input_model, \
                        input_variable, input_level, input_hour)
                        
#clevels = ['200']
#for input_stream in ['1999','2003','2007','2011']:
#    for input_model in ['reanalysis']:
#        for input_variable in ['T']:
#            for input_hour in ['006']:
#                for input_level in clevels:
#                    print (input_stream,input_model,input_variable,input_hour,input_level)
#                    istat = fcst_error_stats(input_stream, input_model, \
#                        input_variable, input_level, input_hour)                                                  
                        
print ('done')