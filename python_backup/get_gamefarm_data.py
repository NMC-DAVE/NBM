""" get_gamefarm_data.py : from Art DeGaetano, an example script to access Ithaca Game farm 
    max, min, and precipitation data. """
import requests
import json
import numpy as np

input_dict = {"sid":"304174","sdate":"2000-01-01","edate":"2019-12-31","elems":"maxt,mint,pcpn"}
req = requests.post("http://data.rcc-acis.org/StnData", json = input_dict)
data = req.json()
data_vals = data['data']
print (np.shape(data_vals))
ndates, nvalues = np.shape(data_vals)
for idate in range(ndates):
    print (data_vals[idate][0:4])


