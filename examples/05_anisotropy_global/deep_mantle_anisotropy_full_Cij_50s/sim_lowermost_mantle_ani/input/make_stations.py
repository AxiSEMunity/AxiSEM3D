#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 13:17:42 2021

@author: jw2449
"""

import numpy as np
import os
import obspy.taup
from obspy.taup import TauPyModel
import obspy.taup.taup_geo 
import obspy.signal
import os
import shutil
from obspy.clients.iris import Client
client = Client()
model = TauPyModel(model="prem")

lons = [90.0,75.0,60.0,45.0,30.0,15.0,0.0,-15.0,-30.0,-45.0,-60.0,-75.0]
lonsm = ['90','75','60','45','30','15','0','m15','m30','m45','m60','m75']

lats = np.arange(-30,46) + 0.0
latsm = []

for lat in lats:
    if lat > 0:
        s = str(int(lat))
    else:
        s = 'm'+str(int(np.abs(lat)))
    latsm.append(s)
    
    
f_write = open("stations.txt", 'a')
for i in range(len(lats)):
    for j in range(len(lons)):
        f_write.write('L'+latsm[i]+'L'+lonsm[j]+" \t"+"JW \t"+str(lats[i])+" \t"+str(lons[j])+" \t"+"0.0"+" \t"+"0.0"+"\n")
                
f_write.close()
