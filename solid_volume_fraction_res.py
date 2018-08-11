#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 15:38:15 2018

@author: zack
"""

import csv
import matplotlib.pyplot as plt
import numpy as np

timestep = '800'
filelocation = '/wd2/csv_data_files/'
resfile = 'box512_ht_20_res'

with open(filelocation + resfile + '.' + timestep + '.csv', 'r') as csvres:
    res_data = list(csv.reader(csvres, delimiter=","))
resarray = np.array(res_data[1:], dtype=np.float)

reduced_res_array = []
for index in range(0, len(resarray)):
    if 1 - resarray[index, 0] < 0.46 and 1 - resarray[index, 0] > 0:
        reduced_res_array.append(resarray[index, :])

reduced_res_array = np.array(reduced_res_array)

solid_volume_fraction = 1 - np.mean(reduced_res_array[:, 0])