#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 16:13:44 2018

@author: zack
"""
import matplotlib.pyplot as plt
import numpy as np

r9osm = [1.29, 1.21, 1.22]
r9sm = [1.29, 1.21, 1.07]
r1osm = [1.04, 0.86, 0.33]
r20osm = [1.52, 1.43, 1.46]

phi9osm = [9.3, 9.5, 10.5]
phi9sm = [9.3, 9.5, 9.9]
phi1osm = [1.2, 1.4, 9.5]
phi20osm = [20, 20, 20.8]

m = 0.66 / 60
y = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
x = []
for i in y:
    x.append((i - 1.05) / m)
x = np.array(x)

fig = plt.figure()
plt.scatter(phi9osm, r9osm, s=160, marker='*', c='r', label='9% OSM')
plt.scatter(phi9sm, r9sm, s=80, marker='*', c='b', label='9% SM')
plt.scatter(phi1osm, r1osm, s=160, marker='*', c='g', label='1% OSM')
plt.scatter(phi20osm, r20osm, s=160, marker='*', c='y', label='20% OSM')
plt.plot(x, y, 'k-', label='RSDL')
plt.legend(loc = 'best')
plt.xlim(0, 25)
plt.ylim(0, 1.6) 
plt.show()