#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 10:40:56 2018

@author: zack
"""
#import numpy as np

time = [0, 27*8, 400, 800]
stepsize = 8
d = 0.4/100
r = d/2
nu = .2
rhof = 2650
rhos = 3300
delta_rho = rhos - rhof
g = 9.81

set_v = (2/9)*(delta_rho*g*r**2)/(nu)
dist = 2.56
tau = []

for i in range(0, len(time)):
    tau.append((time[i]/stepsize)/(dist/set_v))
    
print (tau)