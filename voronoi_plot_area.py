# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
###############################################################################
# voronoi_plot_area
#
# This code takes the particle positions and creates a Voronoi diagram.
# The area of the diagram is found and Eigen values and vectors are 
# used as proxies for the aspect ratio.
#
# date created:
# 01/02/2018
#
# Author:
# Zack McIntire
#
###############################################################################
# importing the necessary python modules

import csv
import matplotlib.pyplot as plt
import numpy as np
#import matplotlib as mpl
#import matplotlib.cm as cm
#import mplstereonet
from scipy.spatial import Voronoi, voronoi_plot_2d

###############################################################################
# Variables and file information

csvarray = []
data_array = []
filelocation = '/Users/zack/Documents/csv_data_files/'
#filename = 'test_case_20'
filename = 'box512_ordered_loc'
timestep = '0'

# defining the bounding region:
x_min = 0
x_max = 512
y_min = 0
y_max = 256

###############################################################################
# Opening the necessary files

#with open(filelocation + filename + '.csv', 'r') as csvfile:
#    csv_data = list(csv.reader(csvfile, delimiter=","))

with open(filelocation + filename + '.' + timestep + '.csv', 'r') as csvfile:
    csv_data = list(csv.reader(csvfile, delimiter=","))
#print (csv_data[:3])
csvarray = np.array(csv_data[1:], dtype=np.float) # Remove headers and convert to floats
#print (csvarray[:3])

data_array = np.array(csvarray)

#for i in range(0, len(csvarray)):
#    if csvarray[i, 7] < x_max and csvarray[i, 7] > x_min and \
#    csvarray[i, 8] < y_max and csvarray[i, 8] > y_min:
#        data_array.append(csvarray[i])
#data_array = np.array(data_array)

###############################################################################
# Finding the Voronoi points and plotting the diagram
 
points = []
for i in range(0, len(data_array)):
    points.append([data_array[i, 7], data_array[i, 8]])
points = np.array(points)

vor = Voronoi(points, qhull_options='Qbb Qc Qx')
voronoi_plot_2d(vor, show_vertices = False)
plt.show()

vor_vert = vor.vertices
position = [-1]
for index in range(0, len(vor_vert)):
    if vor_vert[index, 0] > x_max or vor_vert[index, 0] < x_min \
    or vor_vert[index, 1] > y_max or vor_vert[index, 1] < y_min:
        position.append(index)
position = np.array(position)


vor_reg = vor.regions

###############################################################################
# making sure the Voronoi vertices are within the domain boundaries.

reg = []
for i in range(0, len(vor_reg)):
    temp_reg = vor_reg[i]
    ind = 0
    for n in range(0, len(temp_reg)):
        for index in position:
            if index == temp_reg[n]:
                ind = 1
    if ind != 1:
        reg.append(temp_reg)


#x = [vor_vert[13, 0], vor_vert[6, 0], vor_vert[5, 0], vor_vert[4, 0],\
#     vor_vert[12, 0], data_array[0, 7]]
#y = [vor_vert[13, 1], vor_vert[6, 1], vor_vert[5, 1], vor_vert[4, 1],\
#     vor_vert[12, 1], data_array[0, 8]]
#fig = plt.figure()
#plt.scatter(x, y)

###############################################################################
# Finding the area of the Voronoi polygons

def PolygonSort(corners):
    n = len(corners)
    cx = float(sum(x for x, y in corners)) / n
    cy = float(sum(y for x, y in corners)) / n
    cornersWithAngles = []
    for x, y in corners:
        an = (np.arctan2(y - cy, x - cx) + 2.0 * np.pi) % (2.0 * np.pi)
        cornersWithAngles.append((x, y, an))
    cornersWithAngles.sort(key = lambda tup: tup[2])
    return cornersWithAngles

def PolygonArea(corners):
    n = len(corners)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

area = []
eigvalarray = []
eigvecarray = []
#max_x = 512
#max_y = 256
for i in range(0, len(reg)):
    corners = []
    vor_reg_temp = np.array(reg[i])
    for n in vor_reg_temp:
        corners.append(vor_vert[n])
    corners = np.array(corners)
#    for k in range(0, len(corners)):
#        if corners[k, 0] < max_x:
#            max_x = corners[k, 0]
#        if corners[k, 1] < max_y:
#            max_y = corners[k, 1]
    xy = corners.T
    eigvals, eigvecs = np.linalg.eig(np.cov(xy))
    eigvalarray.append(eigvals)
    eigvecarray.append(eigvecs)
    corners_sorted = PolygonSort(corners)
    area.append(PolygonArea(corners_sorted))
    
data_mean_area = np.mean(area)
    
###############################################################################
# creating a random array of particles for an average area of Voronoi cells to 
# normalize to.

num_of_part = len(csvarray)
y_mid = np.mean(csvarray[:, 8])
y_top = max(csvarray[:, 8])
y_height = (y_top - y_mid) * 2

rand_array = np.random.random((num_of_part, 2))
rand_array[:, 0] = rand_array[:, 0] * 512
rand_array[:, 1] = rand_array[:, 1] * y_height

rand_vor = Voronoi(rand_array, qhull_options='Qbb Qc Qx')

rand_vor_vert = rand_vor.vertices

rand_position = [-1]

for indexr in range(0, len(rand_vor_vert)):
    if rand_vor_vert[indexr, 0] > x_max or rand_vor_vert[indexr, 0] < x_min \
    or rand_vor_vert[indexr, 1] > y_max or rand_vor_vert[indexr, 1] < y_min:
        rand_position.append(indexr)
        
rand_position = np.array(rand_position)


rand_vor_reg = rand_vor.regions

rand_reg = []
for i in range(0, len(rand_vor_reg)):
    temp_reg = rand_vor_reg[i]
    ind = 0
    for n in range(0, len(temp_reg)):
        for index in rand_position:
            if index == temp_reg[n]:
                ind = 1
    if ind != 1:
        rand_reg.append(temp_reg)

rand_area = []
rand_eigvalarray = []
rand_eigvecarray = []
#max_x = 512
#max_y = 256
for i in range(0, len(rand_reg)):
    rand_corners = []
    vor_reg_temp = np.array(rand_reg[i])
    for n in vor_reg_temp:
        rand_corners.append(rand_vor_vert[n])
    rand_corners = np.array(rand_corners)
#    for k in range(0, len(corners)):
#        if corners[k, 0] < max_x:
#            max_x = corners[k, 0]
#        if corners[k, 1] < max_y:
#            max_y = corners[k, 1]
    xy = rand_corners.T
    eigvals, eigvecs = np.linalg.eig(np.cov(xy))
    rand_eigvalarray.append(eigvals)
    rand_eigvecarray.append(eigvecs)
    rand_corners_sorted = PolygonSort(rand_corners)
    rand_area.append(PolygonArea(rand_corners_sorted))
    
rand_mean_area = np.mean(rand_area)

norm_area = []
for i in range(0, len(area)):
    norm_area.append(area[i]/rand_mean_area)

###############################################################################    
# plotting the normalized cell area as a pdf plot
    
fig = plt.figure()
plt.hist(norm_area, bins = 'auto', density=True ) 
plt.title('pdf of normalized cell area')
plt.xlabel('normalized Voronoi cell area')
plt.ylabel('P.D.F')
plt.show()   

eigvalarray = np.array(eigvalarray)
eigvecarray = np.array(eigvecarray)

###############################################################################
# finding the aspect ratio of the eigen vectors

a1 = []
a2 = []
ar = []
for i in range(0, len(eigvalarray)):
    a1.append(eigvalarray[i, 0] * np.sqrt(eigvecarray[i, 0, 0]**2 + \
              eigvecarray[i, 0, 1]**2))
    a2.append(eigvalarray[i, 1] * np.sqrt(eigvecarray[i, 1, 0]**2 + \
              eigvecarray[i, 1, 1]**2))
    
a1 = np.array(a1)
a2 = np.array(a2)

al = []
ashort = []
for i in range(0, len(a1)): #a1 and j in a2:
    if a1[i] >= a2[i]:
        al.append([1, a1[i]])
        ashort.append([1, a2[i]])
    else:
        al.append([2, a2[i]])
        ashort.append([2, a1[i]])
        
al = np.array(al)
ashort = np.array(ashort)
        
for i in range(0, len(al)):
    ar.append(al[i, 1] / ashort[i, 1])
    
ar = np.array(ar)

###############################################################################
# attempted to color the Voronoi cell by aspect ratio

#minima = min(ar)
#maxima = max(ar)
#
#
#norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
#mapper = cm.ScalarMappable(norm=norm, cmap=cm.Blues_r)
#
#
#voronoi_plot_2d(vor, show_points=True, show_vertices=False, s=1)
#for r in range(len(reg)):
#    polygon = [vor_vert[i] for i in reg]
#    plt.fill(*zip(*polygon), color=mapper.to_rgba(ar[r]))
#plt.show()

###############################################################################
# plotting the pdf of the aspect ratio.

fig = plt.figure()
plt.hist(ar, bins = 'auto', density=True)
plt.title('pdf of aspect ratio')
plt.xlabel('Voronoi cell aspect ratio')
plt.ylabel('P.D.F')
plt.show()

###############################################################################
# calculating the angle of the longest eigen vector and creating a rose diagram
# of the orientation of that vector.

theta = []
        
for i in range(0, len(al)):
    if al[i, 0] == 1:
        theta.append(np.arctan(eigvecarray[i, 0, 1] / eigvecarray[i, 0, 0]))
    elif al[i, 0] == 2:
        theta.append(np.arctan(eigvecarray[i, 1, 1] / eigvecarray[i, 1, 0]))
theta = np.array(theta)

# making sure the theta values can be plotted on the rose diagram. increasing 
# values below the rose diagram threshould by 2pi
for i in range(0, len(theta)):
    if theta[i] <= np.deg2rad(-5):
        theta[i] = theta[i] + (2 * np.pi)
    else:
        theta[i] = theta[i]

# creating a theta array that is the opposite of the original data so the 
# diagram is symmetric       
theta_reverse = []
for i in range(0, len(theta)):
    theta_reverse.append(theta[i] + np.pi)
theta_reverse = np.array(theta_reverse)

# making sure there are no values above the maximum threshold of the rose plot
for i in range(0, len(theta_reverse)):
    if theta_reverse[i] > 2*np.pi:
        theta_reverse[i] = theta_reverse[i] - (2 * np.pi)
    else:
        theta_reverse[i] = theta_reverse[i]
        
# combining the two arrays        
new_theta = np.concatenate([theta, theta_reverse])

# creating the bins and bin values
bin_edge = np.deg2rad(np.arange(-5, 360, 10))
number_of_theta, bin_edge = np.histogram(new_theta, bin_edge)
number_of_theta = np.array(number_of_theta)

#for i in range(0, len(number_of_theta)):
#    number_of_theta[i] = number_of_theta[i] / 2.0

# plotting the rose diagram
fig = plt.figure()

ax = plt.subplot(111, projection='polar')

ax.bar(np.deg2rad(np.arange(0, 360, 10)), number_of_theta/2, width=np.deg2rad(10),\
       bottom=0.0, color='g', edgecolor='k')
ax.set_theta_zero_location('E')
ax.set_theta_direction(-1)
ax.set_title('Rose Diagram of the "polygon orientation"', y=1.10, fontsize=15)
plt.show()

