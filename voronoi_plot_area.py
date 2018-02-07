# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
#import mplstereonet
from scipy.spatial import Voronoi, voronoi_plot_2d

csvarray = []
data_array = []
filelocation = '/Users/zack/Documents/csv_data_files/'
#filename = 'test_case_20'
filename = 'box512_ht_1_loc'
timestep = '0'

# defining the bounding region:
x_min = 0
x_max = 512
y_min = 0
y_max = 256

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
for i in range(0, len(reg)):
    corners = []
    vor_reg_temp = np.array(reg[i])
    for n in vor_reg_temp:
        corners.append(vor_vert[n])
    corners = np.array(corners)
    xy = corners.T
    eigvals, eigvecs = np.linalg.eig(np.cov(xy))
    eigvalarray.append(eigvals)
    eigvecarray.append(eigvecs)
    corners_sorted = PolygonSort(corners)
    area.append(PolygonArea(corners_sorted))
    
mean_area = np.mean(area)

norm_area = []
for i in range(0, len(area)):
    norm_area.append(area[i]/mean_area)
fig = plt.figure()
plt.hist(norm_area, bins = 'auto', density=True ) 
plt.title('pdf of normalized cell area')
plt.xlabel('normalized Voronoi cell area')
plt.ylabel('P.D.F')
plt.show()   

eigvalarray = np.array(eigvalarray)
eigvecarray = np.array(eigvecarray)

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
    ar.append(ashort[i, 1] / al[i, 1])
    
ar = np.array(ar)

fig = plt.figure()
plt.hist(ar, bins = 'auto', density=True)
plt.title('pdf of aspect ratio')
plt.xlabel('Voronoi cell aspect ratio')
plt.ylabel('P.D.F')
plt.show()

theta = []
        
for i in range(0, len(al)):
    if al[i, 0] == 1:
        theta.append(np.arctan(eigvecarray[i, 0, 1] / eigvecarray[i, 0, 0]))
    elif al[i, 0] == 2:
        theta.append(np.arctan(eigvecarray[i, 1, 1] / eigvecarray[i, 1, 0]))
theta = np.array(theta)

for i in range(0, len(theta)):
    if theta[i] <= np.deg2rad(-5):
        theta[i] = theta[i] + (2 * np.pi)
    else:
        theta[i] = theta[i]
        
theta_reverse = []
for i in range(0, len(theta)):
    theta_reverse.append(theta[i] + np.pi)
theta_reverse = np.array(theta_reverse)

for i in range(0, len(theta_reverse)):
    if theta_reverse[i] > 2*np.pi:
        theta_reverse[i] = theta_reverse[i] - (2 * np.pi)
    else:
        theta_reverse[i] = theta_reverse[i]
        
new_theta = np.concatenate([theta, theta_reverse])

bin_edge = np.deg2rad(np.arange(-5, 360, 10))
number_of_theta, bin_edge = np.histogram(new_theta, bin_edge)
number_of_theta = np.array(number_of_theta)

#for i in range(0, len(number_of_theta)):
#    number_of_theta[i] = number_of_theta[i] / 2.0

fig = plt.figure()

ax = plt.subplot(111, projection='polar')

ax.bar(np.deg2rad(np.arange(0, 360, 10)), number_of_theta/2, width=np.deg2rad(10),\
       bottom=0.0, color='.8', edgecolor='k')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title('Rose Diagram of the "polygon orientation"', y=1.10, fontsize=15)
plt.show()

