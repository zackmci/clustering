#######################################################################################################################
#  clustering.py
#
#  Calculating R, L(r), and g(r)
#
#  Date 12/20/2017 
#
#  Last Modified 1/25/2018
#
#
#######################################################################################################################
# importing the needed modules

import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata

#######################################################################################################################
# variables for reading and writing files

csvarray = []
<<<<<<< HEAD
timestep = '160'
filelocation = '/Users/zack/Documents/csv_data_files/'
filename = 'box512_ht_1_loc'
resfile = 'box512_ht_1_res'
=======
timestep = '80'
filelocation = '/wd2/csv_data_files/'
savelocation = '/home/zack/Documents/csv_data_files/'
filename = 'box512_ht_loc'
resfile = 'box512_ht_res'
>>>>>>> a95313cb3d6a7e89d11e40a8ffbe2df22f894479
tablename = '_minimum_distance'
lrdata = '_lr_kr_'
grdata = '_gr_'

#######################################################################################################################
# Opening and reading the csv files from the simulations

# particle data
with open(filelocation + filename + '.' + timestep + '.csv', 'r') as csvfile:
    csv_data = list(csv.reader(csvfile, delimiter=","))
#print (csv_data[:3])
csvarray = np.array(csv_data[1:], dtype=np.float) # Remove headers and convert to floats
#print (csvarray[:3])

original_part_num = len(csvarray)

# fluid data
with open(filelocation + resfile + '.' + timestep + '.csv', 'r') as csvres:
    res_data = list(csv.reader(csvres, delimiter=","))
resarray = np.array(res_data[1:], dtype=np.float)

#######################################################################################################################
# Seperating the particles into two arrays based on volume fraction.  One array is the 
# settled bed, and the other is of the particles still settling.

vol_frac = 0.46
box_dist = 5
part_high_vol_array = []
part_low_vol_array = []

for posindex in range(0, len(csvarray)):
   # if posindex % 1000 == 0:
    #    print (posindex)
    volarray = []
    for vfindex in range(0, len(resarray)):
        # creating a smaller fluid array based on how close the corner is to the particle.
        if resarray[vfindex, 6] < csvarray[posindex, 7] + box_dist and \
        resarray[vfindex, 6] > csvarray[posindex, 7] - box_dist and \
        resarray[vfindex, 7] < csvarray[posindex, 8] + box_dist and \
        resarray[vfindex, 7] > csvarray[posindex, 8] - box_dist:
            volarray.append(resarray[vfindex])
    volarray = np.array(volarray)
    distarray = []
    ref_dist = 2
    for volindex in range(0, len(volarray)):
        dist_to_corner = np.sqrt((csvarray[posindex, 7] - volarray[volindex, 6])**2 +\
                                (csvarray[posindex, 8] - volarray[volindex, 7])**2)
        if dist_to_corner < ref_dist:
            ref_dist = dist_to_corner
            distarray = volarray[volindex]
    distarray = np.array(distarray)
    if 1-distarray[0] >= vol_frac:
        part_high_vol_array.append(csvarray[posindex]) # Settled bed
    else:
        part_low_vol_array.append(csvarray[posindex]) # Still settling

csvarray = np.array(part_low_vol_array)
part_high_vol_array = np.array(part_high_vol_array)

x_len = 512 # Length of x axis
y_len = max(csvarray[:, 8]) # height of the heighest particle.

# Some more variables to be taken from the new array
new_part_num = len(csvarray)
difference = original_part_num - new_part_num
x = csvarray[:,7]
y = csvarray[:,8]
rad = np.mean(csvarray[:, 1])/2

#######################################################################################################################
# Calculating the area of the settled bed
if len(part_high_vol_array) > 0:
    part_west = []
    part_east = []
for pos_index in range(0, len(part_high_vol_array)):
    if part_high_vol_array[pos_index, 7] < 250:
        part_west.append(part_high_vol_array[pos_index])
    else:
        part_east.append(part_high_vol_array[pos_index])
if len(part_high_vol_array) > 0:
    if len(part_west) > 0:
        part_west = np.array(part_west)
        min_x_west = min(part_west[:, 7])
        max_x_west = max(part_west[:, 7])
        max_y_west = max(part_west[:, 8])
        area_west = (max_x_west - min_x_west) * max_y_west
    if len(part_east) > 0:
        part_east = np.array(part_east)
        min_x_east = min(part_east[:, 7])
        max_x_east = max(part_east[:, 7])
        max_y_east = max(part_east[:, 8])
        area_east = (max_x_east - min_x_east) * max_y_east
    if len(part_west) > 0 and len(part_east) > 0:
        total_area = area_west + area_east
    elif len(part_west) == 0 and len(part_east) > 0:
        total_area = area_east
    elif len(part_west) > 0 and len(part_east) ==  0:
        total_area = area_west
else:
    total_area = 0

#######################################################################################################################
# A function that creates an array of the particles within a set distance from the selected particle

def domain_reduction(arr, j, x, y, radius):
    arr_red = []
    x_pos = arr[j, x]
    y_pos = arr[j, y]
    #print (x_pos)
    for i in range(0, len(arr)):
        if arr[i, x] - arr[i, 1]/2 < x_pos + arr[j, 1]/2 + radius and \
        arr[i, x] + arr[i, 1]/2 > x_pos - arr[j, 1]/2 - radius:
            if arr[i, y] < y_pos + radius and arr[i, y] > y_pos - radius:
#             if arr[i, y] - arr[i, 1]/2 < y_pos + arr[j, 1]/2 + radius and arr[i, y] + arr[i, 1]/2 > \
#             :y_pos - arr[j, y]/2 - radius:
                arr_red.append(arr[i])
    arr_red = np.array(arr_red)
    #print (arr_red)
    return arr_red

#######################################################################################################################
# Creating an array of the nearest neighbors and which particles are the nearest.

#print (num_of_it)
min_array = []

for i in range(0, new_part_num):

#    if i % 1000 == 0:
#        print (i)
#    print (i)
    near_neigh = []
    #print (near_neigh)

    # Creating an array that is reduced to particles within a set radius.
    red_arr = domain_reduction(csvarray, i, 7, 8, 50)
    #print (len(red_arr))
    #print (red_arr)

    for n in range(0, len(red_arr)):
        if csvarray[i, 0] != red_arr[n, 0]:
            distance = np.sqrt((csvarray[i, 7] - red_arr[n, 7])**2 + (csvarray[i, 8] - red_arr[n, 8])**2)
#             distance = np.sqrt((csvarray[i, 7] - red_arr[n, 7])**2 + (csvarray[i, 8] - red_arr[n, 8])**2) - \
#             ((csvarray[i, 1]/2) + (red_arr[n, 1]/2))
            arr = np.array([csvarray[i, 0], red_arr[n, 0], distance])
            #print (arr)
            near_neigh.append(arr)

    near_neigh = np.array(near_neigh) # Converts the list into an array
    #print (near_neigh)
    min_dist = np.amin(near_neigh[:, 2])
    min_index = np.argmin(near_neigh[:, 2]) # Gives the index of the minimum value in the near_neigh array
    #print (min_index)
    min_array.append(near_neigh[min_index])
min_array = np.array(min_array)
#print (min_array)

#######################################################################################################################
# Writing nearest neighbor data to csv

with open(savelocation + filename + tablename + '.' + timestep + '.csv', 'w', newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['original particle', 'nearest neighbor', 'distance'])
    writer.writerows(min_array)

#######################################################################################################################
# Finding the mean distance rA

#print(min_array[:, 2])
rA = np.mean(min_array[:, 2])

#######################################################################################################################
# Calculating rE the expected value

# Intensity with the area of the settled particles taken into consideration
lam = new_part_num / ((x_len * y_len) - total_area) 

rE = 1 / (2 * np.sqrt(lam))

#######################################################################################################################
# Calculating the R value rA/rE

R = rA / rE

#######################################################################################################################
# Writing rA, rE, and R to a csv file

with open(savelocation + filename + '_R.' + timestep + '.csv', 'w', newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['rA', 'rE', 'R'])
    writer.writerow([rA, rE, R])

#######################################################################################################################
# The purpose of this function is to calculate the distance between the selected particle and all the other particles.
# To do this an array of particles is created where the particles are all within a distance slightly greater than the
# radius of interest.  Anoter array is created that has only the particles within the defined radius, and the number
# of particles within that array is returned to the calling script.

def part_in_rad(arr, part_num, x_val, y_val, radius):

    num_part = []
    smallarr = []

    # creating the box around the selected particle and only keeping the particles within that box.
    for k in range(0, len(arr)):
        if arr[part_num, x_val] - radius - 0.5 < arr[k, x_val] and \
        arr[part_num, x_val] + radius + 0.5 > arr[k, x_val]  and \
        arr[part_num, y_val] - radius - 0.5 < arr[k, y_val] and \
        arr[part_num, y_val] + radius + 0.5 > arr[k, y_val] and part_num != k:
            smallarr.append(arr[k])
            #print (smallarr)
    smallarr = np.array(smallarr)
    #print (smallarr)

    # Finding the particles within the radius of the selected particle.
    for i in range(0, len(smallarr)):
        #print (i)
        dist = np.sqrt((arr[part_num, x_val] - smallarr[i, x_val])**2 + \
        (arr[part_num, y_val] - smallarr[i, y_val])**2)
        if dist <= radius and dist > 0:
            #print (arr[i, :])
            num_part.append(smallarr[i])
            #print(num_part)
        #print (num_part)
    num_part = np.array(num_part)
    #print (num_part)
    #print (len(num_part))
    Np = len(num_part)
    #print (Np)
#    return num_part

    # Returning the number of particles within the radius.
    if len(num_part) == 0:
        return 0
    else:
        return Np

#######################################################################################################################
# Calculating Ripley's K and L values.

Kr = []

radius = [rad]
index_rad = 0
while (radius[index_rad] <(y_len - (total_area / (x_len - 32)))/2  - (rad * 4) and radius[index_rad] < 7):
    index_rad = index_rad + 1
    radius.append(rad * (index_rad + 1))
radius = np.array(radius)
#radius = [rad, rad*2, rad*3, rad*4, rad*5, rad*6, rad*7, rad*8, rad*9, rad*10, rad*11, rad*12, rad*13, rad*14, \
#         rad*15, rad*16, rad*17, rad*18, rad*19, rad*20]
#print (len(radius))
for i in radius:
    kr = []
    for n in range(0, new_part_num):
        #print (n)
        # This if statement is to correct for edge effects.  If the particle surface is within range of the edge
        # of the domain the particle is not used as a source, but can be counted as within the radial distance
        # of another particle.
        if len(part_high_vol_array) == 0:
            if csvarray[n, 7] >= i and csvarray[n, 7] <= x_len - i and \
            csvarray[n, 8] >= i and csvarray[n, 8] <= y_len - i:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                #arr.append(num_part)
                #np.append(arr, num_part)
                #print (arr)
                #particle = len(arr)
                kr.append(num_part / lam)
                #print (kr)
        elif len(part_west) == 0 and len(part_east) > 0:
            if csvarray[n, 7] >= min_x_east and csvarray[n, 7] <= x_len - i and \
            csvarray[n, 8] >= max_y_east + i and csvarray[n, 8] <= y_len - i:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                kr.append(num_part / lam)
            if csvarray[n, 7] >= i  and csvarray[n, 7] < min_x_east  and \
            csvarray[n, 8] >= max_y_east and csvarray[n, 8] <= y_len - i:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                kr.append(num_part / lam)
            if csvarray[n, 7] >= i and csvarray[n, 7] <= min_x_east - i and \
            csvarray[n, 8] >= i and csvarray[n, 8] < max_y_east:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                kr.append(num_part / lam)
        elif len(part_west) > 0 and len(part_east) == 0:
            if csvarray[n, 7] >= i and csvarray[n, 7] < max_x_west and \
            csvarray[n, 8] >= max_y_west + i and csvarray[n, 8] <= y_len - i:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                kr.append(num_part / lam)
            if csvarray[n, 7] >= max_x_west and csvarray[n, 7] < x_len - i  and \
            csvarray[n, 8] >= i and csvarray[n, 8] <= y_len - i:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                kr.append(num_part / lam)
            if csvarray[n, 7] >= max_x_west + i and csvarray[n, 7] <= x_len -  i and \
            csvarray[n, 8] >= i and csvarray[n, 8] < max_y_west:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                kr.append(num_part / lam)
        else:
            if csvarray[n, 7] >= i and csvarray[n, 7] <= max_x_west and \
            csvarray[n, 8] >= max_y_west + i and csvarray[n, 8] <= y_len - i:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                kr.append(num_part / lam)
            if csvarray[n, 7] > max_x_west and csvarray[n, 7] <= x_len - i and \
            csvarray[n, 8] >= max_y_east + i and csvarray[n, 8] <= y_len - i:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                kr.append(num_part / lam)
            if csvarray[n, 7] >= max_x_west + i and csvarray[n, 7] <= min_x_east - i and \
            csvarray[n, 8] >= i and csvarray[n, 8] < max_y_east:
                num_part = part_in_rad(csvarray, n, 7, 8, i)
                kr.append(num_part / lam)
    Kr.append(np.mean(kr))
Kr = np.array(Kr)

# This is to correct for overlapping particles where the distance between
# particles in the simulation could be less than 0.
Kr[0] = 0
#print (Kr)

Lr = np.sqrt(Kr / np.pi)

#######################################################################################################################
# Writing L(r) and K(r) to a csv file

with open(savelocation + filename + lrdata + '.' + timestep + '.csv', 'w', newline='') as f:
    writer=csv.writer(f)
    writer.writerows([Lr, Kr])

#######################################################################################################################
# Plotting L(r) vs radius

fig=plt.figure(figsize=(9,9))
plt.plot(radius, Lr, 'r*-')
plt.plot(radius, radius, 'k-')
plt.title("Ripley's L value", fontsize = 18 )
plt.xlabel('radius (r)', fontsize = 18)
plt.ylabel('L(r)', fontsize = 18)
plt.savefig(savelocation + filename + '_lr.' + timestep + '.svg', format='svg')

#######################################################################################################################
# Calulating Ripley's g value

gr = []

for i in range(1, len(radius)):
    dkr = (Kr[i] - Kr[i-1]) / (radius[i] - radius[i-1])
    #print (dkr)
    gr.append(dkr / (2 * np.pi * radius[i]))
gr = np.array(gr)

#######################################################################################################################
# writing g(r) to a csv file

with open(savelocation + filename + grdata + '.' + timestep + '.csv', 'w', newline='') as f:
    writer=csv.writer(f)
    writer.writerow(gr)

#######################################################################################################################
# Plotting g(r) vs radius

fig=plt.figure(figsize=(9,9))
plt.plot(radius[1:len(radius)], gr, 'bo-')
plt.plot(radius[1:len(radius)], np.ones(len(gr)), 'k-')
plt.title("Ripley's g value (pair correlation)", fontsize = 18)
plt.xlabel('radius (r)', fontsize = 18)
plt.ylabel('g(r)', fontsize = 18)
plt.savefig(savelocation + filename + '_gr.' + timestep + '.svg', format='svg')
