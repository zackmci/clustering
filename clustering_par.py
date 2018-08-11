#######################################################################################################################
#  clustering_par.py
#
#  Calculating R, L(r), and g(r) using parallel processing
#
#  Date 2/1/2018 
#
#  Last Modified 1/30/2018
#
# Author:
# Zack McIntire
#
# Imputs (CSV)
# particle location file
# fluid file
#
# Outputs:
#
# --arrays-- (CSV)
# - An array of the particle data for particles in the low particle volume 
#   fraction portions of the domain.
# - An array of the particled, nearest neighbor, and distance
# - The rA, rE, and R value for the timestep
# - Ripley's k value and L value
# - Ripley's g value 
#
# --Figures-- (SVG)
# - A plot of Ripley's L value over the radius away from the particle.
# -A plot of Ripley's g value over the radius away from the particle.
#
#
###############################################################################
# importing the needed modules

import csv
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
#import os
import time

###############################################################################
# variables for reading and writing files

csv_array = []
timestep = '0'
filelocation = '/wd2/csv_data_files/'
#filelocation = '/home/zack/Documents/csv_data_files/'
savelocation = '/home/zack/Documents/csv_data_files/'
figurelocation = '/home/zack/Documents/clustering_figures/'
filename = 'box512_ht_20_loc'
resfile = 'box512_ht_20_res'
tablename = '_minimum_distance'
lrdata = '_lr_kr_'
grdata = '_gr_'

###############################################################################
# Opening and reading the csv files from the simulations

# particle data
with open(filelocation + filename + '.' + timestep + '.csv', 'r') as csvfile:
    csv_data = list(csv.reader(csvfile, delimiter=","))
#print (csv_data[:3])

# Remove headers and convert to floats
csv_array = np.array(csv_data[1:], dtype=np.float) 
#print (csvarray[:3])

original_part_num = len(csv_array)

# fluid data
with open(filelocation + resfile + '.' + timestep + '.csv', 'r') as csvres:
    res_data = list(csv.reader(csvres, delimiter=","))
resarray = np.array(res_data[1:], dtype=np.float)

###############################################################################
# Seperating the particles into two arrays based on volume fraction. 
# One array is the settled bed, and the other is of the particles 
# still settling.

def settling_array(posindex):

    vol_frac = 0.46
    box_dist = 5

   # print (posindex)

#    for posindex in range(0, len(csvarray)):
       # if posindex % 1000 == 0:
        #    print (posindex)
    volarray = []
    for vfindex in range(0, len(resarray)):
        # creating a smaller fluid array based on how close the corner is 
        # to the particle.
        if resarray[vfindex, 6] < csv_array[posindex, 7] + box_dist and \
        resarray[vfindex, 6] > csv_array[posindex, 7] - box_dist and \
        resarray[vfindex, 7] < csv_array[posindex, 8] + box_dist and \
        resarray[vfindex, 7] > csv_array[posindex, 8] - box_dist:
            volarray.append(resarray[vfindex])
    volarray = np.array(volarray)    
    distarray = []
    ref_dist = 2
    for volindex in range(0, len(volarray)):
        dist_to_corner = np.sqrt((csv_array[posindex, 7] - \
                                  volarray[volindex, 6])**2 + \
                                  (csv_array[posindex, 8] - \
                                  volarray[volindex, 7])**2)
        if dist_to_corner < ref_dist:
            ref_dist = dist_to_corner
            distarray = volarray[volindex]
    distarray = np.array(distarray)
    if 1-distarray[0] >= vol_frac:
        return [csv_array[posindex, 0], 1]
    else:
        return [csv_array[posindex, 0], 0]
    
###############################################################################
# calling settling_array()

start = time.time()
print ('Looking for settled beds at: ', start)
pool = mp.Pool(processes=5)
settled_result = pool.map(settling_array, range(0, original_part_num)) 
settled_result = np.array(settled_result)
#print (result)
pool.close()
pool.join()
csvarray = []
part_high_vol_array = []
for i in range(0, len(csv_array)):
    for n in range(0, len(settled_result)):
        if csv_array[i, 0] == settled_result[n, 0]:
            if settled_result[n, 1] == 1:
                part_high_vol_array.append(csv_array[i, :])
            else:
                csvarray.append(csv_array[i, :])
csvarray = np.array(csvarray)
part_high_vol_array = np.array(part_high_vol_array)

end = time.time()
#os.system('spd-say "first annoying noise"')
print ('particles seperated into high and low volume fraction arrays: ',\
       end - start)
###############################################################################
# writing the array of particles with low volume fractions to file only if
# there is an array of high volume fraction particles.

if len(part_high_vol_array) > 0:
    with open(savelocation + filename + '_low_vol_frac_2.' + timestep + '.csv', \
              'w', newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['ID', 'Diameter', 'Density', 'Velocity:0', \
                         'Velocity:1', 'Velocity:2', 'Temperature', \
                         'coordinates:0', 'coordinates:1', 'coordinates:2' ])
        writer.writerows(csvarray)

###############################################################################
# variables

x_len = 512 # Length of x axis
y_len = max(csvarray[:, 8]) # height of the heighest particle.

# Some more variables to be taken from the new array
new_part_num = len(csvarray)
difference = original_part_num - new_part_num
x = csvarray[:,7]
y = csvarray[:,8]
rad = np.mean(csvarray[:, 1])/2

###############################################################################
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

###############################################################################
# A function that creates an array of the particles within a set distance from 
# the selected particle

def domain_reduction(arr, j, x, y, radius):
    arr_red = []
    x_pos = arr[j, x]
    y_pos = arr[j, y]
    #print (x_pos)
    while len(arr_red) == 0:
        for point in range(0, len(arr)):
            #print(point)
            if point != j:
                if arr[point, x] < x_pos + radius and arr[point, x] > x_pos - radius \
                and arr[point, y] < y_pos + radius and arr[point, y] > y_pos - radius:
                    arr_red.append(arr[point])
        if len(arr_red) == 0:
            radius = radius + 5
                    
    arr_red = np.array(arr_red)            
                    
    
                
    
    #print (arr_red)
    return arr_red

###############################################################################
# Creating an array of the nearest neighbors and which particles are
# the nearest.

def near_neigh(i):
    #print (i)

    #for i in range(0, new_part_num):

    #    if i % 1000 == 0:
    #        print (i)
    #    print (i)
    near_neighb = []
    #print (near_neigh)

    # Creating an array that is reduced to particles within a set radius.
    red_arr = domain_reduction(csvarray, i, 7, 8, 5)
    #print (len(red_arr))
    #print (red_arr)

    for n in range(0, len(red_arr)):
        if csvarray[i, 0] != red_arr[n, 0]:
            distance = np.sqrt((csvarray[i, 7] - red_arr[n, 7])**2 + \
                               (csvarray[i, 8] - red_arr[n, 8])**2)
            arr = np.array([csvarray[i, 0], red_arr[n, 0], distance])
            #print (arr)
            near_neighb.append(arr)

    near_neighb = np.array(near_neighb) # Converts the list into an array
    #print (near_neighb)
    
    # Gives the index of the minimum value in the near_neigh array
    min_index = np.argmin(near_neighb[:, 2]) 
    #print (min_index)
    neigh_array = near_neighb[min_index]
    #print (min_array)
    return neigh_array

###############################################################################
# calling near_neigh

start = time.time()
print ('Finding the nearest neighbor starting at: ', start)
neigh_pool = mp.Pool(processes=5)
neigh_result = neigh_pool.map(near_neigh, range(0, new_part_num))
neigh_result = np.array(neigh_result) 
neigh_pool.close()
neigh_pool.join()
min_array = np.array(neigh_result)
end = time.time()
#os.system('spd-say "annoying noise two"')
print ('Nearest neighbor found in: ', end - start)



###############################################################################
# Writing nearest neighbor data to csv

with open(savelocation + filename + tablename + '_2.' + timestep + '.csv', 'w', newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['original particle', 'nearest neighbor', 'distance'])
    writer.writerows(min_array)

###############################################################################
# Finding the mean distance rA

#print(min_array[:, 2])
rA = np.mean(min_array[:, 2])

###############################################################################
# Calculating rE the expected value

# Intensity with the area of the settled particles taken into consideration
lam = new_part_num / ((x_len * y_len) - total_area)  

rE = 1 / (2 * np.sqrt(lam))

###############################################################################
# Calculating the R value rA/rE

R = rA / rE
print ('R actual = ', rA, '\n', 'R estimated = ', rE, '\n', 'R = ', R)
###############################################################################
# Writing rA, rE, and R to a csv file

with open(savelocation + filename + '_R_2.' + timestep + '.csv', 'w', newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['rA', 'rE', 'R'])
    writer.writerow([rA, rE, R])

###############################################################################
# The purpose of this function is to calculate the distance between the 
# selected particle and all the other particles. To do this an array of
# particles is created where the particles are all within a distance slightly
# greater than the radius of interest.  Anoter array is created that has only
# the particles within the defined radius, and the number of particles within
# that array is returned to the calling script.

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
# defining the distance from particle to search

radius = [rad]
index_rad = 0
while radius[index_rad] < ((y_len - min(csvarray[:, 8])) / 2 - (rad * 4)) and radius[index_rad] < 7:
    index_rad = index_rad + 1
    radius.append(rad * (index_rad + 1))
radius = np.array(radius)
print ('The radial distances for number of particle search is: ', radius)

#######################################################################################################################
# Calculating Ripley's K and L values.

def ripley_k(i):

    kr = []
    print (radius[i])

    for n in range(0, new_part_num):
        # This if statement is to correct for edge effects.  If the particle 
        # surface is within range of the edge of the domain the particle is 
        # not used as a source, but can be counted as within the radial 
        # distance of another particle. To do this the domain is divided by 
        # the injection location, and looking at the array of particles that
        # had a solid volume fraction of greater than 0.46.
        if len(part_high_vol_array) == 0:
            if csvarray[n, 7] >= radius[i] and csvarray[n, 7] <= x_len - radius[i] and \
            csvarray[n, 8] >= radius[i] and csvarray[n, 8] <= y_len - radius[i]:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)
        elif len(part_west) == 0 and len(part_east) > 0:
            if csvarray[n, 7] >= min_x_east and csvarray[n, 7] <= x_len - radius[i] and \
            csvarray[n, 8] >= max_y_east + radius[i] and csvarray[n, 8] <= y_len - radius[i]:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)
            elif csvarray[n, 7] >= radius[i]  and csvarray[n, 7] < min_x_east  and \
            csvarray[n, 8] >= max_y_east and csvarray[n, 8] <= y_len - radius[i]:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)
            elif csvarray[n, 7] >= radius[i] and csvarray[n, 7] <= min_x_east - radius[i] and \
            csvarray[n, 8] >= radius[i] and csvarray[n, 8] < max_y_east:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)
        elif len(part_west) > 0 and len(part_east) == 0:
            if csvarray[n, 7] >= radius[i] and csvarray[n, 7] < max_x_west and \
            csvarray[n, 8] >= max_y_west + radius[i] and csvarray[n, 8] <= y_len - radius[i]:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)
            elif csvarray[n, 7] >= max_x_west and csvarray[n, 7] < x_len - radius[i]  and \
            csvarray[n, 8] >= radius[i] and csvarray[n, 8] <= y_len - radius[i]:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)
            elif csvarray[n, 7] >= max_x_west + radius[i] and csvarray[n, 7] <= x_len -  radius[i] and \
            csvarray[n, 8] >= radius[i] and csvarray[n, 8] < max_y_west:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)
        else:
            if csvarray[n, 7] >= radius[i] and csvarray[n, 7] <= max_x_west and \
            csvarray[n, 8] >= max_y_west + radius[i] and csvarray[n, 8] <= y_len - radius[i]:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)
            elif csvarray[n, 7] > max_x_west and csvarray[n, 7] <= x_len - radius[i] and \
            csvarray[n, 8] >= max_y_east + radius[i] and csvarray[n, 8] <= y_len - radius[i]:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)
            elif csvarray[n, 7] >= max_x_west + radius[i] and csvarray[n, 7] <= min_x_east - radius[i] and \
            csvarray[n, 8] >= radius[i] and csvarray[n, 8] < max_y_east:
                num_part = part_in_rad(csvarray, n, 7, 8, radius[i])
                kr.append(num_part / lam)

    Kr = np.mean(kr)

    return Kr

#######################################################################################################################
# calling ripleys_k and calculating Lr

start = time.time()
print ("Ripley's K value started at: ", start)
pool = mp.Pool(processes=5)
kr_result = pool.map(ripley_k, range(0, len(radius)))
kr_result = np.array(kr_result) 
pool.close()
pool.join()
Kr = np.array(kr_result)
end = time.time()
# os.system('spd-say "annoying noise again "')
print ("ripley's k value found in: ", end - start)

# This is to correct for overlapping particles where the distance between
# particles in the simulation could be less than 0.
Kr[0] = 0 # This is the distance from the particle center to the 
          # particle boundary.
Kr[1] = 0 # This is the distance to the next particle center.

Lr = np.sqrt(Kr / np.pi)
print ('k(r) = ', Kr, '\n', 'L(r) = ', Lr)

#######################################################################################################################
# Writing L(r) and K(r) to a csv file

with open(savelocation + filename + lrdata + '_2.' + timestep + '.csv', 'w', newline='') as f:
    writer=csv.writer(f)
    writer.writerows([Lr, Kr])

#######################################################################################################################
# Plotting L(r) vs radius

fig=plt.figure(figsize=(9,9))
plt.plot(radius, Lr, 'r*-')
plt.plot(radius, radius, 'k-')
plt.title("Ripley's L value", fontsize = 18 )
plt.xlabel('radius (cm)', fontsize = 18)
plt.ylabel('L(r)', fontsize = 18)
plt.savefig(figurelocation + filename + '_lr_2.' + timestep + '.svg', format='svg')

print ('Look for the L(r) plot.')

#######################################################################################################################
# Calulating Ripley's g value

gr = []

for i in range(1, len(radius)):
    dkr = (Kr[i] - Kr[i-1]) / (radius[i] - radius[i-1])
    #print (dkr)
    gr.append(dkr / (2 * np.pi * radius[i]))
gr = np.array(gr)

print ('g(r) = ', gr)
#######################################################################################################################
# writing g(r) to a csv file

with open(savelocation + filename + grdata + '_2.' + timestep + '.csv', 'w', newline='') as f:
    writer=csv.writer(f)
    writer.writerow(gr)

#######################################################################################################################
# Plotting g(r) vs radius

fig=plt.figure(figsize=(9,9))
plt.plot(radius[1:len(radius)], gr, 'bo-')
plt.plot(radius[1:len(radius)], np.ones(len(gr)), 'k-')
plt.title("Ripley's g value (pair correlation)", fontsize = 18)
plt.xlabel('radius (cm)', fontsize = 18)
plt.ylabel('g(r)', fontsize = 18)
plt.savefig(figurelocation + filename + '_gr_2.' + timestep + '.svg', format='svg')

print ('Look for the g(r) plot.')
