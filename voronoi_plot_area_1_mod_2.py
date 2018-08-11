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
# last modified:
# 02/22/2018
#
# Author:
# Zack McIntire
#
# Imputs: (CSV)
# - particle location file --Note: for timesteps greater than 0 the saved 
#   particle file from the clustering code should be used. (low_vol_frac) 
# - fluid file
# - normalized area for timestep zero of non mixing bowl simulation
# 
# Outputs:
#
# --arrays-- (CSV)
# - An array that contains the Voronoi cell area, and the normalized area.
# - An array of the aspect ratio data.
#
# --figures-- (SVG)
# - Voronoi diagram of the settling particles
# - Histogram of the normalized area compared to the initial non mixing bowl
#   and the eq from Ferenc and Neda
# - Voronoi diagram colored by aspect ratio and the color bar (seperate figure)
# - Histogram of the aspect ratio compared to the first timestep of the non
#   mixing bowl run.
# - Voronoi diagram colored to normalized area with colorbar (seperate figure)
# - rose diagram of long eigen vector direction.
# 
###############################################################################
# importing the necessary python modules

import csv
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
#import matplotlib.mlab as mlab
#from scipy.stats import norm
import matplotlib as mpl
import matplotlib.cm as cm
#import mplstereonet
from scipy.spatial import Voronoi, voronoi_plot_2d

###############################################################################
# Variables and file information

csvarray = []
data_array = []
filelocation2 = '/wd2/csv_data_files/'
filelocation = '/home/zack/Documents/csv_data_files/'
figurelocation = '/home/zack/Documents/clustering_figures/'
savelocation = '/home/zack/Documents/csv_data_files/'
#randfilename = 'box512_ht_loc'
#filename = 'test_case_20'
filename = 'box512_ht_loc_low_vol_frac'#_low_vol_frac'
resfile = 'box512_ht_res'
timestep = '400'

###############################################################################
# Opening the necessary files

#with open(filelocation + filename + '.csv', 'r') as csvfile:
#    csv_data = list(csv.reader(csvfile, delimiter=","))

with open(filelocation + filename + '.' + timestep + '.csv', 'r') as csvfile:
    csv_data = list(csv.reader(csvfile, delimiter=","))

csvarray = np.array(csv_data[1:], dtype=np.float) # Remove headers and convert to floats

data_array = np.array(csvarray)

with open(filelocation2 + resfile + '.' + timestep + '.csv', 'r') as csvres:
    res_data = list(csv.reader(csvres, delimiter=","))
resarray = np.array(res_data[1:], dtype=np.float)

###############################################################################
# Finding the Voronoi points and plotting the diagram

# defining the bounding region:
x_min = 0
x_max = 512
y_min = min(data_array[:, 8])
y_max = max(data_array[:, 8])

# finding the total area of space opccupied by the settling particles.
res_array = []
for volume_index in range(0, len(resarray)):
    if 1 - resarray[volume_index, 0] < 0.46 and 1 - resarray[volume_index, 0] \
    > 0:
        res_array.append(resarray[volume_index, :])
        
res_array = np.array(res_array)

# finding the height to make the random domain
phis = 1 - np.mean(res_array[:, 0]) # solid volume fraction
number_of_particles = len(data_array)
particle_volume = 4 / 3 * np.pi *0.2**3 * number_of_particles # volume of particles
total_volume = particle_volume / phis # volume of the domain
total_area = total_volume / 0.4 # divide by the depth of the domain (particle radius)
domain_height = total_area / 512 # The length of the domain should not change
 
points = []
for i in range(0, len(data_array)):
    points.append([data_array[i, 7], data_array[i, 8]])
points = np.array(points)

vor = Voronoi(points, qhull_options='Qbb Qc Qx')
voronoi_plot_2d(vor, show_vertices = False)
plt.savefig(figurelocation + filename + '_Voronoi_regions_2.' + timestep + \
            '.svg', format = 'svg')

vor_vert = np.array(vor.vertices)

position = np.array(-1)

###############################################################################      
def vert_remover(posindex):
    
    vol_frac = 0.46
    box_dist = 13

    if vor_vert[posindex, 0] < x_max and vor_vert[posindex, 0] > x_min \
    and vor_vert[posindex, 1] < y_max and vor_vert[posindex, 1] > y_min:

        volarray = []
        for vfindex in range(0, len(resarray)):
            # creating a smaller fluid array based on how close the corner is 
            # to the particle.
            if resarray[vfindex, 6] < vor_vert[posindex, 0] + box_dist and \
            resarray[vfindex, 6] > vor_vert[posindex, 0] - box_dist and \
            resarray[vfindex, 7] < vor_vert[posindex, 1] + box_dist and \
            resarray[vfindex, 7] > vor_vert[posindex, 1] - box_dist:
                volarray.append(resarray[vfindex])
        volarray = np.array(volarray)  
        #print(volarray[:3])
        
        # Making sure the 8 closest fluid cell corners have a volume fraction 
        # greater than zero.
        distarray = []
        distarray_2 = []
        distarray_3 = []
        distarray_4 = []
        distarray_5 = []
        distarray_6 = []
        distarray_7 = []
        distarray_8 = []
        ref_dist = 2
        holder = ref_dist
        ref_dist_2 = 4
        holder2 = ref_dist_2
        ref_dist_3 = 6
        holder_3 = ref_dist_3
        ref_dist_4 = 8
        holder_4 = ref_dist_4
        ref_dist_5 = 10
        holder_5 = ref_dist_5
        ref_dist_6 = 11
        holder_6 = ref_dist_6
        ref_dist_7 = 12
        holder_7 = ref_dist_7
        ref_dist_8 = 13
        for volindex in range(0, len(volarray)):
            dist_to_corner = np.sqrt((vor_vert[posindex, 0] - \
                                      volarray[volindex, 6])**2 + \
                                      (vor_vert[posindex, 1] - \
                                      volarray[volindex, 7])**2)
            if dist_to_corner < ref_dist:
                holder = ref_dist
                ref_dist = dist_to_corner
                distarray = volarray[volindex]
            elif dist_to_corner < ref_dist_2 and dist_to_corner > holder:
                holder2 = ref_dist_2
                ref_dist_2 = dist_to_corner
                distarray_2 = volarray[volindex]
            elif dist_to_corner < ref_dist_3 and dist_to_corner > holder2:
                holder_3 = ref_dist_3
                ref_dist_3 = dist_to_corner
                distarray_3 = volarray[volindex]
            elif dist_to_corner < ref_dist_4 and dist_to_corner > holder_3:
                holder_4 = ref_dist_4
                ref_dist_4 = dist_to_corner
                distarray_4 = volarray[volindex]
            elif dist_to_corner < ref_dist_5 and dist_to_corner > holder_4:
                holder_5 = ref_dist_5
                ref_dist_5 = dist_to_corner
                distarray_5 = volarray[volindex]
            elif dist_to_corner < ref_dist_6 and dist_to_corner > holder_5:
                holder_6 = ref_dist_6
                ref_dist_6 = dist_to_corner
                distarray_6 = volarray[volindex]
            elif dist_to_corner < ref_dist_7 and dist_to_corner > holder_6:
                holder_7 = ref_dist_7
                ref_dist_7 = dist_to_corner
                distarray_7 = volarray[volindex]
            elif dist_to_corner < ref_dist_8 and dist_to_corner > holder_7:
                ref_dist_8 = dist_to_corner
                distarray_8 = volarray[volindex]
        distarray = np.array(distarray)
        distarray_2 = np.array(distarray_2)
        distarray_3 = np.array(distarray_3)
        distarray_4 = np.array(distarray_4)
        distarray_5 = np.array(distarray_5)
        distarray_6 = np.array(distarray_6)
        distarray_7 = np.array(distarray_7)
        distarray_8 = np.array(distarray_8)
        if 1-distarray[0] >= vol_frac or (1-distarray[0] == 0 and \
                      1-distarray_2[0] == 0 and 1-distarray_3[0] == 0 \
                      and 1-distarray_4[0] == 0 and 1-distarray_5[0] == 0 \
                      and 1-distarray_6[0] == 0 and 1-distarray_7[0] == 0 \
                      and 1-distarray_8[0] == 0):
            return [posindex, 0]
        else:
            return [posindex, 1]
    else:
        return [posindex, 0]
    
###############################################################################            
print ('finding the position of outlying vertices')
vert_len = len(vor_vert)
pool = mp.Pool(processes=3)
result = pool.map(vert_remover, range(0, vert_len)) 
result = np.array(result)
pool.close()
pool.join()
modresult = []
for vert_index in range(0, len(result)):
    if result[vert_index, 1] == 0:
        position = np.append(position, vert_index)

position = np.array(position) 

print ('outlying vertices found')

vor_reg = vor.regions

###############################################################################
# making sure the Voronoi vertices are within the domain boundaries and within
# the specified volume fraction range.

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

print ("Finding Voronoi cell area.")
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
    
eigvalarray = np.array(eigvalarray)
eigvecarray = np.array(eigvecarray)
    
data_mean_area = np.mean(area)
    
###############################################################################
# creating a random array of particles to find the average area of Voronoi
# polygons to find the ratio.

num_of_part = len(csvarray)
#y_mid = np.mean(csvarray[:, 8])
#y_top = max(csvarray[:, 8])
#y_height = (y_top - y_mid) * 2

rand_array = np.random.random((num_of_part, 2))
rand_array[:, 0] = rand_array[:, 0] * 512
rand_array[:, 1] = rand_array[:, 1] * domain_height
y_bottom = min(rand_array[:, 1])

rand_vor = Voronoi(rand_array, qhull_options='Qbb Qc Qx')

rand_vor_vert = rand_vor.vertices

rand_position = [-1]

for indexr in range(0, len(rand_vor_vert)):
    if rand_vor_vert[indexr, 0] > x_max or rand_vor_vert[indexr, 0] < x_min \
    or rand_vor_vert[indexr, 1] > domain_height or rand_vor_vert[indexr, 1] < \
    y_bottom:
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
for i in range(0, len(rand_reg)):
    rand_corners = []
    vor_reg_temp = np.array(rand_reg[i])
    for n in vor_reg_temp:
        rand_corners.append(rand_vor_vert[n])
    rand_corners = np.array(rand_corners)
    xy = rand_corners.T
    rand_corners_sorted = PolygonSort(rand_corners)
    rand_area.append(PolygonArea(rand_corners_sorted))
 
rand_mean_area = np.mean(rand_area)

print ("Voronoi cell area found")

norm_area = []
for i in range(0, len(area)):
    norm_area.append(area[i]/rand_mean_area)
    
###############################################################################
# saving Voronoi area data
with open(savelocation + filename + '_Voronoi_area_2.' + timestep + \
          '.csv', 'w', newline='') as f_area:
    writer=csv.writer(f_area)
    writer.writerows([area, norm_area])
    
###############################################################################
# Reading in the area of the random just settling simulation

with open(savelocation + 'box512_ht_zero_loc_Voronoi_area_2.0.csv', \
          'r', newline='') as areadata:
    area_data = list(csv.reader(areadata, delimiter=","))
norm_area_c = np.array(area_data[1], dtype=np.float)

###############################################################################    
# plotting the normalized cell area as a pdf plot
    
fig = plt.figure()

bins = np.linspace(0, 4, 100)

n_area_c, bins_area_c, patches_area_c = plt.hist(norm_area_c, bins = bins, \
                                                 normed = 1, histtype = \
                                                 'step', label = \
                                                 'compared area histogram')

n_area, bins_area, patches_area = plt.hist(norm_area, bins = bins, normed = 1,\
                                           histtype = 'step', label = \
                                           'area histogram')

area_x = [0]
area_x_index = 100
while area_x_index > 0:
    area_x.append(4 / area_x_index)
    area_x_index = area_x_index - 1
    
f_area = []
i_area = 0
while i_area <= len(area_x) - 1:
    y_area = area_x[i_area]
    f_area.append(343/15*np.sqrt(7/(2*np.pi))*y_area**(5/2)*np.exp(-\
                  (7/2)*y_area))
    i_area = i_area+1
    
plt.plot(area_x, f_area, 'k-.', label = 'Ferenc and Neda')
plt.legend(loc = 'best')
plt.title('pdf of normalized cell area')
plt.xlabel('normalized Voronoi cell area')
plt.ylabel('P.D.F')
plt.xlim(xmin = 0, xmax = 4)
plt.yscale('log', nonposy='clip')

plt.savefig(figurelocation + filename + '_area_2_pdf.' + timestep + \
            '.svg', format = 'svg')   

print ("Finding aspect ratio.")

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
for i in range(0, len(a1)):
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
    
#ar = np.array(ar)

print ("Aspect ratio found.")

###############################################################################
# saving Voronoi aspect ratio data
with open(savelocation + filename + '_Voronoi_aspect_ratio_2.' + timestep + \
          '.csv', 'w', newline='') as f_ar:
    writer=csv.writer(f_ar)
    writer.writerow(ar)

###############################################################################
# Coloring the Voronoi polygons by aspect ratio

minima = min(ar)
maxima = max(ar)
cmap_ar = cm.autumn

norms = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
mapper = cm.ScalarMappable(norm=norms, cmap=cmap_ar)

ar_vor = voronoi_plot_2d(vor, show_points=True, show_vertices=False, s=1)
for r in range(0, len(reg)):
    reg_seg = reg[r]
    polygon = [vor_vert[i] for i in reg_seg]
    plt.fill(*zip(*polygon), color=mapper.to_rgba(ar[r]))
plt.title('Voronoi regions colored by aspect ratio')
plt.savefig(figurelocation + filename + '_Voronoi_regions_ar_2.' + timestep + \
            '.svg', format = 'svg')


fig_ar = plt.figure(figsize=(8, 3))
ax1_ar = fig_ar.add_axes([0.05, 0.80, 0.9, 0.15])
cb1_ar = mpl.colorbar.ColorbarBase(ax1_ar, cmap=cmap_ar, norm=norms, orientation \
                                = 'horizontal')
cb1_ar.set_label('Voronoi region aspect ratio')
plt.savefig(figurelocation + filename + '_Voronoi_regions_ar_colorbar_2.' + \
            timestep + '.svg', format = 'svg')

###############################################################################
# Reading in the aspect ratio of the random just settling simulation

with open(savelocation + 'box512_ht_zero_loc_Voronoi_aspect_ratio_2.0.csv', \
          'r', newline='') as ardata:
    ar_data = list(csv.reader(ardata, delimiter=","))
ar_c = np.array(ar_data[0],  dtype=np.float)

###############################################################################
# plotting the pdf of the aspect ratio.

# Starting with the plot of just settling time zero.
fig = plt.figure()

n_ar_c, bins_ar_c, patches_ar_c = plt.hist(ar_c, bins = 100, normed = 1, \
                                     histtype = 'step', \
                                     label = 'Compared aspect ratio')

# Adding in the plot of the current runs aspect ratio
n_ar, bins_ar, patches_ar = plt.hist(ar, bins = 100, normed = 1, \
                                     histtype = 'step', \
                                     label = 'aspect ratio histogram')

plt.legend(loc = 'best')
plt.title('histogram of aspect ratio')
plt.xlabel('Voronoi cell aspect ratio')
plt.ylabel('P.D.F')
plt.xlim(xmin = 0.5, xmax = max(ar) + 0.5)
plt.yscale('log', nonposy='clip')

plt.savefig(figurelocation + filename + '_ar_hist_2.' + timestep + '.svg', \
            format = 'svg')

print ('Finding Voronoi cell orientation.')

###############################################################################
# Coloring the Voronoi polygons by normalized area

minima_area = min(norm_area)
maxima_area = max(norm_area)
cmap_area = cm.autumn

norms_area = mpl.colors.Normalize(vmin=minima_area, vmax=maxima_area, clip=True)
mapper_area = cm.ScalarMappable(norm=norms_area, cmap=cmap_area)


vor_area = voronoi_plot_2d(vor, show_points=True, show_vertices=False, s=0.5)
for r in range(0, len(reg)):
    reg_seg_area = reg[r]
    polygon_area = [vor_vert[i] for i in reg_seg_area]
    plt.fill(*zip(*polygon_area), color=mapper_area.to_rgba(norm_area[r]))
plt.title('Voronoi regions colored by normalized area')
plt.savefig(figurelocation + filename + '_Voronoi_regions_area_2.' + timestep + \
            '.svg', format = 'svg')

fig_area = plt.figure(figsize=(8, 3))
ax1_area = fig_area.add_axes([0.05, 0.80, 0.9, 0.15])
cb1_area = mpl.colorbar.ColorbarBase(ax1_area, cmap=cmap_area, norm=norms_area, \
                                     orientation = 'horizontal')
cb1_area.set_label('Voronoi region normalized area')
plt.savefig(figurelocation + filename + '_Voronoi_regions_area_colorbar_2.' + timestep + \
            '.svg', format = 'svg')
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
# values below the rose diagram threshold by 2pi
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

print ('Voronoi cell orinetation found')

# creating the bins and bin values
weight_rose = np.ones_like(new_theta)/float(len(new_theta))

bin_edge = np.deg2rad(np.arange(-5, 360, 10))
number_of_theta, bin_edge = np.histogram(new_theta, bin_edge, weights=weight_rose)
number_of_theta = np.array(number_of_theta)
less_than_mean = np.zeros_like(number_of_theta)

mean_theta = np.mean(number_of_theta)

for thetaindex in range(0, len(number_of_theta)):
    if number_of_theta[thetaindex] <= mean_theta:
        less_than_mean[thetaindex] = number_of_theta[thetaindex]
        number_of_theta[thetaindex] = 0

# plotting the rose diagram
fig = plt.figure()

ax = plt.subplot(111, projection='polar')

above_mean = ax.bar(np.deg2rad(np.arange(0, 360, 10)), number_of_theta/2, width=np.deg2rad(10),\
       bottom=0.0, color='g', edgecolor='k')
below_mean = ax.bar(np.deg2rad(np.arange(0, 360, 10)), less_than_mean/2, width=np.deg2rad(10),\
       bottom=0.0, color='b', edgecolor='k')
ax.legend((above_mean, below_mean), ('above mean', 'below mean'))
ax.set_theta_zero_location('E')
ax.set_theta_direction(1)
ax.set_title('Rose Diagram of the "polygon orientation"', y=1.10, fontsize=15)
plt.savefig(figurelocation + filename + '_ar_rose_2.' + timestep + '.svg', \
            format = 'svg')

