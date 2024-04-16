import numpy as np
import math
import matplotlib.pyplot as plt
from gc import collect
from scipy.optimize import curve_fit

path = "/home/sq-wm/OpenFOAM/sq-wm-v2306/run/calc14/"
# read size of new structured array
with open(path + 'grid_size.txt', 'r') as fp:
    count_of_cells_axis = fp.read().split('\n')
    x_cells = int(count_of_cells_axis[0])
    y_cells = int(count_of_cells_axis[1])
    z_cells = int(count_of_cells_axis[2])
    print('structured grid has ' + str(x_cells) + ' * ' + str(y_cells) + ' * ' + str(z_cells) + ' cells')
cells_number_axis = min([x_cells, y_cells, z_cells])

left_boundary_x =   -1.0*math.pi
right_boundary_x =  math.pi
down_boundary_y =   -1.0*math.pi
up_boundary_y =     math.pi
rear_boundary_z =   -1.0*math.pi
front_boundary_z =  math.pi
nu = 0.01
OMEGA = 50.0

# return: 
# x, y, z - arrays with coordinates of nodes - not returned!
# u, v, w - projection of velocity vector in this nodes
def get_data_from_file(full_path_to_files, size_x, size_y, size_z):
    file = open(full_path_to_files + 'U_x.txt', 'r')
    data = file.read().split('\n')
    slice = []
    u = []
    for ind_z in range(size_z):
        for ind_y in range(size_y):
            slice.append(data[ind_z*(size_y+1) + ind_y].split(' ')[ : size_x])
        u.append(slice)
        slice = []
    file.close()

    file = open(full_path_to_files + 'U_y.txt', 'r')
    data = file.read().split('\n')
    slice = []
    v = []
    for ind_z in range(size_z):
        for ind_y in range(size_y):
            slice.append(data[ind_z*(size_y+1) + ind_y].split(' ')[ : size_x])
        v.append(slice)
        slice = []
    file.close()

    file = open(full_path_to_files + 'U_z.txt', 'r')
    data = file.read().split('\n')
    slice = []
    w = []
    for ind_z in range(size_z):
        for ind_y in range(size_y):
            slice.append(data[ind_z*(size_y+1) + ind_y].split(' ')[ : size_x])
        w.append(slice)
        slice = []
    file.close()

    del slice
    collect()
    print("end reading velocity from file " + full_path_to_files)
    u = np.array(u, dtype=float)
    v = np.array(v, dtype=float)
    w = np.array(w, dtype=float)
    return u, v, w

def get_vorticity_z_from_file(full_path_to_files, size_x, size_y, size_z):
    file = open(full_path_to_files + 'vort_z.txt', 'r')
    data = file.read().split('\n')
    slice = []
    vort = []
    for ind_z in range(size_z):
        for ind_y in range(size_y):
            slice.append(data[ind_z*(size_y+1) + ind_y].split(' ')[ : size_x])
        vort.append(slice)
        slice = []
    file.close()

    del slice
    collect()
    print("end reading vorticity from file " + full_path_to_files)
    vort = np.array(vort, dtype=float)
    return vort

def get_mean_data_from_file(full_path_to_files, size_x, size_y, size_z):
    file = open(full_path_to_files + 'Ux_mean.txt', 'r')
    data = file.read().split('\n')
    slice = []
    u = []
    for ind_z in range(size_z):
        for ind_y in range(size_y):
            slice.append(data[ind_z*(size_y+1) + ind_y].split(' ')[ : size_x])
        u.append(slice)
        slice = []
    file.close()

    file = open(full_path_to_files + 'Uy_mean.txt', 'r')
    data = file.read().split('\n')
    slice = []
    v = []
    for ind_z in range(size_z):
        for ind_y in range(size_y):
            slice.append(data[ind_z*(size_y+1) + ind_y].split(' ')[ : size_x])
        v.append(slice)
        slice = []
    file.close()

    file = open(full_path_to_files + 'Uz_mean.txt', 'r')
    data = file.read().split('\n')
    slice = []
    w = []
    for ind_z in range(size_z):
        for ind_y in range(size_y):
            slice.append(data[ind_z*(size_y+1) + ind_y].split(' ')[ : size_x])
        w.append(slice)
        slice = []
    file.close()

    del slice
    collect()
    print("end reading velocity from file " + full_path_to_files)
    u = np.array(u, dtype=float)
    v = np.array(v, dtype=float)
    w = np.array(w, dtype=float)
    return u, v, w

def get_mean_vorticity_from_file(full_path_to_files, size_x, size_y, size_z):
    file = open(full_path_to_files + 'vort_z_mean.txt', 'r')
    data = file.read().split('\n')
    slice = []
    vort = []
    for ind_z in range(size_z):
        for ind_y in range(size_y):
            slice.append(data[ind_z*(size_y+1) + ind_y].split(' ')[ : size_x])
        vort.append(slice)
        slice = []
    file.close()

    del slice
    collect()
    print("end reading velocity from file " + full_path_to_files)
    vort = np.array(vort, dtype=float)
    return vort

# return: radius and angle (-pi, pi) from (x0, y0) to (x, y)
def get_polar_coordinate(x, x0, y, y0):
    r = math.sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0))
    phi = math.atan2((y - y0), (x - x0))
    return r, phi

u_mean, v_mean, w_mean = get_mean_data_from_file(path, x_cells, y_cells, z_cells)
vort_z_mean = get_mean_vorticity_from_file(path, x_cells, y_cells, z_cells)

u_slice = np.zeros((y_cells, x_cells))
v_slice = np.zeros((y_cells, x_cells))
w_slice = np.zeros((y_cells, x_cells))
vort_slice = np.zeros((y_cells, x_cells))

for z_ind in range(z_cells):
    u_slice += u_mean[z_ind, : , :]
    v_slice += v_mean[z_ind, : , :]
    w_slice += w_mean[z_ind, : , :]
    vort_slice += vort_z_mean[z_ind, : , :]
vort_slice[0 : 3, : ] = 0.0
vort_slice[y_cells-3 : , : ] = 0.0
vort_slice[ : , 0 : 3] = 0.0
vort_slice[ : , x_cells-3 : ] = 0.0

u_slice /= z_cells
v_slice /= z_cells
w_slice /= z_cells
vort_slice /= z_cells
#mod_vel_slice = np.sqrt(u_slice*u_slice + v_slice*v_slice + w_slice*w_slice)

min_vort = np.min(vort_slice)
ind_y_min = np.where(vort_slice == min_vort)[0][0]
ind_x_min = np.where(vort_slice == min_vort)[1][0]
print('min vorticity = ' + str(vort_slice[ind_y_min, ind_x_min]))
print('in point (' + str(ind_y_min) + ', ' + str(ind_x_min) + ')')

#r1, _ = get_polar_coordinate(0, ind_x_min, 0, ind_y_min)
#r2, _ = get_polar_coordinate(x_cells, ind_x_min, 0, ind_y_min)
#r3, _ = get_polar_coordinate(0, ind_x_min, y_cells, ind_y_min)
#r4, _ = get_polar_coordinate(x_cells, ind_x_min, y_cells, ind_y_min)
max_radius = min([ind_x_min, x_cells - ind_x_min, ind_y_min, y_cells - ind_y_min])
print(max_radius)
count_of_radius_point = 100
radius_array_vel_projection = np.zeros(count_of_radius_point)
radius_array_counter = np.zeros(count_of_radius_point)

for ind_y in range(y_cells):
    for ind_x in range(x_cells):
        r, phi = get_polar_coordinate(ind_x, ind_x_min, ind_y, ind_y_min)
        if r <= max_radius:
           azimuthal_vel = v_slice[ind_y, ind_x]*math.cos(phi) - u_slice[ind_y, ind_x]*math.sin(phi)
           ind = int((count_of_radius_point - 1) * r / (max_radius - 0.0))
           radius_array_vel_projection[ind] += azimuthal_vel
           radius_array_counter[ind] += 1

#print(radius_array_vel_projection)
radius_array_vel_projection /= radius_array_counter
#print(radius_array_vel_projection)

def optimisation_function(x, A, B):
    return A * x * np.log(x / B)

#fig, axs = plt.subplots(1, 2)
#axs.set_title("mean z vorticity")
#plt.imshow(vort_slice, cmap='seismic')
#axs[0].imshow(u_slice)
#axs[1].imshow(v_slice)
#plt.imshow(mod_vel_slice, cmap='seismic')
#plt.show()

#print(np.linspace(0.0, max_radius, count_of_radius_point))
x_gr = np.linspace(0.0, max_radius, count_of_radius_point)[1:int(count_of_radius_point)]
params, _ = curve_fit(optimisation_function, x_gr, radius_array_vel_projection[1:int(count_of_radius_point)], bounds=([0.1, 80.0], [2.0, 120.0]))
print(params)

fig, axs = plt.subplots()
axs.plot(np.linspace(0.0, max_radius, count_of_radius_point), radius_array_vel_projection, label='results')
axs.plot(x_gr, optimisation_function(x_gr, *params), linestyle='dashed', label='approximation')
axs.grid(True, linestyle='--')
axs.set_ylabel('mean azimuthal velocity')
axs.set_xlabel('radius')
plt.xlim(0.0, max_radius*1.1)
plt.legend()
plt.show()
