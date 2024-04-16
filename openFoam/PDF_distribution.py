import numpy as np
import math
import matplotlib.pyplot as plt
from gc import collect

path = "/home/sq-wm/OpenFOAM/sq-wm-v2306/run/calc12/"
# read size of new structured array
with open(path + 'grid_size.txt', 'r') as fp:
    count_of_cells_axis = fp.read().split('\n')
    x_cells = int(count_of_cells_axis[0])
    y_cells = int(count_of_cells_axis[1])
    z_cells = int(count_of_cells_axis[2])
    print('structured grid has ' + str(x_cells) + ' * ' + str(y_cells) + ' * ' + str(z_cells) + ' cells')
cells_number_axis = min([x_cells, y_cells, z_cells])

start_time = 15
end_time = 20
time_step = 0.5
time_array = np.round(np.linspace(start_time, end_time, num=int((end_time - start_time)/time_step + 1), endpoint=True), decimals=1)
print('time points being processed:')
print(time_array)

left_boundary_x =   -1.0*math.pi
right_boundary_x =  math.pi
down_boundary_y =   -1.0*math.pi
up_boundary_y =     math.pi
rear_boundary_z =   -1.0*math.pi
front_boundary_z =  math.pi
nu = 0.01
OMEGA = 50.0
Lf = 1.0    # fitting factor

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
    print("end reading mean vorticity from file " + full_path_to_files)
    vort = np.array(vort, dtype=float)
    return vort

def write_3d_field(field, nx, ny, nz, name_with_path):
    file = open(name_with_path, "w+")
    for ind_z in range(nz):
        for ind_y in range(ny):
            for ind_x in range(nx):
                file.write(str(field[ind_z, ind_y, ind_x]) + ' ')
            file.write('\n')
        file.write('\n')
    file.close()

# return: radius and angle (-pi, pi) from (x0, y0) to (x, y)
def get_polar_coordinate(x, x0, y, y0):
    r = math.sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0))
    phi = math.atan2((y - y0), (x - x0))
    return r, phi

# U_RMS = sqrt( sum( ux^2 + uy^2 + uz^2 ) / N )
def find_RMS_velocity(u_x, u_y, u_z, N):
    return  math.sqrt(np.sum(np.square(u_x) + np.square(u_y) + np.square(u_z)) / N)

# Re makro = RMS U * Lf / nu
def find_makro_Reynolds_number(U, Lf, nu):
    return U * Lf / nu

# Ro makro = RMS U / (2*OMEGA*Lf) 
def find_makro_Rossby_number(U, Lf, OMEGA):
    return U / (2.0 * OMEGA * Lf)

def get_z_mean_slice(field, z_cells, y_cells, x_cells):
    new_field = np.zeros((y_cells, x_cells))
    for z_coord in range(z_cells):
        new_field += field[z_coord, : , : ]
    return new_field / z_cells

def find_distance_to_center(vort_z_slice, y_cells, x_cells):
    y_geom = y_cells / 2 - 0.5
    x_geom = x_cells / 2 - 0.5
    min_vort = np.min(vort_z_slice)
    y_real = np.where(vort_z_slice == min_vort)[0][0]
    x_real = np.where(vort_z_slice == min_vort)[1][0]
    r, _ = get_polar_coordinate(x_geom, x_real, y_geom, y_real)
    return r

def find_scalar_distribution(scalar_field, z_cells, y_cells, x_cells, left_bound, right_bound, steps_count):
    distribution = np.zeros(steps_count)
    dx = (right_bound - left_bound) / (steps_count - 1)
    for ind_z in range(z_cells):
        for ind_y in range(y_cells):
            for ind_x in range(x_cells):
                value = scalar_field[ind_z, ind_y, ind_x]
                if value <= left_bound:
                    ind = 0
                elif value >= right_bound:
                    ind = steps_count - 1
                else:
                    ind = int((value - left_bound) // dx)
                distribution[ind] += 1
    distribution[0] = 0.0
    distribution[steps_count - 1] = 0.0
    return distribution / (z_cells * y_cells * x_cells)

def find_rang_distribution(scalar_field, z_cells, y_cells, x_cells, left_bound, right_bound, steps_count):
    distribution_plus = find_scalar_distribution(scalar_field, z_cells, y_cells, x_cells, 0.0, right_bound, steps_count)
    distribution_minus = find_scalar_distribution(scalar_field, z_cells, y_cells, x_cells, left_bound, 0.0, steps_count)
    distribution_plus = np.sort(distribution_plus, kind='mergesort')
    distribution_minus = np.sort(distribution_minus, kind='mergesort')
    return distribution_plus, distribution_minus

def find_rang_in_array(scalar_field):
    arr_sorted = np.sort(scalar_field, kind='mergesort', axis=None)
    arr_minus = arr_sorted[0 : len(arr_sorted) // 2]
    for i in range(len(arr_minus)-1, 0, -1):
        if arr_minus[i] > 0.0:
            arr_minus[i] = 0.0
        elif arr_minus[i] < 0.0:
            i = 0
    arr_minus = np.abs(np.flip(arr_minus))
    arr_plus = arr_sorted[len(arr_sorted) // 2 : len(arr_sorted)]
    for i in range(len(arr_plus)):
        if arr_plus[i] < 0.0:
            arr_plus[i] = 0.0
        elif arr_plus[i] > 0.0:
            i = len(arr_plus)

    del arr_sorted
    collect()
    return arr_minus, arr_plus

dx = (right_boundary_x - left_boundary_x) / x_cells
dy = (up_boundary_y - down_boundary_y) / y_cells
dz = (front_boundary_z - rear_boundary_z) / z_cells

count_of_vort_steps = 250
left_boundary_vort = -150.0
right_boundary_vort = 150.0
#vort_distribution_from_time = []
vort_distribution_from_time_plus = []
vort_distribution_from_time_minus = []
vort_plus_time = []
vort_minus_time = []
needed_count = int(z_cells * y_cells * x_cells // 2)
#vort_z_mean = np.zeros((z_cells, y_cells, x_cells))
#d_vort = (right_boundary_vort - left_boundary_vort) / (count_of_vort_steps - 1)
#vort_scale = np.linspace(left_boundary_vort, right_boundary_vort, count_of_vort_steps) / (2.0 * OMEGA)
#vort_line = np.zeros(z_cells * y_cells * x_cells)

for time in time_array:
    if time % 1 == 0:
        time = int(time)
    print('time = ' + str(time))
    #u_center, v_center, w_center = get_data_from_file(path + str(time) + '/', x_cells, y_cells, z_cells)
    #vort_z = get_vorticity_z_from_file(path + str(time) + '/', x_cells, y_cells, z_cells)
    #vort_z_mean += vort_z
    
    vort_z = get_vorticity_z_from_file(path + str(time) + '/', x_cells, y_cells, z_cells)
    #vort_z = get_mean_vorticity_from_file(path, x_cells, y_cells, z_cells)
    distr_plus, distr_minus = find_rang_distribution(vort_z, z_cells, y_cells, x_cells, left_boundary_vort, right_boundary_vort, count_of_vort_steps)
    vort_distribution_from_time_plus.append(distr_plus)
    vort_distribution_from_time_minus.append(distr_minus)
    vort_minus, vort_plus = find_rang_in_array(vort_z)
    vort_minus_time.append(vort_minus)
    vort_plus_time.append(vort_plus)

    del distr_plus
    del distr_minus
    del vort_z
    del vort_minus
    del vort_plus
    collect()

#vort_z_mean /= len(time_array)
#write_3d_field(vort_z_mean, x_cells, y_cells, z_cells, path + 'vort_z_mean.txt')
#vort_distribution_from_time = np.array(vort_distribution_from_time)
vort_distribution_from_time_plus = np.array(vort_distribution_from_time_plus)
vort_distribution_from_time_minus = np.array(vort_distribution_from_time_minus)
vort_plus_time = np.array(vort_plus_time)
vort_minus_time = np.array(vort_minus_time)
vort_plus_time = np.mean(vort_plus_time, axis=0)
vort_minus_time = np.mean(vort_minus_time, axis=0)

#vort_z_mean = get_mean_vorticity_from_file(path, x_cells, y_cells, z_cells)
#vort_mean_distribution = find_scalar_distribution(vort_z_mean, z_cells, y_cells, x_cells, left_boundary_vort, right_boundary_vort, count_of_vort_steps)

#fig, axs = plt.subplots()
#axs.set_title("mean z vorticity")
#plt.imshow(vort_slice, cmap='seismic')
#for i in range(len(time_array)):
#    plt.plot(vort_scale, vort_distribution_from_time[i], label=str(time_array[i]))
#axs[0].set_title('vorticity')
#axs[1].plot(vort_scale, np.gradient(vort_distribution))
#axs[1].set_title('derivative(vorticity)')
#plt.imshow(mod_vel_slice, cmap='seismic')
#plt.grid(True, linestyle='--')
#plt.legend(loc='upper right')
#plt.show()

#vort_distribution_mean = np.mean(vort_distribution_from_time, axis=0)
#fig, axs = plt.subplots(1, 2)
#axs[0].set_title('meaning after calculate')
#axs[0].plot(vort_scale, vort_distribution_mean)
#axs[0].grid(True, linestyle='--')
#axs[1].set_title('calculating after mean')
#axs[1].plot(vort_scale, vort_mean_distribution)
#axs[1].grid(True, linestyle='--')
#plt.show()

vort_distribution_plus_mean = np.mean(vort_distribution_from_time_plus, axis=0)
vort_distribution_minus_mean = np.mean(vort_distribution_from_time_minus, axis=0)
fig, axs = plt.subplots()
axs.plot(np.arange(count_of_vort_steps), vort_distribution_minus_mean, color='blue', label='vort < 0', linewidth=3)
axs.plot(np.arange(count_of_vort_steps), vort_distribution_plus_mean, color='red', label='vort > 0', linewidth=3)
axs.grid(True, linestyle='--')
axs.set_xlabel('rang')
axs.set_ylabel('frequency')
plt.legend()
plt.show()

fig, axs = plt.subplots()
axs.plot(np.arange(needed_count), vort_minus_time, color='blue', label='vort < 0', linewidth=3)
axs.plot(np.arange(needed_count), vort_plus_time, color='red', label='vort > 0', linewidth=3)
axs.grid(True, linestyle='--')
axs.set_xlabel('rang')
axs.set_ylabel('|vort|')
plt.legend()
plt.show()

#fig, axs = plt.subplots(1, 2)
#axs[0].plot(np.arange(count_of_vort_steps), np.sort(vort_distribution_mean))
#axs[0].grid(True, linestyle='--')
#axs[1].plot(np.arange(count_of_vort_steps), np.sort(vort_mean_distribution))
#axs[1].grid(True, linestyle='--')
#plt.show()

#fig, axs = plt.subplots(1, 2)
#axs[0].plot(np.arange(z_cells * y_cells * x_cells), np.sort(vort_distribution_mean))
#axs[0].grid(True, linestyle='--')
#axs[1].plot(np.arange(count_of_vort_steps), np.sort(vort_mean_distribution))
#axs[1].grid(True, linestyle='--')
#plt.show()
