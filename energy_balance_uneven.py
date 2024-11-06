import os
import numpy as np
from gc import collect
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy
# import math

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)
directory_name = path.split('/')[-1]
print(filename)
print(path)
print(directory_name)
#path = '/home/sq-wm/OpenFOAM/sq-wm-v2306/run'

#path += '/calc13'
path += '/'
print(path)

f = open(path + '/calculation_parameters.txt', 'r')
calc_param = f.read().split('\n')

start_time = float(calc_param[0].split(' ')[1])
end_time = float(calc_param[1].split(' ')[1])
time_step = float(calc_param[2].split(' ')[1])
grid_x = int(calc_param[3].split(' ')[1])
grid_y = int(calc_param[3].split(' ')[2])
grid_z = int(calc_param[3].split(' ')[3])
left_boundary_x = float(calc_param[4].split(' ')[1])
right_boundary_x = float(calc_param[4].split(' ')[2])
down_boundary_y = float(calc_param[5].split(' ')[1])
up_boundary_y = float(calc_param[5].split(' ')[2])
rear_boundary_z = float(calc_param[6].split(' ')[1])
front_boundary_z = float(calc_param[6].split(' ')[2])
nu = float(calc_param[7].split(' ')[1])
omega = float(calc_param[8].split(' ')[1])
f0 = float(calc_param[9].split(' ')[1])
x_cells = int(calc_param[10].split(' ')[1])
y_cells = int(calc_param[10].split(' ')[2])
z_cells = int(calc_param[10].split(' ')[3])
f.close()
print('end reading system info')
kx = 1.5
ky = 2.5
kz = 4.0
Ax = 5.0
Ay = 3.0
gradient_coeff = 10.0

def convert_time_name(x):
    if x%1.0 == 0.0:
        return round(x)
    else:
        return x

time = np.round(np.arange(start_time, end_time + time_step, time_step), decimals=2)
time = list(map(convert_time_name, time))
print(time)

def cells_bound(gradient_coeff):
    multiplyer_x = pow(gradient_coeff, 2.0/(grid_x - 2))
    start_x_size = (right_boundary_x - left_boundary_x) * (1.0 - multiplyer_x) / (2.0 - 2.0 * pow(multiplyer_x, grid_x // 2))
    x_grad_cells = np.empty(grid_x // 2, dtype=float)
    x_bounds = np.empty(grid_x + 1)

    multiplyer_y = pow(gradient_coeff, 2.0/(grid_y - 2))
    start_y_size = (up_boundary_y - down_boundary_y) * (1.0 - multiplyer_y) / (2.0 - 2.0 * gradient_coeff)
    y_grad_cells = np.empty(grid_y // 2, dtype=float)
    y_bounds = np.empty(grid_y + 1)

    multiplyer_z = pow(gradient_coeff, 2.0/(grid_z - 2))
    start_z_size = (front_boundary_z - rear_boundary_z) * (1.0 - multiplyer_z) / (2.0 - 2.0 * gradient_coeff)
    z_grad_cells = np.empty(grid_z // 2, dtype=float)
    z_bounds = np.empty(grid_z + 1)

    for i in range(grid_x // 2):
        x_grad_cells[i] = start_x_size * pow(multiplyer_x, i)
    x_grad_cells = np.concatenate([x_grad_cells, np.flip(x_grad_cells)])
    x_bounds[0] = left_boundary_x
    for i in range(1, grid_x + 1):
        x_bounds[i] = x_bounds[i - 1] + x_grad_cells[i - 1]

    for i in range(grid_y // 2):
        y_grad_cells[i] = start_y_size * pow(multiplyer_y, i)
    y_grad_cells = np.concatenate([y_grad_cells, np.flip(y_grad_cells)])
    y_bounds[0] = down_boundary_y
    for i in range(1, grid_y + 1):
        y_bounds[i] = y_bounds[i - 1] + y_grad_cells[i - 1]

    for i in range(grid_z // 2):
        z_grad_cells[i] = start_z_size * pow(multiplyer_z, i)
    z_grad_cells = np.concatenate([z_grad_cells, np.flip(z_grad_cells)])
    z_bounds[0] = rear_boundary_z
    for i in range(1, grid_z + 1):
        z_bounds[i] = z_bounds[i - 1] + z_grad_cells[i - 1]

    del x_grad_cells, y_grad_cells, z_grad_cells
    collect()
    return x_bounds, y_bounds, z_bounds

def binary_found_bounds(bounds_array, point):
    ind_min = np.where(bounds_array < point)[0][-1]
    return ind_min, ind_min + 1

def cells_volume(bounds_x, bounds_y, bounds_z, point_x, point_y, point_z):
    ind1, ind2 = binary_found_bounds(bounds_x, point_x)
    d_x = bounds_x[ind2] - bounds_x[ind1]
    ind1, ind2 = binary_found_bounds(bounds_y, point_y)
    d_y = bounds_y[ind2] - bounds_y[ind1]
    ind1, ind2 = binary_found_bounds(bounds_z, point_z)
    d_z = bounds_z[ind2] - bounds_z[ind1]
    return d_x * d_y * d_z

#dx = (right_boundary_x - left_boundary_x) / x_cells
#dy = (up_boundary_y - down_boundary_y) / y_cells
#dz = (front_boundary_z - rear_boundary_z) / z_cells
dv_characteric = (right_boundary_x - left_boundary_x) * (up_boundary_y - down_boundary_y) * (front_boundary_z - rear_boundary_z) / (grid_x * grid_y * grid_z)

def get_data_from_file(full_path_to_files, size_x, size_y, size_z):
    file = open(full_path_to_files, 'r')
    data = file.read().split('\n')
    slice = []
    u = []
    for ind_z in range(size_z):
        for ind_y in range(size_y):
            slice.append(data[ind_z*(size_y+1) + ind_y].split(' ')[ : size_x])
        u.append(slice)
        slice = []
    file.close()

    del slice, data
    collect()
    print("end reading velocity from file " + full_path_to_files)
    u = np.array(u, dtype=float)
    return u

def kinetic_energy_in_cells(u, v, w, dv):
    return 0.5 * np.sum(u*u + v*v + w*w) * dv        # rho = 1.0

def write_energy_distr(distr, n, name_with_path):
    file = open(name_with_path, "w")
    file.write(str(n) + '\n')
    for ind in range(n):
        file.write(str(distr[ind]) + '\n')
    file.close()

def get_distibution_from_file(name_with_path):
    file = open(name_with_path, "r")
    data = file.read().split('\n')
    file.close()
    return int(data[0]), np.array(data[1 : -1], dtype=float)

def create_force(x, y, z, f_0):
    force_x = Ax * np.cos(kx * x) * np.sin(ky * y) * np.sin(kz * z) / kx
    force_y = Ay * np.sin(kx * x) * np.cos(ky * y) * np.sin(kz * z) / ky
    force_z = - (Ax + Ay) * np.sin(kx * x) * np.sin(ky * y) * np.cos(kz * z) / kz

    return f_0 * force_x, f_0 * force_y, f_0 * force_z

# Read the source file
def get_arrays_with_vel_and_cells(name_of_file):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(name_of_file)
    reader.Update()  # Needed because of GetScalarRange
    U_field = reader.GetOutput().GetPointData().GetArray('U')
    data_vel = vtk_to_numpy(U_field)        # 2D array with components of velocity in points. data[i] = v[x, y, z] in i-point

    nodes_vtk_array = reader.GetOutput().GetPoints().GetData()
    numpy_nodes = vtk_to_numpy(nodes_vtk_array)
    count_of_cells = len(numpy_nodes)      # = 151^3
    print(count_of_cells)
    del U_field
    del nodes_vtk_array
    collect()
    return data_vel, numpy_nodes

def get_arrays_with_grad_vel(name_of_file):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(name_of_file)
    reader.Update()  # Needed because of GetScalarRange
    grad_field = reader.GetOutput().GetPointData().GetArray('grad(U)')
    data_grad = vtk_to_numpy(grad_field)        # 2D array with components of velocity in points. data[i] = v[x, y, z] in i-point

    del grad_field
    collect()
    return data_grad

def get_arrays_with_vel_grad_vel_and_cells(name_of_file):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(name_of_file)
    reader.Update()  # Needed because of GetScalarRange

    U_field = reader.GetOutput().GetPointData().GetArray('U')
    data_vel = vtk_to_numpy(U_field)        # 2D array with components of velocity in points. data[i] = v[x, y, z] in i-point
    grad_field = reader.GetOutput().GetPointData().GetArray('grad(U)')
    data_grad = vtk_to_numpy(grad_field)
    nodes_vtk_array = reader.GetOutput().GetPoints().GetData()
    numpy_nodes = vtk_to_numpy(nodes_vtk_array)

    count_of_cells = len(numpy_nodes)      # = 151^3
    #print(count_of_cells)
    del U_field, nodes_vtk_array, grad_field
    collect()
    return data_vel, data_grad, numpy_nodes

#for i in range(len(cells)):
cells_count = int(max(x_cells, y_cells, z_cells))
dissipation = []
pumping_power = []
U2 = []
U2XY = []
E_kinetic = []

table_name = '/table_with_statistics.txt'
file_table = open(path + table_name, "w")
file_table.write('time' + '\t' + 'U2' + '\t' + 'U2_2D' + '\t' + 'dissipation' + '\t' + 'pumping_power' + '\t' + 'kinetic energy' + '\n')

times = []
old_names = []

# create bash command for convert foam format to vtk
command_str = "foamToVTK -fields '(U grad(U))' -time '"
for i in time:
    if i != end_time:
        command_str += str(i) + '; '
    else:
        command_str += str(i) + "'"
print(command_str)
os.system(command_str)

with open(path + '/VTK/' + directory_name + '.vtm.series', 'r') as fp:
    vtk_files_info = fp.read().split('\n')
# rename vtk files
times = []
old_names = []
count_of_files = len(vtk_files_info) - 6
for i in range(3, 3+count_of_files):
    old_name = vtk_files_info[i].split(' ')[7].split('"')[1]
    old_name = old_name[ : len(old_name)-4]
    old_names.append(old_name)
    time_name = vtk_files_info[i].split(' ')[10]
    times.append(time_name)
    os.rename(path + '/VTK/' + old_name, path + '/VTK/' + time_name)
    os.rename(path + '/VTK/' + old_name + '.vtm', path + '/VTK/' + time_name + '.vtm')
print('vtk files were renaming')

for time_ind in time:
    U_vector, grad_U, points = get_arrays_with_vel_grad_vel_and_cells(path + '/VTK/' + str(time_ind) + '/internal.vtu')

    x_cells_bounds, y_cells_bounds, z_cells_bounds = cells_bound(gradient_coeff)

    E = 0.0
    U_rms = 0.0
    U_rms_2D = 0.0
    dissipation_ = 0.0
    pumping_power_ = 0.0
    for i in range(len(U_vector)):
        u = U_vector[i, 0]
        v = U_vector[i, 1]
        w = U_vector[i, 2]

        x = points[i, 0]
        y = points[i, 1]
        z = points[i, 2]

        #du/dx, dv/dx, dw/dx, du/dy, dv/dy, dw/dy, du/dz, dv/dz, dw/dz
        #du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz, dw/dx, dw/dy, dw/dz - ?
        grad_u = grad_U[i, 0]
        grad_v = grad_U[i, 4]
        grad_w = grad_U[i, 8]

        dv = cells_volume(x_cells_bounds, y_cells_bounds, z_cells_bounds, x, y, z)
        force_x, force_y, force_z = create_force(x, y, z, f0)

        E += kinetic_energy_in_cells(u, v, w, dv)
        pumping_power_ += (force_x * u + force_y * v + force_z * w) * dv
        U_rms += (u*u + v*v + w*w) * dv
        U_rms_2D += (u*u + v*v) * dv
        dissipation_ += (grad_u*grad_u + grad_v*grad_v + grad_w*grad_w) * dv

    E_kinetic.append(E)                                         # E = sum( V^2 *rho * dv / 2)
    pumping_power.append(pumping_power_ / dv_characteric)       # pumping power = sum( (fx*u + fy*v + fz*w) * dv ) / dv_
    dissipation.append(nu * dissipation_ / dv_characteric)      # dissipation = nu * sum( ((du/dx)^2 + (dv/dy)^2 + (dw/dz)^2) * dv) / dv_
    U2.append(U_rms / dv_characteric)                           # U rms = sum( (Vx^2 + Vy^2 + Vz^2) * dv ) / dv_
    U2XY.append(U_rms_2D / dv_characteric)                      # U rms 2D = sum( (Vx^2 + Vy^2) * dv ) / dv_

    file_table.write(str(time_ind) + '\t' + str(U_rms) + '\t' + str(U_rms_2D) + '\t' + str(nu * dissipation_) + '\t' + str(pumping_power_) + '\t' + str(E) + '\n')

file_table.close()
#print('cells size = ' + str(cells_count))
#print(pumping_power)
#print(dissipation)

size_of_text = 15
size_of_numbers = 15
size_of_legend = 20

fig, axs = plt.subplots(figsize=(10, 6), dpi=600)
axs.plot(time, pumping_power, label = r'$ \langle \bar{f} \cdot \bar{U} \rangle $', linewidth=1, linestyle='-')
axs.plot(time, dissipation, label = r'$ \nu \langle \vert \nabla \bar{U} \vert^{2} \rangle $', linewidth=1, linestyle='--',)

axs.grid(True, linestyle='--')
axs.set_ylabel('dissipation, force pumping', fontsize=size_of_text)
axs.set_xlabel('t', fontsize=size_of_text)
axs.legend(fontsize=size_of_legend)
plt.savefig(path + '/pumping_dissipation.png')

U2 = np.array(U2)
U2XY = np.array(U2XY)

fig, axs = plt.subplots(figsize=(10, 6), dpi=600)
axs.plot(time, U2, label = r'$ U^{2} $', linewidth=3, color='blue')
axs.plot(time, U2XY, label = r'$ U^{2}_{2D} $', linewidth=3, color='orange')
axs.grid(True, linestyle='--')
#axs.set_ylabel(r'$ \Sigma $' + ' , ' + r'$ \omega_{z} $', fontsize=size_of_text)
axs.set_xlabel('t', fontsize=size_of_text)
axs.tick_params(labelsize=size_of_numbers)
#plt.xlim(0.0, max_radius * 1.0 * (right_boundary_x - left_boundary_x) / x_cells * 1.1)
plt.ylim(0.0, np.max(U2) * 1.1)
axs.legend(fontsize=size_of_legend)
plt.savefig(path + '/U2_U2D.png')

fig, axs = plt.subplots(figsize=(10, 6), dpi=600)
axs.plot(time, E_kinetic, label = r'$ \sum \frac{ \rho U^{2} }{2} dv $', linewidth=3, color='blue')
axs.grid(True, linestyle='--')
#axs.set_ylabel(r'$ \Sigma $' + ' , ' + r'$ \omega_{z} $', fontsize=size_of_text)
axs.set_xlabel('t', fontsize=size_of_text)
axs.tick_params(labelsize=size_of_numbers)
#plt.xlim(0.0, max_radius * 1.0 * (right_boundary_x - left_boundary_x) / x_cells * 1.1)
plt.ylim(0.0, np.max(E_kinetic) * 1.1)
axs.legend(fontsize=size_of_legend)
plt.savefig(path + '/kinetic_energy.png')
