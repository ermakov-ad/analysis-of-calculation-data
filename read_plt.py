import matplotlib.pyplot as plt
import numpy as np
import struct
from scipy.fft import fftn, fft, fftfreq
import math

# info from parameters.dat file:
left_boundary_x =   -3.1415
right_boundary_x =  3.1415
down_boundary_y =   -3.1415
up_boundary_y =     3.1415
rear_boundary_z =   -3.1415
front_boundary_z =  3.1415
cells_number_x = 100
cells_number_y = 100
cells_number_z = 400
time_step = 0.001

dx = (right_boundary_x - left_boundary_x) / cells_number_x
dy = (up_boundary_y - down_boundary_y) / cells_number_y
dz = (front_boundary_z - rear_boundary_z) / cells_number_z

def construct_file_path(n):
    path = 'data_numerical_solution\\numerical_solution_'
    l = int(len(str(n)))
    help_string = ''
    for i in range(0, 6-l, 1):
        help_string += '0'
    return path + help_string + str(n) + '.plt'

# return: 
# x, y, z = np.array[z][y][x]; size = (cells_z + 1)*(cells_y + 1)*(cells_x + 1)
# p, u, v, w = np.array[z][y][x]; size = cells_z * cells_y * cells_x
def get_data_from_file(file_number):
    file = open(construct_file_path(file_number), 'rb')
    # skip header
    header = file.read(284)
    # print("header text:")
    # print(header)
    integers = []
    for i in range(0, 28, 1):
        integers.append(int(struct.unpack('<i', file.read(4))[0]))
    # one last int, which not exists in first files
    if integers[0] != 0:
        integers.append(int(struct.unpack('<i', file.read(4))[0]))
    #print("system information from header: ")
    #print(integers)

    # read data in double format
    # check count of elements:
    # count of x, y, z: 101*101*401
    # count of p, u, v, w: 100*100*400
    # it = int(0)
    # try:
    #     while it >= 0:
    #         value = file.read(8)
    #         double_value = struct.unpack('<d', value)[0]
    #         structured_data.append(float(double_value))
    #         it += 1

    # read arrays with coordinates and normal variables
    x = np.empty([cells_number_z + 1, cells_number_y + 1, cells_number_x + 1], dtype=float)
    y = np.empty([cells_number_z + 1, cells_number_y + 1, cells_number_x + 1], dtype=float)
    z = np.empty([cells_number_z + 1, cells_number_y + 1, cells_number_x + 1], dtype=float)
    p = np.empty([cells_number_z, cells_number_y, cells_number_x], dtype=float)
    u = np.empty([cells_number_z, cells_number_y, cells_number_x], dtype=float)
    v = np.empty([cells_number_z, cells_number_y, cells_number_x], dtype=float)
    w = np.empty([cells_number_z, cells_number_y, cells_number_x], dtype=float)

    p_new = np.empty([cells_number_x, cells_number_x, cells_number_x], dtype=float) # change array's size for 100 * 100 * 100
    u_new = np.empty([cells_number_x, cells_number_x, cells_number_x], dtype=float)
    v_new = np.empty([cells_number_x, cells_number_x, cells_number_x], dtype=float)
    w_new = np.empty([cells_number_x, cells_number_x, cells_number_x], dtype=float)

    for z_ind in range(0, cells_number_z + 1):
        for y_ind in range(0, cells_number_y + 1):
            for x_ind in range(0, cells_number_x + 1):
                value = file.read(8)
                x[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])

    for z_ind in range(0, cells_number_z + 1):
        for y_ind in range(0, cells_number_y + 1):
            for x_ind in range(0, cells_number_x + 1):
                value = file.read(8)
                y[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])

    for z_ind in range(0, cells_number_z + 1):
        for y_ind in range(0, cells_number_y + 1):
            for x_ind in range(0, cells_number_x + 1):
                value = file.read(8)
                z[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                p[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    for it in range(0, cells_number_x):
        p_new[it] = (9*p[4*it + 1] + 9*p[4*it + 2] - p[4*it + 0] - p[4*it + 3]) / 16
    p = p_new

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                u[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    for it in range(0, cells_number_x):
        u_new[it] = (9*u[4*it + 1] + 9*u[4*it + 2] - u[4*it + 0] - u[4*it + 3]) / 16
    u = u_new

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                v[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    for it in range(0, cells_number_x):
        v_new[it] = (9*v[4*it + 1] + 9*v[4*it + 2] - v[4*it + 0] - v[4*it + 3]) / 16
    v = v_new

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                w[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    for it in range(0, cells_number_x):
        w_new[it] = (9*w[4*it + 1] + 9*w[4*it + 2] - w[4*it + 0] - w[4*it + 3]) / 16
    w = w_new

    mistake = file.read()
    if len(mistake) > 0:
        print("remained " + str(len(mistake)) + " bytes in the end of the file")
        print(mistake)
        print("something doesn't work")
    file.close()
    return x, y, z, p, u, v, w

# chose time interval
# start_time = 0.0
# end_time = 0.003
start_time = 0.231
end_time = 0.234

print(int(start_time / time_step))
print(int(end_time / time_step + 1))

# convert 3D coordinates to 1D-index for array x/y/z
def coordinates_to_index_xyz_arrays(x_coordinate, y_coordinate, z_coordinate):
    if z_coordinate > 0.0:
        z_number = int((z_coordinate - rear_boundary_z) / dz) + 1 # count of layers by axes z from z_start to z_coordinate
    else:
        z_number = int((z_coordinate - rear_boundary_z) / dz)
    z_index = int((cells_number_x + 1) * (cells_number_y + 1) * z_number)
    if y_coordinate > 0.0:
        y_number = int((y_coordinate - down_boundary_y) / dy) + 1   # count of rows by axes y from y_start to y_coordinate
    else:
        y_number = int((y_coordinate - down_boundary_y) / dy)
    y_index = int((cells_number_x + 1) * y_number)
    if x_coordinate > 0.0:
        x_index = int((x_coordinate - left_boundary_x) / dx) + 1
    else:
        x_index = int((x_coordinate - left_boundary_x) / dx)
    return z_index + y_index + x_index

# convert 3D coordinates to 1D-index for array p/u/v/w
def coordinates_to_index(x_coordinate, y_coordinate, z_coordinate):
    z_number = int((z_coordinate - rear_boundary_z) / dz)       # count of layers by axes z from z_start to z_coordinate
    z_index = int(cells_number_x * cells_number_y * z_number)
    y_number = int((y_coordinate - down_boundary_y) / dy)       # count of rows by axes y from y_start to y_coordinate
    y_index = int(cells_number_x * y_number)
    x_index = int((x_coordinate - left_boundary_x) / dx)
    return z_index + y_index + x_index

# find the energy spectrum at the choosen moment of time
# E(q) = summ(Vk^2) / (2*pi)^3, if |k - q| <= 1/4;
# where q = sqrt(q1*q1 + q2*q2 + q3*q3), k = sqrt(kx*kx + ky*ky + kz*kz)
# Vk^2 = (u^2 + v^2 + w^2); u, v, w - projections of spectum of speed; 
# q - vector of initial conditions
def modulus_of_vector(vec, l):
    sqr = 0.0
    for i in range(0, l):
        sqr += vec[i] * vec[i]
    return math.sqrt(sqr)

def sqr_of_vector(vec, l):
    sqr = 0.0
    for i in range(0, l):
        sqr += vec[i] * vec[i]
    return sqr

# spectrum of energy depend on vector q (E(q))
max_wave_deviation = 0.25     # constant in this calculation
def find_energy_spectrum(u, v, w, size_z, size_y, size_x, count_of_cells, freq_z, freq_y, freq_x, amplitude):
    spectrum = np.zeros(count_of_cells)
    for ind_z in range(0, size_z):
        for ind_y in range(0, size_y):
            for ind_x in range(0, size_x):
                k = modulus_of_vector([freq_z[ind_z], freq_y[ind_y], freq_x[ind_x]], 3)
                i = int((k + max_wave_deviation) * count_of_cells / amplitude)
                E = sqr_of_vector([u[ind_z][ind_y][ind_x], v[ind_z][ind_y][ind_x], w[ind_z][ind_y][ind_x]], 3)
                if i >= count_of_cells:
                    i = count_of_cells - 1
                while abs(k - i * amplitude / count_of_cells) <= max_wave_deviation and i >= 0:
                    spectrum[i] += E
                    i -= 1

    return spectrum / (8.0 * math.pi * math.pi * math.pi)

# find full kinetic energy in part of volume
def find_full_energy(u, v, w, size_z, size_y, size_x):
    E = 0.0
    for ind_z in range(0, size_z):
        for ind_y in range(0, size_y):
            for ind_x in range(0, size_x):
                E += sqr_of_vector([u[ind_z][ind_y][ind_x], v[ind_z][ind_y][ind_x], w[ind_z][ind_y][ind_x]], 3)

    return E / (16.0 * math.pi * math.pi * math.pi)

# calculation of wave vectors arrays along coordinate axes
kx = np.arange(cells_number_x // 2)
ky = np.arange(cells_number_y // 2)
kz = np.arange(cells_number_z // 4 // 2)

energy_spectrum = []
amplitude_wave_vec = modulus_of_vector([cells_number_z // 4 // 2 - 1, cells_number_y // 2 - 1, cells_number_x // 2 - 1], 3)

def find_velocity_vec(file_number):
    x, y, z, p, u_center, v_center, w_center = get_data_from_file(file_number)
    print("end reading " + str(file_number) + " time itteration")

    # fourier image of velocity components
    fourier_u = np.real(fftn(u_center))
    fourier_v = np.real(fftn(v_center))
    fourier_w = np.real(fftn(w_center))
    
    fourier_u = fourier_u[:cells_number_z // 4 // 2][:cells_number_y // 2][:cells_number_x // 2]
    fourier_v = fourier_v[:cells_number_z // 4 // 2][:cells_number_y // 2][:cells_number_x // 2]
    fourier_w = fourier_w[:cells_number_z // 4 // 2][:cells_number_y // 2][:cells_number_x // 2]
    # vector velocity(k) = [fourier_u, fourier_v, fourier_w]
    return fourier_u, fourier_v, fourier_w

def find_spectrum_from_file(file_number, count_of_vec_points):
    fourier_u, fourier_v, fourier_w = find_velocity_vec(file_number)
    # find E(q) in this file:
    E_i = find_energy_spectrum(fourier_u, fourier_v, fourier_w, cells_number_z // 4 // 2, cells_number_y // 2, cells_number_x // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec)
    max_energy = np.max(E_i)
    q_max = np.argmax(E_i) * amplitude_wave_vec / count_of_vec_points
    full_e = find_full_energy(fourier_u, fourier_v, fourier_w, cells_number_z // 4 // 2, cells_number_y // 2, cells_number_x // 2)
    return E_i, max_energy, q_max, full_e

def spherical_coord_from_vec(x, y, z):
    xy = x**2 + y**2
    r = np.sqrt(xy + z**2)
    theta = np.arctan2(np.sqrt(xy), z) # for elevation angle defined from Z-axis down
    phi = np.arctan2(y, x)  # for elevation angle defined from XY-plane up
    return r, theta, phi

u_k = []
v_k = []
w_k = []

count_of_vec_steps = 200
# cycle according to the time we are interested in
for i in range(0, int(end_time / time_step + 1) - int(start_time / time_step), 1):    
    energy_sp, maximum_E, Q_max, full_E = find_spectrum_from_file(i + int(start_time / time_step), count_of_vec_steps)
    print("time moment: " + str((i + int(start_time / time_step)) * time_step))
    print("full kinetic energy = " + str(full_E))
    print("maximum spectrum of energy = " + str(maximum_E) + " in wave vector = " + str(Q_max))
    energy_spectrum.append(energy_sp)

    vel_vec = find_velocity_vec(i + int(start_time / time_step))
    u_k.append(vel_vec[0])
    v_k.append(vel_vec[1])
    w_k.append(vel_vec[2])
energy_spectrum = np.array(energy_spectrum)

l = int(end_time / time_step + 1) - int(start_time / time_step)
u_wk = np.real(fft(u_k))
v_wk = np.real(fft(v_k))
w_wk = np.real(fft(w_k))
u_wk = u_wk[:l//2]
v_wk = v_wk[:l//2]
w_wk = w_wk[:l//2]
freq_ = fftfreq(l, time_step)[:l//2]    # w 
# v (w, k) = (u, v, w)[w][kz, ky, kx]
# A_w_theta = summ (v(w, k)^2), where |theta_k - theta| < delta_theta / 2
delta_theta = 0.15
count_of_theta = int(math.pi / delta_theta) + 1
A_w_theta = []
theta_array = np.linspace(0.0, math.pi, num=count_of_theta)
for i in range(0, l//2):
    A = np.zeros(count_of_theta)
    for ind_z in range(0, cells_number_z // 4 // 2):
        for ind_y in range(0, cells_number_y // 2):
            for ind_x in range(0, cells_number_x // 2):
                for theta_ind in range(0, count_of_theta):
                    mod_k, theta_k, phi = spherical_coord_from_vec(ind_x, ind_y, ind_z)
                    if (abs(theta_array[theta_ind] - theta_k) < delta_theta / 2):
                        A[theta_ind] += sqr_of_vector([u_wk[i][ind_z][ind_y][ind_x], v_wk[i][ind_z][ind_y][ind_x], w_wk[i][ind_z][ind_y][ind_x]], 3)
    A_w_theta.append(A)
A_w_theta = np.array(A_w_theta)

# draw energy spectrum
wave_coord = np.linspace(0, amplitude_wave_vec + max_wave_deviation, count_of_vec_steps)

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 14))
fig.suptitle('energy spectrum', fontsize=14)

axs[0, 0].set_title("time = " + str((0 + int(start_time / time_step)) * time_step))
axs[0, 0].plot(np.log10(wave_coord), np.log10(energy_spectrum[0]))
axs[0, 0].set_ylabel('log(E)')
axs[0, 0].set_xlabel('lengh of wave vector')
axs[0, 0].grid(True)

axs[0, 1].set_title("time = " + str((1 + int(start_time / time_step)) * time_step))
axs[0, 1].plot(np.log10(wave_coord), np.log10(energy_spectrum[1]))
axs[0, 1].set_ylabel('log(E)')
axs[0, 1].set_xlabel('lengh of wave vector')
axs[0, 1].grid(True)

axs[1, 0].set_title("time = " + str((2 + int(start_time / time_step)) * time_step))
axs[1, 0].plot(np.log10(wave_coord), np.log10(energy_spectrum[2]))
axs[1, 0].set_ylabel('log(E)')
axs[1, 0].set_xlabel('lengh of wave vector')
axs[1, 0].grid(True)

axs[1, 1].set_title("time = " + str((3 + int(start_time / time_step)) * time_step))
axs[1, 1].plot(np.log10(wave_coord), np.log10(energy_spectrum[3]))
axs[1, 1].set_ylabel('log(E)')
axs[1, 1].set_xlabel('lengh of wave vector')
axs[1, 1].grid(True)

plt.show()

# draw inertial waves
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10, 14))
fig.suptitle('amplitudes of inertial waves in phase space', fontsize=14)

axs[0].set_title("frequency = " + str(freq_[0]))
axs[0].plot(theta_array, A_w_theta[0])
axs[0].set_ylabel('A(frequency, angle)')
axs[0].set_xlabel('angle (theta)')
axs[0].grid(True)

axs[1].set_title("frequency = " + str(freq_[1]))
axs[1].plot(theta_array, A_w_theta[0])
axs[1].set_ylabel('A(frequency, angle)')
axs[1].set_xlabel('angle (theta)')
axs[1].grid(True)

plt.show()

print("frequency = " + str(freq_[0]))
print("maximum of A = " + str(np.max(A_w_theta[0])) + " in angle "+ str(theta_array[np.argmax(A_w_theta[0])]))
print("frequency = " + str(freq_[1]))
print("maximum of A = " + str(np.max(A_w_theta[1])) + " in angle "+ str(theta_array[np.argmax(A_w_theta[1])]))
