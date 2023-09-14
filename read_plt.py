import matplotlib.pyplot as plt
import numpy as np
import struct
from scipy.fft import fftn, fftfreq
import math
from scipy.signal import resample

# info from parameters.dat file:
left_boundary_x =   -1.0*math.pi
right_boundary_x =  math.pi
down_boundary_y =   -1.0*math.pi
up_boundary_y =     math.pi
rear_boundary_z =   -1.0*math.pi
front_boundary_z =  math.pi
cells_number_x = 100
cells_number_y = 100
cells_number_z = 400
time_step = 0.1
#cells_number_x = 200
#cells_number_y = 200
#cells_number_z = 500
#time_step = 0.0004

# chose time interval
start_number = 195      # numbers of file
end_number = 234        # the number of the last file inclusive

kinematic_viscosity = 0.01  # table value, a constant for a given environment. 
OMEGA = 10.0                # angular rotation speed - parameters of the simulated system. Info from parameters.dat file
max_wave_deviation_in_energy_calc = 0.25     # constants in this calculation
max_wave_deviation_in_dencity_velocity_dissipation_calc = 0.5

min_cells_number = min([cells_number_x, cells_number_y, cells_number_z])
dx = (right_boundary_x - left_boundary_x) / min_cells_number
dy = (up_boundary_y - down_boundary_y) / min_cells_number
dz_new = (front_boundary_z - rear_boundary_z) / min_cells_number

# name of txt file, which save information about Reynolds and Rossby number, characteristic wave vector, energy and dissipation
txt_file_name = "system_parameters.txt"

def construct_file_path(n):
    path = 'C:\\Users\\Admin\\Documents\\analysis_of_calculations_of_column_vortices\\data_numerical_solution\\numerical_solution_'
    #path = 'C:\\Users\\Admin\\Documents\\analysis_of_calculations_of_column_vortices\\data_numerical_solution2\\numerical_solution_'
    l = int(len(str(n)))
    help_string = ''
    for i in range(0, 6-l, 1):
        help_string += '0'
    return path + help_string + str(n) + '.plt'

# return: 
# x, y, z = np.array[z][y][x]; size = (cells_z + 1)*(cells_y + 1)*(cells_x + 1)
# p, u, v, w = np.array[z][y][x]; size = cells_z * cells_y * cells_x
# in practice, programm work with cubic arrays
# so new arrays size = min_cells_number * min_cells_number * min_cells_number, which function returned
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

    p_new = np.empty([cells_number_x, cells_number_x, cells_number_x], dtype=float) # change array's size for 200 * 200 * 200
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
    # for 400 -> 100:
    #for it in range(0, cells_number_x):
    #    p_new[it] = (9*p[4*it + 1] + 9*p[4*it + 2] - p[4*it + 0] - p[4*it + 3]) / 16
    # for 500 -> 200:
    #for it in range(0, cells_number_x, 2):
    #    p_new[it] = 415.0*p[5*it//2 + 0]/256.0 - 635.0*p[5*it//2 + 1]/128.0 + 141.0*p[5*it//2 + 2]/16.0 - 765.0*p[5*it//2 + 3]/128 + 385.0*p[5*it//2 + 4]/256.0
    #    p_new[it + 1] = 385.0*p[5*it//2 + 0]/256.0 - 765.0*p[5*it//2 + 1]/128.0 + 141.0*p[5*it//2 + 2]/16.0 - 635.0*p[5*it//2 + 3]/128 + 415.0*p[5*it//2 + 4]/256.0

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                u[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    # for 400 -> 100:
    #for it in range(0, cells_number_x):
    #    u_new[it] = (9*u[4*it + 1] + 9*u[4*it + 2] - u[4*it + 0] - u[4*it + 3]) / 16
    # for 500 -> 200:
    #for it in range(0, cells_number_x, 2):
    #    u_new[it] = 415.0*u[5*it//2 + 0]/256.0 - 635.0*u[5*it//2 + 1]/128.0 + 141.0*u[5*it//2 + 2]/16.0 - 765.0*u[5*it//2 + 3]/128 + 385.0*u[5*it//2 + 4]/256.0
    #    u_new[it + 1] = 385.0*u[5*it//2 + 0]/256.0 - 765.0*u[5*it//2 + 1]/128.0 + 141.0*u[5*it//2 + 2]/16.0 - 635.0*u[5*it//2 + 3]/128 + 415.0*u[5*it//2 + 4]/256.0

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                v[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    # for 400 -> 100:
    #for it in range(0, cells_number_x):
    #    v_new[it] = (9*v[4*it + 1] + 9*v[4*it + 2] - v[4*it + 0] - v[4*it + 3]) / 16
    # for 500 -> 200:
    #for it in range(0, cells_number_x, 2):
    #    v_new[it] = 415.0*v[5*it//2 + 0]/256.0 - 635.0*v[5*it//2 + 1]/128.0 + 141.0*v[5*it//2 + 2]/16.0 - 765.0*v[5*it//2 + 3]/128 + 385.0*v[5*it//2 + 4]/256.0
    #    v_new[it + 1] = 385.0*v[5*it//2 + 0]/256.0 - 765.0*v[5*it//2 + 1]/128.0 + 141.0*v[5*it//2 + 2]/16.0 - 635.0*v[5*it//2 + 3]/128 + 415.0*v[5*it//2 + 4]/256.0

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                w[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    # for 400 -> 100:
    #for it in range(0, cells_number_x):
    #    w_new[it] = (9*w[4*it + 1] + 9*w[4*it + 2] - w[4*it + 0] - w[4*it + 3]) / 16
    # for 500 -> 200:
    #for it in range(0, cells_number_x, 2):
    #    w_new[it] = 415.0*w[5*it//2 + 0]/256.0 - 635.0*w[5*it//2 + 1]/128.0 + 141.0*w[5*it//2 + 2]/16.0 - 765.0*w[5*it//2 + 3]/128 + 385.0*w[5*it//2 + 4]/256.0
    #    w_new[it + 1] = 385.0*w[5*it//2 + 0]/256.0 - 765.0*w[5*it//2 + 1]/128.0 + 141.0*w[5*it//2 + 2]/16.0 - 635.0*w[5*it//2 + 3]/128 + 415.0*w[5*it//2 + 4]/256.0

    # resampling from z axis
    if cells_number_z > min_cells_number:
        p_new = resample(p, min_cells_number)
        u_new = resample(u, min_cells_number)
        v_new = resample(v, min_cells_number)
        w_new = resample(w, min_cells_number)
    # resampling from y axis, along the z axis, the size of array is already correct
    if cells_number_y > min_cells_number:
        p_new[:] = resample(p[:], min_cells_number)
        u_new[:] = resample(u[:], min_cells_number)
        v_new[:] = resample(v[:], min_cells_number)
        w_new[:] = resample(w[:], min_cells_number)
    # resampling from x axis, along the z and y axes, the size of array is already correct
    if cells_number_x > min_cells_number:
        p_new[: , :] = resample(p[: , :], min_cells_number)
        u_new[: , :] = resample(u[: , :], min_cells_number)
        v_new[: , :] = resample(v[: , :], min_cells_number)
        w_new[: , :] = resample(w[: , :], min_cells_number)
    p = p_new
    u = u_new
    v = v_new
    w = w_new

    mistake = file.read()
    if len(mistake) > 0:
        print("remained " + str(len(mistake)) + " bytes in the end of the file")
        print(mistake)
        print("something doesn't work")
    file.close()
    return x, y, z, p, u, v, w

# find the energy spectrum at the choosen moment of time
# E(q) = summ(Vk^2) / (2*pi)^3, if |k - q| <= 1/4;
# k = sqrt(kx*kx + ky*ky + kz*kz)
# Vk^2 = (u^2 + v^2 + w^2); u, v, w - projections of spectum of speed; 
# q - array of modulus of vectors 
def find_energy_spectrum(u, v, w, size_of_array, count_of_cells, freq_z, freq_y, freq_x, amplitude, max_wave_deviation_):
    spectrum = np.zeros(count_of_cells)
    for ind_z in range(0, size_of_array):
        for ind_y in range(0, size_of_array):
            for ind_x in range(0, size_of_array):
                k = np.linalg.norm(np.array([freq_z[ind_z], freq_y[ind_y], freq_x[ind_x]]))
                i = int((k + max_wave_deviation_) * count_of_cells / amplitude)
                E = np.sum(np.square(np.array([u[ind_z][ind_y][ind_x], v[ind_z][ind_y][ind_x], w[ind_z][ind_y][ind_x]])))
                if i >= count_of_cells:
                    i = count_of_cells - 1
                while abs(k - i * amplitude / count_of_cells) <= max_wave_deviation_ and i >= 0:
                    spectrum[i] += E
                    i -= 1

    return spectrum / ((right_boundary_x - left_boundary_x) * (up_boundary_y - down_boundary_y) * (front_boundary_z - rear_boundary_z))

# find full kinetic energy in part of volume
def find_full_energy(u, v, w):
    E_without_koeff = np.sum(np.square(u) + np.square(v) + np.square(w)) * dx * dy * dz_new
    return E_without_koeff / (2.0 * (right_boundary_x - left_boundary_x) * (up_boundary_y - down_boundary_y) * (front_boundary_z - rear_boundary_z))

# calculation of wave vectors arrays along coordinate axes
kx = np.arange(min_cells_number // 2)
ky = np.arange(min_cells_number // 2)
kz = np.arange(min_cells_number // 2)

kx_3d_array = np.multiply(np.ones((min_cells_number // 2, min_cells_number // 2, min_cells_number // 2)), kz)
ky_3d_array = np.multiply(np.ones((min_cells_number // 2, min_cells_number // 2, min_cells_number // 2)), ky.reshape((min_cells_number // 2, 1)))
kz_3d_array = np.multiply(np.ones((min_cells_number // 2, min_cells_number // 2, min_cells_number // 2)), kx.reshape((min_cells_number // 2, 1, 1)))
# 3d array with square of wave vector in every points
k_sqr_vector_array = np.square(kx_3d_array) + np.square(ky_3d_array) + np.square(kz_3d_array)

amplitude_wave_vec = np.linalg.norm(np.array([min_cells_number // 2 - 1, min_cells_number // 2 - 1, min_cells_number // 2 - 1]))

# real part of fourier transform spatial velocity vector
# function return three 3-D arrays with components of vector vk
def find_fourier_velocity_vec_from_file(file_number):
    x, y, z, p, u_center, v_center, w_center = get_data_from_file(file_number)
    #print("end reading " + str(file_number) + " time itteration")

    # fourier image of velocity components
    fourier_u = np.real(fftn(u_center))
    fourier_v = np.real(fftn(v_center))
    fourier_w = np.real(fftn(w_center))
    
    fourier_u = fourier_u[:min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    fourier_v = fourier_v[:min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    fourier_w = fourier_w[:min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    # vector velocity(k) = [fourier_u, fourier_v, fourier_w]
    return fourier_u / (math.pi**1.5), fourier_v / (math.pi**1.5), fourier_w / (math.pi**1.5)

def find_fourier_velocity_vec(u_center, v_center, w_center):
    fourier_u = np.real(fftn(u_center))
    fourier_v = np.real(fftn(v_center))
    fourier_w = np.real(fftn(w_center))
    
    fourier_u = fourier_u[:min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    fourier_v = fourier_v[:min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    fourier_w = fourier_w[:min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    # vector velocity(k) = [fourier_u, fourier_v, fourier_w]
    return fourier_u / (math.pi**1.5), fourier_v / (math.pi**1.5), fourier_w / (math.pi**1.5)

# temporary fourier transform of real part of fourier transform spatial velocity vector
def find_inertial_wave_from_files(file_number_array):
    u_center_arr = []
    v_center_arr = []
    w_center_arr = []
    for num in file_number_array:
        x, y, z, p, u_center, v_center, w_center = get_data_from_file(num)
        print("end reading " + str(num) + " file")
        u_center_arr.append(u_center)
        v_center_arr.append(v_center)
        w_center_arr.append(w_center)
    fourier_u = np.real(fftn(u_center_arr))
    fourier_v = np.real(fftn(v_center_arr))
    fourier_w = np.real(fftn(w_center_arr))

    l = len(file_number_array)
    time_period = l * time_step
    fourier_u = fourier_u[:l//2, :min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    fourier_v = fourier_v[:l//2, :min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    fourier_w = fourier_w[:l//2, :min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    return fourier_u / time_period, fourier_v / time_period, fourier_w / time_period

def find_inertial_wave_from_arrays(number_of_files, u_center_arr, v_center_arr, w_center_arr):
    fourier_u = np.real(fftn(u_center_arr))
    fourier_v = np.real(fftn(v_center_arr))
    fourier_w = np.real(fftn(w_center_arr))

    fourier_u = fourier_u[:number_of_files//2, :min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    fourier_v = fourier_v[:number_of_files//2, :min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    fourier_w = fourier_w[:number_of_files//2, :min_cells_number//2, :min_cells_number//2, :min_cells_number//2]
    time_period = number_of_files * time_step
    return fourier_u / time_period, fourier_v / time_period, fourier_w / time_period

def find_energy_spectrum_from_file(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz):
    fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec_from_file(file_number)
    # find E(q) in this file:
    E_i = find_energy_spectrum(fourier_u, fourier_v, fourier_w, min_cells_number // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
    q_max = (np.argmax(E_i)) * amplitude_wave_vec / count_of_vec_points
    full_e = find_full_energy(fourier_u, fourier_v, fourier_w)
    return E_i, q_max, full_e

# it's a copy of upper function which added returns: spatial fourier transform of velocity vectors 
def find_energy_spectrum_from_file_with_add_return(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz):
    fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec_from_file(file_number)
    # find E(q) in this file:
    E_i = find_energy_spectrum(fourier_u, fourier_v, fourier_w, min_cells_number // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
    q_max = (np.argmax(E_i)) * amplitude_wave_vec / count_of_vec_points
    full_e = find_full_energy(fourier_u, fourier_v, fourier_w)
    return E_i, q_max, full_e, fourier_u, fourier_v, fourier_w

# find the spectral density of the velocity of viscous kinetic energy dissipation per unit volume
# E(q) = summ(Vk^2 * k^2) / (2*pi)^3, if |k - q| <= 1/2;
# k = sqrt(kx*kx + ky*ky + kz*kz)
# Vk^2 = (u^2 + v^2 + w^2); u, v, w - projections of spectum of speed;
def spectral_density_of_velocity_of_dissipation(u, v, w, size_of_array, count_of_cells, freq_z, freq_y, freq_x, amplitude, max_wave_deviation_):
    spectrum = np.zeros(count_of_cells)
    for ind_z in range(0, size_of_array):
        for ind_y in range(0, size_of_array):
            for ind_x in range(0, size_of_array):
                k = np.linalg.norm(np.array([freq_z[ind_z], freq_y[ind_y], freq_x[ind_x]]))
                i = int((k + max_wave_deviation_) * count_of_cells / amplitude)
                E = np.sum(np.array([u[ind_z][ind_y][ind_x], v[ind_z][ind_y][ind_x], w[ind_z][ind_y][ind_x]]) ** 2) * k * k
                if i >= count_of_cells:
                    i = count_of_cells - 1
                while abs(k - i * amplitude / count_of_cells) <= max_wave_deviation_ and i >= 0:
                    spectrum[i] += E
                    i -= 1

    return spectrum / ((right_boundary_x - left_boundary_x) * (up_boundary_y - down_boundary_y) * (front_boundary_z - rear_boundary_z))

def find_full_rate_of_dissipation(u, v, w):
    v_sqr_vector_array = np.square(u) + np.square(v) + np.square(w)
    E_without_koeff = np.sum(v_sqr_vector_array * k_sqr_vector_array)
    return E_without_koeff * dx * dy * dz_new / ((right_boundary_x - left_boundary_x) * (up_boundary_y - down_boundary_y) * (front_boundary_z - rear_boundary_z))

def find_rate_of_viscous_kinetic_energy_dissipation_from_file(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz, max_wave_deviation_):
    fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec_from_file(file_number)
    spectrum = spectral_density_of_velocity_of_dissipation(fourier_u, fourier_v, fourier_w,
                                                           min_cells_number // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_)
    k_nu = (np.argmax(spectrum)) * amplitude_wave_vec / count_of_vec_points
    full_velocity = find_full_rate_of_dissipation(fourier_u, fourier_v, fourier_w)
    return spectrum, k_nu, full_velocity

# it's a copy of upper function which added arguments: spatial fourier transform of velocity vectors 
def find_rate_of_viscous_kinetic_energy_dissipation_from_file_with_add_arg(count_of_vec_points, amplitude_wave_vec, kx, ky, kz, max_wave_deviation_, fourier_u, fourier_v, fourier_w):
    spectrum = spectral_density_of_velocity_of_dissipation(fourier_u, fourier_v, fourier_w,
                                                           min_cells_number // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_)
    k_nu = (np.argmax(spectrum)) * amplitude_wave_vec / count_of_vec_points
    full_velocity = find_full_rate_of_dissipation(fourier_u, fourier_v, fourier_w)
    return spectrum, k_nu, full_velocity

# Re = sqrt(2E)/(nu*k_energy)
def find_posteriori_Reynolds_number(full_energy, energy_characteristic_wavelength, kinematic_viscosity):
    return math.sqrt(2.0 * full_energy) / (kinematic_viscosity * energy_characteristic_wavelength)

# Ro = ke * sqrt(E/2) / OMEGA 
def find_posteriori_Rossby_number(full_energy, energy_characteristic_wavelength, angular_rotation_speed):
    return energy_characteristic_wavelength * math.sqrt(0.5 * full_energy) / angular_rotation_speed

def spherical_coord_from_vec(x, y, z):
    xy = x*x + y*y
    r = np.sqrt(xy + z*z)
    theta = np.arctan2(np.sqrt(xy), z) # for elevation angle defined from Z-axis down
    phi = np.arctan2(y, x)  # for elevation angle defined from XY-plane up
    return r, theta, phi

def save_energy_spectrum(amplitude_wave_vec, max_wave_deviation, count_of_vec_steps, energy_spectrum_array, time):
    wave_coord = np.linspace(0, amplitude_wave_vec + max_wave_deviation, count_of_vec_steps)
    fig, axs = plt.subplots(figsize=(10, 10), dpi = 100)
    fig.suptitle('energy spectrum', fontsize=20)
    axs.set_title("time = " + str(time))
    axs.plot(np.log10(wave_coord), np.log10(energy_spectrum_array))
    axs.set_ylabel('log(E)')
    axs.set_xlabel('log(lengh of wave vector)')
    axs.grid(True, linestyle='--')
    plt.savefig('spectrum_time=' + str(time) + '.png')

def save_spectral_density_of_the_dissipation_rate(amplitude_wave_vec, max_wave_deviation, count_of_vec_steps, spectrum_array, time):
    wave_coord = np.linspace(0, amplitude_wave_vec + max_wave_deviation, count_of_vec_steps)
    fig, axs = plt.subplots(figsize=(10, 10), dpi = 100)
    fig.suptitle('Spectral density of the viscous dissipation rate', fontsize=20)
    axs.set_title("time = " + str(time))
    axs.plot(np.log10(wave_coord), np.log10(spectrum_array))
    axs.set_ylabel('log(epsilon)')
    axs.set_xlabel('log(lengh of wave vector)')
    axs.grid(True, linestyle='--')
    plt.savefig('spectrum_eps_time=' + str(time) + '.png')

# v (w, k) = (u, v, w)[w][kz, ky, kx]
# A_w_theta = summ (v(w, k)^2), where |theta_k - theta| < delta_theta / 2
def find_amplitude_of_inertial_waves_from_files(start_num, end_num, delta_theta, count_of_theta):
    u_wk, v_wk, w_wk = find_inertial_wave_from_files(np.arange(start_number, end_number+1))    # np.arange return numpy array [start, end)
    l = end_num - start_num + 1
    freq_ = 2.0*math.pi*fftfreq(l, time_step)[:l // 2]    # w 
    A_w_theta = []
    maximum_of_A = []
    theta_array = np.linspace(0.0, math.pi, num=count_of_theta)
    for omega in range(0, l // 2):
        A = np.zeros(count_of_theta)
        for ind_z in range(0, min_cells_number // 2):
            for ind_y in range(0, min_cells_number // 2):
                for ind_x in range(0, min_cells_number // 2):
                    mod_k, theta_k, phi = spherical_coord_from_vec(ind_x, ind_y, ind_z)
                    for theta_ind in range(0, count_of_theta):
                        if (abs(theta_array[theta_ind] - theta_k) < delta_theta / 2):
                            A[theta_ind] += np.sum(np.array([u_wk[omega][ind_z][ind_y][ind_x], v_wk[omega][ind_z][ind_y][ind_x], w_wk[omega][ind_z][ind_y][ind_x]]) ** 2)
        A_w_theta.append(A)
        maximum_of_A.append(np.max(A))
        print("end calculation amplitude in frequency " + str(freq_[omega]))
    A_w_theta = np.array(A_w_theta)
    maximum_of_A = np.array(maximum_of_A)
    # delete maximum in frequency = 0.0, because this maximum occurs due to the average flow velocity
    A_w_theta[0] = np.zeros(count_of_theta)
    A_w_theta[1] = np.zeros(count_of_theta)
    A_w_theta[2] = np.zeros(count_of_theta)
    maximum_of_A[0] = 0.0
    maximum_of_A[1] = 0.0
    maximum_of_A[2] = 0.0
    return freq_, A_w_theta, maximum_of_A, theta_array

# function like upper
def find_amplitude_of_inertial_waves_from_arrays(number_of_files, u_array, v_array, w_array, delta_theta, count_of_theta):
    u_wk, v_wk, w_wk = find_inertial_wave_from_arrays(number_of_files, u_array, v_array, w_array)
    freq_ = 2.0*math.pi*fftfreq(number_of_files, time_step)[:number_of_files // 2]    # frequency
    A_w_theta = []
    maximum_of_A = []
    theta_array = np.linspace(0.0, math.pi, num=count_of_theta)
    for omega in range(0, number_of_files // 2):
        A = np.zeros(count_of_theta)
        for ind_z in range(0, min_cells_number // 2):
            for ind_y in range(0, min_cells_number // 2):
                for ind_x in range(0, min_cells_number // 2):
                    mod_k, theta_k, phi = spherical_coord_from_vec(ind_x, ind_y, ind_z)
                    for theta_ind in range(0, count_of_theta):
                        if (abs(theta_array[theta_ind] - theta_k) < delta_theta / 2):
                            A[theta_ind] += np.sum(np.array([u_wk[omega][ind_z][ind_y][ind_x], v_wk[omega][ind_z][ind_y][ind_x], w_wk[omega][ind_z][ind_y][ind_x]]) ** 2)
        A_w_theta.append(A)
        maximum_of_A.append(np.max(A))
        print("end calculation amplitude in frequency " + str(freq_[omega]))
    A_w_theta = np.array(A_w_theta)
    maximum_of_A = np.array(maximum_of_A)
    # delete maximum in frequency = 0.0, because this maximum occurs due to the average flow velocity
    A_w_theta[0] = np.zeros(count_of_theta)
    A_w_theta[1] = np.zeros(count_of_theta)
    A_w_theta[2] = np.zeros(count_of_theta)
    maximum_of_A[0] = 0.0
    maximum_of_A[1] = 0.0
    maximum_of_A[2] = 0.0
    return freq_, A_w_theta, maximum_of_A, theta_array

# Function which find an energy spectrum and a rate of viscous kinetic energy dissipation and plot the graphs
# then find amplitude of inertial waves. Returns array with amplitude of inertial waves in relation to frequency and angle
def find_spectrums_and_ampl_of_inertial_waves_from_file(start_file_number, end_file_number, bool_will_save_spectrum, count_of_vec_steps, count_of_theta, delta_theta, bool_save_txt_info):
    if bool_save_txt_info == True:
        file = open(txt_file_name, "w+")
        file.close()

    u_array = []
    v_array = []
    w_array = []

    for i in np.arange(start_file_number, end_file_number+1, 1, dtype=int):
        x, y, z, p, u_center, v_center, w_center = get_data_from_file(i)
        fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec(u_center, v_center, w_center)
        energy_sp = find_energy_spectrum(fourier_u, fourier_v, fourier_w, min_cells_number // 2, count_of_vec_steps, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
        q_max = (np.argmax(energy_sp)) * amplitude_wave_vec / count_of_vec_steps
        #full_e = find_full_energy(fourier_u, fourier_v, fourier_w)
        full_e = find_full_energy(u_center, v_center, w_center)

        print("time moment: " + str(i * time_step))
        print("full kinetic energy = " + str(full_e))
        print("maximum spectrum of energy in wave vector = " + str(q_max))

        dencity_velocity_sp, maximum_nu, full_eps = find_rate_of_viscous_kinetic_energy_dissipation_from_file_with_add_arg(
            count_of_vec_steps, amplitude_wave_vec, kx, ky, kz, max_wave_deviation_in_dencity_velocity_dissipation_calc, fourier_u, fourier_v, fourier_w)
        Re = find_posteriori_Reynolds_number(full_e, q_max, kinematic_viscosity)
        Ro = find_posteriori_Rossby_number(full_e, q_max, OMEGA)

        print("maximum spectrum of density of viscous dissipation rate = " + str(maximum_nu))
        print("full velocity of viscous kinetic energy dissipation = " + str(full_eps))
        print("posteriori Reynolds number = " + str(Re))
        print("posteriori Rossby number = " + str(Ro))

        if bool_will_save_spectrum == True:
            save_energy_spectrum(amplitude_wave_vec, max_wave_deviation_in_energy_calc, count_of_vec_steps, energy_sp, i * time_step)
            save_spectral_density_of_the_dissipation_rate(amplitude_wave_vec, max_wave_deviation_in_dencity_velocity_dissipation_calc, count_of_vec_steps, dencity_velocity_sp, i * time_step)
        if bool_save_txt_info == True:
            file = open(txt_file_name, "a+")
            file.write("time moment: " + str(i * time_step) + "\n")
            file.write("full kinetic energy = " + str(full_e) + "\n")
            file.write("maximum spectrum of energy in wave vector = " + str(q_max) + "\n")
            file.write("maximum spectrum of density of viscous dissipation rate = " + str(maximum_nu) + "\n")
            file.write("full velocity of viscous kinetic energy dissipation = " + str(full_eps) + "\n")
            file.write("posteriori Reynolds number = " + str(Re) + "\n")
            file.write("posteriori Rossby number = " + str(Ro) + "\n")
            file.close()

        u_array.append(u_center)
        v_array.append(v_center)
        w_array.append(w_center)
    u_array = np.array(u_array)
    v_array = np.array(v_array)
    w_array = np.array(w_array)

    freq_, A_w_theta, maximum_of_A, theta_array = find_amplitude_of_inertial_waves_from_arrays(end_file_number - start_file_number + 1, u_array, v_array, w_array, delta_theta, count_of_theta)

    return freq_, A_w_theta, maximum_of_A, theta_array

# function which draw and save picture with inertial waves
A_w_theta_file_name = "Amplitude_(theta_frequency)_" + str(end_number - start_number + 1) + "_files.png"
def save_amplitude_of_inertial_waves(amplitude, angle_ticks_array, frequency_ticks_array, file_name):
    fig, ax = plt.subplots()
    fig.suptitle("Amplitude(frequency, angle)", fontsize=20)
    ax.imshow(amplitude)
    ax.set_xlabel("angle (theta), rad")
    ax.set_ylabel("frequency")
    ax.set_xticks(np.arange(len(angle_ticks_array)), labels=np.round(angle_ticks_array, decimals=1))
    ax.set_yticks(np.arange(len(frequency_ticks_array)), labels=np.round(frequency_ticks_array, decimals=1))
    plt.xticks(rotation=45)
    plt.savefig(file_name)

"""
fig, ax = plt.subplots()
ax.set_title("maximum of A(frequency, angle)")
ax.plot(np.arange(len(maximum_of_A)), maximum_of_A)
ax.set_xlabel("frequency")
ax.set_ylabel("Amplitude")
#ax.set_yticks(np.arange(count_of_theta), labels=theta_array)
ax.set_xticks(np.arange(len(freq_)), labels=np.round(freq_, decimals=1))
ax.grid(True, linestyle='--')
plt.xticks(rotation=45)
plt.show()"""

delta_theta = 0.15          # some choosen constant
count_of_theta = 70         # count of points which divided full angle range
count_of_vec_steps = 200    # wave vector discretization, wave vector discretization

freq_, A_w_theta, maximum_of_A, theta_array = find_spectrums_and_ampl_of_inertial_waves_from_file(
    start_number, end_number, False, count_of_vec_steps, count_of_theta, delta_theta, True)

save_amplitude_of_inertial_waves(A_w_theta, theta_array, freq_, A_w_theta_file_name)
