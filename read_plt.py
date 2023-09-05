import matplotlib.pyplot as plt
import numpy as np
import struct
from scipy.fft import fftn, fftfreq
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
time_step = 0.1
#cells_number_x = 200
#cells_number_y = 200
#cells_number_z = 500
#time_step = 0.0004

dx = (right_boundary_x - left_boundary_x) / cells_number_x
dy = (up_boundary_y - down_boundary_y) / cells_number_y
dz = (front_boundary_z - rear_boundary_z) / cells_number_z

kinematic_viscosity = 0.01  # table value, a constant for a given environment. 

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
    for it in range(0, cells_number_x):
        p_new[it] = (9*p[4*it + 1] + 9*p[4*it + 2] - p[4*it + 0] - p[4*it + 3]) / 16
    # for 500 -> 200:
    #for it in range(0, cells_number_x, 2):
    #    p_new[it] = 415.0*p[5*it//2 + 0]/256.0 - 635.0*p[5*it//2 + 1]/128.0 + 141.0*p[5*it//2 + 2]/16.0 - 765.0*p[5*it//2 + 3]/128 + 385.0*p[5*it//2 + 4]/256.0
    #    p_new[it + 1] = 385.0*p[5*it//2 + 0]/256.0 - 765.0*p[5*it//2 + 1]/128.0 + 141.0*p[5*it//2 + 2]/16.0 - 635.0*p[5*it//2 + 3]/128 + 415.0*p[5*it//2 + 4]/256.0
    p = p_new

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                u[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    # for 400 -> 100:
    for it in range(0, cells_number_x):
        u_new[it] = (9*u[4*it + 1] + 9*u[4*it + 2] - u[4*it + 0] - u[4*it + 3]) / 16
    # for 500 -> 200:
    #for it in range(0, cells_number_x, 2):
    #    u_new[it] = 415.0*u[5*it//2 + 0]/256.0 - 635.0*u[5*it//2 + 1]/128.0 + 141.0*u[5*it//2 + 2]/16.0 - 765.0*u[5*it//2 + 3]/128 + 385.0*u[5*it//2 + 4]/256.0
    #    u_new[it + 1] = 385.0*u[5*it//2 + 0]/256.0 - 765.0*u[5*it//2 + 1]/128.0 + 141.0*u[5*it//2 + 2]/16.0 - 635.0*u[5*it//2 + 3]/128 + 415.0*u[5*it//2 + 4]/256.0
    u = u_new

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                v[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    # for 400 -> 100:
    for it in range(0, cells_number_x):
        v_new[it] = (9*v[4*it + 1] + 9*v[4*it + 2] - v[4*it + 0] - v[4*it + 3]) / 16
    # for 500 -> 200:
    #for it in range(0, cells_number_x, 2):
    #    v_new[it] = 415.0*v[5*it//2 + 0]/256.0 - 635.0*v[5*it//2 + 1]/128.0 + 141.0*v[5*it//2 + 2]/16.0 - 765.0*v[5*it//2 + 3]/128 + 385.0*v[5*it//2 + 4]/256.0
    #    v_new[it + 1] = 385.0*v[5*it//2 + 0]/256.0 - 765.0*v[5*it//2 + 1]/128.0 + 141.0*v[5*it//2 + 2]/16.0 - 635.0*v[5*it//2 + 3]/128 + 415.0*v[5*it//2 + 4]/256.0
    v = v_new

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                w[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])
    # for 400 -> 100:
    for it in range(0, cells_number_x):
        w_new[it] = (9*w[4*it + 1] + 9*w[4*it + 2] - w[4*it + 0] - w[4*it + 3]) / 16
    # for 500 -> 200:
    #for it in range(0, cells_number_x, 2):
    #    w_new[it] = 415.0*w[5*it//2 + 0]/256.0 - 635.0*w[5*it//2 + 1]/128.0 + 141.0*w[5*it//2 + 2]/16.0 - 765.0*w[5*it//2 + 3]/128 + 385.0*w[5*it//2 + 4]/256.0
    #    w_new[it + 1] = 385.0*w[5*it//2 + 0]/256.0 - 765.0*w[5*it//2 + 1]/128.0 + 141.0*w[5*it//2 + 2]/16.0 - 635.0*w[5*it//2 + 3]/128 + 415.0*w[5*it//2 + 4]/256.0
    w = w_new

    mistake = file.read()
    if len(mistake) > 0:
        print("remained " + str(len(mistake)) + " bytes in the end of the file")
        print(mistake)
        print("something doesn't work")
    file.close()
    return x, y, z, p, u, v, w
#cells_number_z_new = cells_number_z // 4
cells_number_z_new = cells_number_x

# chose time interval
start_number = 220  # numbers of file
end_number = 234

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
# k = sqrt(kx*kx + ky*ky + kz*kz)
# Vk^2 = (u^2 + v^2 + w^2); u, v, w - projections of spectum of speed; 
# q - array of modulus of vectors 
max_wave_deviation_in_energy_calc = 0.25     # constants in this calculation
max_wave_deviation_in_dencity_velocity_dissipation_calc = 0.5
def find_energy_spectrum(u, v, w, size_z, size_y, size_x, count_of_cells, freq_z, freq_y, freq_x, amplitude, max_wave_deviation_):
    spectrum = np.zeros(count_of_cells)
    for ind_z in range(0, size_z):
        for ind_y in range(0, size_y):
            for ind_x in range(0, size_x):
                k = np.linalg.norm(np.array([freq_z[ind_z], freq_y[ind_y], freq_x[ind_x]]))
                i = int((k + max_wave_deviation_) * count_of_cells / amplitude)
                E = np.sum(np.array([u[ind_z][ind_y][ind_x], v[ind_z][ind_y][ind_x], w[ind_z][ind_y][ind_x]]) ** 2)
                if i >= count_of_cells:
                    i = count_of_cells - 1
                while abs(k - i * amplitude / count_of_cells) <= max_wave_deviation_ and i >= 0:
                    spectrum[i] += E
                    i -= 1

    return spectrum / (8.0 * math.pi * math.pi * math.pi)

# find full kinetic energy in part of volume
def find_full_energy(u, v, w, size_z, size_y, size_x):
    E = 0.0
    for ind_z in range(0, size_z):
        for ind_y in range(0, size_y):
            for ind_x in range(0, size_x):
                E += np.sum(np.array([u[ind_z][ind_y][ind_x], v[ind_z][ind_y][ind_x], w[ind_z][ind_y][ind_x]]) ** 2)

    return E / (16.0 * math.pi * math.pi * math.pi)

# calculation of wave vectors arrays along coordinate axes
kx = np.arange(cells_number_x // 2)
ky = np.arange(cells_number_y // 2)
kz = np.arange(cells_number_z_new // 2)

amplitude_wave_vec = np.linalg.norm(np.array([cells_number_z_new // 2 - 1, cells_number_y // 2 - 1, cells_number_x // 2 - 1]))

# real part of fourier transform spatial velocity vector
def find_fourier_velocity_vec(file_number):
    x, y, z, p, u_center, v_center, w_center = get_data_from_file(file_number)
    #print("end reading " + str(file_number) + " time itteration")

    # fourier image of velocity components
    fourier_u = np.real(fftn(u_center))
    fourier_v = np.real(fftn(v_center))
    fourier_w = np.real(fftn(w_center))
    
    fourier_u = fourier_u[:cells_number_z_new // 2][:cells_number_y // 2][:cells_number_x // 2]
    fourier_v = fourier_v[:cells_number_z_new // 2][:cells_number_y // 2][:cells_number_x // 2]
    fourier_w = fourier_w[:cells_number_z_new // 2][:cells_number_y // 2][:cells_number_x // 2]
    # vector velocity(k) = [fourier_u, fourier_v, fourier_w]
    return fourier_u, fourier_v, fourier_w

# temporary fourier transform of real part of fourier transform spatial velocity vector
def find_inertial_wave(file_number_array):
    u_center_arr = []
    v_center_arr = []
    w_center_arr = []
    for num in file_number_array:
        x, y, z, p, u_center, v_center, w_center = get_data_from_file(num)
        u_center_arr.append(u_center)
        v_center_arr.append(v_center)
        w_center_arr.append(w_center)
    fourier_u = np.real(fftn(u_center_arr))
    fourier_v = np.real(fftn(v_center_arr))
    fourier_w = np.real(fftn(w_center_arr))

    l = len(file_number_array)
    fourier_u = fourier_u[:l // 2][:cells_number_z_new // 2][:cells_number_y // 2][:cells_number_x // 2]
    fourier_v = fourier_v[:l // 2][:cells_number_z_new // 2][:cells_number_y // 2][:cells_number_x // 2]
    fourier_w = fourier_w[:l // 2][:cells_number_z_new // 2][:cells_number_y // 2][:cells_number_x // 2]
    return fourier_u, fourier_v, fourier_w

# maximum search coefficient - part of the array where programm find local maximum. It is searched visually according to the schedule
max_search_coef = 0.375 
def find_energy_spectrum_from_file(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz):
    fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec(file_number)
    # find E(q) in this file:
    E_i = find_energy_spectrum(fourier_u, fourier_v, fourier_w, cells_number_z_new // 2, cells_number_y // 2, cells_number_x // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
    max_energy = np.max(E_i[int(len(E_i)*max_search_coef):-1])
    q_max = (np.argmax(E_i[int(len(E_i)*max_search_coef) :]) + int(len(E_i)*max_search_coef)) * amplitude_wave_vec / count_of_vec_points
    full_e = find_full_energy(fourier_u, fourier_v, fourier_w, cells_number_z_new // 2, cells_number_y // 2, cells_number_x // 2)
    return E_i, max_energy, q_max, full_e

# it's a copy of upper function which added returns: spatial fourier transform of velocity vectors 
def find_energy_spectrum_from_file_with_add_return(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz):
    fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec(file_number)
    # find E(q) in this file:
    E_i = find_energy_spectrum(fourier_u, fourier_v, fourier_w, cells_number_z_new // 2, cells_number_y // 2, cells_number_x // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
    max_energy = np.max(E_i[int(len(E_i)*max_search_coef):-1])
    q_max = (np.argmax(E_i[int(len(E_i)*max_search_coef) :]) + int(len(E_i)*max_search_coef)) * amplitude_wave_vec / count_of_vec_points
    full_e = find_full_energy(fourier_u, fourier_v, fourier_w, cells_number_z_new // 2, cells_number_y // 2, cells_number_x // 2)
    return E_i, max_energy, q_max, full_e, fourier_u, fourier_v, fourier_w

# find the spectral density of the velocity of viscous kinetic energy dissipation per unit volume
# E(q) = summ(Vk^2 * k^2) / (2*pi)^3, if |k - q| <= 1/2;
# k = sqrt(kx*kx + ky*ky + kz*kz)
# Vk^2 = (u^2 + v^2 + w^2); u, v, w - projections of spectum of speed;
def spectral_density_of_velocity_of_dissipation(u, v, w, size_z, size_y, size_x, count_of_cells, freq_x, freq_y, freq_z, amplitude, max_wave_deviation_):
    spectrum = np.zeros(count_of_cells)
    for ind_z in range(0, size_z):
        for ind_y in range(0, size_y):
            for ind_x in range(0, size_x):
                k = np.linalg.norm(np.array([freq_z[ind_z], freq_y[ind_y], freq_x[ind_x]]))
                i = int((k + max_wave_deviation_) * count_of_cells / amplitude)
                E = np.sum(np.array([u[ind_z][ind_y][ind_x], v[ind_z][ind_y][ind_x], w[ind_z][ind_y][ind_x]]) ** 2) * k * k
                if i >= count_of_cells:
                    i = count_of_cells - 1
                while abs(k - i * amplitude / count_of_cells) <= max_wave_deviation_ and i >= 0:
                    spectrum[i] += E
                    i -= 1

    return spectrum / (8.0 * math.pi * math.pi * math.pi)

def find_full_rate_of_dissipation(u, v, w, size_z, size_y, size_x, freq_z, freq_y, freq_x):
    E = 0.0
    for ind_z in range(0, size_z):
        for ind_y in range(0, size_y):
            for ind_x in range(0, size_x):
                k = np.linalg.norm(np.array([freq_z[ind_z], freq_y[ind_y], freq_x[ind_x]]))
                E += np.sum(np.array([u[ind_z][ind_y][ind_x], v[ind_z][ind_y][ind_x], w[ind_z][ind_y][ind_x]]) ** 2) * k * k

    return E / (8.0 * math.pi * math.pi * math.pi)

def find_rate_of_viscous_kinetic_energy_dissipation_from_file(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz, max_wave_deviation_):
    fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec(file_number)
    spectrum = spectral_density_of_velocity_of_dissipation(fourier_u, fourier_v, fourier_w,
                                                           cells_number_z_new // 2, cells_number_y // 2, cells_number_x // 2, 
                                                           count_of_vec_points, kx, ky, kz, amplitude_wave_vec, max_wave_deviation_)
    k_nu = (np.argmax(spectrum[int(len(spectrum)*max_search_coef) :]) + int(len(spectrum)*max_search_coef)) * amplitude_wave_vec / count_of_vec_points
    full_velocity = find_full_rate_of_dissipation(fourier_u, fourier_v, fourier_w, cells_number_z_new // 2, cells_number_y // 2, cells_number_x // 2, kz, ky, kx)
    return spectrum, k_nu, full_velocity

# it's a copy of upper function which added arguments: spatial fourier transform of velocity vectors 
def find_rate_of_viscous_kinetic_energy_dissipation_from_file_with_add_arg(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz, max_wave_deviation_, fourier_u, fourier_v, fourier_w):
    spectrum = spectral_density_of_velocity_of_dissipation(fourier_u, fourier_v, fourier_w,
                                                           cells_number_z_new // 2, cells_number_y // 2, cells_number_x // 2, 
                                                           count_of_vec_points, kx, ky, kz, amplitude_wave_vec, max_wave_deviation_)
    k_nu = (np.argmax(spectrum[int(len(spectrum)*max_search_coef) :]) + int(len(spectrum)*max_search_coef)) * amplitude_wave_vec / count_of_vec_points
    full_velocity = find_full_rate_of_dissipation(fourier_u, fourier_v, fourier_w, cells_number_z_new // 2, cells_number_y // 2, cells_number_x // 2, kz, ky, kx)
    return spectrum, k_nu, full_velocity

# Re = sqrt(2E)/(nu*k_energy)
def find_posteriori_Reynolds_number(full_energy, energy_characteristic_wavelength):
    return math.sqrt(2.0 * full_energy) / (kinematic_viscosity * energy_characteristic_wavelength)

def spherical_coord_from_vec(x, y, z):
    xy = x**2 + y**2
    r = np.sqrt(xy + z**2)
    theta = np.arctan2(np.sqrt(xy), z) # for elevation angle defined from Z-axis down
    phi = np.arctan2(y, x)  # for elevation angle defined from XY-plane up
    return r, theta, phi

count_of_vec_steps = 200
def save_energy_spectrum(amplitude_wave_vec, max_wave_deviation, count_of_vec_steps, energy_spectrum_array, time):
    wave_coord = np.linspace(0, amplitude_wave_vec + max_wave_deviation, count_of_vec_steps)
    fig, axs = plt.subplots(figsize=(10, 10), dpi = 100)
    fig.suptitle('energy spectrum', fontsize=20)
    axs.set_title("time = " + str(time))
    axs.plot(np.log10(wave_coord), np.log10(energy_spectrum_array))
    axs.set_ylabel('log(E)')
    axs.set_xlabel('log(lengh of wave vector)')
    axs.grid(True)
    plt.savefig('spectrum_time=' + str(time) + '.png')

def save_spectral_density_of_the_dissipation_rate(amplitude_wave_vec, max_wave_deviation, count_of_vec_steps, spectrum_array, time):
    wave_coord = np.linspace(0, amplitude_wave_vec + max_wave_deviation, count_of_vec_steps)
    fig, axs = plt.subplots(figsize=(10, 10), dpi = 100)
    fig.suptitle('Spectral density of the viscous dissipation rate', fontsize=20)
    axs.set_title("time = " + str(time))
    axs.plot(np.log10(wave_coord), np.log10(spectrum_array))
    axs.set_ylabel('log(epsilon)')
    axs.set_xlabel('log(lengh of wave vector)')
    axs.grid(True)
    plt.savefig('spectrum_eps_time=' + str(time) + '.png')

# v (w, k) = (u, v, w)[w][kz, ky, kx]
# A_w_theta = summ (v(w, k)^2), where |theta_k - theta| < delta_theta / 2
def find_amplitude_of_inertial_waves(start_num, end_num, delta_theta, count_of_theta):
    u_wk, v_wk, w_wk = find_inertial_wave(np.arange(start_number, end_number+1))    # np.arange return numpy array [start, end)
    l = end_num - start_num + 1
    freq_ = 2.0*math.pi*fftfreq(l, time_step)[:l // 2]    # w 
    A_w_theta = []
    maximum_of_A = []
    theta_array = np.linspace(0.0, math.pi, num=count_of_theta)
    for omega in range(0, l // 2):
        A = np.zeros(count_of_theta)
        for ind_z in range(0, cells_number_z_new // 2):
            for ind_y in range(0, cells_number_y // 2):
                for ind_x in range(0, cells_number_x // 2):
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
    return freq_, A_w_theta, maximum_of_A

# cycle according to the time we are interested in
for i in range(0, end_number - start_number + 1, 1):    
    energy_sp, maximum_E, Q_max, full_E, fourier_u, fourier_v, fourier_w = find_energy_spectrum_from_file_with_add_return(
        i + start_number, count_of_vec_steps, amplitude_wave_vec, kx, ky, kz)
    print("time moment: " + str((i + start_number) * 0.001))
    print("full kinetic energy = " + str(full_E))
    print("maximum spectrum of energy = " + str(maximum_E) + " in wave vector = " + str(Q_max))

    dencity_velocity_sp, maximum_nu, full_eps = find_rate_of_viscous_kinetic_energy_dissipation_from_file_with_add_arg(
        i + start_number, count_of_vec_steps, amplitude_wave_vec, kx, ky, kz, max_wave_deviation_in_dencity_velocity_dissipation_calc, fourier_u, fourier_v, fourier_w)
    Re = find_posteriori_Reynolds_number(full_E, Q_max)
    print("maximum spectrum of density of viscous dissipation rate = " + str(maximum_nu))
    print("full velocity of viscous kinetic energy dissipation = " + str(full_eps))
    print("posteriori Reynolds number = " + str(Re))
    save_energy_spectrum(amplitude_wave_vec, max_wave_deviation_in_energy_calc, count_of_vec_steps, energy_sp, (i + start_number) * 0.001)
    save_spectral_density_of_the_dissipation_rate(amplitude_wave_vec, max_wave_deviation_in_dencity_velocity_dissipation_calc, count_of_vec_steps, dencity_velocity_sp, (i + start_number) * 0.001)

delta_theta = 0.15      # some choosen constant
count_of_theta = 70     # count of points which divided full angle range
freq_, A_w_theta, maximum_of_A = find_amplitude_of_inertial_waves(start_number, end_number, delta_theta, count_of_theta)
theta_array = np.linspace(0.0, math.pi, num=count_of_theta)

# draw inertial waves
fig, ax = plt.subplots()
ax.set_title("A(frequency, angle)")
ax.imshow(A_w_theta)
ax.set_xlabel("angle (theta), rad")
ax.set_ylabel("frequency")
ax.set_xticks(np.arange(count_of_theta), labels=np.round(theta_array, decimals=1))
ax.set_yticks(np.arange(len(freq_)), labels=np.round(freq_, decimals=1))
plt.xticks(rotation=45)
plt.show()

fig, ax = plt.subplots()
ax.set_title("maximum of A(frequency, angle)")
ax.plot(np.arange(len(maximum_of_A)), maximum_of_A)
ax.set_xlabel("frequency")
ax.set_ylabel("Amplitude")
#ax.set_yticks(np.arange(count_of_theta), labels=theta_array)
ax.set_xticks(np.arange(len(freq_)), labels=np.round(freq_, decimals=1))
ax.grid(True, linestyle='--')
plt.xticks(rotation=45)
plt.show()
