import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fftn, fftfreq
import math
#from scipy.interpolate import griddata
#from fluidfoam import readmesh, readvector#, readscalar
from gc import collect

# info from parameters.dat file:
left_boundary_x =   -1.0*math.pi
right_boundary_x =  math.pi
down_boundary_y =   -1.0*math.pi
up_boundary_y =     math.pi
rear_boundary_z =   -1.0*math.pi
front_boundary_z =  math.pi
cells_number_axis = 600

kinematic_viscosity = 0.01  # table value, a constant for a given environment. 
OMEGA = 50.0                # angular rotation speed - parameters of the simulated system. Info from parameters.dat file
max_wave_deviation_in_energy_calc = 0.25     # constants in this calculation
max_wave_deviation_in_dencity_velocity_dissipation_calc = 0.5
count_of_boundary_cells = int(round((cells_number_axis * np.sqrt(kinematic_viscosity / OMEGA) / (front_boundary_z - rear_boundary_z) + 0.5), 0)) * 2

dx = (right_boundary_x - left_boundary_x) / cells_number_axis
dy = (up_boundary_y - down_boundary_y) / cells_number_axis
dz = (front_boundary_z - rear_boundary_z) / cells_number_axis

# name of txt file, which save information about Reynolds and Rossby number, characteristic wave vector, energy and dissipation
txt_file_name = "system_parameters.txt"
txt_file_name_table = "Re_Ro_E_table.txt"
spectrum_file_name = "energy_spectrum.txt"

#path = "/storage/p1/ermakov/calc12/"
path = "/home/sq-wm/OpenFOAM/sq-wm-v2306/run/calc13/"
# read size of new structured array
with open(path + 'grid_size.txt', 'r') as fp:
    count_of_cells_axis = fp.read().split('\n')
    x_cells = int(count_of_cells_axis[0])
    y_cells = int(count_of_cells_axis[1])
    z_cells = int(count_of_cells_axis[2])
    print('structured grid has ' + str(x_cells) + ' * ' + str(y_cells) + ' * ' + str(z_cells) + ' cells')
cells_number_axis = min([x_cells, y_cells, z_cells])

start_time = 19
end_time = 20
time_step = 0.1
time_array = np.round(np.linspace(start_time, end_time, num=int((end_time - start_time)/time_step + 1), endpoint=True), decimals=1)
print('time points being processed:')
print(time_array)

def create_structured_mesh(cells_number_axis):
    xi = np.linspace(left_boundary_x, right_boundary_x, cells_number_axis)
    yi = np.linspace(down_boundary_y, rear_boundary_z, cells_number_axis)
    zi = np.linspace(rear_boundary_z, front_boundary_z, cells_number_axis)
    xinterp, yinterp, zinterp = np.meshgrid(xi, yi, zi)
    return xinterp, yinterp, zinterp

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

power_coeff = 10.0
proportion_of_inner_size = 0.95
# for simmetric 3D array
def make_structured_grid_with_zero_boudary(array, x_boundaries, y_boundaries, z_boundaries, coefficient, proportion):
    cells_z = len(array)
    cells_y = len(array[0])
    cells_x = len(array[0, 0])
    array_returned = np.copy(array)
    for z_ind in range(0, cells_z):
        for y_ind in range(0, cells_y):
            for x_ind in range(0, cells_x):
                x = (2.0 * x_ind / (cells_x - 1) - 1.0) * x_boundaries
                y = (2.0 * y_ind / (cells_y - 1) - 1.0) * y_boundaries
                z = (2.0 * z_ind / (cells_z - 1) - 1.0) * z_boundaries
                x_multiplier = np.exp(- coefficient * (x - proportion*x_boundaries)) / (1.0 + np.exp(- coefficient * (x + proportion*x_boundaries))) / (1.0 + np.exp(- coefficient * (x - proportion*x_boundaries)))
                y_multiplier = np.exp(- coefficient * (y - proportion*y_boundaries)) / (1.0 + np.exp(- coefficient * (y + proportion*y_boundaries))) / (1.0 + np.exp(- coefficient * (y - proportion*y_boundaries)))
                z_multiplier = np.exp(- coefficient * (z - proportion*z_boundaries)) / (1.0 + np.exp(- coefficient * (z + proportion*z_boundaries))) / (1.0 + np.exp(- coefficient * (z - proportion*z_boundaries)))
                array_returned[z_ind, y_ind, x_ind] *= (x_multiplier * y_multiplier * z_multiplier)
    print("boundary layer was deleted")
    return array_returned

# find the energy spectrum at the choosen moment of time
# E(q) = summ(Vk^2) / (2*pi)^3, if |k - q| <= 1/4;
# k = sqrt(kx*kx + ky*ky + kz*kz)
# Vk^2 = (u^2 + v^2 + w^2); u, v, w - projections of spectum of speed; 
# q - array of modulus of vectors; can be more speedy with numba!
def find_energy_spectrum(u, v, w, size_of_array, count_of_cells, freq_z, freq_y, freq_x, amplitude, max_wave_deviation_):
    spectrum = np.zeros(count_of_cells)
    k_array = np.linspace(math.sqrt(3.0)/2.0, amplitude, num=count_of_cells)
    for ind_z in range(0, size_of_array):
        for ind_y in range(0, size_of_array):
            for ind_x in range(0, size_of_array):
                k = np.linalg.norm(np.array([freq_z[ind_z], freq_y[ind_y], freq_x[ind_x]]))
                #i = int((k + max_wave_deviation_) * (count_of_cells - 1) / amplitude)
                E = np.sum(np.square(np.array([u[ind_z, ind_y, ind_x], v[ind_z, ind_y, ind_x], w[ind_z, ind_y, ind_x]])))
                #if i >= count_of_cells:
                #    i = count_of_cells - 1
                #while abs(k - i * amplitude / count_of_cells) <= max_wave_deviation_ and i >= 0:
                #    spectrum[i] += E
                #    i -= 1
                for i in range(0, count_of_cells):
                    if abs(k_array[i] - k) <= max_wave_deviation_:
                        spectrum[i] += E

    maximum_ind = np.where(spectrum == np.max(spectrum))[0][0]

    return spectrum / ((right_boundary_x - left_boundary_x) * (up_boundary_y - down_boundary_y) * (front_boundary_z - rear_boundary_z)), k_array[maximum_ind]

# find full kinetic energy in part of volume
def find_full_energy(u, v, w):
    E_without_koeff = np.sum(np.square(u) + np.square(v) + np.square(w)) * dx * dy * dz
    return E_without_koeff / (2.0 * (right_boundary_x - left_boundary_x) * (up_boundary_y - down_boundary_y) * (front_boundary_z - rear_boundary_z))

# calculation of wave vectors arrays along coordinate axes
kx1 = np.arange(0.0, (x_cells // 2) / 2.0, 0.5, dtype=float) + 0.5
kx = np.concatenate([kx1, np.flip(-1.0*kx1, 0)])
ky1 = np.arange(0.0, (y_cells // 2) / 2.0, 0.5, dtype=float) + 0.5
ky = np.concatenate([ky1, np.flip(-1.0*ky1, 0)])
kz1 = np.arange(0.0, (z_cells // 2) / 2.0, 0.5, dtype=float) + 0.5
kz = np.concatenate([kz1, np.flip(-1.0*kz1, 0)])
#kx1 = np.arange(0.0, (x_cells // 2 - count_of_boundary_cells) / 2.0, 0.5, dtype=float) + 0.5
#kx = np.concatenate([kx1, np.flip(-1.0*kx1, 0)])
#ky1 = np.arange(0.0, (y_cells // 2 - count_of_boundary_cells) / 2.0, 0.5, dtype=float) + 0.5
#ky = np.concatenate([ky1, np.flip(-1.0*ky1, 0)])
#kz1 = np.arange(0.0, (z_cells // 2 - count_of_boundary_cells) / 2.0, 0.5, dtype=float) + 0.5
#kz = np.concatenate([kz1, np.flip(-1.0*kz1, 0)])

kx_3d_array = np.multiply(np.ones((z_cells, y_cells, x_cells)), kz)
ky_3d_array = np.multiply(np.ones((z_cells, y_cells, x_cells)), ky.reshape((y_cells, 1)))
kz_3d_array = np.multiply(np.ones((z_cells, y_cells, x_cells)), kx.reshape((x_cells, 1, 1)))
#kx_3d_array = np.multiply(np.ones((z_cells - 2*count_of_boundary_cells, y_cells - 2*count_of_boundary_cells, x_cells - 2*count_of_boundary_cells)), kz)
#ky_3d_array = np.multiply(np.ones((z_cells - 2*count_of_boundary_cells, y_cells - 2*count_of_boundary_cells, x_cells - 2*count_of_boundary_cells)), ky.reshape((y_cells - 2*count_of_boundary_cells, 1)))
#kz_3d_array = np.multiply(np.ones((z_cells - 2*count_of_boundary_cells, y_cells - 2*count_of_boundary_cells, x_cells - 2*count_of_boundary_cells)), kx.reshape((z_cells - 2*count_of_boundary_cells, 1, 1)))
# 3d array with square of wave vector in every points
k_sqr_vector_array = np.square(kx_3d_array) + np.square(ky_3d_array) + np.square(kz_3d_array)

#amplitude_wave_vec = np.linalg.norm(np.array([cells_number_axis // 2 - 1, cells_number_axis // 2 - 1, cells_number_axis // 2 - 1]))
#amplitude_wave_vec = np.linalg.norm(np.array([(cells_number_axis // 2 - 1) / 2.0, (cells_number_axis // 2 - 1) / 2.0, (cells_number_axis // 2 - 1) / 2.0]))
amplitude_wave_vec = np.linalg.norm(np.array([z_cells / 4.0, y_cells / 4.0, x_cells / 4.0]))

# real part of fourier transform spatial velocity vector
# function return three 3-D arrays with components of vector vk
def find_fourier_velocity_vec_from_file(file_number):
    u_center, v_center, w_center = get_data_from_file(file_number)
    #print("end reading " + str(file_number) + " time itteration")

    # fourier image of velocity components
    fourier_u = np.real(fftn(u_center))
    fourier_v = np.real(fftn(v_center))
    fourier_w = np.real(fftn(w_center))
    
    fourier_u = fourier_u[:z_cells//2, :y_cells//2, :x_cells//2]
    fourier_v = fourier_v[:z_cells//2, :y_cells//2, :x_cells//2]
    fourier_w = fourier_w[:z_cells//2, :y_cells//2, :x_cells//2]
    # vector velocity(k) = [fourier_u, fourier_v, fourier_w]
    return fourier_u / (math.pi**1.5), fourier_v / (math.pi**1.5), fourier_w / (math.pi**1.5)

def find_fourier_velocity_vec(u_center, v_center, w_center):
    fourier_u = np.abs(fftn(u_center))
    fourier_v = np.abs(fftn(v_center))
    fourier_w = np.abs(fftn(w_center))
    
    #fourier_u = fourier_u[ : cells_number_axis//2, : cells_number_axis//2, : cells_number_axis//2]
    #fourier_v = fourier_v[ : cells_number_axis//2, : cells_number_axis//2, : cells_number_axis//2]
    #fourier_w = fourier_w[ : cells_number_axis//2, : cells_number_axis//2, : cells_number_axis//2]
    # vector velocity(k) = [fourier_u, fourier_v, fourier_w]
    return fourier_u / (math.pi**1.5), fourier_v / (math.pi**1.5), fourier_w / (math.pi**1.5)

# temporary fourier transform of real part of fourier transform spatial velocity vector
def find_inertial_wave_from_files(file_number_array):
    u_center_arr = []
    v_center_arr = []
    w_center_arr = []
    for num in file_number_array:
        u_center, v_center, w_center = get_data_from_file(num)
        print("end reading " + str(num) + " file")
        u_center_arr.append(u_center)
        v_center_arr.append(v_center)
        w_center_arr.append(w_center)
    fourier_u = np.real(fftn(u_center_arr))
    fourier_v = np.real(fftn(v_center_arr))
    fourier_w = np.real(fftn(w_center_arr))

    l = len(file_number_array)
    time_period = l * time_step
    fourier_u = fourier_u[: l//2, : z_cells//2, : y_cells//2, : x_cells//2]
    fourier_v = fourier_v[: l//2, : z_cells//2, : y_cells//2, : x_cells//2]
    fourier_w = fourier_w[: l//2, : z_cells//2, : y_cells//2, : x_cells//2]
    return fourier_u / time_period, fourier_v / time_period, fourier_w / time_period

def find_inertial_wave_from_arrays(number_of_files, u_center_arr, v_center_arr, w_center_arr):
    fourier_u = np.real(fftn(u_center_arr))
    fourier_v = np.real(fftn(v_center_arr))
    fourier_w = np.real(fftn(w_center_arr))

    fourier_u = fourier_u[: number_of_files//2, : z_cells//2, : y_cells//2, : x_cells//2]
    fourier_v = fourier_v[: number_of_files//2, : z_cells//2, : y_cells//2, : x_cells//2]
    fourier_w = fourier_w[: number_of_files//2, : z_cells//2, : y_cells//2, : x_cells//2]
    time_period = number_of_files * time_step
    return fourier_u / time_period, fourier_v / time_period, fourier_w / time_period

def find_energy_spectrum_from_file(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz):
    fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec_from_file(file_number)
    # find E(q) in this file:
    E_i, q_max = find_energy_spectrum(fourier_u, fourier_v, fourier_w, cells_number_axis // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
    full_e = find_full_energy(fourier_u, fourier_v, fourier_w)
    return E_i, q_max, full_e

# it's a copy of upper function which added returns: spatial fourier transform of velocity vectors 
def find_energy_spectrum_from_file_with_add_return(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz):
    fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec_from_file(file_number)
    # find E(q) in this file:
    E_i, q_max = find_energy_spectrum(fourier_u, fourier_v, fourier_w, cells_number_axis // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
    full_e = find_full_energy(fourier_u, fourier_v, fourier_w)
    return E_i, q_max, full_e, fourier_u, fourier_v, fourier_w

# find the spectral density of the velocity of viscous kinetic energy dissipation per unit volume
# E(q) = summ(Vk^2 * k^2) / (2*pi)^3, if |k - q| <= 1/2;
# k = sqrt(kx*kx + ky*ky + kz*kz)
# Vk^2 = (u^2 + v^2 + w^2); u, v, w - projections of spectum of speed
def spectral_density_of_velocity_of_dissipation(u, v, w, size_of_array, count_of_cells, freq_z, freq_y, freq_x, amplitude, max_wave_deviation_):
    spectrum = np.zeros(count_of_cells)
    k_array = np.linspace(math.sqrt(3.0)/2.0, amplitude, num=count_of_cells)
    for ind_z in range(0, size_of_array):
        for ind_y in range(0, size_of_array):
            for ind_x in range(0, size_of_array):
                k2 = np.sum(np.square(np.array([freq_z[ind_z], freq_y[ind_y], freq_x[ind_x]])))
                E = np.sum(np.square(np.array([u[ind_z, ind_y, ind_x], v[ind_z, ind_y, ind_x], w[ind_z, ind_y, ind_x]]))) * k2
                for i in range(0, count_of_cells):
                    if abs(k_array[i] - math.sqrt(k2)) <= max_wave_deviation_:
                        spectrum[i] += E

    maximum_ind = np.where(spectrum == np.max(spectrum))[0][0]

    return spectrum / ((right_boundary_x - left_boundary_x) * (up_boundary_y - down_boundary_y) * (front_boundary_z - rear_boundary_z)), k_array[maximum_ind]

def find_full_rate_of_dissipation(u, v, w):
    v_sqr_vector_array = np.square(u) + np.square(v) + np.square(w)
    E_without_koeff = np.sum(v_sqr_vector_array * k_sqr_vector_array)
    return E_without_koeff * dx * dy * dz / ((right_boundary_x - left_boundary_x) * (up_boundary_y - down_boundary_y) * (front_boundary_z - rear_boundary_z))

def find_rate_of_viscous_kinetic_energy_dissipation_from_file(file_number, count_of_vec_points, amplitude_wave_vec, kx, ky, kz, max_wave_deviation_):
    fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec_from_file(file_number)
    spectrum, k_nu = spectral_density_of_velocity_of_dissipation(fourier_u, fourier_v, fourier_w,
                                                           cells_number_axis // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_)
    full_velocity = find_full_rate_of_dissipation(fourier_u, fourier_v, fourier_w)
    return spectrum, k_nu, full_velocity

# it's a copy of upper function which added arguments: spatial fourier transform of velocity vectors 
def find_rate_of_viscous_kinetic_energy_dissipation_from_file_with_add_arg(count_of_vec_points, amplitude_wave_vec, kx, ky, kz, max_wave_deviation_, fourier_u, fourier_v, fourier_w):
    spectrum, k_nu = spectral_density_of_velocity_of_dissipation(fourier_u, fourier_v, fourier_w,
                                                           cells_number_axis // 2, count_of_vec_points, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_)
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

def save_energy_spectrum(amplitude_wave_vec, count_of_vec_steps, energy_spectrum_array, time):
    wave_coord = np.linspace(0, amplitude_wave_vec, count_of_vec_steps)
    for i in range(0, count_of_vec_steps):
        if energy_spectrum_array[i] < 1.0:
            energy_spectrum_array[i] = 1.0

    fig, axs = plt.subplots(figsize=(10, 10), dpi = 100)
    fig.suptitle('energy spectrum', fontsize=20)
    axs.set_title("time = " + str(time))
    #plt.axes(xlim=(0, amplitude_wave_vec), ylim=(1.0, 1.1*np.max(energy_spectrum_array)))
    plt.xscale("log")
    plt.yscale("log")
    axs.plot(wave_coord, energy_spectrum_array)
    axs.set_ylabel('E')
    axs.set_xlabel('lengh of wave vector')
    axs.grid(True, linestyle='--')
    plt.savefig(path + 'spectrum_time=' + str(time) + '.png')

def save_spectral_density_of_the_dissipation_rate(amplitude_wave_vec, count_of_vec_steps, spectrum_array, time):
    wave_coord = np.linspace(0, amplitude_wave_vec, count_of_vec_steps)
    for i in range(0, count_of_vec_steps):
        if spectrum_array[i] < 1.0:
            spectrum_array[i] = 1.0

    fig, axs = plt.subplots(figsize=(10, 10), dpi = 100)
    fig.suptitle('Spectral density of the viscous dissipation rate', fontsize=20)
    axs.set_title("time = " + str(time))
    plt.xscale("log")
    plt.yscale("log")
    axs.plot(wave_coord, spectrum_array)
    axs.set_ylabel('epsilon')
    axs.set_xlabel('lengh of wave vector')
    axs.grid(True, linestyle='--')
    plt.savefig(path + 'spectrum_eps_time=' + str(time) + '.png')

# v (w, k) = (u, v, w)[w][kz, ky, kx]
# A_w_theta = summ (v(w, k)^2), where |theta_k - theta| < delta_theta / 2
def find_amplitude_of_inertial_waves_from_files(start_num, end_num, delta_theta, count_of_theta):
    u_wk, v_wk, w_wk = find_inertial_wave_from_files(np.arange(start_time, end_time + time_step, time_step))    # np.arange return numpy array [start, end, step)
    l = end_num - start_num + 1
    freq_ = 2.0*math.pi*fftfreq(l, time_step)[:l // 2]    # w 
    A_w_theta = []
    maximum_of_A = []
    theta_array = np.linspace(0.0, math.pi, num=count_of_theta)
    for omega in range(0, l // 2):
        A = np.zeros(count_of_theta)
        for ind_z in range(0, cells_number_axis // 2):
            for ind_y in range(0, cells_number_axis // 2):
                for ind_x in range(0, cells_number_axis // 2):
                    mod_k, theta_k, phi = spherical_coord_from_vec(ind_x, ind_y, ind_z)
                    for theta_ind in range(0, count_of_theta):
                        if (abs(theta_array[theta_ind] - theta_k) < delta_theta / 2):
                            A[theta_ind] += np.sum(np.array([u_wk[omega, ind_z, ind_y, ind_x], v_wk[omega, ind_z, ind_y, ind_x], w_wk[omega, ind_z, ind_y, ind_x]]) ** 2)
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
        for ind_z in range(0, cells_number_axis // 2):
            for ind_y in range(0, cells_number_axis // 2):
                for ind_x in range(0, cells_number_axis // 2):
                    mod_k, theta_k, phi = spherical_coord_from_vec(ind_x, ind_y, ind_z)
                    for theta_ind in range(0, count_of_theta):
                        if (abs(theta_array[theta_ind] - theta_k) < delta_theta / 2):
                            A[theta_ind] += np.sum(np.array([u_wk[omega, ind_z, ind_y, ind_x], v_wk[omega, ind_z, ind_y, ind_x], w_wk[omega, ind_z, ind_y, ind_x]]) ** 2)
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
def find_spectrums_and_ampl_of_inertial_waves_from_file(start_time, end_time, time_step, bool_save_spectrum, count_of_vec_steps, count_of_theta, delta_theta, bool_save_txt_info):
    if bool_save_txt_info == True:
        file = open(path + txt_file_name, "w+")
        file.close()

    u_array = []
    v_array = []
    w_array = []

    for i in np.arange(start_time, end_time + time_step, time_step):
        u_center, v_center, w_center = get_data_from_file(i)
        fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec(u_center, v_center, w_center)
        energy_sp, q_max = find_energy_spectrum(fourier_u, fourier_v, fourier_w, cells_number_axis // 2, count_of_vec_steps, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
        full_e = find_full_energy(u_center, v_center, w_center)

        print("time moment: " + str(start_time + i * time_step))
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
        u_2_mean = math.sqrt(np.mean(np.square(u_center) + np.square(v_center) + np.square(w_center)))
        print("U RMS = " + str(u_2_mean))

        if bool_save_spectrum == True:
            save_energy_spectrum(amplitude_wave_vec, count_of_vec_steps, energy_sp, i)
            save_spectral_density_of_the_dissipation_rate(amplitude_wave_vec, count_of_vec_steps, dencity_velocity_sp, start_time + i*time_step)
        if bool_save_txt_info == True:
            file = open(path + txt_file_name, "a+")
            file.write("time moment: " + str(np.arange(start_time, end_time + time_step, time_step)) + "\n")
            file.write("full kinetic energy = " + str(full_e) + "\n")
            file.write("maximum spectrum of energy in wave vector = " + str(q_max) + "\n")
            file.write("maximum spectrum of density of viscous dissipation rate = " + str(maximum_nu) + "\n")
            file.write("full velocity of viscous kinetic energy dissipation = " + str(full_eps) + "\n")
            file.write("U RMS = " + str(u_2_mean) + "\n")
            file.write("posteriori Reynolds number = " + str(Re) + "\n")
            file.write("posteriori Rossby number = " + str(Ro) + "\n" + "\n")
            file.close()

        u_array.append(u_center)
        v_array.append(v_center)
        w_array.append(w_center)
    u_array = np.array(u_array)
    v_array = np.array(v_array)
    w_array = np.array(w_array)

    freq_, A_w_theta, maximum_of_A, theta_array = find_amplitude_of_inertial_waves_from_arrays((end_time - start_time)//time_step, u_array, v_array, w_array, delta_theta, count_of_theta)

    return freq_, A_w_theta, maximum_of_A, theta_array

def find_spectrums_and_characteristic_numbers_from_files(start_time, end_time, time_step, count_of_vec_steps, bool_save_txt_info, bool_save_spectrum):
    if bool_save_txt_info == True:
        file = open(path + txt_file_name, "w+")
        file.close()

        file_table = open(path + txt_file_name_table, "w+")
        file_table.close()

        file_e_spectrum = open(path + spectrum_file_name, "w+")
        file_e_spectrum.close()

    for time in time_array:
        if time % 1 == 0:
            time = int(time)
        print('time = ' + str(time))
        u_center, v_center, w_center = get_data_from_file(path + str(time) + '/', x_cells, y_cells, z_cells)
        #u_center_new = make_structured_grid_with_zero_boudary(u_center, right_boundary_x, up_boundary_y, front_boundary_z, power_coeff, proportion_of_inner_size)
        #v_center_new = make_structured_grid_with_zero_boudary(v_center, right_boundary_x, up_boundary_y, front_boundary_z, power_coeff, proportion_of_inner_size)
        #w_center_new = make_structured_grid_with_zero_boudary(w_center, right_boundary_x, up_boundary_y, front_boundary_z, power_coeff, proportion_of_inner_size)
        #u_center_new = np.copy(u_center[count_of_boundary_cells : z_cells - count_of_boundary_cells, count_of_boundary_cells : y_cells - count_of_boundary_cells, count_of_boundary_cells : x_cells - count_of_boundary_cells])
        #v_center_new = np.copy(v_center[count_of_boundary_cells : z_cells - count_of_boundary_cells, count_of_boundary_cells : y_cells - count_of_boundary_cells, count_of_boundary_cells : x_cells - count_of_boundary_cells])
        #w_center_new = np.copy(w_center[count_of_boundary_cells : z_cells - count_of_boundary_cells, count_of_boundary_cells : y_cells - count_of_boundary_cells, count_of_boundary_cells : x_cells - count_of_boundary_cells])

        #del u_center
        #del v_center
        #del w_center
        #collect()
        fourier_u, fourier_v, fourier_w = find_fourier_velocity_vec(u_center, v_center, w_center)
        
        #energy_sp, q_max  = find_energy_spectrum(fourier_u, fourier_v, fourier_w, cells_number_axis - 2*count_of_boundary_cells, count_of_vec_steps, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
        energy_sp, q_max  = find_energy_spectrum(fourier_u, fourier_v, fourier_w, cells_number_axis, count_of_vec_steps, kz, ky, kx, amplitude_wave_vec, max_wave_deviation_in_energy_calc)
        # delete maximum near wave vector = 0.0, because this maximum occurs due to the average flow velocity
        energy_sp_copy = np.copy(energy_sp)
        energy_sp_copy[0] = 0.0
        energy_sp_copy[1] = 0.0
        full_e = find_full_energy(u_center, v_center, w_center)

        print("time moment: " + str(time))
        print("full kinetic energy = " + str(full_e))
        print("maximum spectrum of energy in wave vector = " + str(q_max))

        dencity_velocity_sp, maximum_nu, full_eps = find_rate_of_viscous_kinetic_energy_dissipation_from_file_with_add_arg(
            count_of_vec_steps, amplitude_wave_vec, kx, ky, kz, max_wave_deviation_in_dencity_velocity_dissipation_calc, fourier_u, fourier_v, fourier_w)
        Re = find_posteriori_Reynolds_number(full_e, q_max, kinematic_viscosity)
        Ro = find_posteriori_Rossby_number(full_e, q_max, OMEGA)

        del energy_sp_copy
        collect()

        print("maximum spectrum of density of viscous dissipation rate = " + str(maximum_nu))
        print("full velocity of viscous kinetic energy dissipation = " + str(full_eps))
        print("posteriori Reynolds number = " + str(Re))
        print("posteriori Rossby number = " + str(Ro))
        u_2_mean = math.sqrt(np.mean(np.square(u_center) + np.square(v_center) + np.square(w_center)))
        print("U RMS = " + str(u_2_mean))

        if bool_save_txt_info == True:
            file = open(path + txt_file_name, "a+")
            file.write("time moment: " + str(time) + "\n")
            file.write("full kinetic energy = " + str(full_e) + "\n")
            file.write("maximum spectrum of energy in wave vector = " + str(q_max) + "\n")
            file.write("maximum spectrum of density of viscous dissipation rate = " + str(maximum_nu) + "\n")
            file.write("full velocity of viscous kinetic energy dissipation = " + str(full_eps) + "\n")
            file.write("posteriori Reynolds number = " + str(Re) + "\n")
            file.write("U RMS = " + str(u_2_mean) + "\n")
            file.write("posteriori Rossby number = " + str(Ro) + "\n" + "\n")
            file.close()

            file_table = open(path + txt_file_name_table, "a+")
            file_table.write(str(time) + "\t")
            file_table.write(str(Re) + "\t")
            file_table.write(str(Ro) + "\t")
            file_table.write(str(u_2_mean) + "\t")
            file_table.write(str(full_e) + "\n")
            file_table.close()
        
        if bool_save_spectrum == True:
            file_e_spectrum = open(path + spectrum_file_name, "a+")
            for ind in range(0, len(energy_sp)):
                file_e_spectrum.write(str(energy_sp[ind]) + " ")
            file_e_spectrum.close()
            save_energy_spectrum(amplitude_wave_vec, count_of_vec_steps, energy_sp, time)
            save_spectral_density_of_the_dissipation_rate(amplitude_wave_vec, count_of_vec_steps, dencity_velocity_sp, time)

        del u_center
        del v_center
        del w_center
        del fourier_u
        del fourier_v
        del fourier_w
        del energy_sp
        del dencity_velocity_sp
        collect()

    return 'end calculation characteristic numbers in files from ' + str(start_time) + ' to ' + str(end_time)

# function which draw and save picture with inertial waves
A_w_theta_file_name = "Amplitude_(theta_frequency)_" + str(end_time - start_time + time_step) + "_files.png"
def save_amplitude_of_inertial_waves(amplitude, angle_ticks_array, frequency_ticks_array, file_name):
    len_x = len(angle_ticks_array)
    len_y = len(frequency_ticks_array)
    max_len = max(len_x, len_y)
    len_x_ticks = 10    # count of ticks along axis
    len_y_ticks = 10
    fig, ax = plt.subplots(figsize=(20 * len_x / max_len, 20 * len_y / max_len))    # max size of picture = 20 * 20; all pixels are square
    fig.suptitle("Amplitude(frequency, angle)", fontsize=15)
    plt.xticks(np.linspace(0, len_x, num=len_x_ticks), labels=np.round(np.linspace(np.min(angle_ticks_array), np.max(angle_ticks_array), num=len_x_ticks), decimals=1))
    plt.yticks(np.linspace(0, len_y, num=len_y_ticks), labels=np.round(np.linspace(np.min(frequency_ticks_array), np.max(frequency_ticks_array), num=len_y_ticks), decimals=1))
    ax.set_xlabel("angle (theta), rad")
    ax.set_ylabel("frequency")
    ax.imshow(amplitude)
    plt.savefig(file_name)

delta_theta = 0.15          # some choosen constant
count_of_theta = 70         # count of points which divided full angle range
count_of_vec_steps = 200    # wave vector discretization, wave vector discretization

#freq_, A_w_theta, maximum_of_A, theta_array = find_spectrums_and_ampl_of_inertial_waves_from_file(
#    start_time, end_time, time_step, False, count_of_vec_steps, count_of_theta, delta_theta, True)

#save_amplitude_of_inertial_waves(A_w_theta, theta_array, freq_, A_w_theta_file_name)
ans = find_spectrums_and_characteristic_numbers_from_files(start_time, end_time, time_step, count_of_vec_steps, True, True)
print(ans)
