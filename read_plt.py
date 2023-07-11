import matplotlib.pyplot as plt
import numpy as np
import struct
from scipy.fft import fftn, fftfreq, fftshift
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

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                u[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                v[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])

    for z_ind in range(0, cells_number_z):
        for y_ind in range(0, cells_number_y):
            for x_ind in range(0, cells_number_x):
                value = file.read(8)
                w[z_ind][y_ind][x_ind] = float(struct.unpack('<d', value)[0])

    mistake = file.read()
    if len(mistake) > 0:
        print("remained " + str(len(mistake)) + " bytes in the end of the file")
        print(mistake)
        print("something doesn't work")
    file.close()
    return x, y, z, p, u, v, w

x = []
y = []
z = []
p = []
u_center = []
v_center = []
w_center = []
# chose time interval
# start_time = 0.0
# end_time = 0.003
start_time = 0.231
end_time = 0.234

print(int(start_time / time_step))
print(int(end_time / time_step + 1))

for i in range(int(start_time / time_step), int(end_time / time_step + 1), 1):
    ans = get_data_from_file(i)
    x.append(ans[0])
    y.append(ans[1])
    z.append(ans[2])
    p.append(ans[3])
    u_center.append(ans[4])
    v_center.append(ans[5])
    w_center.append(ans[6])
    print("end reading " + str(i) + " time itteration")

# summary arrays with primitive variables for further processing
x = np.array(x)
y = np.array(y)
z = np.array(z)
p = np.array(p)
u_center = np.array(u_center)
v_center = np.array(v_center)
w_center = np.array(w_center)

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
max_fr_deviation = 0.25     # constant in this calculation
def find_energy_spectrum(u, v, w, size_z, size_y, size_x, count_of_cells, freq_z, freq_y, freq_x, amplitude):
    spectrum = np.zeros(count_of_cells)
    for ind_z in range(0, size_z):
        for ind_y in range(0, size_y):
            for ind_x in range(0, size_x):
                k = modulus_of_vector([freq_z[ind_z], freq_y[ind_y], freq_x[ind_x]], 3)
                i = int((k + max_fr_deviation) * count_of_cells / amplitude)
                E = sqr_of_vector([u[ind_z][ind_y][ind_x], v[ind_z][ind_y][ind_x], w[ind_z][ind_y][ind_x]], 3)
                if (abs(k - i * amplitude / count_of_cells) >= max_fr_deviation):
                    i -= 1
                if (i >= count_of_cells):
                    i = count_of_cells - 1
                while (abs(k - i * amplitude / count_of_cells) < max_fr_deviation and i >= 0):
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

# calculation of frequency arrays along coordinate axes
kx = fftfreq(cells_number_x * 2, dx)
ky = fftfreq(cells_number_y * 2, dy)
kz = fftfreq(cells_number_z * 2, dz)
delta_fr = (kx[1] - kx[0]) / 2.0
kx += delta_fr
delta_fr = (ky[1] - ky[0]) / 2.0
ky += delta_fr
delta_fr = (kz[1] - kz[0]) / 2.0
kz += delta_fr

kx = kx[0:len(kx) // 2]
ky = ky[0:len(ky) // 2]
kz = kz[0:len(kz) // 2]

energy_spectrum = []
# cycle according to the time we are interested in
for i in range(0, int(end_time / time_step + 1) - int(start_time / time_step), 1):
    # fourier image of velocity components
    fourier_u = np.real(fftn(u_center[i]))
    fourier_v = np.real(fftn(v_center[i]))
    fourier_w = np.real(fftn(w_center[i]))

    # find E(q):
    amplitude_freq = max(max(kz), max(ky), max(kx))
    count_of_fr_points = 50
    
    E_i = find_energy_spectrum(fourier_u, fourier_v, fourier_w, cells_number_z, cells_number_y, cells_number_x, count_of_fr_points, kz, ky, kx, amplitude_freq)
        
    energy_spectrum.append(E_i)
    max_energy = np.max(E_i)
    q_max = np.argmax(E_i) * amplitude_freq / count_of_fr_points

    print("time moment: " + str((i + int(start_time / time_step)) * time_step))
    print("full kinetic energy = ")
    print(find_full_energy(fourier_u, fourier_v, fourier_w, cells_number_z, cells_number_y, cells_number_x))
    print("maximum spectrum of energy = " + str(max_energy) + " in wave vector = " + str(q_max))

# draw pressure
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

axs[0, 0].set_title("p, time = " + str((0 + int(start_time / time_step)) * time_step) + " , z = 0")
axs[0, 0].imshow(p[0][200][:-1][:-1])

axs[0, 1].set_title("p, time = " + str((1 + int(start_time / time_step)) * time_step) + " , z = 0")
axs[0, 1].imshow(p[1][200][:-1][:-1])

axs[1, 0].set_title("p, time = " + str((2 + int(start_time / time_step)) * time_step) + " , z = 0")
axs[1, 0].imshow(p[2][200][:-1][:-1])

axs[1, 1].set_title("p, time = " + str((3 + int(start_time / time_step)) * time_step) + " , z = 0")
axs[1, 1].imshow(p[3][200][:-1][:-1])

plt.show()

# draw energy spectrum
energy_spectrum = np.array(energy_spectrum)
freq_coord = np.linspace(0, amplitude_freq + max_fr_deviation, 50)

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
fig.suptitle('energy spectrum', fontsize=14)

axs[0, 0].set_title("time = " + str((0 + int(start_time / time_step)) * time_step))
axs[0, 0].plot(np.log10(freq_coord), np.log10(energy_spectrum[0]))
axs[0, 0].set_ylabel('log(E)')
axs[0, 0].set_xlabel('wave vector')

axs[0, 1].set_title("time = " + str((1 + int(start_time / time_step)) * time_step))
axs[0, 1].plot(np.log10(freq_coord), np.log10(energy_spectrum[1]))
axs[0, 1].set_ylabel('log(E)')
axs[0, 1].set_xlabel('wave vector')

axs[1, 0].set_title("time = " + str((2 + int(start_time / time_step)) * time_step))
axs[1, 0].plot(np.log10(freq_coord), np.log10(energy_spectrum[2]))
axs[1, 0].set_ylabel('log(E)')
axs[1, 0].set_xlabel('wave vector')

axs[1, 1].set_title("time = " + str((3 + int(start_time / time_step)) * time_step))
axs[1, 1].plot(np.log10(freq_coord), np.log10(energy_spectrum[3]))
axs[1, 1].set_ylabel('log(E)')
axs[1, 1].set_xlabel('wave vector')

plt.show()
