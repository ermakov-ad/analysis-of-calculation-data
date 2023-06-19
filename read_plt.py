import matplotlib.pyplot as plt
import numpy as np
import struct
from scipy.fft import fft, fftfreq, ifftn

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
# u, v, w = np.array[z][y][x]; size = cells_z * cells_y * cells_x
def get_data_from_file(file_number):
    file = open(construct_file_path(file_number), 'rb')
    # skip header
    header = file.read(284)
    #print("header text:")
    #print(header)
    integers = []
    for i in range(0, 28, 1):
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

    if len(file.read()) > 0:
        print("remained " + str(len(file.read())) + " bytes in the end of the file")
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
start_time = 0.0
end_time = 0.002

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
    print("end " + str(i) + " time itteration")

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

print("zero: ")
print(u_center[0][int(cells_number_z/2)][0:5][0:5])
print(v_center[0][int(cells_number_z/2)][0:5][0:5])
print(w_center[0][int(cells_number_z/2)][0:5][0:5])
print("one: ")
print(u_center[1][int(cells_number_z/2)][0:5][0:5])
print(v_center[1][int(cells_number_z/2)][0:5][0:5])
print(w_center[1][int(cells_number_z/2)][0:5][0:5])
print("two: ")
print(u_center[2][int(cells_number_z/2)][0:5][0:5])
print(v_center[2][int(cells_number_z/2)][0:5][0:5])
print(w_center[2][int(cells_number_z/2)][0:5][0:5])

after_f = ifftn(u_center[2])
print(len(after_f[0][0]))
print(after_f[0:3][0:3][0:3])
plt.imshow(np.real(after_f[0:-1][0:-1][40]))
plt.grid()
plt.show()

""" picture = []
for i in range(0, 4):
    picture.append(np.empty([cells_number_x, cells_number_y]))
picture = np.array(picture)

for i in range(0, 4):
    for j in range(0, cells_number_y):
        for k in range(0, cells_number_x):
            picture[i][k][j] = p[i][ind_data + k + j*cells_number_x]

fig = plt.figure()

plt.subplot(2, 2, 1)
plt.imshow(picture[0])

plt.subplot(2, 2, 2)
plt.imshow(picture[1])

plt.subplot(2, 2, 3)
plt.imshow(picture[2])

plt.subplot(2, 2, 4)
plt.imshow(picture[3])

plt.show() """
