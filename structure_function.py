import numpy as np
import math
import matplotlib.pyplot as plt
from gc import collect
from scipy.optimize import curve_fit
from scipy.special import gamma
#import os
#import pandas as pd

path = "/home/sq-wm/OpenFOAM/sq-wm-v2306/run/calc18/"

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
kx, ky, kz = 1.5, 2.5, 4.0
Ax, Ay = 5.0, 3.0
k = math.sqrt(kx**2 + ky**2 + kz**2)
Lf = 1.0 / k
L = 2 * math.pi
V0 = math.sqrt(f0 / k)
gradient_coeff = 10.0

dl = (right_boundary_x - left_boundary_x) / x_cells
bound_dev = x_cells//3
max_dx = x_cells//5

size_of_text = 25
size_of_numbers = 25
size_of_legend = 25
marker_size = 110

def convert_time_name(x):
    if x%1.0 == 0.0:
        return round(x)
    else:
        return x

start_time, end_time, time_step = 7, 20, 0.5

time = np.round(np.arange(start_time, end_time + time_step, time_step), decimals=2)
time = list(map(convert_time_name, time))
print(time)

# -------------------------------------------------------------- isotropy

# uu_time, uv_time, uw_time = [], [], []
# vu_time, vv_time, vw_time = [], [], []
# wu_time, wv_time, ww_time = [], [], []

# # first line: [:, y_cells//2, z_cells//2]
# # second line: [:, y_cells//2, z_cells//4]
# # third line: [x_cells//2, y_cells//2, :]
# for i in time:
#     uu_, uv_, uw_ = [], [], []
#     vu_, vv_, vw_ = [], [], []
#     wu_, wv_, ww_ = [], [], []

#     u, v, w = get_data_from_file(path + str(i) + '/', x_cells, y_cells, z_cells)
#     for dx in range(1, max_dx):
#         uu, uv, uw = 0.0, 0.0, 0.0
#         vu, vv, vw = 0.0, 0.0, 0.0
#         wu, wv, ww = 0.0, 0.0, 0.0
#         for x in range(bound_dev, x_cells - bound_dev - dx):
#             uu += u[x_cells//2, y_cells//2, x+dx] * u[x_cells//2, y_cells//2, x]
#             uv += u[x_cells//2, y_cells//2, x+dx] * v[x_cells//2, y_cells//2, x]
#             uw += u[x_cells//2, y_cells//2, x+dx] * w[x_cells//2, y_cells//2, x]

#             vu += v[x_cells//2, y_cells//2, x+dx] * u[x_cells//2, y_cells//2, x]
#             vv += v[x_cells//2, y_cells//2, x+dx] * v[x_cells//2, y_cells//2, x]
#             vw += v[x_cells//2, y_cells//2, x+dx] * w[x_cells//2, y_cells//2, x]

#             wu += w[x_cells//2, y_cells//2, x+dx] * u[x_cells//2, y_cells//2, x]
#             wv += w[x_cells//2, y_cells//2, x+dx] * v[x_cells//2, y_cells//2, x]
#             ww += w[x_cells//2, y_cells//2, x+dx] * w[x_cells//2, y_cells//2, x]

#         uu_.append(uu / (x_cells - 2*bound_dev - dx))
#         uv_.append(uv / (x_cells - 2*bound_dev - dx))
#         uw_.append(uw / (x_cells - 2*bound_dev - dx))

#         vu_.append(vu / (x_cells - 2*bound_dev - dx))
#         vv_.append(vv / (x_cells - 2*bound_dev - dx))
#         vw_.append(vw / (x_cells - 2*bound_dev - dx))

#         wu_.append(wu / (x_cells - 2*bound_dev - dx))
#         wv_.append(wv / (x_cells - 2*bound_dev - dx))
#         ww_.append(ww / (x_cells - 2*bound_dev - dx))

#     uu_time.append(uu_)
#     uv_time.append(uv_)
#     uw_time.append(uw_)

#     vu_time.append(vu_)
#     vv_time.append(vv_)
#     vw_time.append(vw_)

#     wu_time.append(wu_)
#     wv_time.append(wv_)
#     ww_time.append(ww_)


# fig, axs = plt.subplots(figsize=(10, 10))

# axs.plot(np.arange(1, max_dx)*dl, np.mean(np.array(uu_time), axis=0), linewidth=3, label= r'$ \langle u \cdot u \rangle $')
# axs.plot(np.arange(1, max_dx)*dl, np.mean(np.array(uv_time), axis=0), linewidth=3, label= r'$ \langle u \cdot v \rangle $')
# axs.plot(np.arange(1, max_dx)*dl, np.mean(np.array(uw_time), axis=0), linewidth=3, label= r'$ \langle u \cdot w \rangle $')

# axs.plot(np.arange(1, max_dx)*dl, np.mean(np.array(vu_time), axis=0), linewidth=3, label= r'$ \langle v \cdot u \rangle $')
# axs.plot(np.arange(1, max_dx)*dl, np.mean(np.array(vv_time), axis=0), linewidth=3, label= r'$ \langle v \cdot v \rangle $')
# axs.plot(np.arange(1, max_dx)*dl, np.mean(np.array(vw_time), axis=0), linewidth=3, label= r'$ \langle v \cdot w \rangle $')

# axs.plot(np.arange(1, max_dx)*dl, np.mean(np.array(wu_time), axis=0), linewidth=3, label= r'$ \langle w \cdot u \rangle $')
# axs.plot(np.arange(1, max_dx)*dl, np.mean(np.array(wv_time), axis=0), linewidth=3, label= r'$ \langle w \cdot v \rangle $')
# axs.plot(np.arange(1, max_dx)*dl, np.mean(np.array(ww_time), axis=0), linewidth=3, label= r'$ \langle w \cdot w \rangle $')

# axs.grid(True, linestyle='--')
# axs.set_xlabel(r'$ dx $', fontsize=size_of_text)
# axs.set_ylabel(r'$ B(dx) = \langle U(x+dx) \cdot U(x) \rangle $', fontsize=size_of_text)
# axs.legend(fontsize=size_of_legend)
# axs.tick_params(labelsize=size_of_numbers)
# fig.tight_layout()
# plt.show()

# -------------------------------------------------------------- structure function

u_time = []
v_time = []
w_time = []

max_space_dev = 5
max_tau = len(time)//2
min_x, max_x = x_cells//2 - max_space_dev, x_cells//2 + max_space_dev
min_y, max_y = y_cells//2 - max_space_dev, y_cells//2 + max_space_dev
min_z, max_z = z_cells//2 - max_space_dev, z_cells//2 + max_space_dev

for i in range(len(time)):
    u, v, w = get_data_from_file(path + str(time[i]) + '/', x_cells, y_cells, z_cells)
    # u_time.append(u[position])
    # v_time.append(v[position])
    # w_time.append(w[position])

    u_time.append(u[min_z:max_z, min_y:max_y, min_x:max_x])
    v_time.append(v[min_z:max_z, min_y:max_y, min_x:max_x])
    w_time.append(w[min_z:max_z, min_y:max_y, min_x:max_x])

# structure_function = np.zeros(max_tau - 1)
# for cell in range(cells_count):
#     structure_function_u = []
#     structure_function_v = []
#     structure_function_w = []

#     for tau in range(1, len(time)//2):
#         D_u = 0.0
#         D_v = 0.0
#         D_w = 0.0
#         for t in range(0, len(time) - tau):
#             # D_u += (u_array[cell, t+tau] - u_array[cell, t])**2
#             # D_v += (v_array[cell, t+tau] - v_array[cell, t])**2
#             # D_w += (w_array[cell, t+tau] - w_array[cell, t])**2
#             D_u += (u_time[t+tau] - u_time[t])**2
#             D_v += (v_time[t+tau] - v_time[t])**2
#             D_w += (w_time[t+tau] - w_time[t])**2

#         structure_function_u.append(D_u / (len(time) - tau))
#         structure_function_v.append(D_v / (len(time) - tau))
#         structure_function_w.append(D_w / (len(time) - tau))

#     structure_function += np.array(np.array(structure_function_u) + np.array(structure_function_v) + np.array(structure_function_w), dtype=float)

# structure_function /= cells_count


structure_function_u = []
structure_function_v = []
structure_function_w = []

for tau in range(1, len(time)//2):
    D_u = np.zeros((max_z - min_z, max_y - min_y, max_x - min_x))
    D_v = np.zeros((max_z - min_z, max_y - min_y, max_x - min_x))
    D_w = np.zeros((max_z - min_z, max_y - min_y, max_x - min_x))
    for t in range(0, len(time) - tau):
        D_u += np.power(np.array(u_time[t+tau] - u_time[t]), 2.0)
        D_v += np.power(np.array(v_time[t+tau] - v_time[t]), 2.0)
        D_w += np.power(np.array(w_time[t+tau] - w_time[t]), 2.0)

    structure_function_u.append(D_u / (len(time) - tau))
    structure_function_v.append(D_v / (len(time) - tau))
    structure_function_w.append(D_w / (len(time) - tau))

structure_function = np.array(np.array(structure_function_u) + np.array(structure_function_v) + np.array(structure_function_w), dtype=float)
structure_function_mean = np.mean(np.mean(np.mean(structure_function, axis=1), axis=1), axis=1)

print(structure_function.shape, structure_function_mean.shape)

def power_function(x, A, g):
    return A * np.power(x, g)


def spectrum(omega, A, g):
    return (A * math.sin(0.5 * math.pi * g) * gamma(g + 1.0)) * np.power(omega, -g - 1.0) / math.pi


params, _ = curve_fit(power_function, np.arange(1, len(time)//2)*time_step, structure_function_mean)
print(params)

fig, axs = plt.subplots(figsize=(10, 10))
axs.plot(np.arange(1, len(time)//2)*time_step, structure_function_mean, linewidth=3, label= r'$ D (\tau) $')
axs.plot(np.linspace(time_step, len(time)//2*time_step, 100, endpoint=True), power_function(np.linspace(time_step, len(time)//2*time_step, 100, endpoint=True), *params), linewidth=3, label=str(round(params[0], 2)) + r'$ \, \tau $' + f"^ {params[1]:.2f}")

axs.grid(True, linestyle='--')
axs.set_xlabel(r'$ \tau $', fontsize=size_of_text)
axs.set_ylabel(r'$ D_U$ ', fontsize=size_of_text)
axs.legend(fontsize=size_of_legend)
axs.tick_params(labelsize=size_of_numbers)
fig.tight_layout()
plt.show()

fig, axs = plt.subplots(figsize=(10, 10))
axs.plot(np.linspace(0.1, 100, 500, endpoint=True), spectrum(np.linspace(0.1, 100, 500, endpoint=True), *params), linewidth=3, label=r'$ C / \, \omega $' + f"^ {params[1]:.2f}")

axs.grid(True, linestyle='--')
axs.set_xlabel(r'$ \omega $', fontsize=size_of_text)
axs.set_ylabel(r'$ E (\omega) $ ', fontsize=size_of_text)
axs.legend(fontsize=size_of_legend)
axs.tick_params(labelsize=size_of_numbers)
fig.tight_layout()
plt.show()
