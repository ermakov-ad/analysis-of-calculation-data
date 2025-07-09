import numpy as np
import math
import matplotlib.pyplot as plt
from gc import collect
from scipy.optimize import curve_fit
import os
import pandas as pd

size_of_text = 25
size_of_numbers = 25
size_of_legend = 25
marker_size = 110

full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)
print(filename)
print(path)

# df = pd.read_excel(path + '/выход_на_стац_энергия.xlsx')

# time = np.array(df['0'])
# chaostic_energy = np.array(df['0.2'])       # == <u^2 (t)>
# print(chaostic_energy)

# epsilon = nu * <(du/dy + dv/dx)> /2
# epsilon(t) * t =? <u^2 (t)>
nu = 0.0001
rho = 1000.0
file_count = 400

bound_dev = 5
max_dx = 250
cells = 512
min_stacionary_t = 370
max_t = file_count+1
t_count = (max_t - min_stacionary_t)
time_step = 0.25


def read_vtk(path_to_file, step):  #сначала y, потом x
    file = open(path_to_file + '/P_E_V_Vec_' + str(step) + '.vtk')
    data = file.read().split('\n')
    file.close()
    cells = data[4].split(' ')
    cells_x, cells_y, cells_z = int(cells[1]) - 1, int(cells[2]) - 1, int(cells[3])
    if cells_z > 1:
        cells_z -= 1
    strings_count = cells_x * cells_y * cells_z

    u, v, w, vorticity, p = [], [], [], [], []
    position = 14
    for ind_y in range(cells_y):
        line_x, line_y, line_z =[], [], []
        for ind_x in range(cells_x):
            vel = data[position + ind_y*cells_x + ind_x].split(' ')
            line_x.append(vel[0])
            line_y.append(vel[1])
            line_z.append(vel[2])
        u.append(line_x)
        v.append(line_y)
        w.append(line_z)
    del line_x, line_y, line_z, vel
    collect()

    position += strings_count + 4

    for ind_y in range(cells_y):
        line_p = []
        for ind_x in range(cells_x):
            line_p.append(data[position + ind_y*cells_x + ind_x])
        p.append(line_p)
    del line_p
    collect()
    
    position += strings_count + 4

    for ind_y in range(cells_y):
        line_vort = []
        for ind_x in range(cells_x):
            line_vort.append(data[position + ind_y*cells_x + ind_x])
        vorticity.append(line_vort)
    del line_vort
    collect()
    return np.array(u, dtype=float), np.array(v, dtype=float), np.array(w, dtype=float), np.array(vorticity, dtype=float), np.array(p, dtype=float), cells_x


# max_file = 50
# et = np.empty(max_file)
# time_ = np.arange(max_file, dtype=int)
# for i in range(time_):
#     u, v, w, vort, p = read_vtk(path, time_[i])
#     et[i] = 0.5*nu*np.sum(np.gradient(u, axis=1, order=2) + np.gradient(v, axis=0, order=2))        #?
# # сравнить et и chaostic_energy

def kinetic_energy(u, v, w):
    return 0.5 * rho * np.sum(np.square(u) + np.square(v) + np.square(w))


# -------------------------------------- energy spectrum --------------------------------------------------
# E = []
# for i in range(0, file_count+1):
#     u, v, w, vort, p, cells = read_vtk(path + '/G_0_05_mu_0_1_Alpha_0_03_periodic_0-400/', i)
#     E.append(kinetic_energy(u, v, w))

# fig, axs = plt.subplots(figsize=(10, 10))
# axs.plot(np.arange(file_count+1)*0.25, E, linewidth=3, label=r'$ E = \rho U^{2} / 2 $')

# axs.grid(True, linestyle='--')
# axs.set_xlabel(r'$ t, с $', fontsize=size_of_text)
# axs.set_ylabel(r'$ E (t) $ ', fontsize=size_of_text)
# axs.legend(fontsize=size_of_legend)
# axs.tick_params(labelsize=size_of_numbers)
# fig.tight_layout()
# plt.show()


# -------------------------------------- vorticity animation (gif) --------------------------------------------------
# import matplotlib.animation as animation
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# vort_z_slice_time_array = []
# for i in range(file_count+1):
#     u, v, w, vort, p, cells = read_vtk(path + '/G_0_05_mu_0_1_Alpha_0_03_periodic_0-400/', i)
#     vort_z_slice_time_array.append(vort)

# vort_z_slice_time_array = np.array(vort_z_slice_time_array)
# fig, ax = plt.subplots(dpi=500)
# div = make_axes_locatable(ax)
# cax = div.append_axes('right', '5%', '5%')
# ax.set_title('vorticity z')      #tx =         , time = ' + str(time_array[0]))

# ims = []
# for i in range(file_count+1):
#     im = ax.imshow(vort_z_slice_time_array[i], animated=True, vmin=-5.0, vmax=5.0)      # , vmin=-120.0, vmax=120.0, cmap='seismic'
#     if i == 0:
#         ax.imshow(vort_z_slice_time_array[i], vmin=-5.0, vmax=5.0)  # show an initial one first , vmin=-120.0, vmax=120.0, cmap='seismic'
#     #tx.set_text('vorticity z mean, time = {0}'.format(time_array[i]))
#     fig.colorbar(im, cax=cax)
#     ims.append([im])

# ani = animation.ArtistAnimation(fig, ims, interval=30, blit=True, repeat_delay=100)
# ani.save(path + '/2d_processing/' + 'vort_z_chaostic.gif' , fps=3)          
# #plt.show()


# -------------------------------------- isotropy along line --------------------------------------------------
# uu_time, uv_time = [], []
# vu_time, vv_time = [], []

# # first line: [y_cells//2, : ]
# # second line: [y_cells//4, : ]
# # third line: [ : , x_cells//2]
# # fourth line: [ : , x_cells//4]

# for i in range(min_stacionary_t, max_t):
#     uu_, uv_ = [], []
#     vu_, vv_ = [], []

#     u, v, w, vort, p, cells = read_vtk(path + '/G_0_05_mu_0_1_Alpha_0_03_periodic_0-400/', i)

#     for dx in range(0, max_dx):
#         uu, uv = 0.0, 0.0
#         vu, vv = 0.0, 0.0
#         for x in range(bound_dev, cells - bound_dev - dx):
#             uu += u[x+dx, cells//4] * u[x, cells//4]
#             uv += u[x+dx, cells//4] * v[x, cells//4]

#             vu += v[x+dx, cells//4] * u[x, cells//4]
#             vv += v[x+dx, cells//4] * v[x, cells//4]

#         uu_.append(uu / (cells - 2*bound_dev - dx))
#         uv_.append(uv / (cells - 2*bound_dev - dx))

#         vu_.append(vu / (cells - 2*bound_dev - dx))
#         vv_.append(vv / (cells - 2*bound_dev - dx))

#     uu_time.append(uu_)
#     uv_time.append(uv_)

#     vu_time.append(vu_)
#     vv_time.append(vv_)


# fig, axs = plt.subplots(figsize=(10, 10))

# axs.plot(np.arange(0, max_dx), np.mean(np.array(uu_time), axis=0), linewidth=3, label= r'$ \langle u \cdot u \rangle $')
# axs.plot(np.arange(0, max_dx), np.mean(np.array(uv_time), axis=0), linewidth=3, label= r'$ \langle u \cdot v \rangle $')

# axs.plot(np.arange(0, max_dx), np.mean(np.array(vu_time), axis=0), linewidth=3, label= r'$ \langle v \cdot u \rangle $')
# axs.plot(np.arange(0, max_dx), np.mean(np.array(vv_time), axis=0), linewidth=3, label= r'$ \langle v \cdot v \rangle $')

# axs.grid(True, linestyle='--')
# axs.set_xlabel(r'$ dy $', fontsize=size_of_text)
# axs.set_ylabel(r'$ B(dy) = \langle U(y+dy, x/4) \cdot U(y, x/4) \rangle_{x, t} $', fontsize=size_of_text)
# axs.legend(fontsize=size_of_legend)
# axs.tick_params(labelsize=size_of_numbers)
# fig.tight_layout()
# plt.show()


# -------------------------------------- structure function --------------------------------------------------
u_time = []
v_time = []
w_time = []

max_space_dev = 250
max_tau = t_count//2
min_x, max_x = cells//2 - max_space_dev, cells//2 + max_space_dev
min_y, max_y = cells//2 - max_space_dev, cells//2 + max_space_dev

for i in range(min_stacionary_t, max_t):
    u, v, w, vort, p, cells = read_vtk(path + '/G_0_05_mu_0_1_Alpha_0_03_periodic_0-400/', i)

    u_time.append(u[min_y:max_y, min_x:max_x])
    v_time.append(v[min_y:max_y, min_x:max_x])

structure_function_u = []
structure_function_v = []

for tau in range(1, max_tau):
    D_u = np.zeros((max_y - min_y, max_x - min_x))
    D_v = np.zeros((max_y - min_y, max_x - min_x))
    for t in range(0, t_count - tau):
        D_u += np.power(np.array(u_time[t+tau] - u_time[t]), 2.0)
        D_v += np.power(np.array(v_time[t+tau] - v_time[t]), 2.0)

    structure_function_u.append(D_u / (t_count - tau))
    structure_function_v.append(D_v / (t_count - tau))

structure_function = np.array(np.array(structure_function_u) + np.array(structure_function_v), dtype=float)
structure_function_mean = np.mean(np.mean(structure_function, axis=1), axis=1)

print(structure_function.shape, structure_function_mean.shape)


def power_function(x, A, g):
    return A * np.power(x, g)


params, _ = curve_fit(power_function, np.arange(1, t_count//2)*time_step, structure_function_mean)
print(params)

fig, axs = plt.subplots(figsize=(10, 10))
axs.plot(np.arange(1, t_count//2)*time_step, structure_function_mean, linewidth=3, label= r'$ D (\tau) $')
axs.plot(np.linspace(time_step, t_count//2*time_step, 100, endpoint=True), power_function(np.linspace(time_step, t_count//2*time_step, 100, endpoint=True), *params), linewidth=3, label=str(round(params[0], 2)) + r'$ \, \tau $' + f"^ {params[1]:.2f}")

axs.grid(True, linestyle='--')
axs.set_xlabel(r'$ \tau $', fontsize=size_of_text)
axs.set_ylabel(r'$ D_U$ ', fontsize=size_of_text)
axs.legend(fontsize=size_of_legend)
axs.tick_params(labelsize=size_of_numbers)
fig.tight_layout()
plt.show()
