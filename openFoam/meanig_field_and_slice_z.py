import numpy as np
import math
import matplotlib.pyplot as plt
from gc import collect
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from matplotlib.colors import Normalize

path = "/home/sq-wm/OpenFOAM/sq-wm-v2306/run/calc13/"
# read size of new structured array
with open(path + 'grid_size.txt', 'r') as fp:
    count_of_cells_axis = fp.read().split('\n')
    x_cells = int(count_of_cells_axis[0])
    y_cells = int(count_of_cells_axis[1])
    z_cells = int(count_of_cells_axis[2])
    print('structured grid has ' + str(x_cells) + ' * ' + str(y_cells) + ' * ' + str(z_cells) + ' cells')
cells_number_axis = min([x_cells, y_cells, z_cells])

txt_file_name_table = 'table_data.txt'
file_table = open(path + txt_file_name_table, "w+")
#file_table.write('time' + '\t' + 'U_RMS' + '\t' + 'vort_RMS' + '\t' + 'Re_makro' + '\t' + 'Ro_makro' + '\t' + 'Re_mikro' + '\t' + 'Re_mikro' + '\n')
file_table.write('time' + '\t' + '\t' + 'U_RMS' + '\t' + '\t' + 'Re_makro' + '\t' + 'Ro_makro' + '\t' + 'center_dist' + '\t' + 'energy' + '\n')
file_table.close()

start_time = 18
end_time = 18.2
time_step = 0.1
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
    print("end reading velocity from file " + full_path_to_files)
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
    max_rad = 0
    rad = 0
    for ind_y in range(int(y_cells * 0.25), int(y_cells * 0.75)):
        for ind_x in range(int(x_cells * 0.25), int(x_cells * 0.75)):
            r_cyc = find_max_anticyclone_radius(vort_z_slice, ind_x, ind_y, y_cells, x_cells)
            rad_new, _ = get_polar_coordinate(x_geom, ind_x, y_geom, ind_y)
            if r_cyc > max_rad and rad_new < rad:
                rad = rad_new
                max_rad = r_cyc
    return rad

def find_full_energy(ux, uy, uz, dx, dy, dz):
    return np.sum(np.square(ux) + np.square(uy) + np.square(uz)) * dx * dy * dz * 0.5

def bounds_of_part_of_array(ind, size, deviation):
    left = int(ind - deviation)
    right = int(ind + deviation)
    if left < 0:
        left = 0
    if right >= size:
        right = size - 1
    return left, right

def find_max_anticyclone_radius(slice, x0, y0, size_x, size_y):
    radius = 0
    if slice[y0, x0] < 0:
        radius = 1
        sum = 0
        while sum == 0:
            left_y, right_y = bounds_of_part_of_array(y0, size_y, radius)
            for ind_y in range(left_y, right_y + 1):
                left_x, right_x = bounds_of_part_of_array(x0, size_x, radius)
                for ind_x in range(left_x, right_x + 1):
                    r, _ = get_polar_coordinate(ind_x, x0, ind_y, y0)
                    if r <= radius:
                        if slice[ind_y, ind_x] >= 0:
                            sum = 1
                            ind_y = y0 + radius + 2
                            ind_x = x0 + radius + 2
            radius += 1
        radius -= 1
    return radius

dx = (right_boundary_x - left_boundary_x) / x_cells
dy = (up_boundary_y - down_boundary_y) / y_cells
dz = (front_boundary_z - rear_boundary_z) / z_cells

#u_mean = np.zeros((z_cells, y_cells, x_cells))
#v_mean = np.zeros((z_cells, y_cells, x_cells))
#w_mean = np.zeros((z_cells, y_cells, x_cells))
#vort_z_mean = np.zeros((z_cells, y_cells, x_cells))
vort_z_slice_time_array = []

file_table = open(path + txt_file_name_table, "a+")
for time in time_array:
    if time % 1 == 0:
        time = int(time)
    print('time = ' + str(time))
    u_center, v_center, w_center = get_data_from_file(path + str(time) + '/', x_cells, y_cells, z_cells)
    vort_z = get_vorticity_z_from_file(path + str(time) + '/', x_cells, y_cells, z_cells)
    vort_z_slice = get_z_mean_slice(vort_z, z_cells, y_cells, x_cells)
    vort_z_slice[0 : 3, : ] = 0.0
    vort_z_slice[y_cells-3 : , : ] = 0.0
    vort_z_slice[ : , 0 : 3] = 0.0
    vort_z_slice[ : , x_cells-3 : ] = 0.0
    vort_z_slice_time_array.append(vort_z_slice)

    #fig, axs = plt.subplots()
    #axs.set_title("mean z vorticity")
    #plt.imshow(vort_z_slice, cmap='seismic')
    #plt.show()

    u_rms = find_RMS_velocity(u_center, v_center, w_center, x_cells*y_cells*z_cells)
    Re = find_makro_Reynolds_number(u_rms, Lf, nu)
    Ro = find_makro_Rossby_number(u_rms, Lf, OMEGA)
    dist = find_distance_to_center(vort_z_slice, y_cells, x_cells)
    E = find_full_energy(u_center, v_center, w_center, dx, dy, dz)

    file_table.write(str(time) + '\t' + '\t')
    file_table.write(str(round(u_rms, 1)) + '\t' + '\t')
    file_table.write(str(int(round(Re, 0))) + '\t' + '\t')
    file_table.write(str(round(Ro, 3)) + '\t' + '\t')
    file_table.write(str(round(dist, 2)) + '\t' + '\t')
    file_table.write(str(int(round(E, 0))) + '\n')

#    u_mean += u_center
#    v_mean += v_center
#    w_mean += w_center
#    vort_z_mean += vort_z

    del u_center
    del v_center
    del w_center
    del vort_z
    del vort_z_slice
    collect()

file_table.close()
vort_z_slice_time_array = np.array(vort_z_slice_time_array)

fig, ax = plt.subplots()
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
ax.set_title('vorticity z mean')      #tx =         , time = ' + str(time_array[0]))

ims = []
for i in range(len(time_array)):
    im = ax.imshow(vort_z_slice_time_array[i], animated=True)      # , vmin=-120.0, vmax=120.0, cmap='seismic'
    if i == 0:
        ax.imshow(vort_z_slice_time_array[0])  # show an initial one first , vmin=-120.0, vmax=120.0, cmap='seismic'
    #tx.set_text('vorticity z mean, time = {0}'.format(time_array[i]))
    fig.colorbar(im, cax=cax)
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=300, blit=True, repeat_delay=1000)

ani.save(path + 'vorticity_slice' + '.gif', fps=3)
#plt.show()

# meaning
#u_mean /= ((end_time - start_time)/time_step + 1)
#v_mean /= ((end_time - start_time)/time_step + 1)
#w_mean /= ((end_time - start_time)/time_step + 1)
#vort_z_mean /= ((end_time - start_time)/time_step + 1)

#write_3d_field(u_mean, x_cells, y_cells, z_cells, path + 'Ux_mean.txt')
#write_3d_field(v_mean, x_cells, y_cells, z_cells, path + 'Uy_mean.txt')
#write_3d_field(w_mean, x_cells, y_cells, z_cells, path + 'Uz_mean.txt')
#write_3d_field(vort_z_mean, x_cells, y_cells, z_cells, path + 'vort_z_mean.txt')

#u_mean, v_mean, w_mean = get_mean_data_from_file(path, x_cells, y_cells, z_cells)
#vort_z_mean = get_mean_vorticity_from_file(path, x_cells, y_cells, z_cells)

#print('Ux mean value = ' + str(np.mean(u_mean)))
#print('Uy mean value = ' + str(np.mean(v_mean)))
#print('Uz mean value = ' + str(np.mean(w_mean)))

#print('Ux RMS value = ' + str(np.mean(np.square(u_mean)) ** 0.5))
#print('Uy RMS value = ' + str(np.mean(np.square(v_mean)) ** 0.5))
#print('Uz RMS value = ' + str(np.mean(np.square(w_mean)) ** 0.5))

#u_slice = np.zeros((y_cells, x_cells))
#v_slice = np.zeros((y_cells, x_cells))
#w_slice = np.zeros((y_cells, x_cells))
#vort_slice = np.zeros((y_cells, x_cells))

#for z_ind in range(z_cells):
#    u_slice += u_mean[z_ind, : , :]
#    v_slice += v_mean[z_ind, : , :]
#    w_slice += w_mean[z_ind, : , :]
#    vort_slice += vort_z_mean[z_ind, : , :]
#vort_slice[0 : 3, : ] = 0.0
#vort_slice[y_cells-3 : , : ] = 0.0
#vort_slice[ : , 0 : 3] = 0.0
#vort_slice[ : , x_cells-3 : ] = 0.0

#u_slice /= z_cells
#v_slice /= z_cells
#w_slice /= z_cells
#vort_slice /= z_cells
#mod_vel_slice = np.sqrt(u_slice*u_slice + v_slice*v_slice + w_slice*w_slice)

#min_vort = np.min(vort_slice)
#ind_y_min = np.where(vort_slice == min_vort)[0][0]
#ind_x_min = np.where(vort_slice == min_vort)[1][0]
#print('min vorticity = ' + str(vort_slice[ind_y_min, ind_x_min]))
#print('in point (' + str(ind_y_min) + ', ' + str(ind_x_min) + ')')

#r1, _ = get_polar_coordinate(0, ind_x_min, 0, ind_y_min)
#r2, _ = get_polar_coordinate(x_cells, ind_x_min, 0, ind_y_min)
#r3, _ = get_polar_coordinate(0, ind_x_min, y_cells, ind_y_min)
#r4, _ = get_polar_coordinate(x_cells, ind_x_min, y_cells, ind_y_min)
#max_radius = min([r1, r2, r3, r4])
#print(max_radius)
#count_of_radius_point = 200
#radius_array_vel_projection = np.zeros(count_of_radius_point)
#radius_array_counter = np.zeros(count_of_radius_point)

#for ind_y in range(y_cells):
#    for ind_x in range(x_cells):
#        r, phi = get_polar_coordinate(ind_x, ind_x_min, ind_y, ind_y_min)
#        if r <= max_radius:
#           azimuthal_vel = v_slice[ind_y, ind_x]*math.cos(phi) - u_slice[ind_y, ind_x]*math.sin(phi)
#           ind = int((count_of_radius_point - 1) * r / (max_radius - 0.0))
#           radius_array_vel_projection[ind] += azimuthal_vel
#           radius_array_counter[ind] += 1

#print(radius_array_vel_projection)
#radius_array_vel_projection /= radius_array_counter
#print(radius_array_vel_projection)

#def optimisation_function(x, A, B):
#    return A * x * np.log(x / B)

#fig, axs = plt.subplots(1, 2)
#axs.set_title("mean z vorticity")
#plt.imshow(vort_slice, cmap='seismic')
#axs[0].imshow(u_slice)
#axs[1].imshow(v_slice)
#plt.imshow(mod_vel_slice, cmap='seismic')
#plt.show()

#print(np.linspace(0.0, max_radius, count_of_radius_point))
#x_gr = np.linspace(0.0, max_radius, count_of_radius_point)[1:int(count_of_radius_point*2/3)]
#params, _ = curve_fit(optimisation_function, x_gr, radius_array_vel_projection[1:int(count_of_radius_point*2/3)], bounds=([0.1, 80.0], [2.0, 120.0]))
#print(params)

#fig, axs = plt.subplots()
#axs.plot(np.linspace(0.0, max_radius, count_of_radius_point), radius_array_vel_projection, label='results')
#axs.plot(x_gr, optimisation_function(x_gr, *params), linestyle='dashed', label='approximation')
#axs.grid(True, linestyle='--')
#axs.set_ylabel('mean azimuthal velocity')
#axs.set_xlabel('radius')
#plt.xlim(0.0, max_radius*1.1)
#plt.legend()
#plt.show()
