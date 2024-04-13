import numpy as np
import math
import matplotlib.pyplot as plt
from gc import collect
import vtk
from vtk.util.numpy_support import vtk_to_numpy

path = '/home/sq-wm/OpenFOAM/sq-wm-v2306/run/calc11/'
#path = '/home/sq-wm/OpenFOAM/sq-wm-v2306/run/calc13/'
#path = '/home/sq-wm/OpenFOAM/sq-wm-v2306/run/calc14/'

left_boundary_x =   -1.0*math.pi
right_boundary_x =  math.pi
down_boundary_y =   -1.0*math.pi
up_boundary_y =     math.pi
rear_boundary_z =   -1.0*math.pi
front_boundary_z =  math.pi

# read size of new structured array
with open(path + 'grid_size.txt', 'r') as fp:
    count_of_cells_axis = fp.read().split('\n')
    x_cells = int(count_of_cells_axis[0])
    y_cells = int(count_of_cells_axis[1])
    z_cells = int(count_of_cells_axis[2])
    print(x_cells, y_cells, z_cells)

dx = (right_boundary_x - left_boundary_x) / x_cells
dy = (up_boundary_y - down_boundary_y) / y_cells
dz = (front_boundary_z - rear_boundary_z) / z_cells

def write_3d_field(field, nx, ny, nz, name_with_path):
    file = open(name_with_path, "w+")
    for ind_z in range(nz):
        for ind_y in range(ny):
            for ind_x in range(nx):
                file.write(str(field[ind_x, ind_y, ind_z]) + ' ')
            file.write('\n')
        file.write('\n')
    file.close()

# Read the source file
def get_arrays_with_data_and_cells(name_of_file):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(name_of_file)
    reader.Update()  # Needed because of GetScalarRange
    #U_field = reader.GetOutput().GetPointData().GetArray('U')
    vort_field = reader.GetOutput().GetPointData().GetArray('vorticity')
    #data_vel = vtk_to_numpy(U_field)        # 2D array with components of velocity in points. data[i] = v[x, y, z] in i-point
    data_vort = vtk_to_numpy(vort_field)
    count_of_cells = len(data_vort)      # = 151^3
    print(count_of_cells)

    nodes_vtk_array = reader.GetOutput().GetPoints().GetData()
    numpy_nodes = vtk_to_numpy(nodes_vtk_array)
    #del U_field
    del vort_field
    del nodes_vtk_array
    collect()
    return data_vort, numpy_nodes
    #return data_vel, numpy_nodes

# return array with markers in informative cells and new arrays with data
def distribute_data_on_structured_grid(size_x, size_y, size_z, dx, dy, dz, vector_field, coord_field, count_of_points):
    counter_arr = np.zeros((x_cells, y_cells, z_cells), dtype=int)
    counter_mist = 0
    field_x = np.zeros((size_x, size_y, size_z), dtype=float)
    field_y = np.zeros((size_x, size_y, size_z), dtype=float)
    field_z = np.zeros((size_x, size_y, size_z), dtype=float)
    start_vector = [-left_boundary_x, -down_boundary_y, -rear_boundary_z]
    for i in range(0, count_of_points):
        coord_field[i] += start_vector

    # distribute the values by cells in structured array
    for i in range(0, count_of_points):
        ind_x = int(coord_field[i, 0] // dx)
        ind_y = int(coord_field[i, 1] // dy)
        ind_z = int(coord_field[i, 2] // dz)
        if (ind_x > size_x or ind_y > size_y or ind_z > size_z):
            counter_mist += 1
        if ind_x >= size_x:
            ind_x = size_x - 1 
        
        if ind_y >= size_y:
            ind_y = size_y - 1
        
        if ind_z >= size_z:
            ind_z = size_z - 1
        
        counter_arr[ind_x, ind_y, ind_z] += 1
        field_x[ind_x, ind_y, ind_z] += vector_field[i, 0]
        field_y[ind_x, ind_y, ind_z] += vector_field[i, 1]
        field_z[ind_x, ind_y, ind_z] += vector_field[i, 2]

        # handling the case of a hit on bound of cell
        if (coord_field[i, 0] % dx == 0) and (ind_x > 0) and (int(coord_field[i, 0] // dx) < size_x):
            counter_arr[ind_x-1, ind_y, ind_z] += 1
            field_x[ind_x-1, ind_y, ind_z] += vector_field[i, 0]
            field_y[ind_x-1, ind_y, ind_z] += vector_field[i, 1]
            field_z[ind_x-1, ind_y, ind_z] += vector_field[i, 2]
    
        if (coord_field[i, 1] % dy == 0) and (ind_y > 0) and (int(coord_field[i, 1] // dy) < size_y):
            counter_arr[ind_x, ind_y-1, ind_z] += 1
            field_x[ind_x, ind_y-1, ind_z] += vector_field[i, 0]
            field_y[ind_x, ind_y-1, ind_z] += vector_field[i, 1]
            field_z[ind_x, ind_y-1, ind_z] += vector_field[i, 2]

        if (coord_field[i, 2] % dz == 0) and (ind_z > 0) and (int(coord_field[i, 2] // dz) < size_z):
            counter_arr[ind_x, ind_y, ind_z-1] += 1
            field_x[ind_x, ind_y, ind_z-1] += vector_field[i, 0]
            field_y[ind_x, ind_y, ind_z-1] += vector_field[i, 1]
            field_z[ind_x, ind_y, ind_z-1] += vector_field[i, 2]

    print('count of mistakes = ' + str(counter_mist))

    for ind_x in range(0, size_x):
        for ind_y in range(0, size_y):
            for ind_z in range(0, size_z):
                if counter_arr[ind_x, ind_y, ind_z] > 0:
                    field_x[ind_x, ind_y, ind_z] /= counter_arr[ind_x, ind_y, ind_z]
                    field_y[ind_x, ind_y, ind_z] /= counter_arr[ind_x, ind_y, ind_z]
                    field_z[ind_x, ind_y, ind_z] /= counter_arr[ind_x, ind_y, ind_z]
                    counter_arr[ind_x, ind_y, ind_z] = 1
    
    return counter_arr, field_x, field_y, field_z

def bounds_of_part_of_array(ind, size, deviation):
    left = int(ind - deviation)
    right = int(ind + deviation)
    if left < 0:
        left = 0
    if right >= size:
        right = size - 1
    return left, right

def interpolate_field_to_smooth(size_x, size_y, size_z, field, markers):
    new_field = np.empty((size_x, size_y, size_z))
    for ind_x in range(0, size_x):
        for ind_y in range(0, size_y):
            for ind_z in range(0, size_z):
                if markers[ind_x, ind_y, ind_z] > 0:
                    new_field[ind_x, ind_y, ind_z] = field[ind_x, ind_y, ind_z]
                else:
                    radius = 1
                    lx, rx = bounds_of_part_of_array(ind_x, size_x, radius)
                    ly, ry = bounds_of_part_of_array(ind_y, size_y, radius)
                    lz, rz = bounds_of_part_of_array(ind_z, size_z, radius)
                    n = np.sum(markers[lx : rx, ly : ry, lz : rz])
                    while n == 0:
                        radius += 1
                        lx, rx = bounds_of_part_of_array(ind_x, size_x, radius)
                        ly, ry = bounds_of_part_of_array(ind_y, size_y, radius)
                        lz, rz = bounds_of_part_of_array(ind_z, size_z, radius)
                        n = np.sum(markers[lx : rx, ly : ry, lz : rz])
                    new_field[ind_x, ind_y, ind_z] = np.sum(field[lx : rx, ly : ry, lz : rz] * markers[lx : rx, ly : ry, lz : rz]) / n

    return new_field

def interpolate_and_write_from_file(path_to_data, path_to_saving_dir, size_x, size_y, size_z, dx, dy, dz):
    #velocity_arr, nodes_arr = get_arrays_with_data_and_cells(path_to_data)
    vort_arr, nodes_arr = get_arrays_with_data_and_cells(path_to_data)
    data_len = len(vort_arr)
    print('end reading file from ' + path_to_data)
    #markers, vel_x, vel_y, vel_z = distribute_data_on_structured_grid(size_x, size_y, size_z, dx, dy, dz, velocity_arr, nodes_arr, data_len)
    markers, vort_x, vort_y, vort_z = distribute_data_on_structured_grid(size_x, size_y, size_z, dx, dy, dz, vort_arr, nodes_arr, data_len)

    #vel_x_smooth = interpolate_field_to_smooth(size_x, size_y, size_z, vel_x, markers)
    #vel_y_smooth = interpolate_field_to_smooth(size_x, size_y, size_z, vel_y, markers)
    #vel_z_smooth = interpolate_field_to_smooth(size_x, size_y, size_z, vel_z, markers)
    vort_z_smooth = interpolate_field_to_smooth(size_x, size_y, size_z, vort_z, markers)
    print('end interpolating data')
    #write_3d_field(vel_x_smooth, size_x, size_y, size_z, path_to_saving_dir + 'U_x.txt')
    #write_3d_field(vel_y_smooth, size_x, size_y, size_z, path_to_saving_dir + 'U_y.txt')
    #write_3d_field(vel_z_smooth, size_x, size_y, size_z, path_to_saving_dir + 'U_z.txt')
    write_3d_field(vort_z_smooth, size_x, size_y, size_z, path_to_saving_dir + 'vort_z.txt')
    print('end saving data to ' + path_to_saving_dir)
    #del velocity_arr
    del vort_arr
    del nodes_arr
    del markers
    #del vel_x
    #del vel_y
    #del vel_z
    del vort_x
    del vort_y
    del vort_z
    #del vel_x_smooth
    #del vel_y_smooth
    #del vel_z_smooth
    del vort_z_smooth
    collect()

interpolate_and_write_from_file(path + '/VTK/calc11_89446/internal.vtu', path + '15/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_93387/internal.vtu', path + '15.5/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_97403/internal.vtu', path + '16/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_101443/internal.vtu', path + '16.5/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_105858/internal.vtu', path + '17/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_109910/internal.vtu', path + '17.5/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_114386/internal.vtu', path + '18/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_118141/internal.vtu', path + '18.5/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_122173/internal.vtu', path + '19/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_126184/internal.vtu', path + '19.5/', x_cells, y_cells, z_cells, dx, dy, dz)
interpolate_and_write_from_file(path + '/VTK/calc11_130401/internal.vtu', path + '20/', x_cells, y_cells, z_cells, dx, dy, dz)

#interpolate_and_write_from_file(path + '/VTK/calc13_106003/internal.vtu', path + '20/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_105412/internal.vtu', path + '19.9/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_104770/internal.vtu', path + '19.8/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_104099/internal.vtu', path + '19.7/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_103439/internal.vtu', path + '19.6/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_102815/internal.vtu', path + '19.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_102267/internal.vtu', path + '19.4/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_101715/internal.vtu', path + '19.3/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_101205/internal.vtu', path + '19.2/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_100698/internal.vtu', path + '19.1/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_100135/internal.vtu', path + '19/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_99606/internal.vtu', path + '18.9/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_99088/internal.vtu', path + '18.8/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_98593/internal.vtu', path + '18.7/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_98039/internal.vtu', path + '18.6/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_97538/internal.vtu', path + '18.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_97081/internal.vtu', path + '18.4/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_96615/internal.vtu', path + '18.3/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_96165/internal.vtu', path + '18.2/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_95717/internal.vtu', path + '18.1/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_95310/internal.vtu', path + '18/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_94876/internal.vtu', path + '17.9/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_94474/internal.vtu', path + '17.8/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_94067/internal.vtu', path + '17.7/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_93626/internal.vtu', path + '17.6/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_93092/internal.vtu', path + '17.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_92471/internal.vtu', path + '17.4/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_91785/internal.vtu', path + '17.3/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_91070/internal.vtu', path + '17.2/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_90484/internal.vtu', path + '17.1/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_89993/internal.vtu', path + '17/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_89440/internal.vtu', path + '16.9/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_88878/internal.vtu', path + '16.8/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_88450/internal.vtu', path + '16.7/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_88027/internal.vtu', path + '16.6/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_87523/internal.vtu', path + '16.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_86935/internal.vtu', path + '16.4/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_86297/internal.vtu', path + '16.3/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_85664/internal.vtu', path + '16.2/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_85009/internal.vtu', path + '16.1/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_84234/internal.vtu', path + '16/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_83499/internal.vtu', path + '15.9/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_82722/internal.vtu', path + '15.8/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_81888/internal.vtu', path + '15.7/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_81210/internal.vtu', path + '15.6/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_80545/internal.vtu', path + '15.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_79875/internal.vtu', path + '15.4/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_79237/internal.vtu', path + '15.3/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_78645/internal.vtu', path + '15.2/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_78123/internal.vtu', path + '15.1/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc13_77528/internal.vtu', path + '15/', x_cells, y_cells, z_cells, dx, dy, dz)

#interpolate_and_write_from_file(path + '/VTK/calc14_97752/internal.vtu', path + '19/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_99905/internal.vtu', path + '19.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_101859/internal.vtu', path + '20/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_104170/internal.vtu', path + '20.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_107032/internal.vtu', path + '21/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_109487/internal.vtu', path + '21.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_111531/internal.vtu', path + '22/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_114246/internal.vtu', path + '22.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_116548/internal.vtu', path + '23/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_119134/internal.vtu', path + '23.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_121475/internal.vtu', path + '24/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_124655/internal.vtu', path + '24.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_127109/internal.vtu', path + '25/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_129177/internal.vtu', path + '25.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_131865/internal.vtu', path + '26/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_133769/internal.vtu', path + '26.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_136272/internal.vtu', path + '27/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_138092/internal.vtu', path + '27.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_140009/internal.vtu', path + '28/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_142363/internal.vtu', path + '28.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_145364/internal.vtu', path + '29/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_148351/internal.vtu', path + '29.5/', x_cells, y_cells, z_cells, dx, dy, dz)
#interpolate_and_write_from_file(path + '/VTK/calc14_151339/internal.vtu', path + '30/', x_cells, y_cells, z_cells, dx, dy, dz)
