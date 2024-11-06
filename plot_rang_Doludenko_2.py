import numpy as np
import math
import matplotlib.pyplot as plt
from gc import collect
from scipy.optimize import curve_fit

path = "/home/sq-wm/OpenFOAM/sq-wm-v2306/run/rang_data/512_sort_rangs/"
filename_for_answers = 'approximate_parameters_large.txt'
# G = 0.05
names = ['1', '0_3', '0_25', '0_2', '0_15', '0_1', '0_05', '0_03', '0_01', '0_003', '0_002', '0_001', '0_0003', '0_0001']
colors = ['blue', 'orange', 'sienna', 'green', 'red', 'purple', 'yellow', 'black', 'cyan', 'lawngreen', 'magenta', 'darkgreen', 'slategray', 'darkviolet']
labels = ['1', '0.3', '0.25', '0.2', '0.15', '0.1', '0.05', '0.03', '0.01', '0.003', '0.002', '0.001', '0.0003', '0.0001']
# only 1 file for each regimes:
#names = ['1', '0_15', '0_01', '0_001_vort', '0_0001']
#colors = ['blue', 'red', 'black', 'darkgreen', 'darkviolet']
#labels = ['laminar ' + r'$(\alpha = 1)$', 'laminar-turbulent ' + r'$(\alpha = 0.15)$', 'turbulent ' + r'$(\alpha = 0.01)$', 'turbulent-vortex ' + r'$(\alpha = 0.001)$', 'vortex ' + r'$(\alpha = 0.0001)$']
# alpha = ['1', '0.3', '0.25', '0.2', '0.15', '0.1', '0.05', '0.03', '0.01', '0.003', '0.002', '0.001', '0.0003', '0.0001']
# labels = ['l', 'l', 'l', 'l-t', 'l-t', 'l-t', 't', 't', 't', 't', 't', 't-v', 'v', 'v']
# l = laminar; l-t = laminar-turbulent; t = turbulent; t-v = turbulent-vortex; v = vortex
markers = ['o', '*', 'D', 's', '^']
marker_colors = ['red', 'lawngreen', 'green', 'purple', 'gray']
alpha_labels = ['laminar', 'laminar-turbulent', 'turbulent', 'turbulent-vortex', 'vortex']
index_bounds = [0, 3, 6, 11, 12, len(names)]

# G = 0.15
#names = ['1', '0_3', '0_1', '0_03', '0_01', '0_003', '0_002', '0_001', '0_0003', '0_0001']
#colors = ['blue', 'orange', 'green', 'red', 'purple', 'yellow', 'black', 'cyan', 'lawngreen', 'magenta']
#labels = ['1', '0.3', '0.1', '0.03', '0.01', '0.003', '0.002', '0.001', '0.0003', '0.0001']
# labels = ['l', 'l-t', 't', 't', 't', 't-v', 't-v', 't-v', 'v', 'v']
# l = laminar; l-t = laminar-turbulent; t = turbulent; t-v = turbulent-vortex; v = vortex
# only 1 file for each regimes:
#names = ['1', '0_3', '0_03', '0_002', '0_0001']
#colors = ['blue', 'orange', 'red', 'black', 'magenta']
#labels = ['laminar ' + r'$(\alpha = 1)$', 'laminar-turbulent ' + r'$(\alpha = 0.3)$', 'turbulent ' + r'$(\alpha = 0.03)$', 'turbulent-vortex ' + r'$(\alpha = 0.002)$', 'vortex ' + r'$(\alpha = 0.0001)$']
#index_bounds = [0, 1, 2, 5, 8, len(names)]

# G = 0.05, periodic
#names = ['1', '0_3', '0_1', '0_03', '0_01', '0_003', '0_002', '0_001', '0_0003', '0_0001']
#colors = ['blue', 'orange', 'green', 'red', 'purple', 'yellow', 'black', 'cyan', 'lawngreen', 'magenta']
#labels = ['1', '0.3', '0.1', '0.03', '0.01', '0.003', '0.002', '0.001', '0.0003', '0.0001']
# only 1 file for each regimes:
#names = ['1', '0_03', '0_0001']
#colors = ['blue', 'red', 'magenta']
#labels = ['laminar ' + r'$(\alpha = 1)$', 'turbulent ' + r'$(\alpha = 0.03)$', 'vortex ' + r'$(\alpha = 0.0001)$']
#labels = ['l', 'l', 't', 't', 't', 'v', 'v', 'v', 'v', 'v']
#l = laminar; l-t = laminar-turbulent; t = turbulent; t-v = turbulent-vortex; v = vortex
#markers = ['o', 'D', '^']
#marker_colors = ['red', 'green', 'gray']
#alpha_labels = ['laminar', 'turbulent', 'vortex']
#index_bounds = [0, 2, 5, len(names)]

def rang_vort_from_file(path_to_file):
    file = open(path_to_file, 'r')
    data = file.read().split('\n')
    # data[0] - names of columns
    rang = []
    vort = []

    for ind in range(1, len(data) - 1):
        line = data[ind].split(' ')
        rang.append(line[0])
        vort.append(line[1])
    file.close()

    del line
    collect()
    print("end reading rang distribution from file " + path_to_file)
    rang = np.array(rang, dtype=int)
    vort = np.array(vort, dtype=float)
    return rang, vort

rang_arr_plus = []
rang_arr_minus = []
vort_plus = []
vort_minus = []

for i in range(len(names)):
    rang, vort = rang_vort_from_file(path + 'G=0.05/' + 'Vort_minus_sort_G_0_05_Alpha_' + names[i] + '_389-399.dat')
    rang_, vort_ = rang_vort_from_file(path + 'G=0.05/' + 'Vort_plus_sort_G_0_05_Alpha_' + names[i] + '_389-399.dat')

    #rang, vort = rang_vort_from_file(path + 'G=0.15/' + 'Vort_minus_sort_G_0_15_Alpha_' + names[i] + '.dat')
    #rang_, vort_ = rang_vort_from_file(path + 'G=0.15/' + 'Vort_plus_sort_G_0_15_Alpha_' + names[i] + '.dat')

    #rang, vort = rang_vort_from_file(path + 'periodic/' + 'Vort_minus_sort_G_0_05_Alpha_' + names[i] + '_periodic.dat')
    #rang_, vort_ = rang_vort_from_file(path + 'periodic/' + 'Vort_plus_sort_G_0_05_Alpha_' + names[i] + '_periodic.dat')

    rang_arr_minus.append(np.sort(rang, kind='mergesort'))
    rang_arr_plus.append(np.sort(rang_, kind='mergesort'))
    vort_minus.append(np.sort(np.abs(vort), kind='mergesort'))
    vort_plus.append(np.sort(vort_, kind='mergesort'))
    #print(len(rang), len(rang_))
    del rang
    del rang_
    del vort
    del vort_

collect()

def find_max_deviation(experimental_point, approximation_point):
    dev = np.abs(approximation_point - experimental_point) / experimental_point
    return np.max(dev)

def find_scalar_distribution_2D(scalar_field, y_cells, x_cells, left_bound, right_bound, steps_count):
    distribution = np.zeros(steps_count)
    dx = (right_bound - left_bound) / (steps_count - 1)
    for ind_y in range(y_cells):
        for ind_x in range(x_cells):
            value = scalar_field[ind_y, ind_x]
            if (value >= left_bound) and (value < right_bound):
                distribution[int((value - left_bound) // dx)] += 1
    return distribution

def find_analitical_solution(wave_vector, count_of_cells_x, count_of_cells_y, size_x, size_y, amplitude):
    vorticity_field = np.empty([count_of_cells_x, count_of_cells_y])
    x_array = np.linspace(0.0, size_x, count_of_cells_x)
    y_array = np.linspace(0.0, size_y, count_of_cells_y)
    for ind_x in range(count_of_cells_x):
        for ind_y in range(count_of_cells_y):
            vorticity_field[ind_x, ind_y] = math.cos(wave_vector*x_array[ind_x]) + math.cos(wave_vector*y_array[ind_y])

    distribution_plus = find_scalar_distribution_2D(vorticity_field, count_of_cells_y, count_of_cells_x, 0.0, np.max(vorticity_field), 101)
    distribution_minus = find_scalar_distribution_2D(vorticity_field, count_of_cells_y, count_of_cells_x, np.min(vorticity_field), 0.0, 101)
    distribution_plus = np.sort(distribution_plus, kind='mergesort')
    distribution_minus = np.sort(distribution_minus, kind='mergesort')
    ampl = count_of_cells_x * count_of_cells_y
    return distribution_plus[1:] * amplitude / ampl, distribution_minus[1:] * amplitude / ampl

size_of_text = 15
size_of_numbers = 15
size_of_legend = 20
marker_size = 100

max_amplitude = 0
for i in range(len(names)):
    ampl = max([np.max(vort_minus[i]), np.max(vort_plus[i])])
    if ampl > max_amplitude:
        max_amplitude = ampl

approximation_ampl = 0.0
count_of_coeffs = 6
def log_poly_func(x, a, b, c, d, e, f):
    return approximation_ampl*np.log10(a*np.square(x) + b*x + c) - approximation_ampl*np.log10(d*np.square(x) + e*x + f)

#variable_names = ['a1 = ', '\t' + 'b1 = ', '\t' + 'c1 = ', '\n' + '\t' + 'a2 = ', '\t' + 'b2 = ', '\t' + 'c2 = ']

# file = open(path + filename_for_answers, "w+")
# file.write('vorticity for rang, func = k * log( (ax^2 + bx + c) / (dx^2 + ex + f) )' + '\n')

multiplier = 10000.0
max_rang_minus = 0
max_rang_plus = 0
for i in range(len(names)):
    rang_arr_minus[i] = np.array(rang_arr_minus[i], dtype=float)
    rang_arr_plus[i] = np.array(rang_arr_plus[i], dtype=float)
    rang_arr_minus[i] /= multiplier
    rang_arr_plus[i] /= multiplier

    if rang_arr_minus[i][-1] > max_rang_minus:
        max_rang_minus = rang_arr_minus[i][-1]
    if rang_arr_plus[i][-1] > max_rang_plus:
        max_rang_plus = rang_arr_plus[i][-1]

# vorticity for rang
fig, axs = plt.subplots(1, 2)       #, figsize=[15, 9], dpi=1000
count_of_points = 10

for i in range(len(names)):
    axs[0].plot(rang_arr_minus[i], vort_minus[i], linewidth=1, label=labels[i], color=colors[i])
    axs[1].plot(rang_arr_plus[i], vort_plus[i], linewidth=1, label=labels[i], color=colors[i])

    # index_array_minus = np.linspace(0, len(rang_arr_minus_G_0_05[i]) - 1, count_of_points, dtype = int)
    # index_array_plus = np.linspace(0, len(rang_arr_plus_G_0_05[i]) - 1, count_of_points, dtype = int)
    # markers_on_minus = np.empty(count_of_points)
    # markers_minus = np.empty(count_of_points)
    # markers_on_plus = np.empty(count_of_points)
    # markers_plus = np.empty(count_of_points)
    # for j in range(count_of_points):
    #     markers_on_minus[j] = rang_arr_minus_G_0_05[i][index_array_minus[j]]
    #     markers_minus[j] = vort_minus_G_0_05[i][index_array_minus[j]]
    #     markers_on_plus[j] = rang_arr_plus_G_0_05[i][index_array_plus[j]]
    #     markers_plus[j] = vort_plus_G_0_05[i][index_array_plus[j]]

    # axs[0].plot(markers_on_minus, markers_minus, label=labels[i], markersize=marker_size) #, markevery=markers_on_minus; str(markers[i]), color=colors[i],
    # axs[1].plot(markers_on_plus, markers_plus, label=labels[i], markersize=marker_size) #, markevery=markers_on_plus; str(markers[i]), color=colors[i],

    approximation_ampl = vort_minus[i][len(vort_minus[i]) - 1]
    params, errors = curve_fit(log_poly_func, rang_arr_minus[i], vort_minus[i])
    #file.write(str(names[i]) + '\n' + 'vorticity < 0' + '\n')
    #file.write('k = ' + str(approximation_ampl) + '\t')
    #for ind in range(count_of_coeffs):
    #    file.write(variable_names[ind] + str(params[ind]))
    #file.write('\n' + '\n')
    axs[0].plot(rang_arr_minus[i], log_poly_func(rang_arr_minus[i], *params), linewidth=2, linestyle='--', color=colors[i])
    #file.write('the max deviation of the approximation = ' + str(find_max_deviation(vort_minus[i][len(vort_minus[i]) // 2 : ], log_poly_func(rang_arr_minus[i], *params)[len(vort_minus[i]) // 2 : ])) + '\n')

    approximation_ampl = vort_plus[i][len(vort_plus[i]) - 1]
    params, errors = curve_fit(log_poly_func, rang_arr_plus[i], vort_plus[i])
    #file.write('vorticity > 0' + '\n')
    #file.write('k = ' + str(approximation_ampl) + '\t')
    #for ind in range(count_of_coeffs):
    #    file.write(variable_names[ind] + str(params[ind]))
    #file.write('\n' + '\n')
    axs[1].plot(rang_arr_plus[i], log_poly_func(rang_arr_plus[i], *params), linewidth=2, linestyle='--', color=colors[i])
    #file.write('the max deviation of the approximation = ' + str(find_max_deviation(vort_plus[i][len(vort_plus[i]) // 2 : ], log_poly_func(rang_arr_plus[i], *params)[len(vort_plus[i]) // 2 : ])) + '\n')

axs[0].grid(True, linestyle='--')
axs[0].set_xlabel('R', fontsize=size_of_text)
axs[0].set_ylabel('|' + r'$ \omega $' + '|' + r'$ (\alpha)$', fontsize=size_of_text)
axs[0].set_title(r'$ \omega < 0 $', fontsize=size_of_text)
axs[0].set_ylim(-0.05*max_amplitude, 1.05*max_amplitude)
axs[0].legend(fontsize=size_of_legend)
axs[0].tick_params(labelsize=size_of_numbers)
axs[0].text(max_rang_minus*0.35, max_amplitude*0.9, r'$ \vert \omega \vert = k \, log \left( \frac{ a_1 R^2 + b_1 R + c_1 } { a_2 R^2 + b_2 R + c_2 } \right) $', fontsize=size_of_legend)
#axs[0].text(0, max_amplitude*0.3, r'$ \vert \omega \vert = k \, log \left( \frac{ a_1 R^2 + b_1 R + c_1 } { a_2 R^2 + b_2 R + c_2 } \right) $', fontsize=size_of_legend)

axs[1].grid(True, linestyle='--')
axs[1].set_xlabel('R', fontsize=size_of_text)
axs[1].set_ylabel(r'$ \omega (\alpha)$', fontsize=size_of_text)
axs[1].set_title(r'$ \omega > 0 $', fontsize=size_of_text)
axs[1].set_ylim(-0.05*max_amplitude, 1.05*max_amplitude)
axs[1].legend(fontsize=size_of_legend)
axs[1].tick_params(labelsize=size_of_numbers)
axs[1].text(max_rang_plus*0.35, max_amplitude*0.9, r'$ \omega = k \, log \left( \frac{ a_1 R^2 + b_1 R + c_1 } { a_2 R^2 + b_2 R + c_2 } \right) $', fontsize=size_of_legend)
#axs[1].text(0, max_amplitude*0.3, r'$ \omega = k \, log \left( \frac{ a_1 R^2 + b_1 R + c_1 } { a_2 R^2 + b_2 R + c_2 } \right) $', fontsize=size_of_legend)

#file.close()
plt.show()
#plt.savefig(path + 'vorticity_for_rangs_small.eps', format='eps')

def poly_func(x, a, b, c, d, e, f, g, h, i, j, k):
    return a*x + b*np.square(x) + c*np.power(x, 3) + d*np.power(x, 4) + e*np.power(x, 5) + f*np.power(x, 6) + g*np.power(x, 7) + h*np.power(x, 8) + i*np.power(x, 9) + j*np.power(x, 10) + k*np.power(x, 11)

count_of_freq_rang = 100
def find_rang_distribution(scalar_field, rang_count, cells_count):
    distribution = np.zeros(rang_count)
    right_bound = np.max(scalar_field)
    left_bound = np.min(scalar_field)
    dx = (right_bound - left_bound) / (rang_count - 1)
    for ind in range(len(scalar_field)):
        distribution[int((scalar_field[ind] - left_bound) // dx)] += 1                 
    return np.sort(distribution, kind='mergesort') / cells_count

def find_PDF_distribution(scalar_field, rang_count, left, right):
    distribution = np.zeros(rang_count)
    dx = (right - left) / (rang_count - 1)
    for ind in range(len(scalar_field)):
        #if int((scalar_field[ind] - left) // dx) >= rang_count:
        #    distribution[rang_count - 1] += 1
        #elif int((scalar_field[ind] - left) // dx) < 0:
        #    distribution[0] += 1
        #else:
        #    distribution[int((scalar_field[ind] - left) // dx)] += 1
        if int((scalar_field[ind] - left) // dx) >= 0 and int((scalar_field[ind] - left) // dx) < rang_count:
            distribution[int((scalar_field[ind] - left) // dx)] += 1
    return distribution / len(scalar_field)

left_bound = 0.0
right_bound = 0.0
for i in range(len(names)):
    if np.max(vort_minus[i]) > left_bound:
        left_bound = np.max(vort_minus[i])
    if np.max(vort_plus[i]) > right_bound:
        right_bound = np.max(vort_plus[i])
#left_bound *= -1.0
left_bound = -6.0
right_bound = 6.0

distribution_minus = []
distribution_plus = []
PDF_distribution = []
for i in range(len(names)):
    distribution_minus.append(find_rang_distribution(vort_minus[i], count_of_freq_rang, len(rang_arr_minus[i]) + len(rang_arr_plus[i])))
    distribution_plus.append(find_rang_distribution(vort_plus[i], count_of_freq_rang, len(rang_arr_minus[i]) + len(rang_arr_plus[i])))
    PDF_distribution.append(find_PDF_distribution(np.concatenate([-1.0 * vort_minus[i], vort_plus[i]]), count_of_freq_rang * 2, left_bound, right_bound))

distribution_minus = np.array(distribution_minus)
distribution_plus = np.array(distribution_plus)
PDF_distribution = np.array(PDF_distribution)
max_amplitude = max([np.max(distribution_minus), np.max(distribution_plus)])

approximation_minus = []
approximation_plus = []

approximation_ampl = 0.0
count_of_coeffs = 6
def log_poly_func2(x, a, b, c, d, e, f):
    #return e*np.log10(a + b*x + c*np.square(x) + d*np.power(x, 3))# - approximation_ampl*np.log10(d*np.power(x, 2) + e)
    return a + b*x + c*np.square(x) + d*np.power(x, 3) + e*np.power(x, 4) + f*np.power(x, 5)

# frequency of vorticity distribution
# fig, axs = plt.subplots(1, 2, figsize=(16, 10), dpi=600)
fig, axs = plt.subplots(1, 2)
for i in range(len(names)):
    axs[0].plot(np.arange(count_of_freq_rang), distribution_minus[i], linewidth=1, label=labels[i], color=colors[i])
    approximation_ampl = np.max(distribution_minus[i])
    params, errors = curve_fit(log_poly_func, np.arange(count_of_freq_rang), distribution_minus[i])
    approximation_minus.append(log_poly_func(np.linspace(0, count_of_freq_rang, count_of_freq_rang*10), *params))
    axs[0].plot(np.linspace(0, count_of_freq_rang, count_of_freq_rang*10), log_poly_func(np.linspace(0, count_of_freq_rang, count_of_freq_rang*10), *params), linewidth=2, linestyle='--', color=colors[i])

    axs[1].plot(np.arange(count_of_freq_rang), distribution_plus[i], linewidth=1, label=labels[i], color=colors[i])
    approximation_ampl = np.max(distribution_plus[i])
    params, errors = curve_fit(log_poly_func, np.arange(count_of_freq_rang), distribution_plus[i])
    approximation_plus.append(log_poly_func(np.linspace(0, count_of_freq_rang, count_of_freq_rang*10), *params))
    axs[1].plot(np.linspace(0, count_of_freq_rang, count_of_freq_rang*10), log_poly_func(np.linspace(0, count_of_freq_rang, count_of_freq_rang*10), *params), linewidth=2, linestyle='--', color=colors[i])

axs[0].grid(True, linestyle='--')
axs[0].set_xlabel('r', fontsize=size_of_text)
axs[0].set_ylabel(r'$ E (\alpha)$', fontsize=size_of_text)
axs[0].set_title(r'$ \omega < 0 $', fontsize=size_of_text)
axs[0].set_ylim(-0.05*max_amplitude, 1.05*max_amplitude)
axs[0].legend(fontsize=size_of_legend)
axs[0].tick_params(labelsize=size_of_numbers)
axs[0].text(count_of_freq_rang*0.35, max_amplitude*0.9, r'$ E = k \, log \left( \frac{ a_1 r^2 + b_1 r + c_1 } { a_2 r^2 + b_2 r + c_2 } \right) $', fontsize=size_of_legend)
#axs[0].text(0, max_amplitude*0.3, r'$ E = k \, log \left( \frac{ a_1 r^2 + b_1 r + c_1 } { a_2 r^2 + b_2 r + c_2 } \right) $', fontsize=size_of_legend)

axs[1].grid(True, linestyle='--')
axs[1].set_xlabel('r', fontsize=size_of_text)
axs[1].set_ylabel(r'$ E (\alpha)$', fontsize=size_of_text)
axs[1].set_title(r'$ \omega > 0 $', fontsize=size_of_text)
axs[1].set_ylim(-0.05*max_amplitude, 1.05*max_amplitude)
axs[1].legend(fontsize=size_of_legend)
axs[1].tick_params(labelsize=size_of_numbers)
axs[1].text(count_of_freq_rang*0.35, max_amplitude*0.9, r'$ E = k \, log \left( \frac{ a_1 r^2 + b_1 r + c_1 } { a_2 r^2 + b_2 r + c_2 } \right) $', fontsize=size_of_legend)
#axs[1].text(0, max_amplitude*0.3, r'$ E = k \, log \left( \frac{ a_1 r^2 + b_1 r + c_1 } { a_2 r^2 + b_2 r + c_2 } \right) $', fontsize=size_of_legend)

plt.show()
#plt.savefig(path + 'frequency_for_rangs_' + str(count_of_freq_rang) + '.png')

# PDF distribution
max_amplitude = np.max(PDF_distribution)
# fig, axs = plt.subplots(figsize=(16, 10), dpi=600)
fig, axs = plt.subplots()
for i in range(len(names)):
    axs.plot(np.linspace(left_bound, right_bound, count_of_freq_rang * 2), PDF_distribution[i], linewidth=1, label=labels[i], color=colors[i])

axs.grid(True, linestyle='--')
axs.set_xlabel(r'$ \omega $', fontsize=size_of_text)
axs.set_ylabel(r'$ E (\omega)$', fontsize=size_of_text)
axs.set_title('PDF distribution ' + r'$ (\omega) $', fontsize=size_of_text)
axs.set_ylim(0.0, 1.05*max_amplitude)
axs.legend(fontsize=size_of_legend)
axs.tick_params(labelsize=size_of_numbers)
#axs.text(0, np.max(max_amplitude) / 5.0, r'$ E = k \, log \left( \frac{ a_1 r^2 + b_1 r + c_1 } { a_2 r^2 + b_2 r + c_2 } \right) $', fontsize=size_of_legend)

plt.show()
#plt.savefig(path + 'PDF_distr_' + str(count_of_freq_rang) + '.png')

rang_max_minus = []
rang_max_plus = []
rang_max_minus_appr = []
rang_max_plus_appr = []
# derivation of approximative function
fig, axs = plt.subplots(1, 2)
for i in range(len(names)):
    axs[0].plot(np.linspace(0, count_of_freq_rang, count_of_freq_rang*10), np.gradient(approximation_minus[i], edge_order=2), linewidth=1, label=labels[i], color=colors[i])
    axs[1].plot(np.linspace(0, count_of_freq_rang, count_of_freq_rang*10), np.gradient(approximation_plus[i], edge_order=2), linewidth=1, label=labels[i], color=colors[i])
    axs[0].plot(np.linspace(0, count_of_freq_rang, count_of_freq_rang), np.gradient(distribution_minus[i], edge_order=2), linewidth=1, label=labels[i], color=colors[i])
    axs[1].plot(np.linspace(0, count_of_freq_rang, count_of_freq_rang), np.gradient(distribution_plus[i], edge_order=2), linewidth=1, label=labels[i], color=colors[i])

    rang_max_minus_appr.append(round(np.argmax(np.gradient(approximation_minus[i], edge_order=2)) / 10.0))
    rang_max_plus_appr.append(round(np.argmax(np.gradient(approximation_plus[i], edge_order=2)) / 10.0))
    rang_max_minus.append(np.argmax(np.gradient(distribution_minus[i], edge_order=2)))
    rang_max_plus.append(np.argmax(np.gradient(distribution_plus[i], edge_order=2)))

axs[0].grid(True, linestyle='--')
axs[0].set_xlabel('r', fontsize=size_of_text)
axs[0].set_ylabel(r'$ \nabla E (\alpha)$', fontsize=size_of_text)
axs[0].set_title(r'$ \omega < 0 $', fontsize=size_of_text)
#axs[0].set_ylim(-0.05*max_amplitude, 1.05*max_amplitude)
axs[0].legend(fontsize=size_of_legend)
axs[0].tick_params(labelsize=size_of_numbers)

axs[1].grid(True, linestyle='--')
axs[1].set_xlabel('r', fontsize=size_of_text)
axs[1].set_ylabel(r'$ \nabla E (\alpha)$', fontsize=size_of_text)
axs[1].set_title(r'$ \omega > 0 $', fontsize=size_of_text)
#axs[1].set_ylim(-0.05*max_amplitude, 1.05*max_amplitude)
axs[1].legend(fontsize=size_of_legend)
axs[1].tick_params(labelsize=size_of_numbers)
plt.show()

alpha = np.array(labels, dtype=float)

# the inflection point for alpha
#fig, axs = plt.subplots(1, 2, figsize=(16, 10), dpi=600)
fig, axs = plt.subplots(1, 2)

axs[0].plot(alpha, rang_max_minus, linewidth=2, label='initial distribution', linestyle='-')
axs[1].plot(alpha, rang_max_plus, linewidth=2, label='initial distribution', linestyle='-')
axs[0].plot(alpha, rang_max_minus_appr, linewidth=2, label='approximation', linestyle='--')
axs[1].plot(alpha, rang_max_plus_appr, linewidth=2, label='approximation', linestyle='--')

for i in range(len(markers)):
    axs[0].scatter(alpha[index_bounds[i]:index_bounds[i+1]], rang_max_minus[index_bounds[i]:index_bounds[i+1]], color=marker_colors[i], marker=markers[i], s=marker_size, label=alpha_labels[i]) 
    axs[1].scatter(alpha[index_bounds[i]:index_bounds[i+1]], rang_max_plus[index_bounds[i]:index_bounds[i+1]], color=marker_colors[i], marker=markers[i], s=marker_size, label=alpha_labels[i])
    axs[0].scatter(alpha[index_bounds[i]:index_bounds[i+1]], rang_max_minus_appr[index_bounds[i]:index_bounds[i+1]], color=marker_colors[i], marker=markers[i], s=marker_size)
    axs[1].scatter(alpha[index_bounds[i]:index_bounds[i+1]], rang_max_plus_appr[index_bounds[i]:index_bounds[i+1]], color=marker_colors[i], marker=markers[i], s=marker_size)

fig.suptitle(r'coordinates of the inflection point $(\alpha)$', fontsize=size_of_legend)
axs[0].grid(True, linestyle='--')
axs[0].set_xlabel(r'$ \alpha $', fontsize=size_of_text)
axs[0].set_ylabel('r', fontsize=size_of_text)
axs[0].set_title(r'$ \omega < 0 $', fontsize=size_of_text)
#axs[0].set_ylim(-0.05*max_amplitude, 1.05*max_amplitude)
axs[0].legend(fontsize=size_of_legend)
axs[0].tick_params(labelsize=size_of_numbers)
axs[0].set_xscale('log')

axs[1].grid(True, linestyle='--')
axs[1].set_xlabel(r'$ \alpha $', fontsize=size_of_text)
axs[1].set_ylabel('r', fontsize=size_of_text)
axs[1].set_title(r'$ \omega > 0 $', fontsize=size_of_text)
#axs[1].set_ylim(-0.05*max_amplitude, 1.05*max_amplitude)
axs[1].legend(fontsize=size_of_legend)
axs[1].tick_params(labelsize=size_of_numbers)
axs[1].set_xscale('log')
plt.show()
#plt.savefig(path + 'inflection_point_' + str(count_of_freq_rang) + '.png')
