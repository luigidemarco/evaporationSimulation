import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import argparse
import glob
import sys

parser = argparse.ArgumentParser(description='Process output files.')
parser.add_argument("files", type=str, help="List of file names", nargs='*')
parser.add_argument("-f", "--fit", type=str, help="Specify exponential [e] or 2-body [2b] fit.", choices=['e', '2b'])
parser.add_argument("-e", "--exclude", type=float, help="Specify exclusion time in [ms].")
parser.add_argument("-b", "--bin", type=int, help="Specify the number of points over which to bin.")
args = parser.parse_args()

# Constants & Parameters
kB = 1.38E-23
omega = 2 * np.pi * 35
m = 127 * 1.6E-27
bin_width = 300


def two_body(p, t, y):
    return p[0] / (1 + p[0] * p[1] * t) - y


def exp(p, t, y):
    return p[0] * np.exp(-t / p[1]) + p[2] - y


def vector_bin(x, bin_size):
    n_vector = len(x)
    n_binned = int(n_vector / bin_size)
    result = np.zeros(n_binned)
    for k in range(n_binned):
        result[k] = np.sum(x[bin_size * k:bin_size * (k + 1)])
    return result/float(bin_size)


files = []
for f in args.files:
    for g in glob.glob(f):
        files.append(g)

if args.exclude:
    excludeTime = args.exclude
else:
    excludeTime = 0.01

if args.bin:
    bin_width = args.bin
else:
    bin_width = 300

figures = [0]*len(files)
n = 0
for k in files:
    print("File: " + k)
    data = np.loadtxt(k)

    if data.shape[1] == 5:
        nonEq = False
        time = data[:, 0] / 1000.
        temperature = data[:, 1] / kB * 1E9
        potential = data[:, 2] / kB * 1E9
        number = data[:, 3]
        collisions = data[:,4]
        N0 = number[0]
        total_collisions = data[-1, -1]

    elif data.shape[1] == 6:
        nonEq = True
        time = data[:, 0] / 1000.
        temperature = (data[:, 1] + data[:, 2]) / kB * 1E9
        tX = data[:, 1] / kB * 1E9 * 2.0
        tY = data[:, 2] / kB * 1E9 * 2.0
        potential = data[:, 3] / kB * 1E9
        number = data[:, 4]
        collisions = data[:, 5]
        N0 = number[0]
        total_collisions = data[-1, -1]

    time = vector_bin(time, bin_width)
    temperature = vector_bin(temperature, bin_width)
    if nonEq:
        tX = vector_bin(tX, bin_width)
        tY = vector_bin(tY, bin_width)
    potential = vector_bin(potential, bin_width)
    number = vector_bin(number, bin_width)
    collisions = vector_bin(collisions, bin_width)

    volume = 4 * np.pi * (np.sqrt(kB * temperature * 1E-9 / m) / omega) ** 2
    density = number / volume * 1E-4

    if nonEq:
        p0 = [200, 0.5, 300]
        pUpper = [1E6, np.inf, np.inf]
        pLower = [-1E6, 0.0, 0.0]
        resLSQX = least_squares(exp, p0, args=(time[time > excludeTime], tX[time > excludeTime]),
                                bounds=(pLower, pUpper))
        resLSQY = least_squares(exp, p0, args=(time[time > excludeTime], tY[time > excludeTime]),
                                bounds=(pLower, pUpper))
        T0X, tauX, TfX = resLSQX.x
        T0Y, tauY, TfY = resLSQY.x
        col2thermX = np.interp(tauX, time, collisions / number)
        col2thermY = np.interp(tauY, time, collisions / number)
        Xfit = exp([T0X, tauX, TfX], time, 0)
        Yfit = exp([T0Y, tauY, TfY], time, 0)
    else:
        p_linear_temperature = np.polyfit(time[time > excludeTime], temperature[time > excludeTime], 1)

    p_psd = np.polyfit(np.log(temperature[time > excludeTime]), np.log(number[time > excludeTime]), 1)
    if not nonEq:
        print('T0 = {:.1f} nK,  h = {:.1f} nK/s'.format(p_linear_temperature[1], p_linear_temperature[0]))

    if args.fit == '2b' or not args.fit:
        p0 = [4E7, 1E-7]
        pUpper = [np.inf, 1.0]
        pLower = [0, 0.0]
        resLSQ = least_squares(two_body, p0, args=(time, density), bounds=(pLower, pUpper))
        n0, beta = resLSQ.x
        print('N0 = {2}, n0 = {0:.2f}E7 cm^-2,  beta = {1:.2f}E7 cm^2 s^-1'.format(n0 / 1E7, beta * 1E7, int(N0)))
        print('Slope: {:.2f}'.format(p_psd[0]))
        density_plot_title = r'n$_0$ = {0:.2f} x 10$^7$ cm$^{{-2}}$,' \
                             r'  $\beta$ = {1:.2f} x 10$^{{-7}}$ cm$^{{-2}}$ s$^{{-1}}$'.format(n0 / 1E7, beta * 1E7)
        density_fit = two_body([n0, beta], time[time > excludeTime], 0)
    elif args.fit == 'e':
        p0 = [4E7, 0.5]
        pUpper = [np.inf, np.inf]
        pLower = [0, 0.0]
        resLSQ = least_squares(exp, p0, args=(time, density), bounds=(pLower, pUpper))
        n0, tau = resLSQ.x
        print('N0 = {2}, n0 = {0:.2f}E7 cm^-2,  tau = {1:.2f} s'.format(n0 / 1E7, tau, int(N0)))
        print('Slope: {:.2f}'.format(p_psd[0]))
        density_plot_title = r'n$_0$ = {0:.2f} x 10$^7$ cm$^{{-2}}$,' \
                             r'  $\tau$ = {1:.2f} s'.format(n0 / 1E7, tau)
        density_fit = exp([n0, tau], time[time > excludeTime], 0)

    print("Total number of collisions: {}".format(int(total_collisions)))
    if nonEq:
        print("Number of collisions to thermalize: ({:.1f}, {:.1f})".format(col2thermX, col2thermY))

    colRate = np.diff(collisions/number)/np.mean(np.diff(time))

    print("Mean Collision Rate: {:.1f}\n".format(np.mean(colRate)))


    figures[n] = plt.figure(figsize=(11.25, 7))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    ax0 = figures[n].add_subplot(231)
    ax1 = figures[n].add_subplot(232)
    ax2 = figures[n].add_subplot(234)
    ax3 = figures[n].add_subplot(235)

    ax4 = figures[n].add_subplot(233)
    ax5 = figures[n].add_subplot(236)

    all_axes = [ax0, ax1, ax2, ax3, ax4, ax5]
    all_ylabels = ['Number', 'Temperature (nK)', r'Density (10$^7$ cm$^{-2}$)', 'log(N)',
                   'Collisions per Particle', r'Collision Rate per particle (s$^{-1}$)']
    all_xlabels = ['Time (s)', 'Time (s)', 'Time (s)', 'log(T/nK)', 'Time (s)', 'Time (s)']
    all_titles = ['',
                  '',
                  density_plot_title,
                  'Slope: {:.2f}'.format(p_psd[0]),
                  'Collisions',
                  'Collision Rate']
    if nonEq:
        all_titles[1] = "Number of collisions to thermalize: {:.1f}\n".format(min(col2thermX, col2thermY))
    else:
        all_titles[1] = r'T$_0$ = {:.1f} nK,  h = {:.1f} nK/s'.format(p_linear_temperature[1], p_linear_temperature[0])

    ax0.plot(time, number, '.', markersize=15)

    if nonEq:
        ax1.plot(time, tX, '.b', markersize=15)
        ax1.plot(time, tY, '.r', markersize=15)
        ax1.plot(time, Xfit, color='blue')
        ax1.plot(time, Yfit, color='red')
        ax1.legend([r'T$_X$', r'T$_Y$'])
    else:
        ax1.plot(time, temperature, '.', markersize=15)
        ax1.plot(time[time > excludeTime], p_linear_temperature[0] * time[time > excludeTime] + p_linear_temperature[1], 'r')

    ax2.plot(time, density / 1E7, '.', markersize=15)
    ax2.plot(time[time > excludeTime], density_fit/1E7, 'r--')

    ax3.plot(np.log(temperature), np.log(number), '.', markersize=15)
    ax3.plot(np.log(temperature[time > excludeTime]), p_psd[1] + p_psd[0] * np.log(temperature[time > excludeTime]), 'r')

    ax4.plot(time, collisions/number, '.', markersize=15)
    if nonEq:
        ax4.annotate("Collision to thermalize: {:.1f}".format(min(col2thermX, col2thermY)),
                     xy=(0.05, 0.85), xycoords='axes fraction')

    ax5.plot(time[0:-1], colRate, '.', markersize=15)

    j = 0
    for a in all_axes:
        a.tick_params(axis='both', which='both', right=True, top=True, direction='in')
        a.set_ylabel(all_ylabels[j])
        a.set_xlabel(all_xlabels[j])
        a.set_title(all_titles[j])

        j += 1

    figures[n].suptitle(k)
    n += 1

plt.show()
