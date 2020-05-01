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
    return p[0] * np.exp(-t / p[1]) - y


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
    excludeTime = 0.1

if args.bin:
    bin_width = args.bin
else:
    bin_width = 300

figures = [0]*len(files)
n = 0
for k in files:
    print("File: " + k)
    data = np.loadtxt(k)

    time = data[:, 0] / 1000.
    kinetic = data[:, 1] / kB * 1E9
    potential = data[:, 2] / kB * 1E9
    number = data[:, 3]
    N0 = number[0]
    total_collisions = data[-1, -1]

    time = vector_bin(time, bin_width)
    kinetic = vector_bin(kinetic, bin_width)
    potential = vector_bin(potential, bin_width)
    number = vector_bin(number, bin_width)

    temperature = (kinetic + potential) / 2.0
    volume = 4 * np.pi * (np.sqrt(kB * temperature * 1E-9 / m) / omega) ** 2
    density = number / volume * 1E-4

    p_linear_temperature = np.polyfit(time[time > excludeTime], kinetic[time > excludeTime], 1)
    p_psd = np.polyfit(np.log(temperature[time > excludeTime]), np.log(number[time > excludeTime]), 1)
    print('T0 = {:.1f} nK,  h = {:.1f} nK/s'.format(p_linear_temperature[1], p_linear_temperature[0]))

    if args.fit == '2b' or not args.fit:
        p0 = [4E7, 1E-7]
        pUpper = [np.inf, 1.0]
        pLower = [0, 0.0]
        resLSQ = least_squares(two_body, p0, args=(time, density), bounds=(pLower, pUpper))
        n0, beta = resLSQ.x
        print('N0 = {2}, n0 = {0:.2f}E7 cm^-2,  beta = {1:.2f}E7 cm^2 s^-1'.format(n0 / 1E7, beta * 1E7, int(N0)))
        print('Slope: {:.2f}\n'.format(p_psd[0]))
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

    print("Total number of collisions: {}\n".format(int(total_collisions)))

    figures[n] = plt.figure(figsize=(7.5, 6))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    ax0 = figures[n].add_subplot(221)
    ax1 = figures[n].add_subplot(222)
    ax2 = figures[n].add_subplot(223)
    ax3 = figures[n].add_subplot(224)

    all_axes = [ax0, ax1, ax2, ax3]
    all_ylabels = ['Number', 'Temperature (nK)', r'Density (10$^7$ cm$^{-2}$)', 'log(N)']
    all_xlabels = ['Time (s)', 'Time (s)', 'Time (s)', 'log(T/nK)']
    all_titles = ['',
                  r'T$_0$ = {:.1f} nK,  h = {:.1f} nK/s'.format(p_linear_temperature[1], p_linear_temperature[0]),
                  density_plot_title,
                  'Slope: {:.2f}'.format(p_psd[0])]

    ax0.plot(time, number, '.', markersize=15)

    ax1.plot(time, temperature, '.', markersize=15)
    ax1.plot(time, kinetic, 'g.', markersize=15)
    ax1.plot(time[time > excludeTime], p_linear_temperature[0] * time[time > excludeTime] + p_linear_temperature[1], 'r')
    ax1.legend(['Temperature', 'Kinetic Temperature'])

    ax2.plot(time, density / 1E7, '.', markersize=15)
    ax2.plot(time[time > excludeTime], density_fit/1E7, 'r--')

    ax3.plot(np.log(temperature), np.log(number), '.', markersize=15)
    ax3.plot(np.log(temperature[time > excludeTime]), p_psd[1] + p_psd[0] * np.log(temperature[time > excludeTime]), 'r')

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
