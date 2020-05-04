import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import argparse
import glob
import sys

parser = argparse.ArgumentParser(description='Process output files.')
parser.add_argument("files", type=str, help="List of file names", nargs='*')
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

    time = data[:, 0] / 1000.
    kineticX = data[:, 1] / kB * 1E9 * 2
    kineticY = data[:, 2] / kB * 1E9 * 2
    collisions = data[:,-1]
    N = np.mean(data[:,4])
    total_collisions = data[-1, -1]

    time = vector_bin(time, bin_width)
    kineticX = vector_bin(kineticX, bin_width)
    kineticY = vector_bin(kineticY, bin_width)
    collisions = vector_bin(collisions, bin_width)


    p0 = [200, 0.5, 300]
    pUpper = [1E6, np.inf, np.inf]
    pLower = [-1E6, 0.0, 0.0]
    resLSQX = least_squares(exp, p0, args=(time[time > excludeTime], kineticX[time > excludeTime]), bounds=(pLower, pUpper))
    resLSQY = least_squares(exp, p0, args=(time[time > excludeTime], kineticY[time > excludeTime]), bounds=(pLower, pUpper))
    T0X, tauX, TfX = resLSQX.x
    T0Y, tauY, TfY = resLSQY.x

    print('TauX = {0:.2f} s, TauY = {1:.2f} s'.format(tauX, tauY))

    density_plot_title = r'$\tau_X$ = {0:.2f} s, $\tau_Y$ = {0:.2f} s'.format(tauX, tauY)
    Xfit = exp([T0X, tauX, TfX], time, 0)
    Yfit = exp([T0Y, tauY, TfY], time, 0)

    print("Total number of collisions: {}\n".format(int(total_collisions)))

    figures[n] = plt.figure(figsize=(8, 4))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    ax0 = figures[n].add_subplot(121)
    ax1 = figures[n].add_subplot(122)


    all_axes = [ax0, ax1]
    all_ylabels = ['Temperature', 'Collisions per particle']
    all_xlabels = ['Time (s)', 'Time (s)']
    all_titles = ['Cross-dimensional Thermalization', 'Elastic Collisions']

    ax0.plot(time, kineticX, '.', markersize=15, color='blue')
    ax0.plot(time, kineticY, '.', markersize=15, color='red')
    ax0.plot(time, Xfit, color='blue')
    ax0.plot(time, Yfit, color='red')

    ax1.plot(time, collisions/N, '.', markersize=15)

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
