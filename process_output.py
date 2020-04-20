import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares


def two_body(p, t, y):
    return p[0] / (1 + p[0] * p[1] * t) - y


kB = 1.38E-23

omega = 2 * np.pi * 35
m = 127 * 1.6E-27

filePath = './results/'
files = ['NoEvaporation_20200417_00.out']
excludeTime = 0.04

for k in files:
    data = np.loadtxt(filePath + k)
    time = data[:, 0] / 1000.
    kinetic = data[:, 1] / kB * 1E9
    potential = data[:, 2] / kB * 1E9
    number = data[:, 3]

    temperature = (kinetic + potential) / 2.0
    volume = 4 * np.pi * (np.sqrt(kB * temperature * 1E-9 / m) / omega) ** 2
    density = number / volume * 1E-4

    p_linear_temperature = np.polyfit(time[time > excludeTime], kinetic[time > excludeTime], 1)
    p_psd = np.polyfit(np.log(temperature[time > excludeTime]), np.log(number[time > excludeTime]), 1)

    p0 = [4E7, 1E-7]
    pUpper = [np.inf, 1.0]
    pLower = [0, 0.0]
    resLSQ = least_squares(two_body, p0, args=(time, density), bounds=(pLower, pUpper))
    n0, beta = resLSQ.x

    fig0 = plt.figure(figsize=(10, 8))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    ax0 = fig0.add_subplot(221)
    ax1 = fig0.add_subplot(222)
    ax2 = fig0.add_subplot(223)
    ax3 = fig0.add_subplot(224)

    all_axes = [ax0, ax1, ax2, ax3]
    all_ylabels = ['Number', 'Temperature (nK)', r'Density (10$^7$ cm$^{-2}$)', 'log(N)']
    all_xlabels = ['Time (s)', 'Time (s)', 'Time (s)', 'log(T/nK)']
    all_titles = ['',
                  r'T$_0$ = {:.1f} nK,  h = {:.1f} nK/s'.format(p_linear_temperature[1], p_linear_temperature[0]),
                  r'n$_0$ = {0:.2f} x 10$^7$ cm$^{{-2}}$,  $\beta$ = {1:.2f} x 10$^{{-7}}$ cm$^{{-2}}$ s$^{{-1}}$'
                      .format(n0 / 1E7, beta * 1E7),
                  'Slope: {:.2f}'.format(p_psd[0])]

    ax0.plot(time, number)

    ax1.plot(time, temperature)
    ax1.plot(time, kinetic, 'g')
    ax1.plot(time[time > excludeTime], p_linear_temperature[0] * time[time > excludeTime] + p_linear_temperature[1], 'r')
    ax1.legend(['Temperature', 'Kinetic Temperature'])

    ax2.plot(time, density / 1E7)
    ax2.plot(time[time > excludeTime], two_body([n0, beta], time[time > excludeTime], 0) / 1E7, 'r--')

    ax3.plot(np.log(temperature), np.log(number))
    ax3.plot(np.log(temperature[time > excludeTime]), p_psd[1] + p_psd[0] * np.log(temperature[time > excludeTime]), 'r')

    j = 0
    for a in all_axes:
        a.tick_params(axis='both', which='both', right=True, top=True, direction='in')
        a.set_ylabel(all_ylabels[j])
        a.set_xlabel(all_xlabels[j])
        a.set_title(all_titles[j])

        j += 1

    fig0.suptitle(k)

plt.show()
