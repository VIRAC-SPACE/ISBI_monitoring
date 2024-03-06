import sys
import os

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from astropy.time import Time


def read_spectrum_file(spectrum_file):
    data = []
    with open(spectrum_file) as f:
        lines = f.readlines()
        for index in range(0, len(lines)):
            if index < 14:
                continue
            else:
                line = lines[index]
                data_tmp = line.split()
                data_tmp = [d.strip() for d in data_tmp]
                data.append(data_tmp)

    data = np.array(data)
    velocity = np.array(data[:, 4], dtype=float)
    amplitude = np.array(data[:, 5], dtype=float)
    amplitude -= np.mean(amplitude[0:10])
    return velocity, amplitude


def main():
    font = {'family': 'normal',

            'size': 16}

    matplotlib.rc('font', **font)

    data_path = "/home/janis/PycharmProjects/test/ToJanis/ToJanis/"
    spectrum_files = [file for file in os.listdir(data_path) if not file.startswith("contcurve_")]
    monitoring_files = [file for file in os.listdir(data_path) if file.startswith("contcurve_")]
    monitoring_data = [np.loadtxt(data_path + file, dtype=object) for file in monitoring_files]
    dates = monitoring_data[1][:, 0]
    dates = [str(d) for d in dates]
    dates = Time(dates, format='isot', scale='utc')
    mjd = dates.mjd
    contcurve_3C345_flux = np.array(monitoring_data[0][:, 2], dtype=np.float)
    contcurve_G85_flux = np.array(monitoring_data[1][:, 2], dtype=np.float)

    components = [-29.3992, -31.5065]
    c_29_amplitudes = []
    c_31_amplitudes = []
    for file in spectrum_files:
        spectrum = read_spectrum_file(data_path + file)
        max_indexs = []
        for component in components:
            max_index = (np.abs(spectrum[0] - component).argmin())
            max_indexs = range(max_index-5, max_index+5)
            max_amplitudes = []
            for index in max_indexs:
                max_amplitudes.append(spectrum[1][index])
            if components.index(component) == 0:
                c_29_amplitudes.append(np.max(max_amplitudes))
            else:
                c_31_amplitudes.append(np.max(max_amplitudes))

    fig = plt.figure(figsize=(16, 16), dpi=100)
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(mjd[1:len(mjd)], c_29_amplitudes[1:len(mjd)])
    ax.scatter(mjd[1:len(mjd)], c_29_amplitudes[1:len(mjd)], label="Maser flux [Jy] vel: -29.3992 [KM/s]")

    ax.plot(mjd[1:len(mjd)], c_31_amplitudes[1:len(mjd)])
    plt.scatter(mjd[1:len(mjd)], c_31_amplitudes[1:len(mjd)], label="Maser flux [Jy] vel: -31.5065 [KM/s]")

    ax.plot(mjd[1:len(mjd)], contcurve_G85_flux[1:len(mjd)] *1000)
    ax.scatter(mjd[1:len(mjd)], contcurve_G85_flux[1:len(mjd)] *1000 , label="Contcurve G85 flux [Jy]")

    ax.legend()

    datestmp = monitoring_data[1][:, 0]
    datestmp = [str(d) for d in datestmp]
    datestmp = datestmp[1:len(datestmp)]
    ticks = ax.get_xticks()
    ax.set_xticks(ticks, datestmp)
    ax.set_xticklabels(datestmp, rotation=30)

    ax.set_xlabel("Observation Date")
    ax.set_ylabel("Flux [Jy]")
    ax.grid(True)

    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    spectrum_files = spectrum_files[1:len(spectrum_files)]
    mjd = mjd[1:len(mjd)]
    for file in spectrum_files:
        spectrum = read_spectrum_file(data_path + file)
        velocity = spectrum[0]
        amplitude = spectrum[1]
        time = [mjd[spectrum_files.index(file)]] * len(velocity)
        ax2.plot(time, velocity, amplitude)

    ticks = ax.get_xticks()
    ax2.set_xticks(ticks, datestmp)
    ax2.set_xticklabels(datestmp)
    #print(dir(ax2))
    ax2.set_xlabel("Observation dates", labelpad=14)
    ax2.set_ylabel("Velocity [km/s]", labelpad=14)
    ax2.set_zlabel("Flux [Jy]", labelpad=10)
    #ax2.xaxis._axinfo['label']['space_factor'] = 1000
    #ax2.yaxis._axinfo['label']['space_factor'] = 1000
    #ax2.zaxis._axinfo['label']['space_factor'] = 1000
    #ax2.dist = 10

    ax2.grid(True)
    plt.show()


if __name__ == "__main__":
    main()
    sys.exit()
