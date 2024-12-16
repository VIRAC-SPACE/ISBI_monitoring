#! /usr/bin/python3
# -*- coding: utf-8 -*-
import sys
import os

import argparse
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from astropy.time import Time

from parsers.configparser_ import ConfigParser


def read_spectrum_file(spectrum_file):
    vel = []
    amp = []

    with open(spectrum_file) as f:
        lines = f.readlines()
        for index in range(14, len(lines)):
            line = lines[index]
            if "FLAGGED" in line:
                continue
            else:
                line = lines[index]
                if "I" in line:
                    data_tmp = line.split()
                    data_tmp = [d.strip() for d in data_tmp]
                    vel.append(data_tmp[4])
                    amp.append(data_tmp[5])

    return np.array(vel, dtype=float), np.array(amp, dtype=float)


def get_date(spectrum_file):
    date = 0
    with open(spectrum_file) as f:
        lines = f.readlines()
        for index in range(0, len(lines)):
            line = lines[index]
            if "OBS. DATE:" in line:
                date = line.split()[2]
                break

    year = date[0:4]
    month = date[4:6]
    day = date[6:8]
    return year + "-" + month + "-" + day


def get_configs(section, key, config_file_path):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def main(source, config, config_plot):
    plt.style.use(config_plot)
    monitoring_path = get_configs("paths", "monitoring_path", config) + "/" + source + "/"
    monitoring_files = [file for file in os.listdir(monitoring_path) ]
    print("monitoring files:", monitoring_files)

    components = [float(component) for component in
                  get_configs("velocities", source.lower() + "_" + "6668", config).split(",")]
    amp_for_component = {component: [] for component in components}

    dates = [get_date(monitoring_path + file) for file in monitoring_files]
    dates = Time(dates, format='isot', scale='utc')
    mjd = dates.mjd

    fig = plt.figure(figsize=(16, 16), dpi=100)
    ax = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    for file in monitoring_files:
        vel, amp = read_spectrum_file(monitoring_path + file)

        time = [mjd[monitoring_files.index(file)]] * len(vel)
        ax2.plot(time, vel, amp)
        for component in components:
            max_index = (np.abs(amp - component).argmin())
            max_indexs = range(max_index-3, max_index+3)
            max_amplitudes = []
            for index in max_indexs:
                max_amplitudes.append(amp[index])
            amp_for_component[component].append(np.max(max_amplitudes))

    for component in components:
        ax.scatter(mjd, amp_for_component[component], label="Maser flux [Jy] vel: " + str(component) + " [KM/s]")

    '''
    #ax2.xaxis._axinfo['label']['space_factor'] = 1000
    #ax2.yaxis._axinfo['label']['space_factor'] = 1000
    #ax2.zaxis._axinfo['label']['space_factor'] = 1000
    #ax2.dist = 10
    '''
    #ticks = ax.get_xticks()
    #print("ticks", ticks[6])
    #print("dates", dates)
    #print(len(dates), len(ticks))
    #ax2.set_xticks(ticks, dates)
    #ax2.set_xticklabels(datestmp)
    #ax.set_xticks(ticks, dates)
    #ax.set_xticklabels(dates, rotation=30)

    ax2.set_xlabel("Observation dates", labelpad=14)
    ax2.set_ylabel("Velocity [km/s]", labelpad=14)
    ax2.set_zlabel("Flux [Jy]", labelpad=10)
    ax2.grid(True)
    ax.set_xlabel("Observation Date")
    ax.set_ylabel("Flux [Jy]")
    ax.grid(True)
    ax.legend()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Display Monitoring. ''', epilog="""plot monitoring.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str,
                        default="config/config.cfg")
    parser.add_argument("-cp", "--config_plot", help="Matplotlib configuration cfg file", type=str,
                        default="config/plot.style")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    main(args.source, args.config, args.config_plot)
    sys.exit(0)
