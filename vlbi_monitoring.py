#! /usr/bin/python3
# -*- coding: utf-8 -*-
import sys
import os
from datetime import datetime

import argparse
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy.io import ascii

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
    result_path = get_configs("paths", "monitoring_path", config) + "/" + source.upper() + "/results/"
    monitoring_path = get_configs("paths", "monitoring_path", config) + "/" + source.upper() + "/line/"

    print("monitoring_path", monitoring_path)

    if not os.path.isdir(result_path):
        os.system("mkdir -p " + result_path)

    contiuum_data = get_configs("paths", "monitoring_path", config) + "/" + source.upper() + "/UVFIT_" + source + "_.txt"
    dtype = np.dtype([("date", "S12"), ("amp", float), ("error", float)])

    print("contiuum_data", contiuum_data)

    if os.path.isfile(contiuum_data):
        contiuum_date, contiuum_amp, contiuum_amp_error = np.loadtxt(contiuum_data,
                                                                     usecols=(0, 2, 3), unpack=True, dtype=dtype)
        contiuum_date = [str(date).replace("b", "").replace("'", "") for date in contiuum_date]
        format = "%d-%b-%Y"
        mjd_cont = Time([datetime.strptime(date, format) for date in contiuum_date]).mjd

    if os.path.isdir(monitoring_path):
        fig = plt.figure(figsize=(16, 16), dpi=100)
        ax = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')

        monitoring_files = [file for file in os.listdir(monitoring_path)]
        print("monitoring files:", monitoring_files)

        components = [float(component) for component in
                      get_configs("velocities", source.lower() + "_" + "6668", config).split(",")]
        amp_for_component = {component: [] for component in components}

        dates = [get_date(monitoring_path + file) for file in monitoring_files]
        dates = Time(dates, format='isot', scale='utc')
        mjd_line = dates.mjd

        if os.path.isfile(contiuum_data):
            ax.scatter(mjd_cont, contiuum_amp*1000, label="contiuum * 1000")
            ax2.scatter(mjd_cont, [0] * len(mjd_cont), contiuum_amp * 1000, label="contiuum * 1000")

        sources_vrange = ascii.read('DB_vrange.csv')
        source_vrange_index = sources_vrange['name'].tolist().index(source.upper())
        vmin = dict(sources_vrange)["vmin"][source_vrange_index]
        vmax = dict(sources_vrange)["vmax"][source_vrange_index]
        integrate_flux = []

        for file in monitoring_files:
            vel, amp = read_spectrum_file(monitoring_path + file)

            vmin_index = (np.abs(vel - vmin)).argmin()
            vmax_index = (np.abs(vel - vmax)).argmin()

            integrate_flux.append(np.trapezoid(amp[vmin_index:vmax_index], vel[vmin_index:vmax_index]))

            time = [mjd_line[monitoring_files.index(file)]] * len(vel)
            ax2.plot(time, vel, amp)
            for component in components:
                max_index = (np.abs(vel - component).argmin())
                max_indexs = range(max_index-1, max_index+1)
                max_amplitudes = []
                for index in max_indexs:
                    max_amplitudes.append(amp[index])
                amp_for_component[component].append(np.max(max_amplitudes))

        for component in components:
            ax.scatter(mjd_line, amp_for_component[component],
                       label="Maser flux [Jy] vel: " + str(component) + " [KM/s]", alpha=0.9)

        fig5, ax5 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
        ax5.scatter(mjd_line, integrate_flux)

        ax2.set_xlabel("Observation dates", labelpad=14)
        ax2.set_ylabel("Velocity [km/s]", labelpad=14)
        ax2.set_zlabel("Flux [Jy]", labelpad=10)
        ax2.grid(True)
        ax.set_xlabel("Observation Date")
        ax.set_ylabel("Flux [Jy]")
        ax.grid(True)
        ax.legend()

        top = 0.998
        bottom = 0.075
        left = 0.064
        right = 0.98
        hspace = 0.0
        wspace = 0.0
        fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace)
        fig.savefig(result_path + source.lower() + "_timeseries_and_3D_plot", format="png")

    if os.path.isfile(contiuum_data):
        fig3, ax3= plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
        ax3.scatter(mjd_cont, contiuum_amp)

        for i in range(0, len(mjd_cont)):
            ax3.errorbar(mjd_cont[i], contiuum_amp[i], yerr=contiuum_amp_error[i], fmt='', ecolor="r")

        ax3.set_xlabel("MJD")
        ax3.set_ylabel(r'$Flux~(\mathrm{Jy})$')

        top = 0.993
        bottom = 0.12
        left = 0.069
        right = 0.995
        hspace = 0.0
        wspace = 0.0
        fig3.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace)
        fig3.savefig(result_path + source.lower() + "_cont_flux_vs_time", format="png")

        fig4, ax4 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
        ax4.scatter(contiuum_amp, contiuum_amp_error)
        ax4.set_xlabel("contiuum amp " + r'$Flux~(\mathrm{Jy})$')
        ax4.set_ylabel('contiuum amp error ' + r'$Flux~(\mathrm{Jy})$')

        top = 0.993
        bottom = 0.12
        left = 0.084
        right = 0.995
        hspace = 0.0
        wspace = 0.0
        fig4.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace)
        fig4.savefig(result_path + source.lower() + "_cont_error_vs_flux", format="png")

    #plt.show()


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
