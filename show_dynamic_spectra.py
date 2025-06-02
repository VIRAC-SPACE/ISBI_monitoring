import argparse
from datetime import datetime
import os
import sys

import numpy as np
import matplotlib.tri as mtri
from matplotlib import pyplot as plt
from matplotlib import ticker
from astropy.io import ascii
from astropy.time import Time

from parsers.configparser_ import ConfigParser


def get_configs(section, key, config_file_path):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


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


def get_max_min_velocity(vmin, vmax, pipeline_output_files, source, config):
    max_velocitys = []
    min_velocitys = []
    velocity_ = []
    observed_flux = []
    observed_time = []

    monitoring_path = get_configs("paths", "monitoring_path", config) + "/" + source + "/line/"
    for file in pipeline_output_files:
        velocity, amplitude = read_spectrum_file(monitoring_path + file)
        max_velocitys.append(max(velocity))
        min_velocitys.append(min(velocity))

        date = get_date(monitoring_path + file)
        date = Time(date, format='isot', scale='utc')
        mjd = date.mjd

        for i in range(0, len(velocity)):
            if vmin or vmax:
                if vmin <= velocity[i] <= vmax:
                    velocity_.append(velocity[i])
                    observed_flux.append(amplitude[i])
                    observed_time.append(mjd)
                else:
                    velocity_.append(velocity[i])
                    observed_flux.append(amplitude[i])
                    observed_time.append(mjd)

    return max(max_velocitys), min(min_velocitys), velocity_, observed_flux, observed_time


def main(source, config, config_plot):
    plt.style.use(config_plot)
    result_path = get_configs("paths", "monitoring_path", config) + "/" + source + "/results/"
    monitoring_path = get_configs("paths", "monitoring_path", config) + "/" + source + "/line/"
    monitoring_files = [file for file in os.listdir(monitoring_path)]
    print("monitoring files:", monitoring_files)

    if not os.path.isdir(result_path):
        os.mkdir(result_path)

    sources_vrange = ascii.read('DB_vrange.csv')
    source_vrange_index = sources_vrange['name'].tolist().index(source)
    vmin = dict(sources_vrange)["vmin"][source_vrange_index]
    vmax = dict(sources_vrange)["vmax"][source_vrange_index]

    if vmin is None or vmax is None:
        vmin, vmax, velocity, observed_flux, observed_time = get_max_min_velocity(vmin, vmax, monitoring_files, source, config)

    else:
        _, _, velocity, observed_flux, observed_time = get_max_min_velocity(vmin, vmax, monitoring_files, source, config)

    observed_flux = list((np.array(observed_flux).clip(min=1.5)))
    triang = mtri.Triangulation(observed_time, velocity)

    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    lvls = np.linspace(int(np.min(observed_flux)), int(np.max(observed_flux)), 1000)

    cs = ax1.tricontourf(triang, observed_flux, levels=lvls, antialiased=False, locator=ticker.LogLocator, cmap="jet")

    cs.set_clim(vmin=1.5)
    cbar = plt.colorbar(cs, spacing="proportional", label=r'$Flux~(\mathrm{Jy})$', extendrect=False)
    cbar.locator = ticker.LogLocator()

    ax1.set_ylabel('Velocity (km sec$^{-1}$)')
    ax1.set_xlabel("MJD")

    contiuum_data = (get_configs("paths", "monitoring_path", config)
                     + "/" + source + "/UVFIT_" + source + "_.txt")
    dtype=np.dtype([("date", "S12"), ("amp", float), ("error", float)])
    contiuum_date, contiuum_amp, contiuum_amp_error = np.loadtxt(contiuum_data,
                                                                 usecols=(0, 2, 3), unpack=True, dtype=dtype)
    contiuum_date = [str(date).replace("b", "").replace("'", "") for date in contiuum_date]
    format = "%d-%b-%Y"
    mjd = Time([datetime.strptime(date, format) for date in contiuum_date]).mjd

    ax2 = ax1.twinx()
    ax2.scatter(mjd, contiuum_amp, c="r")

    top = 0.993
    bottom = 0.115
    left = 0.084
    right = 1.0
    hspace = 0.0
    wspace = 0.0
    plt.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace)
    #plt.show()
    plt.savefig(result_path + source.lower() + "_dynamic_spectra")


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