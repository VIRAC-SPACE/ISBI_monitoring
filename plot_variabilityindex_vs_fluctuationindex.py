import argparse
import os
import sys
from functools import reduce
from datetime import datetime

import numpy as np
from matplotlib import pyplot as plt
from astropy.time import Time
from scipy.stats import linregress
import pprint

from parsers.configparser_ import ConfigParser


def get_configs(section, key, config_file_path):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def get_configs_items(section, config_file_path):
    """

    :return: None
    """
    config = ConfigParser(config_file_path)
    return config.get_items(section)


def find_nearest_index(array, value):
    """

    :param array: array
    :param value: value to search
    :return: nearest index of value for array
    """
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index

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

def compute_rms(amp):
    N = len(amp)
    ms = 0
    for i in range(0, N):
        ms = ms + amp[i] ** 2

    ms = ms / N
    rms = np.sqrt(ms)
    return rms

def main(config, config_plot):
    plt.style.use(config_plot)
    dtype = np.dtype([("date", "S12"), ("amp", float), ("error", float)])
    format = "%d-%b-%Y"
    fluctuation_indexes = []
    variability_indexes = []

    titles = []

    all_sources = [s.split("_")[0] for s in list(get_configs_items("velocities", config).keys())]
    print(all_sources)

    for source in all_sources:
        source_dir = get_configs("paths", "monitoring_path", config)+ "/" + source.upper() + "/"

        if  not os.path.isdir(source_dir):
            continue

        cont_data = source_dir + "UVFIT_" + source.upper() + "_.txt"
        line_data = source_dir + "line/"

        if os.path.isfile(cont_data):
            print("Process cont data for source " + source.upper())
            titles.append(source + "_cont")

            contiuum_date, contiuum_amp, contiuum_amp_error = np.loadtxt(cont_data,
                                                                         usecols=(0, 2, 3), unpack=True, dtype=dtype)
            contiuum_date = [str(date).replace("b", "").replace("'", "") for date in contiuum_date]
            format = "%d-%b-%Y"
            mjd = Time([datetime.strptime(date, format) for date in contiuum_date]).mjd

            x = mjd
            y = contiuum_amp
            N = len(contiuum_amp)

            error = contiuum_amp_error

            variability_index = ((np.max(y) - error[list(y).index(np.max(y))]) - (
                    np.min(y) + error[list(y).index(np.min(y))])) \
                                / ((np.max(y) - error[list(y).index(np.max(y))]) + (
                    np.min(y) + error[list(y).index(np.min(y))]))

            fluctuation_index = np.sqrt(
                np.abs((N / reduce(lambda x_, y_: x_ + y_, [error[list(y).index(i)] ** 2 for i in y])) *
                       ((reduce(lambda x_, y_: x_ + y_,
                                [i ** 2 * (error[list(y).index(i)]) ** 2 for i in y]) -
                         np.mean(y) * reduce(lambda x_, y_: x_ + y_,
                                             [i * error[list(y).index(i)] ** 2 for i in y]))
                        / (N - 1)) - 1)) / np.mean(y)

            variability_indexes.append(np.float64(variability_index))
            fluctuation_indexes.append(np.float64(fluctuation_index))

            print(source, variability_index, fluctuation_index)

        if len(os.listdir(line_data)) > 0:
            print("Process line data for source " + source.upper())
            monitoring_path = get_configs("paths", "monitoring_path", config) + "/" + source.upper() + "/line/"
            monitoring_files = [file for file in os.listdir(monitoring_path)]

            components = sorted([float(component) for component in
                          get_configs("velocities", source.lower() + "_" + "6668", config).split(",")])

            amp_for_component = {component: [] for component in components}
            error_for_component = {component: [] for component in components}
            for file in monitoring_files:
                vel, amp = read_spectrum_file(monitoring_path + file)

                max_index = (np.abs(vel - np.min(components)).argmin()) +10
                min_index = (np.abs(vel - np.max(components)).argmin()) -10

                amp_tmp = []
                amp_tmp.extend(amp[0:min_index])
                amp_tmp.extend(amp[max_index:-1])
                rms = compute_rms(amp_tmp)

                error = []
                for i in range(0, len(amp)):
                    error.append(2 * rms + amp[i] * 0.05)

                for component in components:
                    max_index = (np.abs(vel - component).argmin())
                    max_indexs = range(max_index - 1, max_index + 1)
                    max_amplitudes = []
                    for index in max_indexs:
                        max_amplitudes.append(amp[index])

                    amp_for_component[component].append(np.max(max_amplitudes))
                    error_for_component[component].append(error[max_index])

            for component in components:
                titles.append(source + "_" + str(component))
                y = amp_for_component[component]
                N = len(y)

                error = error_for_component[component]

                variability_index = ((np.max(y) - error[list(y).index(np.max(y))]) - (
                        np.min(y) + error[list(y).index(np.min(y))])) \
                                    / ((np.max(y) - error[list(y).index(np.max(y))]) + (
                        np.min(y) + error[list(y).index(np.min(y))]))

                fluctuation_index = np.sqrt(
                    np.abs((N / reduce(lambda x_, y_: x_ + y_, [error[list(y).index(i)] ** 2 for i in y])) *
                           ((reduce(lambda x_, y_: x_ + y_,
                                    [i ** 2 * (error[list(y).index(i)]) ** 2 for i in y]) -
                             np.mean(y) * reduce(lambda x_, y_: x_ + y_,
                                                 [i * error[list(y).index(i)] ** 2 for i in y]))
                            / (N - 1)) - 1)) / np.mean(y)

                variability_indexes.append(np.float64(variability_index))
                fluctuation_indexes.append(np.float64(fluctuation_index))

                print(source, component, variability_index, fluctuation_index)


    scatter = plt.scatter(variability_indexes, fluctuation_indexes, alpha=0.3)
    for plot in range(0, len(titles)):
        if np.abs(variability_indexes[plot]) > 5:
            plt.text(variability_indexes[plot], fluctuation_indexes[plot], titles[plot])

    coef = linregress(variability_indexes, fluctuation_indexes)
    plt.plot(variability_indexes, np.array(variability_indexes) * coef.slope + coef.intercept, '--k')
    pprint.pprint(coef, indent=4)

    #handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
    #ranges = ["0.5 < Jy <= 20", "20 < Jy <= 200", "200 < Jy <= 800", "800 < Jy <= 3000"]
    #labels = [labels[l] + "  " + ranges[l] for l in range(0, len(ranges))]

    #plt.legend(handles, labels)
    plt.xlabel("Variability index")
    plt.ylabel("Fluctuation index")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''variabilityindex and  fluctuationindex''',
                                     epilog="plot variabilityindex vs fluctuationindex")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str,
                        default="config/config.cfg")
    parser.add_argument("-cp", "--config_plot", help="Matplotlib configuration cfg file", type=str,
                        default="config/plot.style")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    main(args.config, args.config_plot)
    sys.exit(0)
