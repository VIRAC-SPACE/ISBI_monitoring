#! /usr/bin/python3
# -*- coding: utf-8 -*-
import os
import sys

import argparse
import numpy as np
import matplotlib.pyplot as plt

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


def main(source_name, obs_name, config, config_plot):
    plt.style.use(config_plot)
    print("obs_name:", obs_name)
    monitoring_path = get_configs("paths", "monitoring_path", config) + "/"

    monitoring_file_name = [f for f in os.listdir(monitoring_path + "/" + source_name + "/line/") if f.startswith(obs_name.upper()) and f.endswith("X.txt")]
    if len(monitoring_file_name) == 0:
        print("No such observation")
        sys.exit()

    file = monitoring_path + "/" + source_name + "/line/" + monitoring_file_name[0]
    print("file:", file)

    spectrum = read_spectrum_file(file)
    plt.plot(spectrum[0], spectrum[1], label=source_name + "_" + obs_name)

    plt.xlabel("Velocity [km/s]", labelpad=14)
    plt.ylabel("Flux [Jy]", labelpad=10)
    plt.legend()

    top = 0.979
    bottom = 0.099
    left = 0.066
    right = 0.975
    hspace = 0.2
    wspace = 0.2

    plt.tight_layout()
    plt.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Plot a single file. ''', epilog="""plot single file.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("exper", help="Experiment name", type=str)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str,
                        default="config/config.cfg")
    parser.add_argument("-cp", "--config_plot", help="Matplotlib configuration cfg file", type=str,
                        default="config/plot.style")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    main(args.source, args.exper, args.config, args.config_plot)
    sys.exit(0)
