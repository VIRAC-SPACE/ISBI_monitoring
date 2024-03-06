#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys

import argparse
import numpy as np
import matplotlib.pyplot as plt

from parsers.configparser_ import ConfigParser


def get_configs(section, key):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)

def read_spectrum_file(spectrum_file):
    vel_rr = []
    vel_ll = []
    amp_rr = []
    amp_ll = []

    with open(spectrum_file) as f:
        lines = f.readlines()
        for index in range(0, len(lines)):
            line = lines[index]
            if "FLAGGED" in line:
                continue
            else:
                line = lines[index]
                if "RR" in line or "LL" in line:
                    data_tmp = line.split()
                    data_tmp = [d.strip() for d in data_tmp]

                    pol = data_tmp[2]
                    if pol == "RR":
                        vel_rr.append(data_tmp[4])
                        amp_rr.append(data_tmp[5])
                    elif pol == "LL":
                        vel_ll.append(data_tmp[4])
                        amp_ll.append(data_tmp[5])

    return (np.array(vel_rr, dtype=float), np.array(vel_ll, dtype=float),
            np.array(amp_rr, dtype=float), np.array(amp_ll, dtype=float))


def main(source_name, obs_name):
    monitoring_path = get_configs("paths", "monitoring_path") + "/"
    file = monitoring_path + source_name + "_LINE" + "_" + obs_name + ".TXT"

    spectrum = read_spectrum_file(file)
    plt.plot(spectrum[0], spectrum[2], label="RCP")
    plt.plot(spectrum[1], spectrum[3], label="LRP")

    plt.xlabel("Velocity [km/s]", labelpad=14)
    plt.ylabel("Flux [Jy]", labelpad=10)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Plot a single file. ''', epilog="""plot monitoring.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("exper", help="Experiment name", type=str)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str,
                        default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    main(args.source, args.exper)
    sys.exit()
