import argparse
import sys
from datetime import datetime

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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


def main(sourcea, sourceb, config, config_plot):
    plt.style.use(config_plot)
    dtype = np.dtype([("date", "S12"), ("obs", "S6"), ("amp", float), ("error", float)])
    format = "%d-%b-%Y"

    contiuum_data_a = (get_configs("paths", "monitoring_path", config)
                     + "/" + sourcea + "/UVFIT_" + sourcea + "_.txt")
    contiuum_date_a, contiuum_obs_a, contiuum_amp_a, contiuum_amp_error_a,  = np.loadtxt(contiuum_data_a,
                                                                 usecols=(0, 1, 2, 3), unpack=True, dtype=dtype)
    contiuum_obs_a = [str(obs).replace("b", "").replace("'", "") for obs in contiuum_obs_a]
    contiuum_date_a = [str(date).replace("b", "").replace("'", "") for date in contiuum_date_a]
    mjd_a = Time([datetime.strptime(date, format) for date in contiuum_date_a]).mjd

    contiuum_data_b = (get_configs("paths", "monitoring_path", config)
                       + "/" + sourceb + "/UVFIT_" + sourceb + "_.txt")
    contiuum_date_b, contiuum_obs_b, contiuum_amp_b, contiuum_amp_error_b = np.loadtxt(contiuum_data_b,
                                                                       usecols=(0, 1, 2, 3), unpack=True, dtype=dtype)
    contiuum_obs_b = [str(obs).replace("b", "").replace("'", "") for obs in contiuum_obs_b]
    contiuum_date_b = [str(date).replace("b", "").replace("'", "") for date in contiuum_date_b]
    mjd_b = Time([datetime.strptime(date, format) for date in contiuum_date_b]).mjd

    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)

    for obs_a in contiuum_obs_a:
        if obs_a in contiuum_obs_b:
            obs_a_index = contiuum_obs_a.index(obs_a)
            obs_b_index = contiuum_obs_b.index(obs_a)

            amp_a = contiuum_amp_a[obs_a_index]
            amp_b = contiuum_amp_b[obs_b_index]

            amp_a_error = contiuum_amp_error_a[obs_a_index]
            amp_b_error = contiuum_amp_error_b[obs_b_index]

            ax1.scatter(amp_a, amp_b, c="b", s=mjd_a[obs_a_index]/100)
            ax1.errorbar(amp_a, amp_b, xerr=amp_a_error, yerr=amp_b_error, fmt='', ecolor="r")

    ax1.set_xlabel(sourcea + " " + r'$Flux~(\mathrm{Jy})$')
    ax1.set_ylabel(sourceb + " " + r'$Flux~(\mathrm{Jy})$')

    fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    ax2.scatter(mjd_a, contiuum_amp_a, label=sourcea)
    ax2.scatter(mjd_b, contiuum_amp_b, label=sourceb)

    for i in range(0, len(mjd_a)):
        ax2.errorbar(mjd_a[i], contiuum_amp_a[i], yerr=contiuum_amp_error_a[i], fmt='', ecolor="r")

    for j in range(0, len(mjd_b)):
        ax2.errorbar(mjd_b[j], contiuum_amp_b[j], yerr=contiuum_amp_error_b[j], fmt='', ecolor="r")

    ax2.set_xlabel("MJD")
    ax2.set_ylabel(r'$Flux~(\mathrm{Jy})$')
    ax2.legend()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Compare two cont. sources. ''', epilog="""plot monitoring two cont..""")
    parser.add_argument("sourcea", help="Experiment source A", type=str, default="")
    parser.add_argument("sourceb", help="Experiment source B", type=str, default="")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str,
                        default="config/config.cfg")
    parser.add_argument("-cp", "--config_plot", help="Matplotlib configuration cfg file", type=str,
                        default="config/plot.style")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    main(args.sourcea, args.sourceb, args.config, args.config_plot)
    sys.exit(0)
