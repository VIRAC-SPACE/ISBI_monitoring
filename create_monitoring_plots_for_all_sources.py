import argparse
import os
import sys

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


def main(config):
    all_sources = [s.split("_")[0] for s in list(get_configs_items("velocities", config).keys())]

    for source in all_sources:
        print("Processing source " + source)
        print("python3.12 vlbi_monitoring.py " + source)
        os.system("python3.12 vlbi_monitoring.py " + source)
        print("python3.12 show_dynamic_spectra.py " + source)
        os.system("python3.12 show_dynamic_spectra.py " + source)

        print("\n\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''all monitoring results''',
                                     epilog="create all monitoring plots")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str,
                        default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    main(args.config)
    sys.exit(0)
