import sys

import numpy as np
import matplotlib.pyplot as plt


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


def main():
    file = "W3OH_LINE.TXT"
    spectrum = read_spectrum_file(file)
    plt.plot(spectrum[0], spectrum[2], label="RCP")
    plt.plot(spectrum[1], spectrum[3], label="LRP")

    plt.xlabel("Velocity [km/s]", labelpad=14)
    plt.ylabel("Flux [Jy]", labelpad=10)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
    sys.exit()
