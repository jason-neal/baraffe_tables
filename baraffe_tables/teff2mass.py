#!/usr/bin/env python
"""Stellar mass lookup.

Looks up the stellar mass given its temperature, logg and age.

Inputs
------
Temperature: str
    Temperature.
logg: float
    Logg.
age: float
    Stellar Age. (Closest model is used)
model: str
    Choose between the 2003 and 2015 Baraffe modelling.

"""
# TODO: Interpolate between tables?

import argparse
import sys
from typing import Union
import numpy as np
from astropy.constants import M_jup, M_sun

from baraffe_tables.table_search import baraffe_table_search
from baraffe_tables.table_search import model_ages_15, model_ages_03
import matplotlib.pyplot as plt


def _parser() -> object:
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description='Lookup companion mass from temperature.')
    parser.add_argument('temp', help='Temperature of companion.', type=float)
    parser.add_argument('logg', help='Logg of companion.', type=float)
    parser.add_argument('-p', '--plot', action="store_true",
                        help='Plot the age-logg line.', default='2003')
    parser.add_argument("-f", "--full_table", default=False, action="store_true",
                        help="Print all parameters for found companion.")
    return parser.parse_args()


def main(temp: Union[float, int], logg: float, plot: bool = False) -> int:
    """Compute flux/contrast ratio between a stellar host and companion.

    Parameters
    ----------
    temp: str
        Temperature of companion
    logg: float
        Logg of companion
    plot: bool
        Plot teff vs logg.

    """
    ages_03 = np.asarray(model_ages_03, dtype=float)
    ages_15 = np.asarray(model_ages_15, dtype=float)

    loggs_03 = []
    loggs_15 = []
    masses_03 = []
    masses_15 = []
    radii_03 = []
    radii_15 = []

    for age in ages_03:
        result = baraffe_table_search("Teff", temp, age, model="03")
        masses_03.append(result["M/Ms"])
        loggs_03.append(result["g"])
        radii_03.append(result["R"])

    for age in ages_15:
        result = baraffe_table_search("Teff", temp, age, model="15")
        masses_15.append(result["M/Ms"])
        loggs_15.append(result["g"])
        radii_15.append(result["R/Rs"])

    loggs_03 = np.asarray(loggs_03)
    loggs_15 = np.asarray(loggs_15)
    masses_03 = np.asarray(masses_03)
    masses_15 = np.asarray(masses_15)
    radii_03 = np.asarray(radii_03)
    radii_15 = np.asarray(radii_15)

    sim_age_15 = (ages_15[np.argmin(abs(loggs_15 - logg))])
    sim_age_03 = (ages_03[np.argmin(abs(loggs_03 - logg))])

    if plot:
        plt.subplot(111)
        plt.axhline(logg, alpha=0.5)
        plt.semilogx(ages_15, loggs_15, ".-", label="Baraffe 2015")
        plt.semilogx(ages_03, loggs_03, ".-", label="Baraffe 2003")
        plt.plot(sim_age_03, logg, "h", label="closest 03")
        plt.plot(sim_age_15, logg, "+", label="closest 15")
        plt.xlabel("Age (Gyr)")
        plt.ylabel("logg (dex)")
        plt.legend()
        plt.title("Temperature {} K".format(temp))
        plt.annotate("Exploring evolutionary models to \nfind the table/age that matches \nTemp & logg", xy=(.5, .35), xycoords="axes fraction")
        plt.tight_layout()
        plt.show()

        plt.semilogx(ages_03, masses_03, "+-")
        plt.semilogx(ages_15, masses_15, ".-")
        plt.xlabel("Age (Gyr)")
        plt.ylabel("Mass (M/Ms)")
        plt.title("Temp= {}, Logg = {}".format(temp, logg))
        plt.show()

        plt.semilogx(ages_03, radii_03, "+-")
        plt.semilogx(ages_15, radii_15, ".-")
        plt.xlabel("Age (Gyr)")
        plt.ylabel("Radii (R/Rs)")
        plt.title("Temp= {}, Logg = {}".format(temp, logg))
        plt.show()

        plt.plot(radii_03, masses_03, "+-")
        plt.plot(radii_15, masses_15, ".-")
        plt.xlabel("Radii (R/Rs)")
        plt.ylabel("Mass (M/Ms)")
        plt.title("Temp= {}, Logg = {}".format(temp, logg))
        plt.show()

    # Since my cutoffs are around 3500K I will use the 2015 models
    # TODO Add temperature range checks that the value is inside the temperature range

    result15 = baraffe_table_search("Teff", temp, age=sim_age_15, model="15")
    result03 = baraffe_table_search("Teff", temp, age=sim_age_03, model="03")

    # Add age to table
    result15["age"] = sim_age_15
    result03["age"] = sim_age_03

    # Add jupyter mass
    result15["M/Mjup"] = result15["M/Ms"] * (M_sun / M_jup).value
    result03["M/Mjup"] = result03["M/Ms"] * (M_sun / M_jup).value

    if temp == result15["Teff"]:
        # Preference for 2015 models due to being newer and they have a higher age resolution
        return result15
    elif temp == result03["Teff"]:
        return result03
    else:
        raise ValueError("Temp {} is not within either model grid".format(temp))


if __name__ == '__main__':
    args = vars(_parser())
    full_table = args.pop("full_table", False)

    opts = {k: args[k] for k in args}
    result = main(**opts)

    if full_table:
        print(result)

    print("Teff-logg companion mass")
    print("temp\t= {0:5.0f} K\nlogg\t= {1:4.02}\nage \t= {2:06.4f} Gyr\nmass\t= {3:5.01f} Mjup\n".format(*[
        result[key] for key in ("Teff", "g", "age", "M/Mjup")]))

    sys.exit(0)
