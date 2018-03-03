#!/usr/bin/env python
"""Brown Dwarf Flux ratio calculator.

Calculates the flux/contrast ratio between a host star and a brown dwarf of a specified mass.

This script uses the SIMBAD database to obtain the host star parameters, such as magnitude and age.
The companion/brown dwarf mass is a given input (in M_Jup) and  is used to obtain the band magnitudes
of the companion from the Baraffe tables.

The magnitude difference between the host and companion are used to calculate the flux/contrast ratio.

Inputs
------
Star name: str
    Stellar identification number. eg. HD30501
companion_mass: float
    Mass of companion in Jupiter masses
age: float
    Stellar Age. (Closest model is used)
bands: list of str
    Spectral bands to obtain ratio.
model: str
    Choose between the 2003 and 2015 Baraffe modeling.

"""
# TODO: Interpolate between tables?
from __future__ import division, print_function
import logging
import argparse
import sys
from typing import List, Optional

import numpy as np
from astropy.constants import M_jup, M_sun

try:
    from db_queries import get_stellar_params
    from table_search import mass_table_search
    from calculations import calculate_flux_ratio, calculate_stellar_radius, flux_mag_ratio, absolute_magnitude
except ImportError:
    from baraffe_tables.db_queries import get_stellar_params
    from baraffe_tables.table_search import mass_table_search
    from baraffe_tables.calculations import calculate_stellar_radius, flux_mag_ratio, absolute_magnitude
    #from baraffe_tables.old_code import calculate_flux_ratio


def _parser() -> object:
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description='Determine flux ratio of stellar companion')
    parser.add_argument('star_name', help='Input fits file to calibrate')
    parser.add_argument('companion_mass', help='Mass of companion (M_Jup)', type=float)
    parser.add_argument('stellar_age', help='Star age (Gyr)', type=float)
    parser.add_argument('-b', '--bands', choices=["All", "J", "H", "K"],
                        default=["All"], nargs="+", type=str,
                        help='Spectral Band to measure. Options=["All", "K", ""]')
    parser.add_argument('-m', '--model', choices=['03', '15', '2003', '2015'],
                        help='Baraffe model to use [2003, 2015]', default='2003', type=str)
    parser.add_argument("-a", "--area_ratio", default=False, action="store_true",
                        help="Calculate the area ratio.")
    parser.add_argument("-p", "--paper", default=False, action="store_true",
                        help="Print more parameters for paper.")
    parser.add_argument("-s", "--star_pars", default=False, action="store_true",
                        help="Print star parameters for paper.")
    return parser.parse_args()


def main(star_name: str, companion_mass: float, stellar_age: float, bands: Optional[List[str]] = None,
         model: str = "2003", area_ratio: bool = False, paper: bool = False, star_pars: bool = False) -> int:
    """Compute flux/contrast ratio between a stellar host and companion.

    Parameters
    ----------
    star_name: str
        Stellar identification number. eg. HD30501.
    companion_mass: float
        Mass of companion in Jupiter masses.
    stellar_age: float
        Stellar Age. (Closest model is used).
    bands: list of str
        Spectral bands to obtain ratio.
    model: str (optional)
        Year of Baraffe model to use [2003 (default), 2015].
    area_ratio: bool default=False
        Perform simple radius and area comparisons calculations.
    paper: bool
        Print other parameters need for paper table.
    star_pars: bool
        Print star parameters also.

    """
    if (bands is None) or ("All" in bands):
        bands = ["J", "H", "K"]

    # Obtain Stellar parameters from astroquery
    star_params = get_stellar_params(star_name)  # returns a astroquery result table

    companion_mass_solar = companion_mass * (M_jup / M_sun).value  # transform to solar mass for table search

    # Get parameters for this mass and age
    companion_params = mass_table_search(companion_mass_solar, stellar_age, model=model)

    # flux_ratios = calculate_flux_ratio(star_params, companion_params, bands)
    # print("old ratios", flux_ratios)

    flux_ratios = {}
    for band in bands:
        try:
            mag_label = "FLUX_{0!s}".format(band)
            companion_mag_label = "M{0!s}".format(band.lower())

            # Convert stellar apparent mag to absolute magnitudes.
            apparent_mag = star_params[mag_label]
            parallax = star_params['PLX_VALUE']
            if parallax.unit != "mas":
                raise ValueError("Parallax unit not correct")
            absolute_mag = absolute_magnitude(parallax.data[0], apparent_mag.data[0])

            # Model magnitudes are absolute
            companion_mag = companion_params[companion_mag_label]

            flux_ratio = flux_mag_ratio(absolute_mag, companion_mag)

            flux_ratios.update({band: flux_ratio})
        except:
            logging.warning("Unable to calculate flux ratio for {} band".format(band))

    print("\nFlux ratios:")
    for key, val in flux_ratios.items():
        if key in bands:
            print(("{0!s} band  companion/star Flux ratio = {2:0.4f}"
                   " >>> star/companion Flux ratio = {1:4.2f}").format(key, 1. / val, val))

    if area_ratio:
        # Compare to area ratio
        Rstar = calculate_stellar_radius(star_params)
        print(Rstar)
        Rcomp_Rstar = companion_params["R"] / Rstar

        print("\nRadius Calculation")
        print("Host radius         = {} R_sun".format(Rstar))
        print("companion radius    = {} R_sun".format(np.round(companion_params["R"], 4)))
        print("Radius Ratio of companion/star    = {}".format(Rcomp_Rstar))
        print("Area Ratio of companion/star      = {}".format(Rcomp_Rstar ** 2))

    if paper:
        print(companion_params)
        print(r"Host name, star M_K, companion_mass, companion Teff, companion_Mk, K_ratio")
        print(r"{0!s} & {1:0.2f} & {2:.1f} & {3:4.0f} & {4:.2f} & {5:.1f}\\\\n".format(star_name,
                                                                                       star_params["FLUX_K"][0],
                                                                                       companion_mass,
                                                                                       companion_params["Teff"],
                                                                                       companion_params["Mk"],
                                                                                       flux_ratios["K"]))
    if star_pars:
        print("\nStellar parameters:")
        star_params.pprint(show_unit=True)

    return 0


if __name__ == '__main__':
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))
