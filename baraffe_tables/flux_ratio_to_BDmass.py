#!/usr/bin/env python
"""Brown Dwarf Mass calculator.

Uses stellar parameter databases to find host star parameters. The
magnitude of the low mass companion from the provided flux ratio and the
corresponding mass is looked up in the Baraffe evolutionary models.

Inputs
------
Star name: str
    Stellar identification number. eg. HD30501
flux_ratio: float
    Flux ratio between host and companion.
age: float
    Stellar Age. (Closest model is used)

"""
from __future__ import division, print_function

import argparse
import sys
from typing import List, Optional

from astropy.constants import M_jup, M_sun

try:
    from db_queries import get_stellar_params
    from table_search import magnitude_table_search
    from calculations import calculate_companion_magnitude, absolute_magnitude
except ImportError:
    from baraffe_tables.db_queries import get_stellar_params
    from baraffe_tables.table_search import magnitude_table_search
    from baraffe_tables.calculations import calculate_companion_magnitude, absolute_magnitude


def _parser() -> object:
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description='Determine mass of stellar companion from a flux ratio')
    parser.add_argument('star_name', help='Name of host star.', type=str)
    parser.add_argument('flux_ratio', type=float,
                        help='Flux ratio between host and companion (F_companion/F_host)')
    parser.add_argument('stellar_age', help='Star age (Gyr)', type=float)
    parser.add_argument("-b", "--bands", choices=["All", "J", "H", "K"], default=["K"],
                        help='Magnitude bands for the flux ratio value', nargs="+", type=str)
    parser.add_argument('-m', '--model', choices=['03', '15', '2003', '2015'],
                        help='Baraffe model to use [2003, 2015]',
                        default='2003', type=str)
    return parser.parse_args()


def main(star_name: str, flux_ratio: float, stellar_age: float,
         bands: Optional[List[str]] = None, model: str = "2003") -> int:
    """Compute companion mass from flux ratio value.

    Parameters
    ----------
    star_name: str
        Stellar identification number. eg. HD30501
    flux_ratio: float
        Flux ratio for the system (F_companion/F_host).
    stellar_age: float
        Age of star/system (Gyr).
    bands: str
        Wavelength band to use. (optional)
    model: int (optional)
       Year of Baraffe model to use [2003 (default), 2015].

    """
    Jup_sol_mass = (M_sun / M_jup).value  # Jupiter's in 1 M_sol

    if (bands is None) or ("All" in bands):
        bands = ["H", "J", "K"]

    # Obtain Stellar parameters from astroquery
    star_params = get_stellar_params(star_name)  # returns a astroquery result table

    for band in bands:
        mag_label = "FLUX_{0!s}".format(band)
        companion_mag_label = "M{0!s}".format(band.lower())

        # Convert stellar apparent mag to absolute magnitude.
        apparent_mag = star_params[mag_label]
        parallax = star_params['PLX_VALUE']
        if parallax.unit != "mas":
            raise ValueError("Parallax unit not correct")
        absolute_mag = absolute_magnitude(parallax.data[0], apparent_mag.data[0])

        # Calculate Absolute companion magnitude for this flux ratio
        companion_mag = calculate_companion_magnitude(absolute_mag, flux_ratio)

        print("Magnitude calculation for companion M{0} = {1}".format(band, companion_mag))

        # Find companion parameters that match these magnitudes
        companion_params = magnitude_table_search(companion_mag, stellar_age,
                                                  band=band, model=model)

        # Print flux ratios using a generator
        print("Estimated Companion Mass from {0} band flux ratio".format(band.upper()))
        print("M/M_S = {0} (M_star)".format(companion_params["M/Ms"]) +
              " = {} (M_Jup)".format(Jup_sol_mass * companion_params["M/Ms"]) +
              ", Temp = {} K".format(companion_params["Teff"]))

    return 0


if __name__ == '__main__':
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))
