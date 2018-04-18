#!/usr/bin/env python
"""Brown Dwarf Flux ratio calculator.

Calculates the flux/contrast ratio between a host star and a brown dwarf of a specified mass.

This script uses the SIMBAD database to obtain the host star parameters, such as magnitude and age.
The companion/brown dwarf mass is a given input (in M_Jup) and  is used to obtain the band magnitudes
of the companion from the Baraffe tables.

The magnitude difference between the host and companion are used to calculate the flux/contrast ratio.

"""
import argparse

from baraffe_tables.table_search import baraffe_table_search


def _parser() -> object:
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Baraffe table Query.')
    parser.add_argument('column', help='Input fits file to calibrate',
                        choices=["M/Ms", "Teff", "L/Ls", "g", "R/Rs", "R", "Li/Li0",
                                 "Mv", "Mr", "Mi", "Mj", " Mh", "Mk", "Ml", "Mll", "Mm"])
    parser.add_argument('value', help='Parameter value', type=float)
    parser.add_argument('age', help='Star age (Gyr)', type=float)
    parser.add_argument('-m', '--model', choices=['03', '15', '2003', '2015'],
                        help='Baraffe model to use [2003, 2015]', default='2003', type=str)
    return parser.parse_args()


if __name__ == '__main__':
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    result = baraffe_table_search(**opts)
    print(result)
