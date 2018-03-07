"""Calculations for flux ratios."""
from typing import Any

import numpy as np

try:
    from db_queries import get_temperature
except ImportError:
    from baraffe_tables.db_queries import get_temperature


def flux_mag_ratio(mag1: float, mag2: float) -> float:
    """Calculate the flux ratio between two magnitudes.

    Using the equation f1/f2 = 10**(-0.4 * (mag1 - mag2)).
    A common approximation is the equation f1/f2 approx 2.512**(m2 - m1), but is not here.

    Parameters
    ----------
    mag1 float
        Magnitude of first object.
    mag2: float
        Magnitude of second object.

    Returns
    -------
    flux_ratio: float
        flux/contrast ratio between the two magnitudes.

    """
    flux_ratio = 10 ** (-0.4 * (mag1 - mag2))
    return flux_ratio


def calculate_stellar_radius(star_params: Any) -> float:
    """Based on R/Rs = (Ts/T)^2(L/Ls)^(1/2) equation.

    Parameters
    ----------
    star_params: votable, dict
        Table of Stellar parameters.

    Returns
    -------
    R_Rs: float
        Estimated Stellar Radius in solar radii.

    """
    star_name = star_params['MAIN_ID'][0].decode("utf-8")
    teff_star = get_temperature(star_name, star_params)

    Ts_T = 5800. / teff_star  # Temperature ratio
    Dm = 4.83 - star_params["FLUX_V"][0]  # Difference of absolute magnitude
    L_Ls = 2.51 ** Dm  # Luminosity ratio
    R_Rs = (Ts_T) ** 2 * np.sqrt(L_Ls)  # Radius of Star in Solar Radii

    return R_Rs  # Radius of star in solar radii


def calculate_companion_magnitude(star_mag: float, flux_ratio: float) -> float:
    """Calculate companion magnitude from flux ratio.

    Using the equation m - n = -2.5 * log_10(F_m / F_n).

    Parameters
    ----------
    star_mag: float
        Host star magnitude.
    flux_ratio: float
        Flux ratio for the system (F_companion/F_host).

    Returns
    -------
    magnitude: float
        Companion magnitude.

    Note
    ----
    This is possibly not quite the correct implementation as we are
    only using a single flux_ratio value.

    """

    magnitude = star_mag - 2.5 * np.log10(flux_ratio)
    return magnitude


def distance_modulus(d: float):
    """Calculate distance modulus.

    Input
    -----
    d: float
        Distance in parsec
    Output
    ------
    mu: float
        m-M distance modulus
    """
    mu = 5 * np.log10(d) - 5
    return mu


def absolute_magnitude(parallax, m):
    """Calculate the absolute magnitude based on distance and apparent mag.
    Inputs
    ------
    parallax : float
      The parallax in mas
    m : float
      The apparent magnitude
    Output
    ------
    M : float
      The absolute magnitude
    """
    d = 1. / (parallax * 1e-3)  # Conversion to arcsecond before deriving distance
    mu = distance_modulus(d)
    M = m - mu
    return M


def apparent_magnitude(parallax, M):
    """Calculate the apparent magnitude based on distance and absolute mag.
    Inputs
    ------
    parallax : float
      The parallax in mas
    M : float
      The absolute magnitude
    ------
    m : float
      The apparent magnitude
    """
    d = 1. / (parallax * 1e-3)  # Conversion to arcsecond before deriving distance
    mu = distance_modulus(d)
    m = M + mu
    return m
