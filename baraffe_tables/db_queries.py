"""Access Databases."""
from typing import Any, Optional, Union

import numpy as np
from PyAstronomy import pyasl
from astroquery.simbad import Simbad


def get_stellar_params(star_name: str) -> Any:
    """"Astroquery SIMBAD search for stellar parameters.

    Parameters
    ----------
    star_name: str
        Stellar name to get parameters for.

    Returns
    -------
    result_table: votable, dict-like

    Notes
    -----
    result_table.colnames = ['MAIN_ID',
    'RA', 'DEC', 'RA_PREC', 'DEC_PREC',
    'COO_ERR_MAJA', 'COO_ERR_MINA', 'COO_ERR_ANGLE', 'COO_QUAL', 'COO_WAVELENGTH', 'COO_BIBCODE',
    'PLX_VALUE', 'PLX_PREC', 'PLX_ERROR', 'PLX_QUAL', 'PLX_BIBCODE',
    'SP_TYPE',
    'FILTER_NAME_B', 'FLUX_B', 'FLUX_ERROR_B', 'FLUX_SYSTEM_B', 'FLUX_BIBCODE_B', 'FLUX_VAR_B', 'FLUX_MULT_B', 'FLUX_QUAL_B', 'FLUX_UNIT_B',
    'FILTER_NAME_V', 'FLUX_V', 'FLUX_ERROR_V', 'FLUX_SYSTEM_V', 'FLUX_BIBCODE_V', 'FLUX_VAR_V', 'FLUX_MULT_V', 'FLUX_QUAL_V', 'FLUX_UNIT_V',
    'FILTER_NAME_J', 'FLUX_J', 'FLUX_ERROR_J', 'FLUX_SYSTEM_J', 'FLUX_BIBCODE_J', 'FLUX_VAR_J', 'FLUX_MULT_J', 'FLUX_QUAL_J', 'FLUX_UNIT_J',
    'FILTER_NAME_H', 'FLUX_H', 'FLUX_ERROR_H', 'FLUX_SYSTEM_H', 'FLUX_BIBCODE_H', 'FLUX_VAR_H', 'FLUX_MULT_H', 'FLUX_QUAL_H', 'FLUX_UNIT_H',
    'FILTER_NAME_K', 'FLUX_K', 'FLUX_ERROR_K', 'FLUX_SYSTEM_K', 'FLUX_BIBCODE_K', 'FLUX_VAR_K', 'FLUX_MULT_K', 'FLUX_QUAL_K', 'FLUX_UNIT_K',
    'Fe_H_Teff', 'Fe_H_log_g', 'Fe_H_Fe_H', 'Fe_H_flag', 'Fe_H_CompStar', 'Fe_H_CatNo', 'Fe_H_bibcode',
    'name']

    """
    # return Magnitudes, parallax, Temp
    customSimbad = Simbad()
    # Can add more fluxes here if need to extend flux ranges. Although K is the SIMBAD limit.
    # if want higher need to search for Wise band in VISIER probably.
    customSimbad.add_votable_fields('parallax', 'sp', 'fluxdata(B)',
                                    'fluxdata(V)', 'fluxdata(J)', 'fluxdata(H)', 'fluxdata(K)',
                                    'fe_h')

    result_table = customSimbad.query_object(star_name)

    return result_table


def get_sweet_cat_temp(star_name: str) -> Union[float, int]:
    """Obtain spectroscopic temperature from SWEET-Cat.

    Parameters
    ----------
    star_name: str
        Star identifier. HD number only accepted currently.

    """
    sc = pyasl.SWEETCat()
    data = sc.data

    if star_name[0:2].lower() != "hd":
        # only accept HD numbers atm
        raise NotImplementedError

    # Assuming given as hd******
    hd_number = star_name[2:]
    # print("hd number ", hd_number)
    if hd_number in sc.data.hd.values:
        hd_entry = data[data.hd == hd_number]

        if hd_entry.empty:
            return 0
        elif (hd_entry.iloc[0]["teff"] != 0) and (not np.isnan(hd_entry.iloc[0]["teff"])):
            # Temp = 0 when doesn't exist
            return hd_entry.iloc[0]["teff"]
        else:
            return 0
    else:
        print("{!s} was not in SWEET-Cat.".format(star_name))
        return 0


def get_temperature(star_name: str, star_params: Optional[Any] = None) -> float:
    """Find temperature of the star multiple ways.

    1st - Try Fe_H_Teff param from SIMBAD.
    2nd - Try SWEETCat but the star might not be there (only planet hosts).
    3rd - Calculate from B-V and interpolation.

    """
    if star_params is None:
        star_params = get_stellar_params(star_name)

    good_temp = False
    # This is not the best way to do this but work atm
    if "Fe_H_Teff" in star_params.keys():  # need table and interpolate to this B-V
        # print("star_params['Fe_H_Teff'] =", star_params["Fe_H_Teff"])
        teff = star_params["Fe_H_Teff"][0]
        if teff == 0 or teff == [0]:
            # No teff given by SIMBAD
            print("SIMBAD Temperature was zero.")
            teff = None
        else:
            good_temp = teff
            print("Temperature obtained from Fe_H_Teff = {0:5.0f} K".format(good_temp))
            return teff

    if not good_temp:
        teff = get_sweet_cat_temp(star_name)

        if (teff == 0) or (np.isnan(teff)):  # temp from sweet-cat
            print("No SWEET-Cat temperature, teff was {0} K".format(teff))
            teff = None
        else:
            print("SWEET-Cat teff = {0:.0f} K".format(teff))
            good_temp = True
            return teff

    if not good_temp:
        print("Using the B-V method as last resort.")
        teff = calculate_bv_temp(star_params["FLUX_B"], star_params["FLUX_V"])
        print("Temperature of star was calculated from b-v ~= {} K".format(teff))
    return teff


def calculate_bv_temp(b_mag: float, v_mag: float) -> float:
    """Calculate Stellar Temperature from B-V magnitudes.

    Parameters
    ----------
    b_mag: float
        Stellar B magnitude.
    v_mag: float
        Stellar V magnitude.
    Returns
    -------
    temp: float
       Temperature value in Kelvin.

    """
    b_v = b_mag - v_mag

    # Interpolate from B-V
    b_minus_vs = np.array([-0.31, -0.24, -0.20, -0.12, 0.0, 0.15, 0.29,
                         0.42, 0.58, 0.69, 0.85, 1.16, 1.42, 1.61])
    temps = np.array([34000, 23000, 18500, 13000, 9500, 8500, 7300,
                      6600, 5900, 5600, 5100, 4200, 3700, 3000])

    return np.interp(b_v, b_minus_vs, temps)