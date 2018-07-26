"""Test companion/flux ratio codes with Baraffe tables."""
import sys

import numpy as np
import pytest
from astropy.constants import M_jup, M_sun
from baraffe_tables.BDmass_to_flux_ratio import _parser as mass_parser
from baraffe_tables.BDmass_to_flux_ratio import main as mass_main
from baraffe_tables.calculations import calculate_companion_magnitude
from baraffe_tables.calculations import flux_mag_ratio
from baraffe_tables.db_queries import (get_stellar_params, get_sweet_cat_temp,
                                       get_temperature)
from baraffe_tables.flux_ratio_to_BDmass import _parser as ratio_parser
from baraffe_tables.flux_ratio_to_BDmass import main as ratio_main
from baraffe_tables.table_search import (age_table, magnitude_table_search,
                                         mass_table_search, baraffe_table_search)

org_sysargv = sys.argv


@pytest.mark.parametrize("area_ratio", [True, False])
def test_BD_to_flux_runs(area_ratio):
    """Check it returns 0 (Runs normally).

    If there is no internet then an Exception is raised.
        Exception: Query failed: HTTPConnectionPool(host='simbad.u-strasbg.fr', port=80): Max retries exceeded with
        url: /simbad/sim-script (Caused by NewConnectionError('<requests.packages.urllib3.connection.HTTPConnection
        object at 0x7fe4e509b438>: Failed to establish a new connection: [Errno -3] Temporary failure in name
        resolution')).
    """
    assert mass_main("HD30501", 90, 5, area_ratio=area_ratio) is 0


def test_ratio_to_BD_runs():
    """Check it returns 0 (Runs normally).

    If there is no internet then an Exception is raised.
        Exception: Query failed: HTTPConnectionPool(host='simbad.u-strasbg.fr', port=80): Max retries exceeded with
        url: /simbad/sim-script (Caused by NewConnectionError('<requests.packages.urllib3.connection.HTTPConnection
        object at 0x7fe4e509b438>: Failed to establish a new connection: [Errno -3] Temporary failure in name
        resolution', )).
    """
    assert ratio_main("HD30501", 0.01, 5) is 0


@pytest.mark.parametrize("hd_number, temp", [
    ("HD107383", 4830),
    ("HD99706", 4891),
    ("HD195689", 10170)])
def test_get_sweet_cat_temp(hd_number, temp):
    """Test getting temperature from Sweet-cat."""
    # HD number in SWEETCat
    assert temp == get_sweet_cat_temp(hd_number)


@pytest.mark.parametrize("hd_number", [
    ("HD10"),
    ("HD1234565830"),
])
def test_get_temp_star_not_in_sweet_cat_temp(hd_number):
    """Test getting from Sweet-cat."""
    # HD number in SWEETCat
    assert get_sweet_cat_temp(hd_number) == 0


def test_get_sweet_cat_temp_without_hd_number():
    """Not hd number call not implemented yet."""
    with pytest.raises(NotImplementedError):
        get_sweet_cat_temp("GJ 422")  # in SWEETCat but not an hd number


@pytest.mark.parametrize("hd_number", [
    ("HD145934"),
    ("HD41004B")
])
def test_get_sweet_cat_temp_no_temp(hd_number):
    """Not hd number call not implemented yet.

    # 2 tries that line in Sweet-Cat that has no temperature
    # These are the only two that have hd numbers
    """
    assert get_sweet_cat_temp(hd_number) == 0  # This may change if SWEETCat is updated


def test_get_temperature_without_params_input():
    """Test it calls for params itself, if needed.

    If there is no internet then an Exception is raised.
        Exception: Query failed: HTTPConnectionPool(host='simbad.u-strasbg.fr', port=80): Max retries exceeded with
        url: /simbad/sim-script (Caused by NewConnectionError('<requests.packages.urllib3.connection.HTTPConnection
        object at 0x7fe4e509b438>: Failed to establish a new connection: [Errno -3] Temporary failure in name
        resolution', )).
    """
    name = "HD30501"
    params = get_stellar_params(name)
    assert get_temperature(name, star_params=None) == get_temperature(name, star_params=params)


@pytest.mark.parametrize("name, temp", [
    ("HD215909", 4328),  # SIMBAD temp
    ("HD343246", 5703),  # SweetCat temp
])
def test_get_temperature_examples(name, temp):
    """Test some temperatures."""
    assert get_temperature(name) == temp


def test_get_temperature_ignores_zero_temp():
    star = "HD12116"  # Has zero temp in SIMBAD (Check this is not used)
    assert get_temperature(star) != 0


@pytest.mark.xfail()
@pytest.mark.parametrize("param_key, expected", [
    ("FLUX_B", 8.68),
    ("FLUX_V", 8.01),
    ("FLUX_K", 6.530),
    ('Fe_H_log_g', 4.19),  # log g
    ('Fe_H_Fe_H', 0.16),  # metallicity
])
def test_get_stellar_params(param_key, expected):
    """Test some values from SIMBAD database result.
    Fails as a result of parameters changing in database
    # """
    name = "HD219828"
    params = get_stellar_params(name)
    assert params[param_key] == expected


@pytest.mark.xfail()
@pytest.mark.parametrize("param_key, expected", [
    ("Fe_H_Teff", 5842),
    ("PLX_VALUE", 14.0),  # parallax
])
def test_get_stellar_params_list(param_key, expected):
    """Test list values from SIMBAD database result.
    May Fail as a result of parameters changing in database
    """
    name = "HD219828"
    params = get_stellar_params(name)

    assert params[param_key][0] == expected


# Test Flux ratio to Mass
@pytest.mark.parametrize("mag_1, mag_2", [(3, 7), (4, 2), (5, 11)])
def test_mag_conversions(mag_1, mag_2):
    """Test converting from flux ratio to magnitude back to flux ratio etc.

    Tests show the conversion goes back and forward.
    """
    ratio = flux_mag_ratio(mag_1, mag_2)
    magnitude = calculate_companion_magnitude(mag_1, 1. / ratio)
    assert magnitude == mag_2
    new_ratio = flux_mag_ratio(mag_1, magnitude)

    assert np.allclose(new_ratio, ratio)


@pytest.mark.parametrize("mass", [80, 100])
@pytest.mark.parametrize("age", [1, 3.5, 5, 10])
@pytest.mark.parametrize("band", ["H", "J", "K"])
def test_table_searches(mass, age, band, age_interp, baraffe_model):
    """That a given mass calculates a magnitude and that mag finds a correct mass.

    Narrow mass overlap window.
    """
    start_sol_mass = mass * (M_jup / M_sun).value

    # Find magnitude from mass
    mass_comp_params = mass_table_search(start_sol_mass, age=age, model=baraffe_model,
                                         age_interp=age_interp)
    found_mag = mass_comp_params["M{0!s}".format(band.lower())]

    # Find mass from magnitude
    mag_comp_params = magnitude_table_search(found_mag, age=age, band=band, model=baraffe_model,
                                             age_interp=age_interp)
    found_mass = mag_comp_params["M/Ms"]

    assert np.allclose(found_mass, start_sol_mass)
    assert np.allclose(found_mass * M_sun / M_jup, mass)  # Back int M_jup


@pytest.mark.parametrize("mass, model", [(50, "03"), (160, "15")])
@pytest.mark.parametrize("age", [1, 3.5, 5, 10])
@pytest.mark.parametrize("band", ["H", "J", "K"])
def test_table_searches_individual(mass, age, band, age_interp, model):
    """That a given mass calculates a magnitude and that mag finds a correct mass.
    Some wider extremes that do not work at the with the opposite model"""
    start_sol_mass = mass * (M_jup / M_sun).value

    # Find magnitude from mass
    mass_comp_params = mass_table_search(start_sol_mass, age=age, model=model,
                                         age_interp=age_interp)
    found_mag = mass_comp_params["M{0!s}".format(band.lower())]

    # Find mass from magnitude
    mag_comp_params = magnitude_table_search(found_mag, age=age, band=band, model=model,
                                             age_interp=age_interp)
    found_mass = mag_comp_params["M/Ms"]

    assert np.allclose(found_mass, start_sol_mass)
    assert np.allclose(found_mass * M_sun / M_jup, mass)  # Back int M_jup


@pytest.mark.parametrize("age", [0.1, 0.65, 1, 4.5, 7.893, 10.])
def test_age_table(age, baraffe_model, age_interp):
    """Select a Baraffe table of certain age and model."""
    model_table, cols, __ = age_table(age, model=baraffe_model, age_interp=age_interp)

    assert isinstance(model_table, dict)
    assert isinstance(cols, list)

    for key in model_table:
        print("type of value", type(model_table[key]))
        assert isinstance(model_table[key], np.ndarray)

    # check common columns are present
    common_cols = ["M/Ms", "Teff", "L/Ls", "g", "Mv",
                   "Mr", "Mi", "Mj", "Mh", "Mk", "Mll", "Mm"]
    for header in common_cols:
        assert header in cols


@pytest.fixture(params=[False, True])
def age_interp(request):
    """Fixture for age_interp flag."""
    return request.param


@pytest.fixture(params=["2003", "2015", "03", "15"])
def baraffe_model(request):
    """Fixture for baraffe model."""
    return request.param


@pytest.mark.parametrize("model", ["2016", "", 2003, 15, "word"])
def test_bad_model_age_table(model, age_interp):
    """Call model with wrong values."""
    with pytest.raises(ValueError):
        age_table(5, model, age_interp=age_interp)


@pytest.mark.parametrize("mag, ratio, result", [
    (1, 100, -4),
    (10, 100, 5),
    (5, 10000, -5),
    (1, 1. / 100, 6),
    (10, 1. / 100, 15),
    (5, 1. / 10000, 15)
])
def test_calculate_comp_magnitude(mag, ratio, result):
    """Test flux ratio to magnitude difference."""
    magnitude = calculate_companion_magnitude(mag, ratio)
    print(magnitude, result)
    assert np.allclose(magnitude, result)


def test_mag_table_search_band_fail(age_interp):
    """One one band value is allowed."""
    magnitudes = {"H": 1, "J": 4, "K": 5}

    with pytest.raises(ValueError):
        magnitude_table_search(magnitudes, age=5, band=["H", "J", "K"], model="2003",
                               age_interp=age_interp)


def test_mass_table_search_03():
    """That a value from the table returns the correct row."""
    comp_params = mass_table_search(0.09, age=5, model="2003")
    assert comp_params["M/Ms"] == 0.09
    assert comp_params["Teff"] == 2622
    assert comp_params["R"] == 0.113
    assert comp_params["Mk"] == 10.04


def test_mass_table_search_15():
    """Manual test of a 2015 table mass search."""
    comp_params_15 = mass_table_search(0.09, age=5, model="2015")
    assert comp_params_15["M/Ms"] == 0.09
    assert comp_params_15["Teff"] == 2644
    assert comp_params_15["R/Rs"] == 0.113
    assert comp_params_15["Mk"] == 9.91


def test_magnitude_table_search_03():
    """That a value from the table returns the correct row."""
    mag_params = magnitude_table_search(10.04, age=5, band="K", model="2003")
    assert mag_params["M/Ms"] == 0.09
    assert mag_params["Teff"] == 2622
    assert mag_params["R"] == 0.113
    assert mag_params["Mk"] == 10.04


def test_magnitude_table_search_15():
    """Manual test of a 2015 table magnitude search."""
    mag_params_15 = magnitude_table_search(9.91, age=5, band="K", model="2015")
    assert mag_params_15["M/Ms"] == 0.09
    assert mag_params_15["Teff"] == 2644
    assert mag_params_15["R/Rs"] == 0.113
    assert mag_params_15["Mk"] == 9.91


@pytest.mark.parametrize("band", [
    (["H", "K"]),
    (("K",)),
    (["K"]),
    (["H", "K"]),
])
def test_magnitude_table_search_errors(band, age_interp):
    """Test only takes one band as a str."""
    with pytest.raises(ValueError):
        # More than one band not allowed
        magnitude_table_search(magnitude=5, age=5, band=band, model="2015", age_interp=age_interp)


# Need sys.argv fixture with teardown
def test_BDmass_parser():
    """Test argparse function using sys.argv."""
    sys.argv = "pytest HD30501 90 5 -a -f -s".split()

    args = mass_parser()
    assert args.star_name == "HD30501"
    assert args.companion_mass == 90
    assert args.stellar_age == 5
    assert args.bands == ["All"]
    assert args.model == "2003"
    assert args.area_ratio is True
    assert args.full_table == True
    assert args.star_pars == True
    assert args.age_interp == False
    sys.argv = org_sysargv


# Need sys.argv fixture with teardown
def test_ratio_parser():
    """Test argparse function using sys.argv."""
    sys.argv = "pytest HD30501 0.001 5 -b H -m 2015".split()

    args = ratio_parser()
    assert args.star_name == "HD30501"
    assert args.flux_ratio == 0.001
    assert args.bands == ["H"]
    assert args.stellar_age == 5
    assert args.model == "2015"
    assert args.full_table == False
    assert args.star_pars == False
    assert args.age_interp == False
    sys.argv = org_sysargv


def test_ratio_parser2():
    """Test argparse function using sys.argv more than one band."""
    sys.argv = "pytest HD30501 0.001 5 -m 2015 -b H K -f -s --age_interp".split()

    args = ratio_parser()
    assert args.star_name == "HD30501"
    assert args.flux_ratio == 0.001
    assert args.bands == ["H", "K"]
    assert args.stellar_age == 5
    assert args.model == "2015"
    assert args.full_table == True
    assert args.star_pars == True
    assert args.age_interp == True
    sys.argv = org_sysargv


@pytest.mark.parametrize("parse_string", [
    "pytest HD30501 ratio 5 -b H K -m 2015",
    "pytest HD30501 0.001 5 -b H K -m 2015 -z 4",
    "pytest HD305010 0.01 5 -b Q",
])
def test_failing_ratio_parsers(parse_string):
    """Test argparse function using sys.argv."""
    sys.argv = parse_string.split()
    with pytest.raises(SystemExit):
        ratio_parser()  # ratio is not a number

    with pytest.raises(SystemExit):
        mass_parser()  # ratio is not a number


@pytest.mark.parametrize("col", ["G", 7, "RRs"])
def test_table_search_invalid_parameter(col, baraffe_model, age_interp):
    with pytest.raises(ValueError):
        baraffe_table_search(col, 1, age=5, model=baraffe_model, age_interp=age_interp)


@pytest.mark.parametrize("age", [0.5, 4.5, 5.0])
@pytest.mark.parametrize("col, value", [
    ("M/Ms", 0.08),
])
def test_table_search_returns_the_value_inputed(col, value, age, age_interp, baraffe_model):
    """Test the input value is returned in the result."""
    result = baraffe_table_search(col, value, age=age, model=baraffe_model, age_interp=age_interp)
    assert result[col] == value


@pytest.mark.parametrize("col, value, age, bound", [
    ("M/Ms", 500, 5, "upper"),
    ("M/Ms", 0.0006, 5, "lower"),
    ("Teff", 100, 1, "lower"),
    ("Teff", 10000, 2, "upper"),
    ("Mk", 38, 0.1, "lower"),
    ("Mk", 1, 2, "upper"),
])
# @pytest.mark.xfail(strict=True)
def test_table_search_outside_bound_produces_error(col, value, age, bound, age_interp,
                                                   baraffe_model):
    """Test the warnings when interpolation value outside the reference column values."""
    # Check that the warning is raised.
    with pytest.warns(UserWarning) as record:
        baraffe_table_search(col, value, age=age, model=baraffe_model, age_interp=age_interp)
    # Check correct warning is raised.
    assert len(record) == 1
    assert str(
        record[0].message) == "Interpolated values are outside the {0!s} bound of {1!s}.".format(
        bound, col)
