"""Code to obtain and find row in Baraffe tables."""
import warnings
from typing import Dict, List, Tuple

import numpy as np
import pkg_resources
from scipy.interpolate import interp1d

# Table model details
model_ages_03 = [
    "0.001",
    "0.005",
    "0.010",
    "0.050",
    "0.100",
    "0.120",
    "0.500",
    "1.000",
    "5.000",
    "10.000",
]
cols_03 = [
    "M/Ms",
    "Teff",
    "L/Ls",
    "g",
    "R",
    "Mv",
    "Mr",
    "Mi",
    "Mj",
    "Mh",
    "Mk",
    "Mll",
    "Mm",
]

model_ages_15 = [
    "0.0005",
    "0.001",
    "0.003",
    "0.004",
    "0.005",
    "0.008",
    "0.010",
    "0.015",
    "0.020",
    "0.025",
    "0.030",
    "0.040",
    "0.050",
    "0.080",
    "0.100",
    "0.120",
    "0.200",
    "0.300",
    "0.400",
    "0.500",
    "0.625",
    "0.800",
    "1.000",
    "2.000",
    "3.000",
    "4.000",
    "5.000",
    "8.000",
    "10.000",
]
cols_15 = [
    "M/Ms",
    "Teff",
    "L/Ls",
    "g",
    "R/Rs",
    "Li/Li0",
    "Mv",
    "Mr",
    "Mi",
    "Mj",
    "Mh",
    "Mk",
    "Mll",
    "Mm",
]


def find_bounding_ages(age: float, model_ages: List[str]) -> Tuple[str, str]:
    """ Find the two bounding model ages to age.

    Uses numpy.seachsorted() to find where the age is located.

    Parameters
    ---------
    age: float
        Age to find position in model ages.
    model_ages: List[str]
        List of model ages.

    Returns
    -------
    Lower_age: str
        The smaller age.
    upper_age: str
        The larger age.

    """
    age_array = np.asarray(model_ages)
    m_ages = age_array.astype(np.float)
    sortargs = np.argsort(m_ages)
    indx = m_ages[sortargs].searchsorted(age)  # Where to put age in sorted numpy array
    sorted_ages = age_array[sortargs]  # sort the array of strings
    return sorted_ages[indx - 1], sorted_ages[indx]


def interp_data_dicts(
    age: float,
    lower_age: str,
    lower_data: Dict[str, List[float]],
    upper_age: str,
    upper_data: Dict[str, List[float]],
) -> Dict[str, List[float]]:
    """Interpolate two data dictionaries to a new age.

    The keys should be the same. The lower age may have extra rows at the start which are removed.

    (lower_age < age) and (upper_age > age)

    """

    lower_age, upper_age = float(lower_age), float(upper_age)
    assert (lower_age < age) and (
        age < upper_age
    ), "age is not between lower_age and upper_age!"
    assert set(lower_data.keys()) == set(
        upper_data.keys()
    ), "Data dicts do not have the same keys."
    interp_data_dict = {}
    for key in lower_data:
        data1, data2 = lower_data[key], upper_data[key]

        while len(data1) != len(data2):
            # Lower age data may have more rows at lower masses so remove leading entries
            data1 = data1[1:]

        if np.all(data1 == data2):
            result = data1
        else:
            ages = np.array([lower_age, upper_age])
            data_in = np.vstack((data1, data2))

            interp_function = interp1d(ages, data_in, axis=0)
            data_out = interp_function(age)
            result = np.round(data_out, 3)

        interp_data_dict[key] = result
    return interp_data_dict


def age_table(
    age: float, model: str = "2003", age_interp=False
) -> Tuple[Dict[str, List[float]], List[str], float]:
    """Determine the correct Baraffe table to load.

    Parameters
    ----------
    age: float
        Stellar age (Gyr).
    model: str
        Baraffe model version to use. options=[03, 15, 2003, 2015].
    age_interp: bool
        Interpolate tables across age. Default=False..

    Returns
    -------
    model_data: numpy.ndarray
        The correct model table data.
    column_names: list of str
        List of the columns in the the table.
    model_age: float
        Age of table returned

    """
    if not isinstance(model, str):
        raise ValueError("Model is not the valid type 'str'.")
    elif model not in ["2003", "03", "2015", "15"]:
        raise ValueError("Model value '{}' is not valid".format(model))

    if model in "2003":
        modelages = model_ages_03
        base_name = "data/Baraffe2003/BaraffeCOND2003-"
        skiprows = 18
        cols = cols_03
    else:
        modelages = model_ages_15
        base_name = "data/Baraffe2015/BaraffeBHAC15-"
        skiprows = 22
        cols = cols_15

    closest_age = min(modelages, key=lambda x: abs(float(x) - age))  # Closest one

    if (
        age_interp
        and (float(closest_age) != float(age))
        and (float(modelages[0]) < age)
        and (age < float(modelages[-1]))
    ):
        # Find two closest tables, interp values to given age.
        lower_age, upper_age = find_bounding_ages(age, modelages)
        print(
            "Interpolating tables {0} Gyr and {1} Gyr to {2} Gyr".format(
                lower_age, upper_age, age
            )
        )

        lower_data = model_age_table(base_name, lower_age, skiprows=skiprows)
        upper_data = model_age_table(base_name, upper_age, skiprows=skiprows)

        lower_data_dict = {}
        upper_data_dict = {}
        for i, col in enumerate(cols):
            lower_data_dict[col] = lower_data[i]
            upper_data_dict[col] = upper_data[i]

        # Interpolate age tables together
        data_dict = interp_data_dicts(
            age, lower_age, lower_data_dict, upper_age, upper_data_dict
        )
        model_age = age
    else:
        # Find closest model age table only.
        model_age = closest_age

        model_data = model_age_table(base_name, model_age, skiprows=skiprows)

        # Turn into Dict of values
        data_dict = {col: model_data[i] for i, col in enumerate(cols)}
    return data_dict, cols, model_age


def model_age_table(base_name, model_age, skiprows):
    """Load in model age table."""
    model_id = "p".join(str(model_age).split("."))  # Replace . with p in number str
    model_name = base_name + model_id + "Gyr.dat"
    model_name = pkg_resources.resource_filename("baraffe_tables", model_name)

    model_data = np.loadtxt(model_name, skiprows=skiprows, unpack=False)
    model_data = model_data.T
    return model_data


def mass_table_search(
    companion_mass: float, age: float, model: str = "2003", age_interp: bool = False
) -> Dict[str, float]:
    """Search Baraffe tables to find the companion entry given a mass value.

    Parameters
    ----------
    companion_mass: float
        Companion Mass (M_sun)
    age: float
        Age of star/system (Gyr).
    model: str
       Year of Baraffe model to use [2003 (default), 2015].
    age_interp: bool
        Interpolate tables across age. Default=False.

    Returns
    -------
    companion_parameters: list
        Companion parameters from Baraffe table, interpolated to the provided mass.

    """
    model_data, __, __ = age_table(age, model=model, age_interp=age_interp)

    ref_val = companion_mass
    ref_col = "M/Ms"
    companion_parameters = table_interpolation(model_data, ref_col, ref_val)
    return companion_parameters  # as a dictionary


def magnitude_table_search(
    magnitude: float,
    age: float,
    band: str = "K",
    model: str = "2003",
    age_interp: bool = False,
) -> Dict[str, float]:
    """Search Baraffe tables to find the companion entry given a band magnitude value.

    Parameters
    ----------
    magnitude: float
        Dictionary of (band: magnitude) pairs.
    age: float
        Age of star/system (Gyr).
    band: str
        Wavelength band to use.
    model: str
       Year of Baraffe model to use [2003 (default), 2015].
    age_interp: bool
        Interpolate tables across age. Default=False.

    Returns
    -------
    companion_parameters: list
        Companion parameters from Baraffe table, interpolated between the
        rows to the provided magnitude.

    """
    if not isinstance(band, str):
        raise ValueError(
            "Band {0} was given, when not given as a single string.".format(band)
        )

    ref_col = "M{}".format(band.lower())
    companion_parameters = baraffe_table_search(
        ref_col, magnitude, age, model, age_interp=age_interp
    )

    return companion_parameters  # as a dictionary


def baraffe_table_search(
    column: str, value: float, age: float, model: str, age_interp: bool = False
) -> Dict[str, float]:
    """Search Baraffe tables to find the companion entry given a column and value.

    Parameters
    ----------
    column: str
        Dictionary of (band: magnitude) pairs.
    age: float
        Age of star/system (Gyr).
    value: float
        Parameter value to find parameters for.
    model: str
        Year of Baraffe model to use [2003 (default), 2015].
    age_interp: bool
        Interpolate tables across age. Default=False.

    Returns
    -------
    companion_parameters: Dict[str, float]
        Companion parameters from Baraffe table, interpolated between the
        rows to the provided magnitude.

    """
    found_table, cols, model_age = age_table(age, model=model, age_interp=age_interp)
    if column not in cols:
        raise ValueError(
            "Column {0} not in Baraffe table (age={1}, model={2})".format(
                column, model_age, model
            )
        )

    return table_interpolation(found_table, column, value)


def table_interpolation(
    data: Dict[str, List[float]], ref_col: str, ref_value: float
) -> Dict[str, float]:
    """Interpolate table data from dictionary to the reference value.

    Parameters
    ----------
    data: dict
        Dictionary of table data. keys are the column headers.
    ref_value: float
        Value of reference parameter interpolating to.
    ref_col: str
        Column name string.
    age: float
        Age of star?system (Gyr).
    model: str
       Year of Baraffe model to use [2003 (default), 2015].

    Returns
    -------
    result_parameters: dict(str, float)
        Result from interpolation of each dict item to the reference.

    """
    column_reversed = False
    result_parameters = {}
    for key, y_data in data.items():
        x_data = data[ref_col]
        if x_data[-1] < x_data[0]:
            # Reverse data if not increasing.
            x_data = x_data[::-1]
            y_data = y_data[::-1]
            column_reversed = True

        result_parameters[key] = np.interp(ref_value, x_data, y_data)

        if isinstance(result_parameters[key], (np.ndarray, list)):
            result_parameters[key] = result_parameters[key][0]

    # Raising warning if value outside bounds of table
    result = np.interp(ref_value, x_data, y_data, left=-99999999, right=99999999)
    indicator = result * (-1) ** (column_reversed)
    if indicator == -99999999:
        warnings.warn(
            "Interpolated values are outside the lower bound of {0!s}.".format(ref_col)
        )
    elif indicator == 99999999:
        warnings.warn(
            "Interpolated values are outside the upper bound of {0!s}.".format(ref_col)
        )

    return result_parameters
