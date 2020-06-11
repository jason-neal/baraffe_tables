import random

import numpy as np
import pytest

"""Testing interpolation,
Specifically interpolation between tables of different ages.

"""
from baraffe_tables.table_search import find_bounding_ages, interp_data_dicts
from baraffe_tables.table_search import model_ages_03, model_ages_15


@pytest.mark.parametrize("model_ages", [model_ages_03, model_ages_15])
@pytest.mark.parametrize("age", [0.012, 0.45, 0.89, 1.5, 4.6, 5.21, 9.6])
def test_find_bounding_ages(age, model_ages):
    # Test on input first
    assert age < float(model_ages[-1])
    assert float(model_ages[0]) < age

    lower_age, upper_age = find_bounding_ages(age, model_ages)

    assert float(lower_age) <= age
    assert age < float(upper_age)
    assert isinstance(lower_age, str)
    assert isinstance(upper_age, str)
    assert lower_age in model_ages
    assert upper_age in model_ages


@pytest.mark.parametrize("model_ages", [model_ages_03, model_ages_15])
@pytest.mark.parametrize("age", [0.89, 1.5, 9.6])
def test_find_bounding_works_when_unsorted_ages_given(age, model_ages):
    random.shuffle(model_ages)
    lower_age, upper_age = find_bounding_ages(age, model_ages)
    assert float(lower_age) <= age
    assert age < float(upper_age)
    assert lower_age in model_ages
    assert upper_age in model_ages


@pytest.mark.parametrize("age, lower_age, upper_age",
                         [(0.1, 0.5, 0.6), (0.6, 0.1, 0.5), (0.5, 0.6, 0.1)])
def test_table_interpolation_fails_with_bad_ages(age, lower_age, upper_age):
    with pytest.raises(AssertionError):
        interp_data_dicts(age, lower_age, None, upper_age, None)


@pytest.mark.parametrize("lower, upper", [(1., 2.), (4., 10.)])
def test_table_interpolation_in_middle(lower, upper):
    """Simple fixed example in middle of interpolation range"""
    age = (upper + lower) / 2.
    dict_lower = {"a": [1., 2., 3., 4., 5.], "b": [20., 30., 40., 50.]}
    dict_upper = {"a": [4., 5., 6., 7., 8.], "b": [40., 50., 60., 70.]}
    expected = {"a": np.array([2.5, 3.5, 4.5, 5.5, 6.5]), "b": np.array([30., 40., 50., 60.])}
    result = interp_data_dicts(age, lower, dict_lower, upper, dict_upper)

    for key, value in expected.items():
        assert np.allclose(result[key], value)
