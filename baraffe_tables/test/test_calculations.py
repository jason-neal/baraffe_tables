import numpy as np
from baraffe_tables.calculations import distance_modulus, flux_mag_ratio, absolute_magnitude, apparent_magnitude
import pytest
from baraffe_tables.db_queries import calculate_bv_temp


@pytest.mark.parametrize("mu, d", [
    (-4, 1.6),
    (-3, 2.5),
    (-2, 4.0),
    (-1, 6.3),
    (0, 10),
    (1, 16),
    (2, 25),
    (3, 40),
    (4, 63),
    (5, 100),
    (10, 1e3),
    (20, 1e5)])
def test_distance_modulus(d, mu):
    assert round(distance_modulus(d)) == mu


def test_flux_mag_ratio():
    """Test flux-magnitude ratio.

    A magnitude difference of 5 should be 100 (by definition).

    """
    # positive 5 magnitude difference
    assert np.allclose(flux_mag_ratio(1, 6), 100)

    # Negative 5 magnitude difference
    assert np.allclose(flux_mag_ratio(7, 2), 1. / 100)


@pytest.mark.parametrize("m", [5, 8.5, ])
@pytest.mark.parametrize("parallax", [10, 50, 200.1])
def test_apparent_absolute_magnitude_reversible(m, parallax):
    M = absolute_magnitude(parallax, m)

    assert np.allclose(apparent_magnitude(parallax, absolute_magnitude(parallax, m)), m)


@pytest.mark.parametrize("v, b_v, teff", [
    (4, -0.2, 18500),
    (4.5, 1.5, 3405.26),
    (6.0, 0.0, 9500),
    (0.0, 0.0, 9500),
    (18.0, 0.0, 9500),
    (5, 0.75, 5412.5)])

def test_calculate_bv_temp(v, b_v, teff):
    """Test calculate_bv_temp returns correct temperatures."""
    assert np.allclose(calculate_bv_temp(b_v + v, v), teff)
