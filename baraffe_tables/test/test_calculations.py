import numpy as np
from baraffe_tables.calculations import  flux_mag_ratio

def test_flux_mag_ratio():
    """Test flux-magnitude ratio.

    A magnitude difference of 5 should be 100 (by definition).

    """
    # positive 5 magnitude difference
    assert np.allclose(flux_mag_ratio(1, 6), 100)

    # Negative 5 magnitude difference
    assert np.allclose(flux_mag_ratio(7, 2), 1. / 100)

