import numpy as np, pytest
from pylidarterrain import assert_north_up
from pylidarterrain.derive import slope_percent


def _synthetic_north_facing(n=40, res=10.0):
    """Ground rising to the south: row index increases southward, so
    elevation must increase with row."""
    z = np.tile(np.arange(n, dtype=float)[:, None] * res * 0.2, (1, n))
    aspect = np.zeros_like(z)          # 0 deg = north-facing, correct here
    return z, aspect


def test_north_up_guard_accepts_correct_orientation():
    z, aspect = _synthetic_north_facing()
    assert assert_north_up(z, aspect, res=10.0)


def test_north_up_guard_catches_the_flip():
    z, aspect = _synthetic_north_facing()
    with pytest.raises(AssertionError, match="south-up"):
        assert_north_up(z, np.full_like(aspect, 180.0), res=10.0)


def test_slope_percent_matches_degrees():
    assert slope_percent(45.0) == pytest.approx(100.0)
    assert slope_percent(np.degrees(np.arctan(0.42))) == pytest.approx(42.0)
