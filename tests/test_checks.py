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


def _plane(dz_east=0.0, dz_north=0.0, n=40, res=10.0):
    """A synthetic plane. Row 0 is the NORTHERN edge, so northing decreases with row."""
    cols = np.arange(n)[None, :] * res
    north = (n - 1 - np.arange(n))[:, None] * res
    return dz_east * cols + dz_north * north


def test_aspect_points_downhill_on_a_known_plane():
    from pylidarterrain.derive import slope_aspect
    # ground falling toward the EAST -> faces east -> aspect 90
    s, a = slope_aspect(_plane(dz_east=-0.2), 10.0)
    assert np.allclose(s[2:-2, 2:-2], np.degrees(np.arctan(0.2)))
    assert np.allclose(a[2:-2, 2:-2], 90.0)
    # ground falling toward the NORTH -> aspect 0
    s, a = slope_aspect(_plane(dz_north=-0.2), 10.0)
    assert np.allclose(a[2:-2, 2:-2], 0.0)
    # ground falling toward the SOUTH -> aspect 180
    s, a = slope_aspect(_plane(dz_north=0.2), 10.0)
    assert np.allclose(a[2:-2, 2:-2], 180.0)
    # ground falling toward the WEST -> aspect 270
    s, a = slope_aspect(_plane(dz_east=0.2), 10.0)
    assert np.allclose(a[2:-2, 2:-2], 270.0)


def test_priority_flood_removes_a_pit_without_moving_the_rim():
    from pylidarterrain.derive import fill_depressions
    z = _plane(dz_east=-0.05, n=30) + 100.0
    z[15, 15] -= 8.0                       # a closed pit
    filled = fill_depressions(z, 10.0)
    assert filled[15, 15] > z[15, 15]
    rim = np.ones_like(z, bool); rim[14:17, 14:17] = False
    assert np.allclose(filled[rim], z[rim], atol=1e-2)


def test_mfd_accumulation_conserves_and_concentrates():
    from pylidarterrain.derive import fill_depressions, flow_accumulation
    z = _plane(dz_east=-0.05, n=30) + 100.0
    acc = flow_accumulation(fill_depressions(z, 10.0), 10.0, method="mfd")
    assert acc.min() >= 1.0
    # on a uniform east-facing plane, accumulation must grow toward the east edge
    assert acc[15, -1] > acc[15, 1]
