"""Synthetic-geometry tests: a plane with a known gradient, a canopy at a known
height, and cover fractions known by construction."""
import numpy as np
import pytest
from pylidarterrain.canopy import (ground_plane, heights_above_plane, fit_cell,
                                   canopy_top, stratum_cover, canopy_skewness,
                                   dieback_score, fit_cloud)

RNG = np.random.default_rng(7)
CX, CY, B0, BX, BY = 1005.0, 2005.0, 1000.0, 0.5, -0.2   # 0.5 = 26.6 deg fall


def _ground(n=400, noise=0.10):
    x = RNG.uniform(CX - 5, CX + 5, n)
    y = RNG.uniform(CY - 5, CY + 5, n)
    z = B0 + BX * (x - CX) + BY * (y - CY) + RNG.normal(0, noise, n)
    return x, y, z


def test_plane_recovers_a_steep_gradient():
    x, y, z = _ground()
    b0, bx, by = ground_plane(x, y, z, CX, CY)
    assert abs(b0 - B0) < 0.05
    assert abs(bx - BX) < 0.02 and abs(by - BY) < 0.02


def test_plane_is_robust_to_canopy_contamination():
    x, y, z = _ground()
    x[:40] += 0.0
    z[:40] += 15.0                      # 10% of points lifted into a crown
    b0, bx, by = ground_plane(x, y, z, CX, CY)
    assert abs(b0 - B0) < 0.10 and abs(bx - BX) < 0.05


def test_heights_are_measured_from_the_local_surface_not_the_cell_mean():
    x, y, z = _ground()
    plane = ground_plane(x, y, z, CX, CY)
    r = heights_above_plane(x, y, z, plane, CX, CY)
    assert abs(np.median(r)) < 0.05
    # the raw elevations span metres across the cell; the residuals do not
    assert np.ptp(z) > 4.0 and np.ptp(r) < 1.0


def _stand(n_g=300, n_c=200, top=20.0, spread=1.5):
    xg, yg, zg = _ground(n_g)
    xc = RNG.uniform(CX - 5, CX + 5, n_c)
    yc = RNG.uniform(CY - 5, CY + 5, n_c)
    zc = (B0 + BX * (xc - CX) + BY * (yc - CY) + top
          - np.abs(RNG.normal(0, spread, n_c)))
    x = np.r_[xg, xc]; y = np.r_[yg, yc]; z = np.r_[zg, zc]
    cl = np.r_[np.full(n_g, 2), np.full(n_c, 5)]
    return x, y, z, cl


def test_canopy_top_tracks_the_real_top_and_beats_the_component_mean():
    x, y, z, cl = _stand()
    plane = ground_plane(x[cl == 2], y[cl == 2], z[cl == 2], CX, CY)
    r = heights_above_plane(x, y, z, plane, CX, CY)
    w = np.ones(len(r))
    f = fit_cell(r, w)
    assert f.tree
    emp95 = np.percentile(r[r > 1.5], 95)
    assert abs(canopy_top(f) - emp95) < 2.0
    mean_c = f.xi + f.omega * (f.alpha / np.sqrt(1 + f.alpha ** 2)) * np.sqrt(2 / np.pi)
    assert mean_c < canopy_top(f)          # centroid sits below the top


def test_stratum_cover_matches_construction():
    r = np.array([0.0, 0.1, 2.0, 3.0, 5.0, 9.0, 12.0])
    w = np.ones(7)
    low, mid, relict = stratum_cover(r, w)
    assert low == pytest.approx(2 / 7) and mid == pytest.approx(1 / 7)
    assert relict == pytest.approx(2 / 7)


def test_closed_canopy_is_negatively_skewed():
    x, y, z, cl = _stand()
    plane = ground_plane(x[cl == 2], y[cl == 2], z[cl == 2], CX, CY)
    r = heights_above_plane(x, y, z, plane, CX, CY)
    assert canopy_skewness(r, np.ones(len(r))) < 0


def test_dieback_score_methods():
    r = np.r_[np.zeros(50), np.full(20, 2.5), np.full(10, 12.0)]
    w = np.ones(len(r))
    assert dieback_score(r, w, method="lowlayer") == pytest.approx(20 / 80)
    conj = dieback_score(r, w, method="conjunction")
    assert 0.0 < conj <= 1.0
    # no relict stems -> the conjunction collapses, the low-layer score does not
    r2 = np.r_[np.zeros(50), np.full(20, 2.5)]
    w2 = np.ones(len(r2))
    assert dieback_score(r2, w2, method="conjunction") == pytest.approx(0.0)
    assert dieback_score(r2, w2, method="lowlayer") > 0.2
    with pytest.raises(ValueError):
        dieback_score(r, w, method="nonsense")


def test_fit_cloud_grids_and_labels_cells():
    x1, y1, z1, cl1 = _stand()
    x2, y2, z2, cl2 = _stand()
    x2 = x2 + 10.0                         # a second cell, one grid step east
    cols = fit_cloud(np.r_[x1, x2], np.r_[y1, y2], np.r_[z1, z2],
                     np.ones(len(x1) + len(x2)), np.r_[cl1, cl2], res=10.0)
    assert len(cols["easting"]) == 2
    assert np.allclose(np.diff(np.sort(cols["easting"])), 10.0)
    assert cols["tree"].all()
    assert np.all(cols["canopy_top"] > 15.0)
    assert np.all(np.abs(cols["ground_z"] - B0) < 2.0)
