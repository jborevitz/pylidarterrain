"""Tests for crown-unit delineation and the hard scoring cap."""
import os
import numpy as np
import pytest
from pylidarterrain import units, stand

DATA = os.path.join(os.path.dirname(__file__), "data", "avonside_window.npz")


def test_chm_takes_the_maximum_per_cell():
    # two returns in cell (0,0), one each in (1,1) and (1,0)
    x = np.array([0.2, 0.8, 1.5, 0.5]); y = np.array([0.2, 0.8, 1.5, 1.5])
    h = np.array([3.0, 7.0, 9.0, 5.0])
    c, frac = units.chm(x, y, h, 0.0, 0.0, res=1.0, shape=(2, 2), median_size=0)
    assert c[0, 0] == pytest.approx(7.0)      # max of 3 and 7, not the last
    assert c[1, 1] == pytest.approx(9.0)
    assert c[1, 0] == pytest.approx(5.0)


def test_chm_reports_the_populated_fraction():
    """The number to check before trusting a CHM resolution."""
    x = np.array([0.5]); y = np.array([0.5]); h = np.array([4.0])
    _, frac = units.chm(x, y, h, 0.0, 0.0, res=1.0, shape=(2, 2), median_size=0)
    assert frac == pytest.approx(0.25)


def test_chm_row_zero_is_the_southern_edge():
    """Convention check: this module is origin='lower', unlike dem.py."""
    x = np.array([0.5, 0.5]); y = np.array([0.5, 1.5]); h = np.array([2.0, 20.0])
    c, _ = units.chm(x, y, h, 0.0, 0.0, res=1.0, shape=(2, 1), median_size=0)
    assert c[0, 0] < c[1, 0]


def test_cover_regime_splits_open_scattered_continuous():
    rng = np.random.default_rng(0)
    xs, ys, hs = [], [], []
    for j, cov in enumerate((0.0, 0.3, 0.9)):     # three cells along x
        n = 200
        x = rng.uniform(j * 10, j * 10 + 10, n)
        y = rng.uniform(0, 10, n)
        h = np.where(rng.random(n) < cov, rng.uniform(8, 15, n), 0.0)
        xs.append(x); ys.append(y); hs.append(h)
    x, y, h = map(np.concatenate, (xs, ys, hs))
    w = np.ones_like(x)
    reg, cover, top = units.cover_regime(x, y, h, w, 0.0, 0.0, cell=10.0,
                                         shape=(1, 3), cover_split=0.6)
    assert reg[0, 0] == 0 and reg[0, 1] == 1 and reg[0, 2] == 2


def test_cap_splits_a_wide_blob_and_bounds_the_diameter():
    """The cap is a rule about the unit of scoring, so it must bind."""
    seg = np.zeros((40, 40), np.int32)
    seg[5:35, 5:35] = 1                      # a 30x30 cell blob at res 1 m
    res = 1.0
    before = 2 * np.sqrt((seg == 1).sum() * res * res / np.pi)
    unit = units.cap_units(seg, res, cap_diameter=10.0, grid=10.0, min_area=5.0)
    areas = np.bincount(unit.ravel())[1:] * res * res
    after = 2 * np.sqrt(areas.max() / np.pi)
    assert before > 10.0
    assert unit.max() > 1
    assert after < before


def test_cap_leaves_small_segments_alone():
    seg = np.zeros((20, 20), np.int32)
    seg[8:12, 8:12] = 1                      # 4x4 m, under the cap
    unit = units.cap_units(seg, 1.0, cap_diameter=10.0, grid=10.0, min_area=2.0)
    assert unit.max() == 1


def test_crowns_separates_two_cones():
    yy, xx = np.mgrid[0:60, 0:60].astype(float)
    def cone(cx, cy):
        d = np.hypot(xx - cx, yy - cy)
        return np.clip(18.0 - 1.2 * d, 0, None)
    chm = np.maximum(cone(18, 30), cone(42, 30))
    seg = units.crowns(chm, np.ones_like(chm, bool), res=1.0,
                       min_separation=5.0, min_area=4.0, canopy_cut=2.0)
    assert seg.max() == 2


def test_binary_closing_bridges_gaps_in_the_canopy_mask():
    """Closing is done on the BINARY mask, never on the height surface."""
    chm = np.zeros((20, 20)); chm[5:15, 5:15] = 12.0
    chm[5:15, 10] = 0.0                       # a one-cell slot through it
    closed = units.crowns(chm, np.ones_like(chm, bool), 1.0,
                          min_separation=3.0, min_area=2.0, closing=True)
    raw = units.crowns(chm, np.ones_like(chm, bool), 1.0,
                       min_separation=3.0, min_area=2.0, closing=False)
    assert (closed > 0).sum() > (raw > 0).sum()


def test_closing_mask_defaults_to_the_segmentation_mask():
    """`closing_mask=None` must reproduce passing `mask` explicitly; passing
    the wider region is what changes the answer at a region boundary."""
    chm = np.zeros((20, 20)); chm[4:16, 4:16] = 12.0
    region = np.zeros_like(chm, bool); region[:, :12] = True
    a = units.crowns(chm, region, 1.0, min_area=2.0)
    b = units.crowns(chm, region, 1.0, min_area=2.0, closing_mask=region)
    assert np.array_equal(a, b)


def test_unit_returns_groups_every_return_once():
    unit = np.array([[1, 1], [2, 0]], np.int32)
    x = np.array([0.5, 1.5, 0.5, 1.5]); y = np.array([0.5, 0.5, 1.5, 1.5])
    h = np.array([1.0, 2.0, 3.0, 4.0]); w = np.ones(4)
    ids, rb, wb = units.unit_returns(unit, x, y, h, w, 0.0, 0.0, 1.0)
    assert list(ids) == [1, 2]
    assert sum(len(r) for r in rb) == 3          # the cell labelled 0 is dropped
    assert sorted(np.concatenate(rb)) == [1.0, 2.0, 3.0]


@pytest.mark.skipif(not os.path.exists(DATA), reason="window fixture absent")
def test_end_to_end_on_the_avonside_window():
    """Full pipeline on 150 m of real 0.8 pulses/m2 lidar at the field site.

    Values recorded from the verified farm-wide run; they pin the pipeline,
    they are not an accuracy claim.
    """
    d = np.load(DATA)
    x, y = d["x"].astype(float), d["y"].astype(float)
    h, w = d["h"].astype(float), d["w"].astype(float)
    x0 = y0 = -75.0
    chm, frac = units.chm(x, y, h, x0, y0, res=1.5, shape=(100, 100))
    assert frac == pytest.approx(0.895, abs=0.01)
    reg, cover, top = units.cover_regime(x, y, h, w, x0, y0, cell=10.0,
                                         shape=(15, 15))
    assert np.bincount(reg.ravel(), minlength=3).tolist() == [3, 210, 12]
    ri = np.clip(((np.arange(100) + .5) * 1.5 // 10).astype(int), 0, 14)
    seg = units.crowns(chm, reg[np.ix_(ri, ri)] == 1, 1.5,
                       closing_mask=np.ones_like(chm, bool))
    assert seg.max() == 124
    unit = units.cap_units(seg, 1.5, cap_diameter=10.0, grid=10.0)
    assert unit.max() == 236
    ids, rb, wb = units.unit_returns(unit, x, y, h, w, x0, y0, 1.5)
    idx = np.array([stand.dieback_index(r, ww) for r, ww in zip(rb, wb)])
    assert int(np.isfinite(idx).sum()) == 60
    assert np.nanmedian(idx) == pytest.approx(0.0716, abs=0.001)
    # the cap trades resolution for unassessed area, and that must stay visible
    assert np.isfinite(idx).mean() < 0.5
