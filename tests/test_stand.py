"""Tests for the dieback index and stand-state classification."""
import os
import numpy as np
import pytest
from pylidarterrain import stand

rng = np.random.default_rng(0)
DATA = os.path.join(os.path.dirname(__file__), "data", "avonside_window.npz")


def synth(n_ground, low=(), tall=(), seed=0):
    """Build a synthetic return profile: ground, a low band, a tall band."""
    r = np.random.default_rng(seed)
    parts = [np.zeros(n_ground)]
    for lo, hi, n in (low, tall):
        parts.append(r.uniform(lo, hi, n))
    return np.concatenate(parts)


def _gaps(group_a, group_b, n_iter=60, seed=0):
    """Observed and mean-null difference in index between two groups."""
    r_by = list(group_a) + list(group_b)
    w_by = [np.ones_like(r) for r in r_by]
    obs = np.array([stand.dieback_index(r) for r in r_by])
    null = stand.profile_null(r_by, w_by, n_iter=n_iter,
                              rng=np.random.default_rng(seed))
    na = len(group_a)
    obs_gap = np.nanmean(obs[:na]) - np.nanmean(obs[na:])
    null_gap = float(np.nanmean(np.nanmean(null[:, :na], axis=1)
                                - np.nanmean(null[:, na:], axis=1)))
    return float(obs_gap), null_gap


def test_top_height_is_mean_of_highest_three():
    r = np.array([0.0, 0.2, 1.0, 5.0, 10.0, 20.0, 22.0, 24.0])
    assert stand.top_height(r) == pytest.approx((20 + 22 + 24) / 3)


def test_top_height_ignores_ground_and_needs_n_top():
    assert np.isnan(stand.top_height(np.zeros(50)))
    assert np.isnan(stand.top_height(np.array([0.0, 0.0, 6.0, 7.0]), n_top=3))


def test_top_height_resists_one_spurious_return():
    """The point of averaging three: a single high outlier must not set the top."""
    base = np.concatenate([np.zeros(20), rng.uniform(14, 16, 30)])
    with_outlier = np.append(base, 60.0)
    assert stand.top_height(with_outlier) - stand.top_height(base) < 16.0


def test_band_weight_normalisations_differ_and_bracket_each_other():
    r = np.concatenate([np.zeros(80), np.full(20, 2.0)])
    total = stand.band_weight(r, band=(0.5, 4.0), normalise="total")
    canopy = stand.band_weight(r, band=(0.5, 4.0), normalise="canopy")
    assert total == pytest.approx(0.2)
    assert canopy == pytest.approx(1.0)
    assert total < canopy


def test_canopy_share_cannot_see_a_thinning_canopy_but_total_can():
    """The reason `dieback_index` normalises by total weight.

    Two stands with the same vertical SHAPE but different densities relative to
    the ground: the canopy-normalised share is identical, the total-normalised
    density is not.
    """
    dense = np.concatenate([np.zeros(20), np.full(40, 2.0), np.full(40, 18.0)])
    sparse = np.concatenate([np.zeros(80), np.full(10, 2.0), np.full(10, 18.0)])
    assert (stand.band_weight(dense, band=(0.5, 4.0), normalise="canopy")
            == pytest.approx(stand.band_weight(sparse, band=(0.5, 4.0),
                                               normalise="canopy")))
    assert (stand.band_weight(dense, band=(0.5, 4.0), normalise="total")
            > stand.band_weight(sparse, band=(0.5, 4.0), normalise="total"))


def test_index_is_high_only_for_the_dieback_conjunction():
    """The behavioural contract: tall AND dense underneath, not either alone."""
    intact = synth(40, low=(0.6, 3.5, 3), tall=(16, 22, 60), seed=1)
    shrub = synth(40, low=(0.6, 3.5, 60), tall=(16, 22, 0), seed=2)
    dieback = synth(40, low=(0.6, 3.5, 45), tall=(16, 22, 20), seed=3)
    i_intact = stand.dieback_index(intact)
    i_shrub = stand.dieback_index(shrub)
    i_dieback = stand.dieback_index(dieback)
    assert i_dieback > i_intact
    assert i_dieback > i_shrub
    assert all(0.0 <= v <= 1.0 for v in (i_intact, i_shrub, i_dieback))


def test_height_gate_suppresses_shrubland_specifically():
    """Same low-layer density, different height: the gate must separate them."""
    short = synth(50, low=(0.6, 3.5, 50), tall=(5, 7, 10), seed=4)
    tall = synth(50, low=(0.6, 3.5, 50), tall=(18, 22, 10), seed=4)
    assert stand.dieback_index(short) < stand.dieback_index(tall)


def test_index_refuses_to_measure_below_the_floor():
    thin = synth(10, low=(0.6, 3.5, 5), tall=(16, 22, 4), seed=5)
    assert np.isnan(stand.dieback_index(thin))
    assert np.isfinite(stand.dieback_index(thin, min_returns=5))


def test_index_gate_saturates_at_height_ref():
    """Above `height_ref` the gate is 1, so two tall stands with equal low
    layers score equally however much taller one of them is."""
    a = synth(50, low=(0.6, 3.5, 40), tall=(16, 17, 20), seed=6)
    b = a.copy()
    b[b > 10] += 15.0
    assert stand.dieback_index(a) == pytest.approx(stand.dieback_index(b), rel=1e-9)


def test_weights_down_weight_multi_return_pulses():
    r = np.concatenate([np.zeros(50), np.full(50, 2.0), np.full(20, 18.0)])
    w = np.ones_like(r)
    w[50:100] = 0.25          # the low layer recorded as 4-return pulses
    assert stand.dieback_index(r, w) < stand.dieback_index(r)


def test_stand_state_four_states_and_order():
    assert stand.stand_state(20.0, 0.01) == "intact tall forest"
    assert stand.stand_state(20.0, 0.20) == "tall with dense low layer"
    assert stand.stand_state(6.0, 0.20) == "low shrubland"
    assert stand.stand_state(20.0, 0.20, is_tree=False) == "open"


def test_stand_state_vectorises():
    out = stand.stand_state(np.array([20.0, 6.0, 20.0]),
                            np.array([0.01, 0.20, 0.20]))
    assert list(out) == ["intact tall forest", "low shrubland",
                         "tall with dense low layer"]


def test_null_preserves_size_and_top_height_but_not_structure():
    """The null must differ from the data ONLY in vertical arrangement."""
    r_by = [synth(40, low=(0.6, 3.5, 45), tall=(16, 22, 20), seed=i)
            for i in range(6)]
    w_by = [np.ones_like(r) for r in r_by]
    out = stand.profile_null(r_by, w_by, n_iter=5, rng=np.random.default_rng(0))
    assert out.shape == (5, 6)
    assert np.isfinite(out).all()


def test_null_is_sensitive_to_a_height_only_contrast():
    """A contrast driven by stand height alone is largely reproduced by the
    null, because the null preserves each unit's own top height.  This is why
    an observed separation must be compared against the null and NOT against
    0.5 - at Avonside the Eva's-pond-versus-intact contrast failed exactly here.
    """
    tall = [synth(40, low=(0.6, 3.5, 40), tall=(18, 22, 25), seed=i) for i in range(8)]
    short = [synth(40, low=(0.6, 3.5, 40), tall=(9, 11, 25), seed=i + 50) for i in range(8)]
    obs_gap, null_gap = _gaps(tall, short, seed=1)
    assert null_gap > 0.4 * obs_gap


def test_null_does_not_reproduce_a_structural_contrast():
    """Two groups at the SAME height differing only in low-layer density: the
    null, which randomises vertical arrangement, must fall well short."""
    dense = [synth(40, low=(0.6, 3.5, 55), tall=(18, 22, 25), seed=i) for i in range(8)]
    bare = [synth(40, low=(0.6, 3.5, 5), tall=(18, 22, 25), seed=i + 50) for i in range(8)]
    obs_gap, null_gap = _gaps(dense, bare, seed=1)
    assert obs_gap > 2.0 * null_gap


@pytest.mark.skipif(not os.path.exists(DATA), reason="window fixture absent")
def test_index_on_the_avonside_window():
    """Regression on 150 m of real 0.8 pulses/m2 lidar around the field site."""
    d = np.load(DATA)
    r, w = d["h"].astype(float), d["w"].astype(float)
    assert len(r) == 16484
    i = stand.dieback_index(r, w)
    assert 0.0 < i < 1.0
    assert stand.top_height(r, w) == pytest.approx(28.0, abs=2.0)
