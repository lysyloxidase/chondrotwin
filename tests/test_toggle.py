import numpy as np
from chondrotwin import ToggleSwitch


def test_bistable_at_zero_bias():
    m = ToggleSwitch()
    fps = m.fixed_points(0.0)
    assert sum(fp.stable for fp in fps) == 2
    assert sum(not fp.stable for fp in fps) == 1  # the saddle
    labels = {fp.label for fp in fps if fp.stable}
    assert labels == {"articular (SOX9-high)", "hypertrophic (RUNX2-high)"}


def test_monostable_under_strong_bias():
    m = ToggleSwitch()
    fps = m.fixed_points(0.9)
    stable = [fp for fp in fps if fp.stable]
    assert len(stable) == 1 and stable[0].R > stable[0].S
    fps = m.fixed_points(-0.9)
    stable = [fp for fp in fps if fp.stable]
    assert len(stable) == 1 and stable[0].S > stable[0].R


def test_bistable_range_is_symmetric_for_symmetric_parameters():
    rng = ToggleSwitch().bistable_range()
    assert rng is not None
    lo, hi = rng
    assert abs(lo + hi) < 0.02 and hi > 0.3


def test_fixed_points_are_true_steady_states():
    m = ToggleSwitch()
    for lam in (-0.5, 0.0, 0.3):
        for fp in m.fixed_points(lam):
            assert np.allclose(m.rhs((fp.S, fp.R), lam), 0.0, atol=1e-6)


def test_transient_pulse_switches_permanently():
    m = ToggleSwitch()
    t, S, R, lam_t = m.pulse(lam_pulse=0.8)
    assert R[0] < S[0]            # starts articular
    assert lam_t[-1] == 0.0       # pulse is over
    assert R[-1] > S[-1]          # ...but the cell stays RUNX2-high (hysteresis)


def test_weak_pulse_does_not_switch():
    m = ToggleSwitch()
    assert not m.switches_permanently(lam_pulse=0.2)


def test_no_cooperativity_no_bistability():
    m = ToggleSwitch(n=1.0)
    assert m.bistable_range() is None
