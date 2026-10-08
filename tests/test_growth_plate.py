import numpy as np
from chondrotwin import GrowthPlateModel


def test_steady_state_profiles_are_monotone_and_bounded():
    sol = GrowthPlateModel().steady_state()
    assert sol.converged
    assert np.all(sol.P >= 0) and np.all(sol.I >= 0)
    assert np.all((sol.H >= 0) & (sol.H <= 1))
    assert np.all(np.diff(sol.P) <= 1e-9)        # PTHrP decreases away from its source
    assert np.all(np.diff(sol.H) >= -1e-9)       # commitment increases toward the bone


def test_boundary_in_interior():
    sol = GrowthPlateModel().steady_state()
    xb = sol.boundary()
    assert 0.05 < xb < 0.6


def test_jansen_like_shifts_boundary_distally():
    base = GrowthPlateModel()
    xb_wt = base.steady_state().boundary()
    xb_j = base.jansen_like().steady_state().boundary()
    assert xb_j > xb_wt + 0.05


def test_blomstrand_like_commits_everywhere():
    sol = GrowthPlateModel().blomstrand_like().steady_state()
    assert np.all(sol.H > 0.99)
    assert sol.proliferative_fraction() == 0.0
