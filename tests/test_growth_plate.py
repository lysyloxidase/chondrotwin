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
    sol_j = base.jansen_like().steady_state()
    assert sol_j.converged
    xb_j = sol_j.boundary()
    assert xb_j > xb_wt + 0.05


def test_boundary_is_grid_converged():
    """Boundary positions quoted in the docs (0.205, 0.37) change by < 0.01 between 101 and 201 nodes (401: scripts/validate.py)."""
    from dataclasses import replace
    for model in (GrowthPlateModel(), GrowthPlateModel().jansen_like()):
        a = model.steady_state().boundary()
        b = replace(model, n_grid=101).steady_state().boundary()
        assert abs(a - b) < 0.01


def test_blomstrand_like_commits_everywhere():
    sol = GrowthPlateModel().blomstrand_like().steady_state()
    assert np.all(sol.H > 0.99)
    assert sol.proliferative_fraction() == 0.0
