import numpy as np
from chondrotwin import ConfinedCompression


def test_creep_reaches_equilibrium():
    tissue = ConfinedCompression()
    sigma0 = 0.1e6
    sol = tissue.solve(sigma0, t_end=6 * tissue.tau)
    eq = tissue.equilibrium_strain(sigma0)
    assert np.isclose(sol.surface_displacement[-1] / tissue.h, eq, rtol=0.02)
    assert np.all(sol.pressure[-1] < 0.02 * sigma0)   # fluid pressure relaxed


def test_fluid_carries_load_initially():
    tissue = ConfinedCompression()
    sigma0 = 0.1e6
    sol = tissue.solve(sigma0, t_end=0.02 * tissue.tau)
    assert np.mean(sol.pressure[-1][:-5]) > 0.8 * sigma0


def test_characteristic_time_scaling():
    a = ConfinedCompression(h=1.5e-3); b = ConfinedCompression(h=3.0e-3)
    assert np.isclose(b.tau / a.tau, 4.0)
