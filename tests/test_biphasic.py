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


def _terzaghi(z, t, sigma0, H_A, h, tau, n_terms=400):
    """Analytical creep solution (Mow et al., 1980): surface displacement and pore pressure."""
    n = np.arange(n_terms); M = (2 * n + 1) * np.pi / 2
    b = 2 * sigma0 * h * (-1) ** n / (H_A * M ** 2)
    decay = np.exp(-M ** 2 * t / tau)
    u = (b * np.sin(np.outer(z / h, M)) * decay).sum(1) - sigma0 * z / H_A
    p = (H_A * b * (M / h) * np.cos(np.outer(z / h, M)) * decay).sum(1)
    return u, p


def test_matches_analytical_creep_solution():
    tissue = ConfinedCompression()
    sigma0 = 0.1e6
    sol = tissue.solve(sigma0, t_end=3 * tissue.tau, n_out=301)
    for frac in (0.02, 0.1, 0.3, 1.0, 3.0):
        i = int(np.argmin(np.abs(sol.t / sol.tau - frac)))
        u, p = _terzaghi(sol.z, sol.t[i], sigma0, tissue.H_A, tissue.h, tissue.tau)
        assert abs(sol.u[i, -1] - u[-1]) < 0.015 * sigma0 * tissue.h / tissue.H_A
        assert np.max(np.abs(sol.pressure[i, :-1] - p[:-1])) < 0.05 * sigma0
