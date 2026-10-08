import numpy as np
from chondrotwin.abm import CartilageABM, CartilageSlab, Scenario, cells_for_volume, memory_gb, tile_plan, benchmark, scaling_table

SMALL = dict(slab=CartilageSlab(Lx=0.5, Ly=0.5, h=2.0, nx=5, ny=5, nz=20), seed=3)


def run(sc):
    return CartilageABM(**SMALL).run(sc)


def test_seeding_density_and_zones():
    sim = CartilageABM(**SMALL)
    ag = sim.agents
    expected = 1.4e4 * 0.5 * 0.5 * 2.0
    assert abs(ag.n_initial - expected) / expected < 0.05
    sup = ag.zone[ag.alive] == 0; deep = ag.zone[ag.alive] == 2
    # superficial zone is thin but dense: more cells per mm of depth than the deep zone
    assert sup.sum() / (0.15 * 2.0) > deep.sum() / (0.4 * 2.0)
    assert np.all(ag.S[ag.alive] > ag.R[ag.alive])


def test_physiological_loading_is_homeostatic():
    f = run(Scenario("phys", sigma0=0.05e6, days=120)).final()
    assert f["articular"] > 0.9 and f["matrix_superficial"] > 0.9 and f["dead_ferroptosis"] + f["dead_apoptosis"] < 0.03


def test_overload_degenerates_surface_first():
    res = run(Scenario("over", sigma0=0.12e6, days=120))
    f = res.final()
    assert f["matrix_superficial"] < 0.5 < f["matrix_deep"]
    assert f["senescent"] + f["dead_ferroptosis"] + f["dead_apoptosis"] > 0.3
    assert f["oarsi_like"] > 2.0
    c = res.cells
    assert (c.loc[c["zone"] == 0, "dead"].mean()) > (c.loc[c["zone"] == 2, "dead"].mean())


def test_focal_overload_is_spatially_localized():
    res = run(Scenario("focal", sigma0=0.05e6, focal_factor=1.5, focal_width=0.1, days=120))
    c = res.cells; s = res.slab
    r = np.hypot(c["x"] - s.Lx / 2, c["y"] - s.Ly / 2)
    sup = c["zone"] == 0
    assert c.loc[sup & (r < 0.1), "M"].mean() < c.loc[sup & (r > 0.2), "M"].mean() - 0.2


def test_inflammation_increases_senescence():
    phys = run(Scenario("phys", sigma0=0.05e6, days=120)).final()
    infl = run(Scenario("infl", sigma0=0.05e6, C_synovial=1.0, days=120)).final()
    assert infl["senescent"] > phys["senescent"] and infl["cytokine_mean"] > phys["cytokine_mean"]


def test_overload_episode_leaves_persistent_damage():
    res = run(Scenario("episode", sigma0=lambda d: 0.12e6 if 20 <= d < 60 else 0.05e6, days=150))
    f = res.frame()
    assert f.loc[149, "matrix_superficial"] < f.loc[60, "matrix_superficial"] + 0.05   # no recovery after the episode
    assert f.loc[149, "oarsi_like"] > 1.0


def test_proliferation_respects_capacity():
    sim = CartilageABM(**SMALL)
    res = sim.run(Scenario("over", sigma0=0.12e6, days=120))
    assert sim.agents.n_used <= sim.agents.capacity
    assert res.final()["n_clusters"] >= 1


def test_scaling_helpers():
    assert cells_for_volume(1.0) == 14000
    assert memory_gb(2.8e8) > 10
    plan = tile_plan()
    assert plan["n_rve"] == 10000 and plan["total_cells"] == 280_000_000
    b = benchmark(n_steps=2, Lx=0.3, Ly=0.3, h=1.0)
    table = scaling_table(b)
    assert "whole knee" in table and b.seconds_per_step > 0
