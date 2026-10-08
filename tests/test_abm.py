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


# ---------------------------------------------------------------- stage-2 mechanisms (defaults = 0.2-0.4 model)
from chondrotwin.abm import Rates


def _final(rates, days=90, sigma=0.12e6):
    return CartilageABM(rates=rates, **SMALL).run(Scenario("over", sigma0=sigma, days=days)).final()


def test_stage2_defaults_reproduce_the_original_model():
    a = _final(Rates(), days=60); b = _final(Rates.stage2(damage_D50=0.0, c_over_max=0.9, k_cat=0.0), days=60)
    assert a == b


def test_cumulative_damage_delays_senescence():
    assert _final(Rates.stage2(100.0, 0.9, 0.0))["senescent"] < 0.5 * _final(Rates())["senescent"]


def test_bounded_mechanical_bias_prevents_switching_without_inflammation():
    # max strain bias 0.3 < saddle-node 0.51: overload alone cannot flip SOX9-high cells
    f = CartilageABM(rates=Rates.stage2(1e6, 0.3, 0.0), **SMALL).run(Scenario("over", sigma0=0.12e6, days=90))
    assert (f.cells.query("~dead")["R"] > f.cells.query("~dead")["S"]).mean() < 0.05


def test_switch_independent_catabolism_degrades_matrix():
    no_cat = _final(Rates.stage2(1e6, 0.3, 0.0))["matrix_superficial"]
    cat = _final(Rates.stage2(1e6, 0.3, 1.0))["matrix_superficial"]
    assert cat < no_cat - 0.1


# ---------------------------------------------------------------- stage-3 heterogeneity (defaults = 0.2-0.5 model)
def test_stage3_defaults_reproduce_the_original_model():
    a = _final(Rates(), days=60)
    b = _final(Rates(switch_bias_sd=0.0, senescence_competent_fraction=1.0), days=60)
    assert a == b


def test_switch_heterogeneity_creates_runx2_high_cells_in_intact_tissue():
    res = CartilageABM(rates=Rates.stage2(0.0, 0.15, 0.5), **SMALL).run(Scenario("phys", sigma0=0.05e6, days=60))
    het = CartilageABM(rates=Rates(c_over_max=0.15, k_cat=0.5, switch_bias_sd=0.6), **SMALL).run(Scenario("phys", sigma0=0.05e6, days=60))
    frac = lambda r: float((r.cells.query("~dead")["R"] > r.cells.query("~dead")["S"]).mean())   # noqa: E731
    assert frac(res) < 0.005 and 0.01 < frac(het) < 0.2


def test_senescence_competence_caps_senescence():
    full = _final(Rates(), days=120)["senescent"]
    capped = _final(Rates(senescence_competent_fraction=0.2), days=120)["senescent"]
    assert capped < 0.25 and capped < full


# ---------------------------------------------------------------- stage-4 senescence logic (defaults = earlier model)
def test_stage4_defaults_reproduce_the_original_model():
    assert _final(Rates(), days=60) == _final(Rates(sen_overload_weight=1.0, senescent_death_rate=0.0), days=60)


def test_inflammation_gating_delays_overload_senescence():
    # without cytokines (no senescent / hypertrophic sources yet, no synovial inflammation) gating removes most of it
    assert _final(Rates(sen_overload_weight=0.0), days=60)["senescent"] < 0.5 * _final(Rates(), days=60)["senescent"]


def test_senescent_cells_can_die():
    f0 = _final(Rates(), days=120); f1 = _final(Rates(senescent_death_rate=0.05), days=120)
    assert f1["senescent"] < f0["senescent"] and f1["dead_apoptosis"] > f0["dead_apoptosis"]


# ---------------------------------------------------------------- stage-5 joint-level inflammation (defaults = earlier model)
def test_stage5_defaults_reproduce_the_original_model():
    a = CartilageABM(**SMALL).run(Scenario("o", sigma0=0.12e6, days=60)).final()
    b = CartilageABM(rates=Rates(k_inf=0.0, sen_infl_gain=2.0), **SMALL).run(Scenario("o", sigma0=0.12e6, days=60, C_joint=0.0)).final()
    assert a == b


def test_joint_level_inflammation_reaches_the_whole_tissue():
    sim = CartilageABM(**SMALL)
    sim.run(Scenario("j", sigma0=0.05e6, days=10, C_joint=0.5))
    assert float(sim.agents.C[sim.agents.alive].min()) >= 0.5


def test_inflammation_driven_catabolism_degrades_matrix_without_overload():
    kw = dict(sigma0=0.05e6, days=120, C_joint=1.0)
    base = CartilageABM(rates=Rates.stage2(0.0, 0.15, 0.5), **SMALL).run(Scenario("j", **kw)).final()
    cat = CartilageABM(rates=Rates(c_over_max=0.15, k_cat=0.5, k_inf=2.0), **SMALL).run(Scenario("j", **kw)).final()
    assert cat["matrix_superficial"] < base["matrix_superficial"] - 0.02
