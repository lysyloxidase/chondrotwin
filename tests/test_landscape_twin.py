import numpy as np
from chondrotwin import quasi_potential_langevin, schematic_landscape, MiniTwin, catabolic_bias
from chondrotwin.landscape import basin_depths


def test_schematic_landscape_has_minima_at_basins():
    X, Y, U = schematic_landscape()
    assert U.min() < 0.5 and U.max() <= 2.5


def test_langevin_landscape_has_two_basins_and_a_ridge():
    res = quasi_potential_langevin(n_walkers=800, n_steps=2500, burn_in=800, seed=1)
    d = basin_depths(res["U"], res["S_edges"], res["R_edges"])
    assert d["barrier_from_SOX9"] > 0.5 and d["barrier_from_RUNX2"] > 0.5


def test_catabolic_bias_sign():
    assert catabolic_bias(0.45) < 0      # physiological: anabolic
    assert catabolic_bias(1.0) > 0.5     # overload: catabolic


def test_twin_physiological_is_stable_and_overload_degrades():
    twin = MiniTwin()
    phys = MiniTwin.to_arrays(twin.run(0.05e6, n_cycles=40))
    over = MiniTwin.to_arrays(twin.run(0.15e6, n_cycles=40))
    assert phys["H_A"][-1] >= 0.95 * phys["H_A"][0]
    assert phys["S"][-1] > phys["R"][-1]
    assert over["H_A"][-1] < 0.6 * over["H_A"][0]
    assert over["R"][-1] > over["S"][-1]


def test_twin_remembers_overload_episode():
    twin = MiniTwin()
    recs = twin.run(lambda c: 0.15e6 if 10 <= c < 22 else 0.05e6, n_cycles=50)
    a = MiniTwin.to_arrays(recs)
    assert a["R"][-1] > a["S"][-1]           # committed to hypertrophy after the episode
    assert a["H_A"][-1] < a["H_A"][5]        # and matrix keeps degrading under normal load
