#!/usr/bin/env python
"""Numerical verification of the chondrotwin models (reproduces the numbers in VALIDATION.md).

Usage:  python scripts/validate.py [--quick] [--out validation_results.md]

--quick skips the multi-seed agent-based runs and ablations (~1 min instead of ~8 min).
"""
import argparse
import os
import time
from dataclasses import replace

import numpy as np
import pandas as pd

from chondrotwin import ConfinedCompression, GrowthPlateModel, MiniTwin, ToggleSwitch, quasi_potential_langevin
from chondrotwin.abm import CartilageABM, CartilageSlab, Rates, Scenario
from chondrotwin.landscape import basin_depths

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
lines = []


def out(s=""):
    print(s)
    lines.append(s)


def toggle():
    out("## 1. SOX9-RUNX2 toggle switch")
    sw = ToggleSwitch()
    res = max(np.abs(sw.rhs((fp.S, fp.R), lam)).max() for lam in np.linspace(-0.5, 0.5, 11) for fp in sw.fixed_points(lam))
    lo, hi = 0.3, 0.9
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if sw.is_bistable(mid) else (lo, mid)
    amps = np.arange(0.40, 0.80, 0.01)
    a_min = next(float(a) for a in amps if sw.switches_permanently(lam_pulse=float(a)))
    out(f"- max |f(x*)| over all fixed points, lam in [-0.5, 0.5]: {res:.1e}")
    out(f"- saddle-node bifurcation (bisection): lam* = +/-{lo:.4f} (docs: 0.51)")
    out(f"- smallest pulse amplitude (20 time units) that switches permanently: {a_min:.2f} (> lam*, as expected)")
    out(f"- bistable range for n = 1, 1.5, 2: {[ToggleSwitch(n=n).bistable_range() for n in (1.0, 1.5, 2.0)]}")
    out()


def growth_plate():
    out("## 2. IHH-PTHrP growth plate: time and grid convergence of the boundary")
    out("| grid | WT boundary | WT interp. | Jansen-like | Jansen interp. | converged (WT, J) |")
    out("|---|---|---|---|---|---|")
    for n in (101, 201, 401):
        row = []
        conv = []
        for m in (replace(GrowthPlateModel(), n_grid=n), replace(GrowthPlateModel(), n_grid=n).jansen_like()):
            s = m.steady_state()
            H = s.H; i = int(np.argmax(H >= 0.5))
            xi = s.x[i - 1] + (0.5 - H[i - 1]) * (s.x[i] - s.x[i - 1]) / (H[i] - H[i - 1])
            row += [f"{s.boundary():.4f}", f"{xi:.4f}"]; conv.append(s.converged)
        out(f"| {n} | " + " | ".join(row) + f" | {conv} |")
    old = GrowthPlateModel().jansen_like().steady_state(t_end=12.0)
    out(f"- with the former default t_end = 12 the Jansen-like run was not converged (converged={old.converged}); "
        f"boundary was already {old.boundary():.3f}")
    out()


def biphasic():
    out("## 3. Biphasic confined-compression creep vs analytical solution (Mow et al., 1980)")
    out("| n_z | max rel. error, surface displacement | max error, pressure / sigma0 |")
    out("|---|---|---|")
    for nz in (51, 101, 201):
        t = ConfinedCompression(n_z=nz); s0 = 0.1e6
        sol = t.solve(s0, t_end=3 * t.tau, n_out=301)
        n = np.arange(400); M = (2 * n + 1) * np.pi / 2; b = 2 * s0 * t.h * (-1) ** n / (t.H_A * M ** 2)
        eu, ep = 0.0, 0.0
        for frac in (0.02, 0.1, 0.3, 1.0, 3.0):
            i = int(np.argmin(np.abs(sol.t / sol.tau - frac))); d = np.exp(-M ** 2 * sol.t[i] / t.tau)
            u = (b * np.sin(np.outer(sol.z / t.h, M)) * d).sum(1) - s0 * sol.z / t.H_A
            p = (t.H_A * b * (M / t.h) * np.cos(np.outer(sol.z / t.h, M)) * d).sum(1)
            eu = max(eu, abs(sol.u[i, -1] - u[-1]) / (s0 * t.h / t.H_A)); ep = max(ep, np.abs(sol.pressure[i, :-1] - p[:-1]).max() / s0)
        out(f"| {nz} | {eu:.4f} | {ep:.4f} |")
    out("- first-order convergence (error halves with dz), from the one-sided traction boundary condition")
    out()


def landscape():
    out("## 4. Langevin quasi-potential: inter-basin equilibration and sensitivity to the initial split")
    out("| lam | noise | start SOX9-high | crossed | barrier from SOX9 | barrier from RUNX2 | U_min(SOX9) - U_min(RUNX2) |")
    out("|---|---|---|---|---|---|---|")
    for lam, noise in ((0.0, 0.18), (0.0, 0.25), (0.35, 0.18), (0.35, 0.25)):
        for f in (0.5, 0.9, 0.1):
            r = quasi_potential_langevin(lam=lam, noise=noise, start_fraction_sox9=f)
            d = basin_depths(r["U"], r["S_edges"], r["R_edges"])
            out(f"| {lam} | {noise} | {f} | {r['crossed']:.3f} | {d['barrier_from_SOX9']:.2f} | {d['barrier_from_RUNX2']:.2f} | "
                f"{d['U_min_SOX9_high'] - d['U_min_RUNX2_high']:+.2f} |")
    out()


def twin():
    out("## 5. Digital-twin loop: bias after an overload episode (cycles 10-21 at 0.15 MPa)")
    a = MiniTwin.to_arrays(MiniTwin().run(lambda c: 0.15e6 if 10 <= c < 22 else 0.05e6, n_cycles=60))
    out(f"- lam at cycles 22-30: {np.round(a['lam'][22:31], 3).tolist()}")
    out(f"- H_A (MPa) at cycles 0, 9, 21, 30, 59: {np.round(a['H_A'][[0, 9, 21, 30, 59]] / 1e6, 3).tolist()}; "
        f"final fate {'hypertrophic' if a['R'][-1] > a['S'][-1] else 'articular'}")
    out()


def abm_reproducibility():
    out("## 6. Agent-based model: reproducibility of the reference time courses (data/reference/abm/*.csv)")
    sc = [Scenario("physiological", sigma0=0.05e6), Scenario("overload", sigma0=0.12e6),
          Scenario("overload episode", sigma0=lambda d: 0.12e6 if 30 <= d < 90 else 0.05e6),
          Scenario("focal overload", sigma0=0.05e6, focal_factor=1.5), Scenario("inflamed joint", sigma0=0.05e6, C_synovial=1.0),
          Scenario("senescent seeding 10%", sigma0=0.05e6, initial_senescent_fraction=0.10)]
    for s in sc:
        path = os.path.join(ROOT, "data", "reference", "abm", f"timecourse_{s.name.replace(' ', '_').replace('%', '')}.csv")
        new = CartilageABM(seed=1).run(s).frame()
        ref = pd.read_csv(path, index_col="day")
        out(f"- {s.name}: max |difference| = {(new - ref).abs().max().max():.1e}")
    out()


def abm_seeds(n_seeds=5):
    out(f"## 7. Agent-based model: variability over {n_seeds} seeds (1 x 1 x 2 mm, 28,000 agents, 365 days; mean +/- SD)")
    scs = {"physiological": dict(sigma0=0.05e6), "overload": dict(sigma0=0.12e6),
           "overload episode": dict(sigma0=lambda d: 0.12e6 if 30 <= d < 90 else 0.05e6),
           "focal overload": dict(sigma0=0.05e6, focal_factor=1.5), "inflamed joint": dict(sigma0=0.05e6, C_synovial=1.0),
           "senescent seeding 10%": dict(sigma0=0.05e6, initial_senescent_fraction=0.10)}
    out("| scenario | articular | senescent | dead | matrix superficial | matrix deep | OARSI-like |")
    out("|---|---|---|---|---|---|---|")
    for name, kw in scs.items():
        rows = []
        for seed in range(n_seeds):
            f = CartilageABM(seed=seed).run(Scenario(name, **kw)).final()
            rows.append([f["articular"], f["senescent"], f["dead_ferroptosis"] + f["dead_apoptosis"], f["matrix_superficial"],
                         f["matrix_deep"], f["oarsi_like"]])
        m, sd = np.mean(rows, 0), np.std(rows, 0, ddof=1)
        out(f"| {name} | " + " | ".join(f"{a:.3f} ± {b:.3f}" for a, b in zip(m, sd)) + " |")
    out()


def abm_ablations(n_seeds=3):
    out(f"## 8. Agent-based model: ablations of the 'emergent' behaviours ({n_seeds} seeds, mean)")
    out("| question | variant | matrix sup. | matrix deep | senescent | dead | OARSI-like | dead sup. / deep |")
    out("|---|---|---|---|---|---|---|---|")
    ep = dict(sigma0=lambda d: 0.12e6 if 30 <= d < 90 else 0.05e6)
    cases = [
        ("senescent seeding", "physiological, no seeding", dict(sigma0=0.05e6), {}),
        ("senescent seeding", "10 % senescent", dict(sigma0=0.05e6, initial_senescent_fraction=0.1), {}),
        ("senescent seeding", "10 % senescent, SASP = 0", dict(sigma0=0.05e6, initial_senescent_fraction=0.1), dict(rates=Rates(sasp_production=0.0))),
        ("senescent seeding", "30 % senescent", dict(sigma0=0.05e6, initial_senescent_fraction=0.3), {}),
        ("senescent seeding", "30 % senescent, SASP = 0", dict(sigma0=0.05e6, initial_senescent_fraction=0.3), dict(rates=Rates(sasp_production=0.0))),
        ("damage after episode", "default", ep, {}),
        ("damage after episode", "no matrix -> stiffness feedback", ep, dict(slab=CartilageSlab(matrix_power=0.0))),
        ("damage after episode", "monostable switch (n = 1)", ep, dict(toggle=ToggleSwitch(n=1.0))),
        ("surface-first", "overload, modulus 0.4 -> 1.2 MPa", dict(sigma0=0.12e6), {}),
        ("surface-first", "overload, uniform modulus 0.8 MPa", dict(sigma0=0.12e6), dict(slab=CartilageSlab(H_A_surface=0.8e6, H_A_deep=0.8e6))),
    ]
    for q, label, skw, kw in cases:
        vals = []
        for seed in range(n_seeds):
            sim = CartilageABM(seed=seed, **kw); f = sim.run(Scenario(label, **skw)).final(); c = sim.snapshot()
            vals.append([f["matrix_superficial"], f["matrix_deep"], f["senescent"], f["dead_ferroptosis"] + f["dead_apoptosis"],
                         f["oarsi_like"], c.loc[c.zone == 0, "dead"].mean(), c.loc[c.zone == 2, "dead"].mean()])
        m = np.mean(vals, 0)
        out(f"| {q} | {label} | {m[0]:.3f} | {m[1]:.3f} | {m[2]:.3f} | {m[3]:.3f} | {m[4]:.2f} | {m[5]:.2f} / {m[6]:.2f} |")
    out()


def abm_memory():
    import tracemalloc
    out("## 9. Agent-based model: memory per agent")
    sim = CartilageABM(); ag = sim.agents; sc = Scenario(days=2)
    sim.step(0, sc)
    tracemalloc.start(); sim.step(1, sc); _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
    out(f"- state: {ag.memory_bytes() / ag.capacity:.0f} B per slot, {ag.memory_bytes() / ag.n_initial:.0f} B per initial agent "
        f"(capacity {ag.capacity / ag.n_initial:.1f}x)")
    out(f"- peak temporaries during one step: {peak / ag.n_initial:.0f} B per agent")
    out()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--quick", action="store_true")
    p.add_argument("--out", default=None)
    a = p.parse_args()
    t0 = time.perf_counter()
    out("# chondrotwin numerical verification")
    out()
    toggle(); growth_plate(); biphasic(); landscape(); twin(); abm_reproducibility()
    if not a.quick:
        abm_seeds(); abm_ablations()
    abm_memory()
    out(f"_total time {time.perf_counter() - t0:.0f} s_")
    if a.out:
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
