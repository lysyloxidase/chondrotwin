#!/usr/bin/env python
"""Run the multi-agent cartilage model for a set of scenarios and write figures + a scaling table.

Usage:  python scripts/run_abm.py [OUTDIR] [--days N] [--cells-per-mm3 D] [--quick]
"""
import argparse
import os
import time

from chondrotwin.abm import CartilageABM, CartilageSlab, Scenario, benchmark, scaling_table, tile_plan
from chondrotwin.abm.plotting import figure_cross_sections, figure_depth_profiles, figure_time_courses


def main():
    p = argparse.ArgumentParser()
    p.add_argument("outdir", nargs="?", default="figures/abm")
    p.add_argument("--days", type=int, default=365)
    p.add_argument("--cells-per-mm3", type=float, default=1.4e4)
    p.add_argument("--quick", action="store_true", help="smaller slab and fewer days")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    slab = CartilageSlab(Lx=0.5, Ly=0.5, h=2.0, nx=5, ny=5, nz=20) if a.quick else CartilageSlab()
    days = 120 if a.quick else a.days
    scenarios = [
        Scenario("physiological", sigma0=0.05e6, days=days),
        Scenario("overload", sigma0=0.12e6, days=days),
        Scenario("overload episode", sigma0=lambda d: 0.12e6 if 30 <= d < 90 else 0.05e6, days=days),
        Scenario("focal overload", sigma0=0.05e6, focal_factor=1.5, days=days),
        Scenario("inflamed joint", sigma0=0.05e6, C_synovial=1.0, days=days),
        Scenario("senescent seeding 10%", sigma0=0.05e6, initial_senescent_fraction=0.10, days=days),
    ]
    results = []
    for sc in scenarios:
        t0 = time.perf_counter()
        sim = CartilageABM(slab=slab, cells_per_mm3=a.cells_per_mm3, seed=1)
        res = sim.run(sc)
        f = res.final()
        print(f"{sc.name:22s} agents {res.n_initial:6d}  {time.perf_counter() - t0:5.1f}s  articular {f['articular']:.2f}  "
              f"hypertrophic {f['hypertrophic']:.2f}  senescent {f['senescent']:.2f}  dead {f['dead_ferroptosis'] + f['dead_apoptosis']:.2f}  "
              f"M_sup {f['matrix_superficial']:.2f}  OARSI-like {f['oarsi_like']:.2f}")
        results.append(res)
        res.frame().to_csv(os.path.join(a.outdir, f"timecourse_{sc.name.replace(' ', '_').replace('%', '')}.csv"))
    figure_cross_sections(results[:4], os.path.join(a.outdir, "abm_cross_sections.png"))
    figure_time_courses(results, os.path.join(a.outdir, "abm_time_courses.png"))
    figure_depth_profiles(results, os.path.join(a.outdir, "abm_depth_profiles.png"))
    b = benchmark(n_steps=20)
    table = scaling_table(b)
    plan = tile_plan()
    text = (f"Benchmark on this machine: {b.n_agents} agents, {b.seconds_per_step * 1e3:.1f} ms per step\n\n{table}\n\n"
            f"Tiling plan for a whole knee: {plan}\n")
    with open(os.path.join(a.outdir, "scaling.md"), "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
