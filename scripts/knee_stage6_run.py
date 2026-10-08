#!/usr/bin/env python
"""Stage-6 computations: loss profile over (c_over_max, k_cat) and the multiverse of analysis choices.

Usage:  python scripts/knee_stage6_run.py          (~15 min with 8 workers)
Writes data/derived/stage6_profile.csv and data/derived/stage6_multiverse.csv.
"""
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import product

import numpy as np
import pandas as pd

from chondrotwin import datasets
from chondrotwin.knee import calibrate as C
from chondrotwin.knee import uq

DERIVED = datasets.data_dir().parent / "derived"
SEEDS = (1, 2, 3)


def profile_point(args):
    params, seed = args
    st = C.evaluate(params, seed=seed)
    return params, seed, st


def mv_point(args):
    sim, model, seed = args
    v = uq.SIM_VARIANTS[sim]
    params = uq.STAGE1 if model == "stage 1" else uq.STAGE2
    return sim, model, seed, uq.cells_for(params, seed, v["slab_mm"], v["years"], v["sigmas"])


def main():
    t0 = time.perf_counter()
    ji = C.ji_cell_table()
    cal4, val4 = C.data_targets(ji, C.CAL_PATIENTS, pairs=C.PAIRS4), C.data_targets(ji, C.VAL_PATIENTS, pairs=C.PAIRS4)
    jobs = [(dict(zip(uq.PROFILE_GRID, vals)), s) for vals in product(*uq.PROFILE_GRID.values()) for s in SEEDS]
    rows = []
    with ProcessPoolExecutor() as ex:
        for params, seed, st in ex.map(profile_point, jobs):
            lr = C.model_log_ratios(st, C.PAIRS4)
            rows.append({**params, "seed": seed, "loss_cal": C.loss(lr, cal4), "loss_val": C.loss(lr, val4),
                         "sox9_retention_stage4": float(st.loc[4, "mean_S"] / st.loc[0, "mean_S"]) if st.loc[4, "n"] > 0 else np.nan,
                         "runx2_high_stage4": float(st.loc[4, "runx2_high"]), "senescent_stage4": float(st.loc[4, "senescent"]),
                         "senescent_stage0": float(st.loc[0, "senescent"]), "n_stage4": float(st.loc[4, "n"] or 0)})
    pd.DataFrame(rows).to_csv(DERIVED / "stage6_profile.csv", index=False)
    print(f"profile done ({time.perf_counter() - t0:.0f} s)", flush=True)

    cells_by = {}
    jobs = [(sim, model, s) for sim in uq.SIM_VARIANTS for model in ("stage 1", "stage 2") for s in SEEDS]
    with ProcessPoolExecutor() as ex:
        for sim, model, seed, cells in ex.map(mv_point, jobs):
            cells_by[(sim, model, seed)] = cells
    mv = uq.multiverse(cells_by, ji)
    mv.to_csv(DERIVED / "stage6_multiverse.csv", index=False)
    print(f"multiverse done: {len(mv)} universes, stage 2 better in {int(mv['stage 2 better (val)'].sum())} "
          f"({time.perf_counter() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
