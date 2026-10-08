#!/usr/bin/env python
"""Robust selection: the 5 best parameter sets of both grid rounds (by calibration loss, seed 1) are re-evaluated
with seeds 2 and 3 (the stage-1 model too, for comparison); the selected set has the lowest mean
calibration loss (OA1-OA5).  Validation loss is recorded
but not used.  Writes data/derived/stage2_selection.csv and stage2_selection_tables.csv.
"""
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from chondrotwin import datasets
from chondrotwin.knee import calibrate as C

DERIVED = datasets.data_dir().parent / "derived"
KEYS = ["damage_D50", "c_over_max", "k_cat", "baseline_runx2"]


def run(args):
    params, seed = args
    return params, seed, C.evaluate(params, seed=seed)


def main():
    g = pd.concat([pd.read_csv(DERIVED / "stage2_grid.csv").assign(baseline_runx2=0.0),
                   pd.read_csv(DERIVED / "stage2_grid_refine.csv")], ignore_index=True)
    top = g.sort_values("loss_cal").head(5)
    cand = [{k: float(r[k]) for k in KEYS} for _, r in top.iterrows()]
    cand.append(dict(damage_D50=0.0, c_over_max=0.9, k_cat=0.0, baseline_runx2=0.0))   # stage-1 model, same seeds
    t = C.ji_cell_table()
    cal, val = C.data_targets(t, C.CAL_PATIENTS), C.data_targets(t, C.VAL_PATIENTS)
    rows, tables = [], []
    with ProcessPoolExecutor() as ex:
        for params, seed, st in ex.map(run, [(p, s) for p in cand for s in (1, 2, 3)]):
            lr = C.model_log_ratios(st)
            rows.append({**params, "seed": seed, "loss_cal": C.loss(lr, cal), "loss_val": C.loss(lr, val), "n_stage4": st.loc[4, "n"]})
            tables.append(st.reset_index().assign(seed=seed, **params))
    d = pd.DataFrame(rows)
    summary = d.groupby(KEYS)[["loss_cal", "loss_val", "n_stage4"]].agg(["mean", "std"]).sort_values(("loss_cal", "mean"))
    print(summary.round(1).to_string())
    d.to_csv(DERIVED / "stage2_selection.csv", index=False)
    pd.concat(tables, ignore_index=True).to_csv(DERIVED / "stage2_selection_tables.csv", index=False)


if __name__ == "__main__":
    main()
