#!/usr/bin/env python
"""Stage-2 grid calibration of the cell-state rules on GSE104782 calibration patients (OA1-OA5), run in parallel.

Usage:  python scripts/knee_stage2_grid.py [--workers N]
Writes data/derived/stage2_grid.csv (every grid point: loss on calibration and validation patients) and
data/derived/stage2_grid_tables.csv (the ABM stage tables).  ~7 min with 4 workers.
"""
import argparse
import os
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

from chondrotwin import datasets
from chondrotwin.knee import calibrate as C

DERIVED = datasets.data_dir().parent / "derived"


def run(params):
    st = C.evaluate(params)
    st = st.assign(**params)
    return params, st


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--refine", action="store_true", help="second-round grid (GRID_REFINE)")
    p.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 2) - 1)))
    a = p.parse_args()
    t = C.ji_cell_table()
    cal, val = C.data_targets(t, C.CAL_PATIENTS), C.data_targets(t, C.VAL_PATIENTS)
    rows, tables = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for params, st in ex.map(run, list(C.grid_points(C.GRID_REFINE if a.refine else C.GRID))):
            lr = C.model_log_ratios(st)
            rows.append({**params, "loss_cal": C.loss(lr, cal), "loss_val": C.loss(lr, val),
                         "n_stage4": int(st.loc[4, "n"]) if pd.notna(st.loc[4, "n"]) else 0})
            tables.append(st.reset_index())
            print(rows[-1], flush=True)
    DERIVED.mkdir(parents=True, exist_ok=True)
    tag = "_refine" if a.refine else ""
    pd.DataFrame(rows).sort_values("loss_cal").to_csv(DERIVED / f"stage2_grid{tag}.csv", index=False)
    pd.concat(tables, ignore_index=True).to_csv(DERIVED / f"stage2_grid{tag}_tables.csv", index=False)


if __name__ == "__main__":
    main()
