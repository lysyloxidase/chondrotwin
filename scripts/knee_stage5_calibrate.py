#!/usr/bin/env python
"""Stage-5 calibration: GSE104782 OA1-OA5 (within-joint stages, SenMayo target) + GSE114007 raw counts (between-joint
OA / normal), multi-seed selection and freezing.

1. grid (chondrotwin.knee.calibrate.GRID_STAGE5, seed 1, in parallel); loss = cell loss (PAIRS4) + between-joint loss;
2. the 5 best sets are re-run with seeds 2 and 3, together with the stage-1 and stage-2 models; the selected set has
   the lowest mean calibration loss (validation patients OA6-OA10 recorded, not used);
3. the selected parameters are written to docs/stage5_frozen_parameters.json and its SHA-256 is appended to
   docs/preregistration_stage5.md, before the GSE55235 / GSE55457 expression values (independent4.py) are read.

Usage:  python scripts/knee_stage5_calibrate.py          (~10 min)
"""
import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from chondrotwin import datasets
from chondrotwin.knee import calibrate as C

ROOT = Path(__file__).resolve().parents[1]
DERIVED = datasets.data_dir().parent / "derived"
STAGE1 = dict(damage_D50=0.0, c_over_max=0.9, k_cat=0.0, switch_bias_sd=0.0, senescence_competent_fraction=1.0,
              sen_overload_weight=1.0, senescent_death_rate=0.0, k_inf=0.0, sen_infl_gain=2.0, C_joint=0.0,
              baseline_runx2=0.0)
STAGE2 = dict(STAGE1, c_over_max=0.15, k_cat=0.5)
KEYS = list(STAGE1)


def run(args):
    params, seed = args
    return params, seed, C.evaluate_joint(params, seed=seed)


def score(res, cal, val, btw):
    st, oa, nrm = res
    lr = C.model_log_ratios(st, C.PAIRS4)
    b = C.between_loss(oa, nrm, btw)
    return {"loss_cal": C.loss(lr, cal) + b, "loss_val": C.loss(lr, val) + b, "loss_between": b,
            "loss_cells_val": C.loss(lr, val), "n_stage4": st.loc[4, "n"],
            **{f"oa_{k}": v for k, v in oa.items()}, **{f"normal_{k}": v for k, v in nrm.items()}}


def main():
    t0 = time.perf_counter()
    cells = C.ji_cell_table()
    cal = C.data_targets(cells, C.CAL_PATIENTS, pairs=C.PAIRS4)
    val = C.data_targets(cells, C.VAL_PATIENTS, pairs=C.PAIRS4)
    btw = C.between_targets()
    grid = [{**STAGE2, **p} for p in C.grid_points(C.GRID_STAGE5)]
    rows = []
    with ProcessPoolExecutor() as ex:
        for params, _seed, res in ex.map(run, [(p, 1) for p in grid]):
            rows.append({**params, **score(res, cal, val, btw)})
            print({k: params[k] for k in ("C_joint", "k_inf", "sen_infl_gain")}, round(rows[-1]["loss_cal"], 1), flush=True)
    g = pd.DataFrame(rows).sort_values("loss_cal")
    g.to_csv(DERIVED / "stage5_grid.csv", index=False)

    cand = [{k: float(r[k]) for k in KEYS} for _, r in g.head(5).iterrows()] + [STAGE2, STAGE1]
    rows, tables = [], []
    with ProcessPoolExecutor() as ex:
        for params, seed, res in ex.map(run, [(p, s) for p in cand for s in (1, 2, 3)]):
            rows.append({**params, "seed": seed, **score(res, cal, val, btw)})
            tables.append(res[0].reset_index().assign(seed=seed, **params))
    sel = pd.DataFrame(rows)
    sel.to_csv(DERIVED / "stage5_selection.csv", index=False)
    pd.concat(tables, ignore_index=True).to_csv(DERIVED / "stage5_selection_tables.csv", index=False)
    summ = sel.groupby(KEYS)[["loss_cal", "loss_val", "loss_between"]].agg(["mean", "std"]).sort_values(("loss_cal", "mean"))
    print(summ.round(1).to_string())
    best = dict(zip(KEYS, map(float, summ.index[0])))

    frozen = {"stage5_parameters": best,
              "selected_by": "lowest mean calibration loss (GSE104782 OA1-OA5 SenMayo cell loss + GSE114007 between-joint loss, seeds 1-3)",
              "calibration_loss_mean": float(summ.iloc[0][("loss_cal", "mean")]),
              "validation_loss_mean": float(summ.iloc[0][("loss_val", "mean")])}
    fp = ROOT / "docs" / "stage5_frozen_parameters.json"
    fp.write_text(json.dumps(frozen, indent=2) + "\n", encoding="utf-8")
    h = hashlib.sha256(fp.read_bytes()).hexdigest()
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(ROOT / "docs" / "preregistration_stage5.md", "a", encoding="utf-8") as fh:
        fh.write(f"\n## Frozen parameters\n\n* Frozen at **{stamp}** (UTC), before any expression value of GSE55235 / GSE55457 "
                 f"was read.\n* SHA-256 of `docs/stage5_frozen_parameters.json`: `{h}`\n")
    print("frozen:", best, "sha256", h, f"({time.perf_counter() - t0:.0f} s)")


if __name__ == "__main__":
    main()
