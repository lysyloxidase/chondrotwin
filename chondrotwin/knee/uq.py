"""Stage 6: identifiability, uncertainty and robustness of the selected (stage-2) model. No new data are used.

* Profile of the calibration loss over the two stage-2 parameters (c_over_max, k_cat), averaged over seeds, with an
  approximate 95 % confidence region.  The loss is a sum of squared standardized errors, but the model is misspecified
  (reduced chi-square well above 1), so the chi-square threshold is scaled by L_min / dof (quasi-likelihood, as for
  over-dispersion): region = { L - L_min <= chi2_{2, 0.95} * L_min / (n_terms - 2) }.
* Multiverse: the stage-1 vs stage-2 comparison on the validation patients is repeated under alternative, equally
  defensible analysis choices (simulation settings, stage bins, pseudo-count, senescence readout).
"""
from __future__ import annotations

from itertools import product
from typing import Dict

import numpy as np
import pandas as pd

from . import calibrate as C

STAGE1 = dict(c_over_max=0.9, k_cat=0.0)
STAGE2 = dict(c_over_max=0.15, k_cat=0.5)
PROFILE_GRID = {"c_over_max": (0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5),
                "k_cat": (0.0, 0.25, 0.4, 0.5, 0.6, 0.75, 1.0)}
BINS = {"0.2/0.4/0.6/0.8 (default)": (-1, 0.2, 0.4, 0.6, 0.8, 2),
        "0.3/0.5/0.7/0.9": (-1, 0.3, 0.5, 0.7, 0.9, 2),
        "0.1/0.3/0.5/0.7": (-1, 0.1, 0.3, 0.5, 0.7, 2)}
SIM_VARIANTS = {"default (0.5 mm slab, 2 y, 0.07-0.16 MPa)": dict(slab_mm=(0.5, 0.5, 2.2), years=2, sigmas=C.SIGMAS),
                "3 simulated years": dict(slab_mm=(0.5, 0.5, 2.2), years=3, sigmas=C.SIGMAS),
                "0.7 mm slab": dict(slab_mm=(0.7, 0.7, 2.2), years=2, sigmas=C.SIGMAS),
                "stresses 0.06-0.12 MPa": dict(slab_mm=(0.5, 0.5, 2.2), years=2,
                                               sigmas=(0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12))}
READOUTS = {"SenMayo": C.PAIRS4, "CDKN2A": C.PAIRS}


def chi2_threshold(n_terms: int, l_min: float, n_params: int = 2, level: float = 0.95) -> float:
    from scipy import stats
    return float(stats.chi2.ppf(level, n_params) * max(l_min / max(n_terms - n_params, 1), 1.0))


def cells_for(params: Dict[str, float], seed: int, slab_mm, years: int, sigmas) -> pd.DataFrame:
    from ..abm import Rates
    from .analysis import abm_stage_composition
    _, cells = abm_stage_composition(sigmas, years=years, seed=seed, rates=Rates(**params), slab_mm=slab_mm,
                                     return_cells=True)
    return cells


def stage_table(cells: pd.DataFrame, bins) -> pd.DataFrame:
    d = cells.copy()
    d["stage"] = pd.cut(d["M"], list(bins), labels=[4, 3, 2, 1, 0]).astype(int)
    g = d.groupby("stage")
    out = g[["runx2_high", "senescent", "clustered", "catabolic"]].mean()
    out["mean_S"], out["mean_R"], out["mean_cat"] = g["S"].mean(), g["R"].mean(), g["cat"].mean()
    out["n"] = g.size()
    return out.reindex(range(5))


def log_ratios(stats: pd.DataFrame, eps: float, fraction_cols=C.FRACTIONS) -> pd.DataFrame:
    out = {}
    for c in stats.columns:
        v = stats[c].astype(float) + (eps if c in fraction_cols else 0.0)
        out[c] = np.log(v.loc[[1, 2, 3, 4]].to_numpy() / v.loc[0])
    return pd.DataFrame(out, index=[1, 2, 3, 4])


def model_lr(st: pd.DataFrame, pairs, eps: float) -> pd.DataFrame:
    return log_ratios(pd.DataFrame({d: st[m].astype(float) for d, m in pairs}), eps)


def data_targets(t: pd.DataFrame, patients, pairs, eps: float, n_boot: int = 500, seed: int = 0):
    cols = tuple(d for d, _ in pairs)
    pats = np.array(list(patients))
    sub = t[t.patient.isin(pats)]
    stats = lambda x: x.groupby("stage")[list(cols)].mean().reindex(range(5))   # noqa: E731
    est = log_ratios(stats(sub), eps)
    rng = np.random.default_rng(seed)
    b = np.array([log_ratios(stats(pd.concat([sub[sub.patient == p] for p in rng.choice(pats, len(pats))])), eps).to_numpy()
                  for _ in range(n_boot)])
    var = pd.DataFrame(np.nanvar(np.where(np.isfinite(b), b, np.nan), axis=0), index=est.index, columns=est.columns)
    return {"log_ratio": est, "var": var.clip(lower=0.02)}


def multiverse(cells_by: Dict[tuple, pd.DataFrame], ji: pd.DataFrame) -> pd.DataFrame:
    """cells_by[(sim_variant, model, seed)] -> per-cell table; returns one row per universe."""
    rows = []
    targets = {}
    for (rname, pairs), eps in product(READOUTS.items(), (0.01, 0.05)):
        targets[(rname, eps)] = (data_targets(ji, C.CAL_PATIENTS, pairs, eps), data_targets(ji, C.VAL_PATIENTS, pairs, eps))
    sims = sorted({k[0] for k in cells_by})
    seeds = sorted({k[2] for k in cells_by})
    for sim, (bname, bins), (rname, pairs), eps in product(sims, BINS.items(), READOUTS.items(), (0.01, 0.05)):
        cal, val = targets[(rname, eps)]
        res = {}
        for model in ("stage 1", "stage 2"):
            lv, lc = [], []
            for s in seeds:
                lr = model_lr(stage_table(cells_by[(sim, model, s)], bins), pairs, eps)
                lv.append(C.loss(lr, val))
                lc.append(C.loss(lr, cal))
            res[model] = (float(np.mean(lc)), float(np.mean(lv)))
        rows.append({"simulation": sim, "stage bins": bname, "senescence readout": rname, "pseudo-count": eps,
                     "stage 1 cal": res["stage 1"][0], "stage 2 cal": res["stage 2"][0],
                     "stage 1 val": res["stage 1"][1], "stage 2 val": res["stage 2"][1],
                     "stage 2 better (val)": res["stage 2"][1] < res["stage 1"][1],
                     "val loss ratio stage1/stage2": res["stage 1"][1] / res["stage 2"][1]})
    return pd.DataFrame(rows)
