"""Stage-2 calibration of the cell-state rules on single-cell data (GSE104782), with a patient split fixed in advance.

Calibration patients OA1-OA5, validation patients OA6-OA10.  Targets: for stages 1-4 relative to stage 0,

    data (Ji et al. 2019)                               model (ABM, live cells, local matrix stage)
    fraction of cells with >= 1 CDKN2A UMI      <->     senescent fraction
    fraction of cells with >= 1 RUNX2 UMI       <->     RUNX2-high fraction (R > S)
    fraction of cells with >= 1 MMP13 UMI       <->     catabolically active fraction
    mean SOX9 UMIs per 10,000                   <->     mean SOX9 activity S

The same correspondence is used in the independent test (:mod:`chondrotwin.knee.independent`).  Fractions get a
pseudo-count of 0.01 before taking ratios.  The loss is the sum of squared log-ratio errors weighted by the inverse
patient-bootstrap variance of the data log-ratios.
"""
from __future__ import annotations

from itertools import product
from typing import Dict, Iterable, Sequence

import numpy as np
import pandas as pd

from .. import datasets

CAL_PATIENTS = tuple(f"OA{i}" for i in range(1, 6))
VAL_PATIENTS = tuple(f"OA{i}" for i in range(6, 11))
PAIRS = (("CDKN2A", "senescent"), ("RUNX2", "runx2_high"), ("MMP13", "catabolic"), ("SOX9", "mean_S"))
# stage 4: senescence read out with the SenMayo score (fraction of senescence-high cells) instead of CDKN2A
PAIRS4 = (("SENMAYO", "senescent"), ("RUNX2", "runx2_high"), ("MMP13", "catabolic"), ("SOX9", "mean_S"))
FRACTIONS = ("CDKN2A", "RUNX2", "MMP13", "SENMAYO")
EPS = 0.01


def ji_cell_table(root=None) -> pd.DataFrame:
    """Per-cell detection of CDKN2A / RUNX2 / MMP13, SOX9 per 10k and the SenMayo senescence-high flag, with patient and stage."""
    from .cells import ji2019_cells
    cells = ji2019_cells(root).set_index("Cell")
    umi = pd.read_csv(datasets.require("gse104782", "GSE104782_allcells_UMI_count.txt.gz", root), sep="\t", index_col=0)
    tot = umi.sum(axis=0)
    t = pd.DataFrame({g: (umi.loc[g] > 0).astype(float) for g in ("CDKN2A", "RUNX2", "MMP13")})
    t["SOX9"] = umi.loc["SOX9"] / tot * 1e4
    from .senescence import ji_scores
    t["SENMAYO"] = ji_scores(root)["sen_high"]
    return t.join(cells[["patient", "stage"]], how="inner")


def _stage_stats(t: pd.DataFrame, cols=("CDKN2A", "RUNX2", "MMP13", "SOX9")) -> pd.DataFrame:
    s = t.groupby("stage")[list(cols)].mean()
    return s.reindex(range(5))


def _log_ratios(stats: pd.DataFrame, fraction_cols: Sequence[str]) -> pd.DataFrame:
    out = {}
    for c in stats.columns:
        v = stats[c].astype(float)
        if c in fraction_cols or c in FRACTIONS:
            v = v + EPS
        out[c] = np.log(v.loc[[1, 2, 3, 4]].to_numpy() / v.loc[0])
    return pd.DataFrame(out, index=[1, 2, 3, 4])


def data_targets(t: pd.DataFrame, patients: Iterable[str], n_boot: int = 1000, seed: int = 0,
                 pairs=PAIRS) -> Dict[str, pd.DataFrame]:
    pats = np.array(list(patients))
    sub = t[t.patient.isin(pats)]
    cols = tuple(d for d, _ in pairs)
    est = _log_ratios(_stage_stats(sub, cols), ("CDKN2A", "RUNX2", "MMP13"))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        pick = rng.choice(pats, len(pats))
        bt = pd.concat([sub[sub.patient == p] for p in pick])
        boots.append(_log_ratios(_stage_stats(bt, cols), ("CDKN2A", "RUNX2", "MMP13")).to_numpy())
    b = np.array(boots)
    var = pd.DataFrame(np.nanvar(np.where(np.isfinite(b), b, np.nan), axis=0), index=est.index, columns=est.columns)
    return {"log_ratio": est, "var": var.clip(lower=0.02)}


def model_log_ratios(stage_table: pd.DataFrame, pairs=PAIRS) -> pd.DataFrame:
    """Model log-ratios (stage s / stage 0) in the data's column names; NaN where a stage has no cells."""
    m = pd.DataFrame({data: stage_table[model].astype(float) for data, model in pairs})
    return _log_ratios(m, ("CDKN2A", "RUNX2", "MMP13"))


def loss(model_lr: pd.DataFrame, targets: Dict[str, pd.DataFrame], missing_penalty: float = 25.0) -> float:
    err = (model_lr - targets["log_ratio"]) ** 2 / targets["var"]
    n_missing = int(err.isna().to_numpy().sum())
    return float(np.nansum(err.to_numpy()) + missing_penalty * n_missing)


GRID = {"damage_D50": (0.0, 30.0, 100.0, 300.0), "c_over_max": (0.9, 0.6, 0.45, 0.3), "k_cat": (0.0, 0.25, 0.5, 1.0)}
# refinement after the first grid put c_over_max on the boundary and left a baseline RUNX2+ population unexplained
GRID_REFINE = {"damage_D50": (0.0, 1000.0), "c_over_max": (0.0, 0.15, 0.3), "k_cat": (0.35, 0.5, 0.7),
               "baseline_runx2": (0.0, 0.05, 0.10)}
# stage 4: senescence gated by inflammation and death of senescent cells (stage-2 switch parameters kept)
GRID_STAGE4 = {"sen_overload_weight": (1.0, 0.3, 0.1, 0.0), "senescent_death_rate": (0.0, 0.002, 0.01, 0.03),
               "c_over_max": (0.15, 0.3), "k_cat": (0.5,)}
# stage 3: cell-to-cell heterogeneity of the switch and of senescence competence (stage-2 mechanisms kept)
GRID_STAGE3 = {"switch_bias_sd": (0.0, 0.3, 0.45, 0.6), "senescence_competent_fraction": (1.0, 0.5, 0.25, 0.1),
               "c_over_max": (0.15, 0.3), "k_cat": (0.5, 0.7)}
SIGMAS = (0.07, 0.08, 0.09, 0.10, 0.12, 0.14, 0.16)


def evaluate(params: Dict[str, float], sigmas=SIGMAS, years: int = 2, seed: int = 1) -> pd.DataFrame:
    """ABM stage table for one parameter set (``baseline_runx2`` = initial RUNX2-high fraction, the rest are Rates)."""
    from ..abm import Rates
    from .analysis import abm_stage_composition
    p = dict(params)
    base = p.pop("baseline_runx2", 0.0)
    return abm_stage_composition(sigmas, years=years, seed=seed, rates=Rates(**p), slab_mm=(0.5, 0.5, 2.2),
                                 initial_hypertrophic_fraction=base)


def grid_points(grid=GRID):
    keys = list(grid)
    for vals in product(*(grid[k] for k in keys)):
        yield dict(zip(keys, vals))


# ---------------------------------------------------------------- stage 5: joint-level inflammation
GRID_STAGE5 = {"C_joint": (0.0, 0.25, 0.5, 1.0), "k_inf": (0.0, 0.5, 1.0, 2.0), "sen_infl_gain": (2.0, 0.0)}
BETWEEN = (("SOX9", "mean_S"), ("RUNX2", "mean_R"), ("MMP13", "mean_cat"), ("CDKN2A", "senescent"))


def between_targets(root=None, n_boot: int = 2000, seed: int = 0) -> pd.DataFrame:
    """OA / normal knee cartilage, GSE114007 raw counts: ln ratio of mean CPM with bootstrap variance."""
    x = pd.ExcelFile(datasets.require("gse114007", "GSE114007_raw_counts.xlsx", root))
    nrm = x.parse("Normal").set_index("symbol").groupby(level=0).sum()
    oa = x.parse("OA").set_index("symbol").groupby(level=0).sum()
    common = nrm.index.intersection(oa.index)
    cn, co = nrm.loc[common] / nrm.loc[common].sum() * 1e6, oa.loc[common] / oa.loc[common].sum() * 1e6
    rng = np.random.default_rng(seed)
    rows = {}
    for g, _ in BETWEEN:
        a, b = co.loc[g].to_numpy(float), cn.loc[g].to_numpy(float)
        bt = [np.log(rng.choice(a, len(a)).mean() / rng.choice(b, len(b)).mean()) for _ in range(n_boot)]
        rows[g] = dict(log_ratio=float(np.log(a.mean() / b.mean())), var=max(float(np.var(bt)), 0.02))
    return pd.DataFrame(rows).T


def evaluate_joint(params: Dict[str, float], seed: int = 1, sigmas=SIGMAS, years: int = 2):
    """(OA-joint stage table, OA-joint bulk means, normal-joint bulk means) for one parameter set.

    OA joint: the overload runs of :func:`evaluate` with the joint-level cytokine background ``C_joint``.
    Normal joint: physiological stress (0.05 MPa) without joint-level inflammation.
    """
    from ..abm import Rates
    from .analysis import abm_stage_composition, bulk_means
    p = dict(params)
    base = p.pop("baseline_runx2", 0.0)
    cj = p.pop("C_joint", 0.0)
    rates = Rates(**p)
    st, cells = abm_stage_composition(sigmas, years=years, seed=seed, rates=rates, slab_mm=(0.5, 0.5, 2.2),
                                      initial_hypertrophic_fraction=base, C_joint=cj, return_cells=True)
    _, normal = abm_stage_composition((0.05,), years=years, seed=seed, rates=rates, slab_mm=(0.5, 0.5, 2.2),
                                      initial_hypertrophic_fraction=base, C_joint=0.0, return_cells=True)
    return st, bulk_means(cells), bulk_means(normal)


def between_loss(oa: Dict[str, float], normal: Dict[str, float], targets: pd.DataFrame) -> float:
    err = 0.0
    for g, col in BETWEEN:
        m = np.log(max(oa[col], 1e-6) / max(normal[col], 1e-6))
        err += (m - targets.loc[g, "log_ratio"]) ** 2 / targets.loc[g, "var"]
    return float(err)
