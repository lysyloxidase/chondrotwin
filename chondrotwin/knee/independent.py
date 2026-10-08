"""Independent out-of-sample test on GSE152805 (Chou et al., Sci Rep 2020).

PRE-REGISTERED PROTOCOL — written and saved before any of the GSE152805 files were opened (2026-10-08,
chondrotwin 0.5.0 development; the stage-2 model was calibrated without these data):

Data      Six 10x Genomics cartilage samples from three OA patients (113, 116, 118): minimally damaged outer
          lateral tibial plateau (oLT) and damaged medial tibial plateau (MT) from the same knee.
QC        Cells with >= 500 detected genes and <= 20 % mitochondrial (MT-) UMIs; no other filtering.
Statistic Per patient, MT / oLT ratio of
            * the fraction of cells with >= 1 UMI of CDKN2A (p16), RUNX2 and MMP13, and
            * the mean SOX9 expression (UMIs per 10,000).
Model     The same ratios for the ABM's damaged tissue (local voxel matrix < 0.4, stages 3-4) versus intact
          tissue (>= 0.8, stage 0): senescent fraction, RUNX2-high fraction, fraction of catabolically active
          cells (R > S or overload-driven catabolism > 0.5) and mean S.
Criterion For each marker the model ratio "agrees" if it lies within [min, max] of the three per-patient
          ratios.  Reported for the stage-1 (0.4.0) and the stage-2 parameterization; no parameter is changed
          after the data are opened.
"""
from __future__ import annotations

import gzip
from typing import Dict

import numpy as np
import pandas as pd

from .. import datasets

SAMPLES = {("113", "oLT"): "GSM4626766_OA_oLT_113", ("116", "oLT"): "GSM4626767_OA_oLT_116",
           ("118", "oLT"): "GSM4626768_OA_oLT_118", ("113", "MT"): "GSM4626769_OA_MT_113",
           ("116", "MT"): "GSM4626770_OA_MT_116", ("118", "MT"): "GSM4626771_OA_MT_118"}
MARKERS_DETECTION = ("CDKN2A", "RUNX2", "MMP13")
MIN_GENES, MAX_MITO = 500, 0.20


def _read_10x(prefix: str, genes_of_interest, root=None):
    from scipy.io import mmread
    with gzip.open(datasets.require("gse152805", f"{prefix}.genes.tsv.gz", root), "rt") as fh:
        genes = pd.read_csv(fh, sep="\t", header=None)
    with gzip.open(datasets.require("gse152805", f"{prefix}.matrix.mtx.gz", root), "rb") as fh:
        m = mmread(fh).tocsc()                       # genes x cells
    symbols = genes.iloc[:, 1].astype(str).to_numpy() if genes.shape[1] > 1 else genes.iloc[:, 0].astype(str).to_numpy()
    n_genes = np.asarray((m > 0).sum(axis=0)).ravel()
    total = np.asarray(m.sum(axis=0)).ravel()
    mito = np.asarray(m[np.flatnonzero(np.char.startswith(symbols.astype(str), "MT-")), :].sum(axis=0)).ravel()
    keep = (n_genes >= MIN_GENES) & (mito <= MAX_MITO * np.maximum(total, 1))
    out = {"n_cells_raw": int(m.shape[1]), "n_cells_qc": int(keep.sum())}
    for g in genes_of_interest:
        idx = np.flatnonzero(symbols == g)
        row = np.asarray(m[idx, :].sum(axis=0)).ravel()[keep] if idx.size else np.zeros(int(keep.sum()))
        out[g] = row
    out["total"] = total[keep]
    return out


def gse152805_ratios(root=None) -> pd.DataFrame:
    """Per-patient MT / oLT ratios as pre-registered above."""
    rows = []
    per = {}
    for (pat, site), prefix in SAMPLES.items():
        d = _read_10x(prefix, list(MARKERS_DETECTION) + ["SOX9"], root)
        rec = {"patient": pat, "site": site, "n_cells_raw": d["n_cells_raw"], "n_cells_qc": d["n_cells_qc"]}
        for g in MARKERS_DETECTION:
            rec[f"{g}_detected"] = float(np.mean(d[g] > 0))
        rec["SOX9_per10k"] = float(np.mean(d["SOX9"] / np.maximum(d["total"], 1) * 1e4))
        per[(pat, site)] = rec
        rows.append(rec)
    table = pd.DataFrame(rows)
    ratios = []
    for pat in ("113", "116", "118"):
        mt, olt = per[(pat, "MT")], per[(pat, "oLT")]
        r = {"patient": pat}
        for k in [f"{g}_detected" for g in MARKERS_DETECTION] + ["SOX9_per10k"]:
            r[k] = mt[k] / olt[k] if olt[k] > 0 else np.inf
        ratios.append(r)
    out = pd.DataFrame(ratios).set_index("patient")
    out.attrs["samples"] = table
    return out


def model_ratios(stage_table: pd.DataFrame) -> Dict[str, float]:
    """Model counterparts from :func:`chondrotwin.knee.analysis.abm_stage_composition` output (damaged = stages 3-4)."""
    w = stage_table["n"]
    dam = stage_table.loc[[3, 4]]
    wd = w.loc[[3, 4]]

    def mean(col, rows, weights):
        return float((rows[col] * weights).sum() / weights.sum())

    intact = stage_table.loc[[0]]
    wi = w.loc[[0]]
    res = {}
    for data_key, col in (("CDKN2A_detected", "senescent"), ("RUNX2_detected", "runx2_high"),
                          ("MMP13_detected", "catabolic"), ("SOX9_per10k", "mean_S")):
        if col not in stage_table:
            res[data_key] = float("nan")
            continue
        a, b = mean(col, dam, wd), mean(col, intact, wi)
        res[data_key] = a / b if b > 0 else float("inf")
    return res


def agreement(data: pd.DataFrame, model: Dict[str, float]) -> pd.DataFrame:
    rows = []
    for k, v in model.items():
        lo, hi = float(data[k].min()), float(data[k].max())
        rows.append(dict(marker=k, data_min=lo, data_max=hi, data_geomean=float(np.exp(np.log(data[k]).mean())),
                         model=v, agrees=bool(lo <= v <= hi)))
    return pd.DataFrame(rows).set_index("marker")
