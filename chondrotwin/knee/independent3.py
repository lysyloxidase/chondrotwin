"""Third independent test (stage 4): OA versus non-OA knee cartilage, GSE117999 (primary) and GSE169077 (secondary),
with a multi-gene senescence readout (SenMayo) in addition to the four markers.

PRE-REGISTERED PROTOCOL — written and saved before any expression value of GSE117999 or GSE169077 was read.  Only
structure was inspected (SenMayo gene list, column names and identifier columns, sample metadata, platform
annotation rows of the four markers); see docs/preregistration_stage4.md.  Stage-4 parameters are calibrated on
GSE104782 patients OA1-5 only and frozen with a SHA-256 before the expression values are opened.

Datasets
    GSE117999 (Rai et al. 2019; Agilent): the authors' normalized matrix (GSE117999_normalized_matrix.xlsx), which
    already applies their exclusions (P4-010, -012, -104, -108).  Columns "P4-0xx" = arthroscopic partial
    meniscectomy without OA (non-OA), "P4-1xx" = total knee arthroplasty for OA.  Rows with ControlType != 0 are
    dropped; genes by the GeneName column.
    GSE169077 (Attur, Abramson; Affymetrix U133A, GPL96): series matrix values (dChip); samples titled "N*" = normal
    pools, "OA*" = late-stage OA pools; genes by GPL96.annot.
Scale     If the largest value of a dataset exceeds 50, values are log2(max(x, 1))-transformed; otherwise used as given.
Probes    Per gene, the probe (row) with the highest mean over all samples.
Readouts  (i) SOX9, RUNX2, MMP13, CDKN2A: difference of group means of log2 expression (OA - non-OA) with a Welch
          95 % interval.  (ii) SenMayo score: for every SenMayo gene present, log2 expression z-scored across the
          dataset's samples; per-sample score = mean z; group difference with a Welch 95 % interval.
Model     Damaged (stages 3-4) vs intact (stage 0) tissue, log2 ratio: SOX9 <-> mean S, RUNX2 <-> mean R, MMP13 <->
          mean catabolic drive, CDKN2A and SenMayo <-> senescent fraction (SenMayo is dominated by SASP genes,
          produced by senescent cells).
Criterion C1 for every readout: if the interval excludes 0 the model must have the same sign; otherwise the model
          must predict |log2 ratio| < 1.  No magnitude criterion (array fold changes are compressed; the SenMayo
          score is in z units).
Reported for the stage-1, stage-2 (= 3) and stage-4 models.
"""
from __future__ import annotations

import gzip
import json
from typing import Dict, List

import numpy as np
import pandas as pd

from .. import datasets

MARKERS = ("SOX9", "RUNX2", "MMP13", "CDKN2A")


def senmayo_genes(root=None) -> List[str]:
    g = json.load(open(datasets.require("senmayo", "SAUL_SEN_MAYO.json", root), encoding="utf-8"))
    return list(g[next(iter(g))]["geneSymbols"])


def _scale(df: pd.DataFrame) -> pd.DataFrame:
    return np.log2(df.clip(lower=1.0)) if float(np.nanmax(df.to_numpy())) > 50 else df


def _collapse(expr: pd.DataFrame, gene_of_row: pd.Series) -> pd.DataFrame:
    """One row per gene: the row with the highest mean over all samples."""
    m = expr.mean(axis=1)
    d = pd.DataFrame({"gene": gene_of_row, "mean": m}).dropna(subset=["gene"])
    best = d.sort_values("mean", ascending=False).groupby("gene").head(1)
    out = expr.loc[best.index]
    out.index = best["gene"].to_numpy()
    return out


def gse117999(root=None):
    """(log2 expression genes x samples, group labels) for GSE117999."""
    x = pd.read_excel(datasets.require("gse117999", "GSE117999_normalized_matrix.xlsx", root))
    x = x[x["ControlType"] == 0]
    samples = [c for c in x.columns if str(c).startswith("P4-")]
    expr = _scale(x[samples].astype(float))
    expr = _collapse(expr, x["GeneName"].astype(str))
    groups = pd.Series(["OA" if s[3] == "1" else "nonOA" for s in samples], index=samples)
    return expr, groups


def gse169077(root=None):
    titles, header, rows = None, None, {}
    with gzip.open(datasets.require("gse169077", "GSE169077_series_matrix.txt.gz", root), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                header = [h.strip('"') for h in next(fh).rstrip("\n").split("\t")]
                for row in fh:
                    if row.startswith("!series_matrix_table_end"):
                        break
                    v = row.rstrip("\n").split("\t")
                    rows[v[0].strip('"')] = [float(t) if t not in ("", "null") else np.nan for t in v[1:]]
                break
    expr = pd.DataFrame(rows, index=titles).T
    ann = {}
    with gzip.open(datasets.require("gse169077", "GPL96.annot.gz", root), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                next(fh)
                break
        for line in fh:
            v = line.rstrip("\n").split("\t")
            if len(v) > 2 and v[2] and "///" not in v[2]:
                ann[v[0]] = v[2]
    expr = _scale(expr)
    expr = _collapse(expr, pd.Series({p: ann.get(p) for p in expr.index}))
    groups = pd.Series(["OA" if t.upper().startswith("OA") else "nonOA" for t in titles], index=titles)
    return expr, groups


def welch(a: np.ndarray, b: np.ndarray, level: float = 0.95):
    from scipy import stats
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    d = a.mean() - b.mean()
    se = np.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    h = stats.t.ppf(0.5 + level / 2, df) * se
    return float(d), float(d - h), float(d + h)


def readouts(expr: pd.DataFrame, groups: pd.Series, senmayo: List[str]) -> pd.DataFrame:
    oa, non = groups[groups == "OA"].index, groups[groups == "nonOA"].index
    rows = {}
    for g in MARKERS:
        if g in expr.index:
            d, lo, hi = welch(expr.loc[g, oa], expr.loc[g, non])
            rows[g] = dict(diff=d, ci_low=lo, ci_high=hi, units="log2")
    present = [g for g in senmayo if g in expr.index]
    z = expr.loc[present]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    score = z.mean(axis=0)
    d, lo, hi = welch(score[oa], score[non])
    rows["SenMayo"] = dict(diff=d, ci_low=lo, ci_high=hi, units="z")
    out = pd.DataFrame(rows).T
    out.attrs.update(n_oa=len(oa), n_non=len(non), senmayo_present=len(present), senmayo_total=len(senmayo))
    return out


def model_log2(stage_table: pd.DataFrame) -> Dict[str, float]:
    from .independent2 import model_log2 as m2
    m = m2(stage_table)
    m["SenMayo"] = m["CDKN2A"]
    return m


def criterion_c1(diff: float, lo: float, hi: float, model: float) -> bool:
    if lo > 0 or hi < 0:
        return bool(np.sign(model) == np.sign(diff))
    return bool(abs(model) < 1.0)
