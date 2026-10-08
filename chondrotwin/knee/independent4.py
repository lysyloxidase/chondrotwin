"""Stage-5 premise test: is the synovium of OA joints inflamed? GSE55235 + GSE55457 (Woetzel et al. 2014).

PRE-REGISTERED PROTOCOL — written and saved before any expression value of GSE55235 / GSE55457 was read (only the
sample metadata lines, table headers and the hallmark gene list were inspected; GPL96 annotation was already in use).
The stage-5 model is calibrated without these data and frozen before they are opened.

The stage-5 model assumes a joint-level inflammatory input (C_joint > 0) in OA joints and none in normal joints.
This test does not check the model's numbers; it checks that premise in independent tissue.

Data      Series matrix values (GSE55235: GCOS; GSE55457: MAS5), Affymetrix HG-U133A (GPL96).  Groups by the
          "disease state" / "clinical status" characteristic: healthy / normal control, osteoarthritis, rheumatoid
          arthritis.  Values > 50 anywhere -> log2(max(x, 1)).  Per gene the probe with the highest mean (GPL96.annot;
          probes mapping to several genes are dropped).
Readouts  (i) Inflammation score = mean over HALLMARK_INFLAMMATORY_RESPONSE genes of the per-gene z-score (z across all
          samples of the dataset); (ii) IL1B, TNF, IL6, CCL2: log2 differences.
Contrast  OA - healthy in each dataset (Welch 95 % interval); the two datasets combined by inverse-variance weighting
          (normal approximation).  RA - healthy is computed the same way as a positive control of the assay.
Criterion The premise is supported if the combined OA - healthy inflammation-score interval lies above 0.
          The positive control is valid if the combined RA - healthy interval lies above 0.
"""
from __future__ import annotations

import gzip
import json
from typing import Dict, List

import numpy as np
import pandas as pd

from .. import datasets
from .independent3 import welch

CYTOKINES = ("IL1B", "TNF", "IL6", "CCL2")
FILES = ("GSE55235_series_matrix.txt.gz", "GSE55457_series_matrix.txt.gz")


def hallmark_genes(root=None) -> List[str]:
    g = json.load(open(datasets.require("hallmark_inflammation", "HALLMARK_INFLAMMATORY_RESPONSE.json", root), encoding="utf-8"))
    return list(g[next(iter(g))]["geneSymbols"])


def _gpl96(root=None) -> Dict[str, str]:
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
    return ann


def _group(label: str) -> str:
    s = label.lower()
    if "rheumatoid" in s:
        return "RA"
    if "osteoarthri" in s:
        return "OA"
    if "healthy" in s or "normal" in s:
        return "healthy"
    return "other"


def read_dataset(name: str, root=None):
    groups, rows = None, {}
    with gzip.open(datasets.require("gse55235", name, root), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_characteristics_ch1") and ("disease state" in line or "clinical status" in line):
                labs = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
                if groups is None or all(g == "other" for g in groups):
                    groups = [_group(t) for t in labs]
            elif line.startswith("!series_matrix_table_begin"):
                next(fh)
                for row in fh:
                    if row.startswith("!series_matrix_table_end"):
                        break
                    v = row.rstrip("\n").split("\t")
                    rows[v[0].strip('"')] = [float(t) if t not in ("", "null") else np.nan for t in v[1:]]
                break
    expr = pd.DataFrame(rows).T
    expr.columns = range(expr.shape[1])
    if float(np.nanmax(expr.to_numpy())) > 50:
        expr = np.log2(expr.clip(lower=1.0))
    ann = _gpl96(root)
    m = expr.mean(axis=1)
    d = pd.DataFrame({"gene": [ann.get(p) for p in expr.index], "mean": m}, index=expr.index).dropna(subset=["gene"])
    best = d.sort_values("mean", ascending=False).groupby("gene").head(1)
    e = expr.loc[best.index]
    e.index = best["gene"].to_numpy()
    return e, pd.Series(groups)


def contrasts(root=None) -> Dict[str, pd.DataFrame]:
    genes = hallmark_genes(root)
    per = {}
    for f in FILES:
        e, g = read_dataset(f, root)
        present = [x for x in genes if x in e.index]
        z = e.loc[present]
        z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
        score = z.mean(axis=0)
        rows = {}
        for case in ("OA", "RA"):
            a, b = g.index[g == case], g.index[g == "healthy"]
            rows[(case, "inflammation score")] = welch(score[a], score[b]) + (len(a), len(b))
            for c in CYTOKINES:
                if c in e.index:
                    rows[(case, c)] = welch(e.loc[c, a], e.loc[c, b]) + (len(a), len(b))
        per[f.split("_")[0]] = pd.DataFrame(rows, index=["diff", "ci_low", "ci_high", "n_case", "n_healthy"]).T
        per[f.split("_")[0]].attrs["hallmark_present"] = len(present)
    return per


def combine(per: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    keys = next(iter(per.values())).index
    out = {}
    for k in keys:
        ds = [p.loc[k] for p in per.values() if k in p.index]
        se = [(r["ci_high"] - r["ci_low"]) / (2 * 1.96) for r in ds]
        w = np.array([1 / s ** 2 for s in se])
        est = float(np.sum(w * np.array([r["diff"] for r in ds])) / w.sum())
        s = float(np.sqrt(1 / w.sum()))
        out[k] = dict(diff=est, ci_low=est - 1.96 * s, ci_high=est + 1.96 * s)
    return pd.DataFrame(out).T
