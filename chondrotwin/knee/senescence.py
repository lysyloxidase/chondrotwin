"""Multi-gene senescence readout (SenMayo; Saul et al., Nat Commun 2022) for the single-cell and bulk datasets.

Single cell: per gene log1p(UMI per 10,000) is z-scored across the cells of a dataset; a cell's score is the mean z
over the SenMayo genes present.  A cell is "senescence-high" if its score exceeds the 90th percentile of the
reference cells (stage-0 cells of GSE104782; the minimally damaged lateral sample of the same patient in GSE152805),
so the reference fraction is 0.10 by construction and stage ratios are comparable with the model's senescent
fraction ratio.
Bulk (GSE57218): per-sample mean z over SenMayo genes (probe with the highest mean per gene), paired OA - preserved.
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd

from .. import datasets
from .independent3 import senmayo_genes

THRESHOLD_Q = 0.90


def _zscore_rows(m: pd.DataFrame) -> pd.DataFrame:
    sd = m.std(axis=1, ddof=1).replace(0, np.nan)
    return m.sub(m.mean(axis=1), axis=0).div(sd, axis=0)


def ji_scores(root=None) -> pd.DataFrame:
    from .cells import ji2019_cells
    cells = ji2019_cells(root).set_index("Cell")
    umi = pd.read_csv(datasets.require("gse104782", "GSE104782_allcells_UMI_count.txt.gz", root), sep="\t", index_col=0)
    genes = [g for g in senmayo_genes(root) if g in umi.index]
    lv = np.log1p(umi.loc[genes] / umi.sum(axis=0) * 1e4)
    score = _zscore_rows(lv).mean(axis=0)
    d = pd.DataFrame({"score": score}).join(cells[["patient", "stage"]], how="inner")
    thr = d.loc[d.stage == 0, "score"].quantile(THRESHOLD_Q)
    d["sen_high"] = (d["score"] > thr).astype(float)
    d.attrs["n_genes"] = len(genes)
    return d


def stage_ratios(d: pd.DataFrame, patients=None, n_boot: int = 2000, seed: int = 0) -> pd.DataFrame:
    sub = d if patients is None else d[d.patient.isin(list(patients))]
    pats = sub.patient.unique()
    f = sub.groupby("stage")["sen_high"].mean()
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        bt = pd.concat([sub[sub.patient == p] for p in rng.choice(pats, len(pats))])
        fb = bt.groupby("stage")["sen_high"].mean()
        boots.append([(fb.get(s, np.nan) + 0.01) / (fb.get(0, np.nan) + 0.01) for s in (1, 2, 3, 4)])
    b = np.array(boots)
    return pd.DataFrame({"fraction": [f.get(s, np.nan) for s in (1, 2, 3, 4)],
                         "ratio_vs_stage0": [(f.get(s) + 0.01) / (f.get(0) + 0.01) for s in (1, 2, 3, 4)],
                         "ci_low": np.nanpercentile(b, 2.5, axis=0), "ci_high": np.nanpercentile(b, 97.5, axis=0)},
                        index=pd.Index([1, 2, 3, 4], name="stage"))


def gse152805_senmayo(root=None) -> pd.DataFrame:
    from scipy.io import mmread
    from .independent import SAMPLES, MIN_GENES, MAX_MITO
    genes = senmayo_genes(root)
    per = {}
    for (pat, site), prefix in SAMPLES.items():
        with gzip.open(datasets.require("gse152805", f"{prefix}.genes.tsv.gz", root), "rt") as fh:
            sym = np.asarray(pd.read_csv(fh, sep="\t", header=None).iloc[:, -1].astype(str).tolist(), dtype=str)
        with gzip.open(datasets.require("gse152805", f"{prefix}.matrix.mtx.gz", root), "rb") as fh:
            m = mmread(fh).tocsc()
        n_genes = np.asarray((m > 0).sum(axis=0)).ravel()
        total = np.asarray(m.sum(axis=0)).ravel()
        mito = np.asarray(m[np.flatnonzero(np.char.startswith(sym, "MT-")), :].sum(axis=0)).ravel()
        keep = (n_genes >= MIN_GENES) & (mito <= MAX_MITO * np.maximum(total, 1))
        idx = [int(np.flatnonzero(sym == g)[0]) for g in genes if np.any(sym == g)]
        sub = np.asarray(m[idx, :][:, keep].todense(), float)
        per[(pat, site)] = np.log1p(sub / np.maximum(total[keep], 1) * 1e4)
    rows = []
    for pat in ("113", "116", "118"):
        both = np.hstack([per[(pat, "oLT")], per[(pat, "MT")]])
        mu, sd = both.mean(axis=1, keepdims=True), both.std(axis=1, ddof=1, keepdims=True)
        sd[sd == 0] = np.nan
        z_olt = np.nanmean((per[(pat, "oLT")] - mu) / sd, axis=0)
        z_mt = np.nanmean((per[(pat, "MT")] - mu) / sd, axis=0)
        thr = np.quantile(z_olt, THRESHOLD_Q)
        f_olt, f_mt = float(np.mean(z_olt > thr)), float(np.mean(z_mt > thr))
        rows.append(dict(patient=pat, frac_oLT=f_olt, frac_MT=f_mt, ratio=f_mt / f_olt, mean_z_diff=float(z_mt.mean() - z_olt.mean())))
    return pd.DataFrame(rows).set_index("patient")


def gse57218_senmayo(root=None):
    from .independent2 import t_interval
    genes = set(senmayo_genes(root))
    probe_gene = {}
    with gzip.open(datasets.require("gse57218", "GPL6947.annot.gz", root), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                next(fh)
                break
        for line in fh:
            v = line.rstrip("\n").split("\t")
            if len(v) > 2 and v[2] in genes:
                probe_gene[v[0]] = v[2]
    titles, rows = None, {}
    with gzip.open(datasets.require("gse57218", "GSE57218_series_matrix.txt.gz", root), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                next(fh)
                for row in fh:
                    if row.startswith("!series_matrix_table_end"):
                        break
                    v = row.rstrip("\n").split("\t")
                    pid = v[0].strip('"')
                    if pid in probe_gene:
                        rows[pid] = [float(t) if t not in ("", "null") else np.nan for t in v[1:]]
                break
    expr = pd.DataFrame(rows, index=titles).T
    best = pd.DataFrame({"g": [probe_gene[p] for p in expr.index], "m": expr.mean(axis=1)}, index=expr.index)
    best = best.sort_values("m", ascending=False).groupby("g").head(1)
    e = expr.loc[best.index]
    score = _zscore_rows(e).mean(axis=0)
    diffs = []
    for t in titles:
        if t.endswith("_OA") and t.replace("_OA", "_Preserved") in titles:
            diffs.append(score[t] - score[t.replace("_OA", "_Preserved")])
    return t_interval(np.array(diffs)), len(best)
