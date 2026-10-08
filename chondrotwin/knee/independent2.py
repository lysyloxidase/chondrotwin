"""Second independent test (stage 3): bulk cartilage transcriptomes, GSE57218 (primary) and GSE114007 (secondary).

PRE-REGISTERED PROTOCOL — written and saved before any expression value of these datasets was read.  Only file
structure was inspected (sample metadata lines, table headers, gene identifiers, platform annotation); this is
disclosed in docs/preregistration_stage3.md.  The stage-3 parameters are calibrated on GSE104782 patients OA1-5
only and frozen (with a SHA-256) before the expression values are opened; no parameter is changed afterwards.

Markers and model counterparts (bulk expression = mean over the live cells of a tissue piece):

    SOX9    <->  mean SOX9 activity S
    RUNX2   <->  mean RUNX2 activity R
    MMP13   <->  mean catabolic drive R + k_cat * overload   (equals R when k_cat = 0)
    CDKN2A  <->  senescent fraction

Model contrast: damaged tissue (local voxel matrix < 0.4, stages 3-4) versus intact tissue (>= 0.8, stage 0),
cell-number weighted and pooled over the three selection seeds; log2 ratio.

Primary — GSE57218 (RAAK; Illumina HumanHT-12, GPL6947; lumi RSN + VST values, used as given on their log2-like
scale).  Pairs: samples "cartilage_RAAK_<k>_OA" and "cartilage_RAAK_<k>_Preserved" of the same joint; "Healthy"
samples are not used.  Probe per gene: the GPL6947 probe with the highest mean value over all 73 samples.
Statistic: per-pair difference OA - Preserved; mean with a two-sided 95 % t interval.
    Criterion C1 (direction): if the interval excludes 0, the model agrees when its log2 ratio has the same sign;
    if the interval contains 0, the model agrees when |model log2 ratio| < 1.  Magnitudes are not compared
    because microarray fold changes are compressed.

Secondary — GSE114007 (bulk RNA-seq; separate files of normalized counts for normal and OA cartilage).  Each
sample is rescaled to counts per million over the genes present in both files; per gene the log2 ratio of the
group means of CPM (OA / normal), with a 95 % bootstrap interval (2000 resamples of samples within groups,
seed 0).
    C1 as above, and C2 (magnitude): the model log2 ratio lies inside the 95 % interval.

Reported for the stage-1 (0.4.0), stage-2 (0.5.0) and stage-3 models.
"""
from __future__ import annotations

import gzip
import re
from typing import Dict, List

import numpy as np
import pandas as pd

from .. import datasets

GENES = ("SOX9", "RUNX2", "MMP13", "CDKN2A")
MODEL_COLUMN = {"SOX9": "mean_S", "RUNX2": "mean_R", "MMP13": "mean_cat", "CDKN2A": "senescent"}


# ---------------------------------------------------------------- GSE57218
def _gpl6947_probes(root=None) -> Dict[str, List[str]]:
    probes: Dict[str, List[str]] = {g: [] for g in GENES}
    with gzip.open(datasets.require("gse57218", "GPL6947.annot.gz", root), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                next(fh)
                break
        for line in fh:
            v = line.rstrip("\n").split("\t")
            if len(v) > 2 and v[2] in probes:
                probes[v[2]].append(v[0])
    return probes


def gse57218_pairs(root=None) -> pd.DataFrame:
    """Per-pair OA - Preserved differences for the four genes (rows: RAAK pair id)."""
    titles, accs = None, None
    rows = {}
    wanted = {p for ps in _gpl6947_probes(root).values() for p in ps}
    with gzip.open(datasets.require("gse57218", "GSE57218_series_matrix.txt.gz", root), "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                header = [h.strip('"') for h in next(fh).rstrip("\n").split("\t")]
                for row in fh:
                    if row.startswith("!series_matrix_table_end"):
                        break
                    v = row.rstrip("\n").split("\t")
                    pid = v[0].strip('"')
                    if pid in wanted:
                        rows[pid] = np.array([float(x) if x not in ("", "null") else np.nan for x in v[1:]])
                break
    if header[1:] != accs:
        raise ValueError("column order differs from sample accessions")
    expr = pd.DataFrame(rows, index=titles).T                       # probes x samples (by title)
    probes = _gpl6947_probes(root)
    chosen = {g: max(ps, key=lambda p: np.nanmean(expr.loc[p])) for g, ps in probes.items() if ps}
    out = {}
    for t in titles:
        m = re.match(r"cartilage_RAAK_(\d+)_OA$", t)
        if m and f"cartilage_RAAK_{m.group(1)}_Preserved" in titles:
            pres = f"cartilage_RAAK_{m.group(1)}_Preserved"
            out[int(m.group(1))] = {g: float(expr.loc[p, t] - expr.loc[p, pres]) for g, p in chosen.items()}
    d = pd.DataFrame(out).T.sort_index()
    d.attrs["probes"] = chosen
    return d


def t_interval(x: np.ndarray, level: float = 0.95):
    from scipy import stats
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    m, se = x.mean(), x.std(ddof=1) / np.sqrt(len(x))
    h = stats.t.ppf(0.5 + level / 2, len(x) - 1) * se
    return float(m), float(m - h), float(m + h), len(x)


# ---------------------------------------------------------------- GSE114007
def _read_counts(name: str, root=None) -> pd.DataFrame:
    with gzip.open(datasets.require("gse114007", name, root), "rt") as fh:
        d = pd.read_csv(fh, sep="\t", index_col=0)
    return d.groupby(level=0).sum()


def gse114007_log2_ratios(root=None, n_boot: int = 2000, seed: int = 0) -> pd.DataFrame:
    normal = _read_counts("GSE114007_normal_normalized.counts.txt.gz", root)
    oa = _read_counts("GSE114007_OA_normalized.counts.txt.gz", root)
    common = normal.index.intersection(oa.index)
    cpm_n = normal.loc[common] / normal.loc[common].sum() * 1e6
    cpm_o = oa.loc[common] / oa.loc[common].sum() * 1e6
    rng = np.random.default_rng(seed)
    rows = {}
    for g in GENES:
        a, b = cpm_o.loc[g].to_numpy(), cpm_n.loc[g].to_numpy()
        est = np.log2(a.mean() / b.mean())
        boots = [np.log2(rng.choice(a, len(a)).mean() / rng.choice(b, len(b)).mean()) for _ in range(n_boot)]
        rows[g] = dict(log2_ratio=float(est), ci_low=float(np.percentile(boots, 2.5)),
                       ci_high=float(np.percentile(boots, 97.5)), n_oa=len(a), n_normal=len(b))
    return pd.DataFrame(rows).T


# ---------------------------------------------------------------- model and criteria
def model_log2(stage_table: pd.DataFrame) -> Dict[str, float]:
    w = stage_table["n"].astype(float)
    dam, intact = [3, 4], [0]

    def wmean(col, rows):
        return float((stage_table.loc[rows, col] * w.loc[rows]).sum() / w.loc[rows].sum())

    out = {}
    for g in GENES:
        a, b = wmean(MODEL_COLUMN[g], dam), wmean(MODEL_COLUMN[g], intact)
        out[g] = float(np.log2(a / b)) if a > 0 and b > 0 else (float("inf") if a > 0 else float("-inf"))
    return out


def criterion_c1(mean: float, lo: float, hi: float, model: float) -> bool:
    if lo > 0 or hi < 0:
        return bool(np.sign(model) == np.sign(mean))
    return bool(abs(model) < 1.0)


def criterion_c2(lo: float, hi: float, model: float) -> bool:
    return bool(lo <= model <= hi)
