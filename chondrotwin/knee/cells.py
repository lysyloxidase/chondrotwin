"""Chondrocyte state composition by degeneration stage from single-cell RNA-seq (Ji et al., Ann Rheum Dis 2019; GEO GSE104782).

Cells are named ``OA<patient>_S<stage>.<n>``: 10 patients undergoing knee arthroplasty, cartilage pieces of
OARSI-type stage 0-4, 32 cells per patient and stage.  The authors' clusters are used as published:
EC (effector), HomC (homeostatic), RegC (regulatory), ProC (proliferative), preHTC (prehypertrophic),
HTC (hypertrophic), FC (fibrocartilage); 136 cells failed quality control (no cluster).

Mapping to the model's fates (an assumption, stated so it can be argued with):

* articular (SOX9-high)     <- EC, HomC, RegC
* hypertrophic (RUNX2-high) <- preHTC, HTC
* proliferative / cluster   <- ProC
* dedifferentiated          <- FC
"""
from __future__ import annotations

import pandas as pd

from .. import datasets

CLUSTER_TO_FATE = {"EC": "articular", "HomC": "articular", "RegC": "articular", "preHTC": "hypertrophic",
                   "HTC": "hypertrophic", "ProC": "proliferative", "FC": "dedifferentiated"}
FATES = ("articular", "hypertrophic", "proliferative", "dedifferentiated")


def ji2019_cells(root=None) -> pd.DataFrame:
    p = datasets.require("gse104782", "GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx", root)
    df = pd.read_excel(p)
    m = df["Cell"].str.extract(r"^(OA\d+)_S(\d)\.(\d+)$")
    if m.isna().any().any():
        raise ValueError("unexpected cell names in GSE104782 table")
    df["patient"] = m[0]
    df["stage"] = m[1].astype(int)
    df["fate"] = df["Cluster"].map(CLUSTER_TO_FATE)
    return df


def composition_by_stage(cells: pd.DataFrame, patients=None) -> pd.DataFrame:
    """Fraction of each model fate among clustered cells, per stage (optionally for a subset of patients)."""
    d = cells.dropna(subset=["fate"])
    if patients is not None:
        d = d[d["patient"].isin(list(patients))]
    tab = pd.crosstab(d["stage"], d["fate"], normalize="index")
    return tab.reindex(columns=list(FATES), fill_value=0.0)
