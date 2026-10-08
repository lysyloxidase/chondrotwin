"""Tibial cartilage geometry from the OAIZIB-CM segmentations (507 OAI knees, medial / lateral tibial cartilage split).

Labels: 1 femur, 2 femoral cartilage, 3 tibia, 4 medial tibial cartilage, 5 lateral tibial cartilage.
Volumes are 3D DESS MRI in RAS orientation (voxel 0.7 x 0.365 x 0.365 mm).  For each tibial plate the
cartilage is collapsed along the superior-inferior axis: the projected (axial-plane) area is the number of
covered (R, A) columns times the column cross-section, thickness is the number of cartilage voxels in a
column times the S voxel size, and the mean thickness is volume / projected area.  This "vertical"
thickness slightly overestimates the surface-normal thickness where the plateau is inclined.
"""
from __future__ import annotations

import gzip
import io
import zipfile
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .. import datasets

LABEL_MEDIAL, LABEL_LATERAL = 4, 5


def _load_nifti_bytes(data: bytes):
    import nibabel as nib
    raw = gzip.decompress(data)
    img = nib.Nifti1Image.from_bytes(raw)
    return img


def plate_metrics(labels: np.ndarray, zooms, axcodes, label: int) -> dict:
    s_axis = axcodes.index("S") if "S" in axcodes else axcodes.index("I")
    other = [a for a in range(3) if a != s_axis]
    from scipy import ndimage
    m = labels == label
    col = m.sum(axis=s_axis)                         # voxels per (column)
    covered = col > 0
    col_area = float(zooms[other[0]] * zooms[other[1]])
    vol = float(m.sum() * np.prod(zooms))
    area = float(covered.sum() * col_area)
    th_map = col * float(zooms[s_axis])
    th = th_map[covered]
    # denuded area: holes inside the plate footprint (full-thickness cartilage loss)
    footprint = ndimage.binary_fill_holes(ndimage.binary_closing(covered, iterations=3))
    denuded = float((footprint & ~covered).sum() * col_area)
    # central load-bearing region: footprint columns within half the equivalent radius of the footprint centroid
    ii, jj = np.nonzero(footprint)
    if ii.size:
        ci, cj = ii.mean(), jj.mean()
        r_eq = np.sqrt(footprint.sum() * col_area / np.pi)
        d = np.hypot((ii - ci) * zooms[other[0]], (jj - cj) * zooms[other[1]])
        central = d < 0.5 * r_eq
        th_c = th_map[ii[central], jj[central]]      # zero where denuded
        central_mean = float(th_c.mean())
    else:
        central_mean = np.nan
    return {"volume_mm3": vol, "area_mm2": area, "thickness_mean_mm": vol / area if area else np.nan,
            "thickness_p10_mm": float(np.percentile(th, 10)) if th.size else np.nan,
            "thickness_p90_mm": float(np.percentile(th, 90)) if th.size else np.nan,
            "thin_fraction_lt1mm": float(np.mean(th < 1.0)) if th.size else np.nan,
            "denuded_mm2": denuded, "central_thickness_mm": central_mean}


def subject_info(root=None) -> pd.DataFrame:
    with zipfile.ZipFile(datasets.require("oaizib_cm", "info.zip", root)) as zf:
        tr = pd.read_excel(io.BytesIO(zf.read("info/subInfo_train.xlsx"))).assign(split="train")
        ts = pd.read_excel(io.BytesIO(zf.read("info/subInfo_test.xlsx"))).assign(split="test")
    info = pd.concat([tr, ts], ignore_index=True)
    info["case"] = info["CMT-ID"].map(lambda i: f"oaizib_{int(i):03d}")
    return info


def tibial_geometry(root=None, cache: Optional[Path] = None, progress: bool = False) -> pd.DataFrame:
    """Per-knee medial / lateral tibial cartilage metrics joined with KL grade, sex, age, BMI and train/test split."""
    if cache is None:
        cache = datasets.data_dir(root).parent / "derived" / "oaizib_cm_tibia.csv"
    if cache.exists():
        return pd.read_csv(cache, index_col="case")
    rows = []
    for zname, folder in (("labelsTr.zip", "labelsTr"), ("labelsTs.zip", "labelsTs")):
        with zipfile.ZipFile(datasets.require("oaizib_cm", zname, root)) as zf:
            for name in sorted(n for n in zf.namelist() if n.endswith(".nii.gz")):
                img = _load_nifti_bytes(zf.read(name))
                import nibabel as nib
                a = np.asanyarray(img.dataobj)
                zooms = img.header.get_zooms()[:3]
                ax = nib.aff2axcodes(img.affine)
                rec = {"case": Path(name).name.replace(".nii.gz", ""), "folder": folder}
                for lab, side in ((LABEL_MEDIAL, "medial"), (LABEL_LATERAL, "lateral")):
                    for k, v in plate_metrics(a, zooms, ax, lab).items():
                        rec[f"{side}_{k}"] = v
                rows.append(rec)
                if progress and len(rows) % 50 == 0:
                    print(f"  {len(rows)} knees")
    geo = pd.DataFrame(rows).set_index("case")
    info = subject_info(root).set_index("case")
    out = geo.join(info[["SubjectID", "KneeSide", "KLGrade", "Gender", "Age", "BMI", "split"]], how="left")
    out["ml_thickness_ratio"] = out["medial_thickness_mean_mm"] / out["lateral_thickness_mean_mm"]
    out["ml_area_ratio"] = out["medial_area_mm2"] / out["lateral_area_mm2"]
    cache.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(cache)
    return out
