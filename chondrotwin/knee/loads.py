"""In vivo knee contact loads from OrthoLoad (Bergmann et al., PLoS ONE 2014) and their medial/lateral split.

The telemetric INNEX implants measure the six load components on the tibial
tray (X lateral, Y anterior, Z superior; origin on the tibial plateau at the
stem axis).  The varus-valgus moment My arises from the axial force acting off
the stem axis, so with two contact points a distance ``l`` apart

    F_med = |Fz| / 2 - My / l,      F_lat = |Fz| - F_med

(Kutzner et al., J Biomech 2011; PLoS ONE 2013).  ``l`` is not given in the open
articles; :func:`calibrate_condylar_distance` chooses it so that the mean early-
and late-stance peaks of F_med across the OrthoLoad subjects match the published
values for the same cohort (176 and 187 %BW; Kutzner et al., PLoS ONE 2013).

Caveat: these are loads in total knee replacements (cruciate-sacrificing,
ultracongruent insert), not in native knees with menisci and cruciate ligaments.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Sequence

import numpy as np
import pandas as pd

from .. import datasets

_trapz = getattr(np, "trapezoid", None) or np.trapz   # NumPy >= 2.0 / < 2.0

SUBJECTS = ("K1L", "K2L", "K3R", "K5R", "K6L", "K7L", "K8L", "K9L")
KUTZNER_2013_PEAKS_BW = (1.76, 1.87)      # mean peak medial force, early / late stance (x body weight)


@dataclass(frozen=True)
class GaitLoads:
    subject: str
    body_weight: float          # N
    varus_deg: float            # static leg alignment (varus > 0)
    t: np.ndarray               # s, one gait cycle
    Fz: np.ndarray              # N, axial force on the tibia (compression positive here)
    My: np.ndarray              # Nm, varus-valgus moment (positive: lateral shift of the axial force)

    @property
    def cycle_time(self) -> float:
        return float(self.t[-1] - self.t[0])

    def split(self, l: float):
        """(F_med, F_lat) in N for condylar distance ``l`` (m); clipped to [0, |Fz|]."""
        f_med = np.clip(self.Fz / 2.0 - self.My / l, 0.0, self.Fz)
        return f_med, self.Fz - f_med

    def medial_ratio(self, l: float) -> float:
        """Impulse-weighted share of the axial force carried by the medial compartment over the cycle."""
        f_med, _ = self.split(l)
        return float(_trapz(f_med, self.t) / _trapz(self.Fz, self.t))


def read_orthoload_walking(subject: str, root=None) -> GaitLoads:
    p = datasets.require("orthoload", f"{subject}_Walking.xlsx", root)
    raw = pd.read_excel(p, sheet_name="Dat_Individ", header=None)
    col0 = raw.iloc[:, 0].astype(str).str.strip()
    body_weight = float(raw.iloc[int(np.flatnonzero(col0.str.startswith("BodyWeight"))[0]), 1])
    varus = float(str(raw.iloc[int(np.flatnonzero(col0.str.startswith("Varus Angle"))[0]), 1]).replace(",", "."))
    n = int(raw.iloc[int(np.flatnonzero(col0.str.startswith("Number of Data Sets"))[0]), 1])
    hdr = int(np.flatnonzero(col0 == "Time")[0])
    if list(raw.iloc[hdr, :9]) != ["Time", "Fx", "Fy", "Fz", "Fres", "Mx", "My", "Mz", "Mres"]:
        raise ValueError(f"{p}: unexpected column layout {list(raw.iloc[hdr, :9])}")
    d = raw.iloc[hdr + 3: hdr + 3 + n, :9].astype(float).to_numpy()
    return GaitLoads(subject, body_weight, varus, d[:, 0], -d[:, 3], d[:, 6])


def load_all(subjects: Iterable[str] = SUBJECTS, root=None) -> Dict[str, GaitLoads]:
    return {s: read_orthoload_walking(s, root) for s in subjects}


def medial_peaks_bw(g: GaitLoads, l: float, split: float = 0.30, stance_end: float = 0.65, min_force: float = 1000.0):
    """Early / late stance peaks of F_med in x body weight (|Fz| > ``min_force``, as in Kutzner et al. 2013)."""
    f_med, _ = g.split(l)
    phase = (g.t - g.t[0]) / g.cycle_time
    ok = g.Fz > min_force
    early = f_med[ok & (phase < split)].max()
    late = f_med[ok & (phase >= split) & (phase < stance_end)].max()
    return early / g.body_weight, late / g.body_weight


def calibrate_condylar_distance(loads: Dict[str, GaitLoads], target=KUTZNER_2013_PEAKS_BW,
                                grid: Sequence[float] = np.arange(0.030, 0.0701, 0.0005)) -> float:
    """Least-squares condylar distance (m) reproducing the published mean medial peaks."""
    def err(l):
        pk = np.array([medial_peaks_bw(g, l) for g in loads.values()]).mean(axis=0)
        return float(np.sum((pk - np.asarray(target)) ** 2))
    return float(min(grid, key=err))


def summary(loads: Dict[str, GaitLoads], l: float) -> pd.DataFrame:
    rows = []
    for s, g in loads.items():
        f_med, f_lat = g.split(l)
        e, lt = medial_peaks_bw(g, l)
        rows.append(dict(subject=s, body_weight_N=g.body_weight, varus_deg=g.varus_deg, cycle_s=g.cycle_time,
                         peak_total_BW=g.Fz.max() / g.body_weight, mean_total_BW=_trapz(g.Fz, g.t) / g.cycle_time / g.body_weight,
                         mean_medial_BW=_trapz(f_med, g.t) / g.cycle_time / g.body_weight,
                         mean_lateral_BW=_trapz(f_lat, g.t) / g.cycle_time / g.body_weight,
                         medial_ratio=g.medial_ratio(l), peak_medial_early_BW=e, peak_medial_late_BW=lt))
    return pd.DataFrame(rows).set_index("subject")
