"""Mechanical bridge: from joint loads and plate geometry to the effective (time-averaged, drained) stress of the ABM.

1. Compartment force.  The OrthoLoad gait cycle gives F_med(t), F_lat(t) in multiples of body weight
   (:mod:`chondrotwin.knee.loads`); their cycle means f_med, f_lat are the walking load of each compartment.
2. Nominal walking stress.  L_c = f_c * BW / A_c, with A_c the tibial plate footprint (OAIZIB-CM).
3. Time averaging.  The consolidation time of cartilage (tau = h^2 / (H_A k) ~ 1 h) is thousands of gait
   cycles, so the solid matrix responds to the time-averaged stress (checked with
   :meth:`ConfinedCompression.solve_history`; see VALIDATION.md).  The ABM is driven by a daily effective
   stress sigma_c = kappa * L_c.  The single factor kappa lumps everything that is common to both
   compartments and not observed here: daily walking duty cycle, other activities, the share of the load
   carried by the menisci, the true contact area, and the ABM's stress scale.  kappa is calibrated on
   data (:mod:`chondrotwin.knee.analysis`), never set by hand.

Assumptions that are *not* absorbed by kappa (and are therefore tested or flagged): the medial share of the
force (from 8 total-knee-replacement subjects, not from the OAI knees), equal meniscal load share in both
compartments, and body height from sex-specific population means (BMI is given, height is not).
"""
from __future__ import annotations

import numpy as np

G = 9.81
MEAN_HEIGHT_M = {1: 1.76, 2: 1.62}   # OAI sex code 1 = male, 2 = female; population means (assumption)


def body_weight_N(bmi, sex) -> np.ndarray:
    bmi = np.asarray(bmi, float)
    h = np.vectorize(lambda s: MEAN_HEIGHT_M.get(int(s), 1.69))(np.asarray(sex))
    return bmi * h ** 2 * G


def nominal_walking_stress(force_bw: float, body_weight: float, area_mm2: float) -> float:
    """Cycle-mean compartment force (x BW) times body weight over plate area, in Pa."""
    return force_bw * body_weight / (area_mm2 * 1e-6)


def stress_ratio_lateral_to_medial(medial_ratio: float, area_medial_mm2: float, area_lateral_mm2: float) -> float:
    """sigma_lat / sigma_med for a given medial share of the axial force."""
    return (1.0 - medial_ratio) / medial_ratio * area_medial_mm2 / area_lateral_mm2
