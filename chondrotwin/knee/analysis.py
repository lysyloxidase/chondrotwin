"""Calibration and out-of-sample tests of the agent-based model against open data.

Tissue test (OAIZIB-CM, cross-sectional, 507 OAI knees, the dataset's own train / test split)
    Observed: central tibial cartilage thickness in KL4 versus KL0 knees, adjusted for sex, BMI and age
    (OLS), separately for the medial and lateral plate.
    Calibration (train split only): the effective medial stress at which the ABM erodes as much medial
    cartilage as observed.  This fixes kappa (see :mod:`chondrotwin.knee.bridge`).
    Prediction: lateral loss, from the lateral/medial stress ratio implied by the in vivo medial force share
    (OrthoLoad) and the plate areas (OAIZIB-CM).  Compared with the observed lateral loss in the train and
    the held-out test split.  Nothing lateral is used for calibration.

Cell test (GSE104782, 10 patients x 5 stages x 32 cells)
    Dropout-robust statistic: ratio of the fraction of cells with detectable RUNX2 / MMP13 / COL10A1 /
    CDKN2A (p16) transcripts at stage 4 versus stage 0, with 95 % patient-bootstrap intervals; compared with
    the ABM's ratio of RUNX2-high (R > S) and senescent cells between degraded and intact tissue.
"""
from __future__ import annotations

from typing import Dict, Sequence

import numpy as np
import pandas as pd

from . import erosion
from .bridge import stress_ratio_lateral_to_medial


# ------------------------------------------------------------------ observed tissue losses
def adjusted_kl_effects(geo: pd.DataFrame, outcome: str) -> pd.DataFrame:
    """OLS of ``outcome`` on sex, BMI, age and KL1-4 dummies; returns KL effects (vs KL0) with SE and n."""
    d = geo.dropna(subset=["KLGrade", "BMI", outcome])
    X = np.column_stack([np.ones(len(d)), (d.Gender == 2).astype(float), d.BMI, d.Age]
                        + [(d.KLGrade == k).astype(float) for k in (1, 2, 3, 4)])
    y = d[outcome].to_numpy()
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    cov = (r @ r) / (len(y) - X.shape[1]) * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    return pd.DataFrame({"effect": b[4:], "se": se[4:], "n": [int((d.KLGrade == k).sum()) for k in (1, 2, 3, 4)]},
                        index=pd.Index([1, 2, 3, 4], name="KL"))


def observed_losses(geo: pd.DataFrame) -> pd.DataFrame:
    """Adjusted central thickness loss (KL4 vs KL0, mm, positive = thinner) per compartment and split."""
    rows = []
    for split in ("train", "test", "all"):
        g = geo if split == "all" else geo[geo.split == split]
        for side in ("medial", "lateral"):
            e = adjusted_kl_effects(g, f"{side}_central_thickness_mm").loc[4]
            rows.append(dict(split=split, compartment=side, loss_mm=-e.effect, se_mm=e.se, n_KL4=int(e.n),
                             n_KL0=int(((g.KLGrade == 0) & g.BMI.notna()).sum())))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ tissue test
def tissue_test(geo: pd.DataFrame, medial_ratios: Dict[str, float], dose: pd.DataFrame,
                h_medial: float, h_lateral: float, year: int = 5) -> pd.DataFrame:
    """Calibrate on train medial loss, predict lateral loss for every medial-force share in ``medial_ratios``."""
    obs = observed_losses(geo).set_index(["split", "compartment"])
    target = float(obs.loc[("train", "medial"), "loss_mm"])
    s_m, e_m = erosion.curve(dose, h_medial, year)
    s_l, e_l = erosion.curve(dose, h_lateral, year)
    sigma_med = erosion.invert(s_m, e_m, target)
    kl4 = geo[(geo.split == "train") & (geo.KLGrade == 4)]
    a_med = float((kl4.medial_area_mm2 + kl4.medial_denuded_mm2).mean())
    a_lat = float((kl4.lateral_area_mm2 + kl4.lateral_denuded_mm2).mean())
    rows = []
    for name, mr in medial_ratios.items():
        r = stress_ratio_lateral_to_medial(mr, a_med, a_lat)
        sig_l = r * sigma_med
        pred = float(np.interp(sig_l, s_l, e_l, left=e_l[0], right=np.nan))
        rows.append(dict(source=name, medial_force_share=mr, stress_ratio_lat_med=r, sigma_medial_MPa=sigma_med,
                         sigma_lateral_MPa=sig_l, predicted_lateral_loss_mm=pred))
    out = pd.DataFrame(rows)
    out.attrs.update(target_medial_loss_mm=target, area_medial_mm2=a_med, area_lateral_mm2=a_lat, year=year,
                     threshold_lateral_MPa=float(s_l[np.argmax(e_l > 0)]) if np.any(e_l > 0) else float("nan"))
    return out


def mixture_test(geo: pd.DataFrame, medial_ratios: Sequence[float], dose: pd.DataFrame,
                 h_medial: float, h_lateral: float, year: int = 5) -> dict:
    """Population version of :func:`tissue_test`.

    KL4 knees differ in alignment; each is assigned one of the in vivo medial force shares with equal
    weight.  The common load level (kappa * BW / area) is calibrated so that the *mean* medial loss over
    the mixture equals the observed train-split medial loss; the prediction is the mean lateral loss.
    """
    obs = observed_losses(geo).set_index(["split", "compartment"])
    target = float(obs.loc[("train", "medial"), "loss_mm"])
    s_m, e_m = erosion.curve(dose, h_medial, year)
    s_l, e_l = erosion.curve(dose, h_lateral, year)
    kl4 = geo[(geo.split == "train") & (geo.KLGrade == 4)]
    a_ratio = float((kl4.medial_area_mm2 + kl4.medial_denuded_mm2).mean() / (kl4.lateral_area_mm2 + kl4.lateral_denuded_mm2).mean())
    mr = np.asarray(list(medial_ratios), float)

    def losses(c):          # c: total walking stress scale (MPa) applied to the medial share
        sm, sl = c * mr, c * (1 - mr) * a_ratio
        em = np.interp(sm, s_m, e_m, left=e_m[0], right=np.nan)
        el = np.interp(sl, s_l, e_l, left=e_l[0], right=np.nan)
        return em, el

    grid = np.linspace(0.05, 0.25, 4001)
    mean_med = np.array([np.nanmean(losses(c)[0]) if not np.any(np.isnan(losses(c)[0])) else np.nan for c in grid])
    ok = np.flatnonzero(mean_med >= target)
    if not ok.size:
        return dict(target_medial_loss_mm=target, calibrated=False)
    c = float(grid[ok[0]])
    em, el = losses(c)
    return dict(target_medial_loss_mm=target, calibrated=True, load_scale_MPa=c, mean_medial_loss_mm=float(np.mean(em)),
                predicted_mean_lateral_loss_mm=float(np.mean(el)), per_knee_medial=em.round(3).tolist(),
                per_knee_lateral=el.round(3).tolist(), area_ratio_med_lat=a_ratio, year=year)


# ------------------------------------------------------------------ cell test
def marker_detection(umi: pd.DataFrame, cells: pd.DataFrame, genes: Sequence[str]) -> pd.DataFrame:
    det = (umi.loc[list(genes)] > 0).T
    return det.join(cells.set_index("Cell")[["patient", "stage"]], how="inner")


def stage_ratio_bootstrap(det: pd.DataFrame, gene: str, stage: int = 4, n_boot: int = 2000, seed: int = 0):
    """Detection fraction at ``stage`` / stage 0 with a 95 % bootstrap interval over patients."""
    pats = det["patient"].unique()
    frac = det.groupby(["patient", "stage"])[gene].agg(["sum", "count"])

    def ratio(sample):
        s = frac.loc[[(p, st) for p in sample for st in (0, stage)]].groupby(level="stage").sum()
        f0, f4 = s.loc[0, "sum"] / s.loc[0, "count"], s.loc[stage, "sum"] / s.loc[stage, "count"]
        return f4 / f0 if f0 > 0 else np.inf, f0, f4

    est, f0, f4 = ratio(pats)
    rng = np.random.default_rng(seed)
    boots = np.array([ratio(rng.choice(pats, len(pats)))[0] for _ in range(n_boot)])
    boots = boots[np.isfinite(boots)]
    return dict(gene=gene, frac_stage0=f0, **{f"frac_stage{stage}": f4}, ratio=est,
                ci_low=float(np.percentile(boots, 2.5)), ci_high=float(np.percentile(boots, 97.5)))


def abm_stage_composition(sigmas_MPa=(0.07, 0.08, 0.09, 0.10, 0.12), years: int = 2, seed: int = 1,
                          bins=(-1, 0.2, 0.4, 0.6, 0.8, 2), rates=None, slab_mm=(0.6, 0.6, 2.2),
                          initial_hypertrophic_fraction: float = 0.0, C_joint: float = 0.0,
                          return_cells: bool = False):
    """Live-cell composition of the ABM by local degeneration stage (voxel matrix bins: M >= 0.8 -> stage 0 ... M < 0.2 -> 4).

    Columns: runx2_high (R > S), senescent, clustered, catabolic (R > S or k_cat * overload > 0.5),
    mean_S, mean_R, mean_cat (mean catabolic drive R + k_cat * overload) and the number of live cells n.
    """
    from ..abm import CartilageABM, CartilageSlab, Rates, Scenario
    rates = rates or Rates()
    lx, ly, h = slab_mm
    parts = []
    for sig in sigmas_MPa:
        sim = CartilageABM(slab=CartilageSlab(Lx=lx, Ly=ly, h=h, nx=int(round(lx / 0.1)), ny=int(round(ly / 0.1)),
                                              nz=int(round(h / 0.1))), seed=seed, rates=rates)
        sim.run(Scenario("cells", sigma0=sig * 1e6, days=365 * years, initial_hypertrophic_fraction=initial_hypertrophic_fraction,
                         C_joint=C_joint))
        ag, s = sim.agents, sim.slab
        ix, iy, iz = s.voxel_index(ag.x, ag.y, ag.z)
        live = ag.alive
        m = ag.strain[live] / rates.strain_scale                      # ag.strain is the cell (PCM-amplified) strain
        overload = 1.0 / (1.0 + np.exp(-(m - rates.overload_threshold) / rates.overload_width))
        infl = ag.C[live] / (rates.K_C + ag.C[live])
        drive = rates.k_cat * overload + rates.k_inf * infl
        parts.append(pd.DataFrame({"M": s.M_grid[ix, iy, iz][live], "runx2_high": (ag.R > ag.S)[live],
                                   "senescent": ag.senescent[live], "clustered": (ag.cluster > 0)[live],
                                   "catabolic": (ag.R > ag.S)[live] | (drive > 0.5),
                                   "S": ag.S[live], "R": ag.R[live], "cat": ag.R[live] + drive}))
    d = pd.concat(parts, ignore_index=True)
    d["stage"] = pd.cut(d["M"], list(bins), labels=[4, 3, 2, 1, 0]).astype(int)
    g = d.groupby("stage")
    out = g[["runx2_high", "senescent", "clustered", "catabolic"]].mean()
    out["mean_S"] = g["S"].mean()
    out["mean_R"] = g["R"].mean()
    out["mean_cat"] = g["cat"].mean()
    out["n"] = g.size()
    out = out.reindex(range(5))
    if return_cells:
        return out, d
    return out


def bulk_means(cells: pd.DataFrame) -> dict:
    """Bulk-tissue counterparts (means over all live cells): S, R, catabolic drive, senescent fraction."""
    return {"mean_S": float(cells["S"].mean()), "mean_R": float(cells["R"].mean()),
            "mean_cat": float(cells["cat"].mean()), "senescent": float(cells["senescent"].mean())}
