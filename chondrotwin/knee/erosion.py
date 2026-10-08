"""Dose-response of the agent-based model: depth of matrix-depleted cartilage versus effective stress.

A representative column of the plate (0.4 x 0.4 mm, full thickness h) is simulated under a constant daily
effective stress.  "Eroded depth" is the depth from the surface over which the voxel matrix content stays below
``M_crit`` (fibrillated / lost tissue); it is the model counterpart of central cartilage thickness loss.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from ..abm import CartilageABM, CartilageSlab, Scenario


def eroded_depth(slab: CartilageSlab, M_crit: float = 0.3) -> float:
    Mz = slab.M_grid.mean(axis=(0, 1))
    if Mz[0] >= M_crit:
        return 0.0
    ok = np.flatnonzero(Mz >= M_crit)
    return float((ok[0] if ok.size else len(Mz)) * slab.dz)


def dose_response(sigmas_MPa: Iterable[float], h_mm: float, years: int = 5, seeds: Iterable[int] = (0, 1),
                  M_crit: float = 0.3, cache: Optional[Path] = None, rates=None,
                  initial_hypertrophic_fraction: float = 0.0) -> pd.DataFrame:
    """Eroded depth (mm) at the end of each simulated year for every stress, thickness and seed."""
    if cache is not None and cache.exists():
        return pd.read_csv(cache)
    rows = []
    nz = max(int(round(h_mm / 0.1)), 5)
    for sig in sigmas_MPa:
        for seed in seeds:
            kw = {"rates": rates} if rates is not None else {}
            sim = CartilageABM(slab=CartilageSlab(Lx=0.4, Ly=0.4, h=h_mm, nx=4, ny=4, nz=nz), seed=seed, **kw)
            sc = Scenario("dose", sigma0=float(sig) * 1e6, days=365 * years,
                          initial_hypertrophic_fraction=initial_hypertrophic_fraction)
            sim.seed_fates(hypertrophic_fraction=initial_hypertrophic_fraction)
            for d in range(sc.days):
                sim.step(d, sc)
                if (d + 1) % 365 == 0:
                    rows.append(dict(sigma_MPa=float(sig), h_mm=h_mm, seed=seed, year=(d + 1) // 365,
                                     eroded_mm=eroded_depth(sim.slab, M_crit)))
    df = pd.DataFrame(rows)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(cache, index=False)
    return df


def curve(df: pd.DataFrame, h_mm: float, year: int) -> tuple:
    """(sigma, mean eroded depth) for one thickness and horizon, monotone in sigma."""
    d = df[(np.isclose(df.h_mm, h_mm)) & (df.year == year)].groupby("sigma_MPa")["eroded_mm"].mean()
    s, e = d.index.to_numpy(), np.maximum.accumulate(d.to_numpy())
    return s, e


def invert(s: np.ndarray, e: np.ndarray, target: float) -> float:
    """Smallest stress whose eroded depth reaches ``target`` (linear interpolation); nan if never reached."""
    if target <= e[0]:
        return float(s[0])
    if target > e[-1]:
        return float("nan")
    k = int(np.argmax(e >= target))
    return float(s[k - 1] + (target - e[k - 1]) * (s[k] - s[k - 1]) / (e[k] - e[k - 1])) if e[k] > e[k - 1] else float(s[k])
