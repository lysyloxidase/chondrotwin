"""
Fate landscapes (paper, Sections 2.1, 5.7 and Figure 1).

Two kinds of landscape are provided:

1. :func:`schematic_landscape` - the drawable surface of Figure 1a, a sum of
   Gaussian wells placed at the fate states of Table 1.  Purely illustrative.
2. :func:`quasi_potential_langevin` - a *computable* landscape: the toggle
   switch of :mod:`chondrotwin.toggle` is simulated as a Langevin (stochastic)
   system with additive noise, the stationary occupancy P(S, R) is estimated by
   histogramming many trajectories, and the quasi-potential is U = -ln P
   (Wang et al., PNAS 2011).  Basin depths and barrier heights can then be read
   off numerically - this is what "computable attractor landscape" means in the
   paper.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import numpy as np

from .toggle import ToggleSwitch

# (name, x, y, depth, width) - schematic basins of Figure 1a
DEFAULT_BASINS: List[Tuple[str, float, float, float, float]] = [
    ("Chondroprogenitor", 1.3, 7.9, 1.4, 1.0),
    ("Proliferative", 3.7, 8.5, 1.5, 0.9),
    ("Prehypertrophic / hypertrophic", 6.3, 7.7, 1.8, 1.0),
    ("Osteoblast", 8.8, 8.4, 1.5, 0.8),
    ("Permanent / articular", 3.0, 5.0, 2.2, 1.25),
    ("Cluster", 1.0, 3.5, 0.9, 0.7),
    ("Dedifferentiated / fibroblastic", 6.1, 4.7, 1.2, 0.9),
    ("Senescent", 4.3, 1.8, 1.9, 1.0),
    ("Ferroptotic", 7.1, 1.5, 1.7, 0.85),
    ("Apoptotic", 8.8, 3.5, 1.6, 0.85),
]


def schematic_landscape(basins: Sequence[Tuple[str, float, float, float, float]] = DEFAULT_BASINS,
                        n: int = 400, baseline: float = 2.4, tilt: Tuple[float, float] = (0.0, 0.0)):
    """Schematic quasi-potential U(x, y) on [0, 10]^2 as a sum of Gaussian wells.

    ``tilt = (gx, gy)`` adds a linear term gx*x + gy*y, the "mechanical tilt" of
    Figure 1.  Returns (X, Y, U).
    """
    X, Y = np.meshgrid(np.linspace(0, 10, n), np.linspace(0, 10, n))
    U = np.full_like(X, baseline)
    for _, cx, cy, depth, width in basins:
        U -= depth * np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / (2 * width ** 2))
    U += tilt[0] * X + tilt[1] * Y
    return X, Y, U


def cross_section(s: np.ndarray, tilt: float = 0.0) -> np.ndarray:
    """Double-well cross-section U(s) = s^4 - 2 s^2 + tilt * s (Figure 1b)."""
    s = np.asarray(s, dtype=float)
    return s ** 4 - 2 * s ** 2 + tilt * s


def quasi_potential_langevin(model: ToggleSwitch | None = None, lam: float = 0.0, noise: float = 0.25,
                             n_walkers: int = 3000, n_steps: int = 6000, burn_in: int = 1500,
                             dt: float = 0.01, s_max: float = 3.0, bins: int = 60,
                             seed: int = 0, start_fraction_sox9: float = 0.5) -> Dict[str, np.ndarray]:
    """Quasi-potential U = -ln P of the toggle switch from Langevin dynamics.

    Many walkers are integrated with Euler-Maruyama (additive Gaussian noise of
    amplitude ``noise``) and their positions after ``burn_in`` steps are
    histogrammed on a (S, R) grid.  Returns a dict with ``S_edges``, ``R_edges``,
    ``P`` (occupancy probability), ``U`` (quasi-potential, min 0, NaN where
    unvisited) and ``crossed`` (fraction of walkers that crossed the diagonal
    S = R at least once).

    A fraction ``start_fraction_sox9`` (default one half) of the walkers starts
    in the SOX9-high basin, the rest in the RUNX2-high basin.  When ``crossed`` is small the
    sampling has not equilibrated *between* basins, so P is quasi-stationary
    within each basin only: the relative depth of the two basins (and the
    barrier seen from the less-populated side) then reflects the initial split
    rather than the dynamics.  Exceptions: at lam = 0 the switch is symmetric and
    50/50 is the exact stationary split; and the escape barrier from the SOX9-high
    basin for lam > 0 is insensitive to the initial split (see VALIDATION.md).
    Barrier heights in units of -ln P scale roughly as 1/noise^2.
    """
    model = model or ToggleSwitch()
    rng = np.random.default_rng(seed)
    # start walkers near both attractors so both basins are sampled
    y = np.empty((n_walkers, 2))
    half = int(start_fraction_sox9 * n_walkers)
    y[:half] = [2.0, 0.1]; y[half:] = [0.1, 2.0]
    y += 0.1 * rng.standard_normal(y.shape)
    side0 = y[:, 0] > y[:, 1]
    crossed = np.zeros(n_walkers, dtype=bool)
    H = np.zeros((bins, bins))
    edges = np.linspace(0.0, s_max, bins + 1)
    sq = noise * np.sqrt(dt)
    for step in range(n_steps):
        S, R = y[:, 0], y[:, 1]
        dS = model.alpha_s * (1 - lam) * model.hill_rep(R) + model.beta_s - model.gamma * S
        dR = model.alpha_r * (1 + lam) * model.hill_rep(S) + model.beta_r - model.gamma * R
        y[:, 0] += dt * dS + sq * rng.standard_normal(n_walkers)
        y[:, 1] += dt * dR + sq * rng.standard_normal(n_walkers)
        np.maximum(y, 0.0, out=y)
        crossed |= (y[:, 0] > y[:, 1]) != side0
        if step >= burn_in:
            h, _, _ = np.histogram2d(y[:, 0], y[:, 1], bins=[edges, edges])
            H += h
    P = H / H.sum()
    with np.errstate(divide="ignore"):
        U = -np.log(P)
    U[~np.isfinite(U)] = np.nan
    U -= np.nanmin(U)
    return {"S_edges": edges, "R_edges": edges, "P": P, "U": U, "crossed": float(crossed.mean())}


def basin_depths(U: np.ndarray, S_edges: np.ndarray, R_edges: np.ndarray) -> Dict[str, float]:
    """Minimum U in the SOX9-high (S > R) and RUNX2-high (R > S) halves, and the lowest ridge value between them.

    The ridge estimate is the minimum of U along the diagonal S = R, which is a
    crude but informative proxy of the barrier for a symmetric switch.  If no
    walker ever visited the diagonal, the ridge is reported as a lower bound
    (``ridge_is_lower_bound`` is True); increase ``noise`` or the sampling.
    """
    Sc = 0.5 * (S_edges[:-1] + S_edges[1:]); Rc = 0.5 * (R_edges[:-1] + R_edges[1:])
    SS, RR = np.meshgrid(Sc, Rc, indexing="ij")
    finite_max = float(np.nanmax(U))
    Uf = np.where(np.isnan(U), np.inf, U)
    sox = float(np.min(Uf[SS > RR])); runx = float(np.min(Uf[RR > SS]))
    diag = np.array([Uf[i, i] for i in range(len(Sc))])
    lower_bound = bool(np.all(~np.isfinite(diag)))
    # unvisited diagonal bins are at least as high as the highest visited bin
    diag = np.where(np.isfinite(diag), diag, finite_max)
    ridge = float(np.min(diag))
    return {"U_min_SOX9_high": sox, "U_min_RUNX2_high": runx, "ridge_on_diagonal": ridge,
            "barrier_from_SOX9": ridge - sox, "barrier_from_RUNX2": ridge - runx,
            "ridge_is_lower_bound": lower_bound}
