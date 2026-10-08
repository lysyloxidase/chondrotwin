"""
The SOX9-RUNX2 toggle switch (paper, Section 3.1 and Figure 4).

Two transcription factors, SOX9 (S) and RUNX2 (R), repress each other through
Hill-type (cooperative) inhibition:

    dS/dt = alpha_s (1 - lam) / (1 + (R/K)^n) + beta_s - gamma S
    dR/dt = alpha_r (1 + lam) / (1 + (S/K)^n) + beta_r - gamma R

The bias ``lam`` in (-1, 1) represents the net catabolic (lam > 0: overload,
inflammation, WNT, HIF-2alpha) or anabolic (lam < 0: TGF-beta, physiological
loading, PTHrP) input integrated by the network.  For the default parameters the
system is bistable around lam = 0, which produces hysteresis: a transient
catabolic pulse can switch a cell permanently from the SOX9-high (articular)
to the RUNX2-high (hypertrophic) attractor.

Units are dimensionless; the parameters are illustrative.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Sequence, Tuple, Union

import numpy as np

LamLike = Union[float, Callable[[float], float]]


@dataclass(frozen=True)
class FixedPoint:
    """A steady state of the toggle switch."""
    S: float
    R: float
    stable: bool
    eigenvalues: Tuple[float, float]

    @property
    def label(self) -> str:
        if not self.stable:
            return "saddle"
        return "articular (SOX9-high)" if self.S > self.R else "hypertrophic (RUNX2-high)"


@dataclass
class ToggleSwitch:
    alpha_s: float = 2.0   # maximal SOX9 production
    alpha_r: float = 2.0   # maximal RUNX2 production
    K: float = 0.5         # half-maximal repression constant
    n: float = 4.0         # Hill coefficient (cooperativity)
    beta_s: float = 0.05   # basal SOX9 production
    beta_r: float = 0.05   # basal RUNX2 production
    gamma: float = 1.0     # first-order decay

    # ------------------------------------------------------------------ core
    def hill_rep(self, x: np.ndarray) -> np.ndarray:
        """Decreasing Hill function 1 / (1 + (x/K)^n)."""
        x = np.asarray(x, dtype=float)
        return 1.0 / (1.0 + (np.maximum(x, 0.0) / self.K) ** self.n)

    def rhs(self, state: Sequence[float], lam: float = 0.0, t: float = 0.0) -> np.ndarray:
        """Right-hand side of the ODE system for state (S, R)."""
        S, R = state
        dS = self.alpha_s * (1.0 - lam) * self.hill_rep(R) + self.beta_s - self.gamma * S
        dR = self.alpha_r * (1.0 + lam) * self.hill_rep(S) + self.beta_r - self.gamma * R
        return np.array([dS, dR], dtype=float)

    def rhs_arrays(self, S: np.ndarray, R: np.ndarray, lam: np.ndarray):
        """Vectorized right-hand side for arrays of cells with per-cell bias (used by the agent-based model)."""
        dS = self.alpha_s * (1.0 - lam) * self.hill_rep(R) + self.beta_s - self.gamma * S
        dR = self.alpha_r * (1.0 + lam) * self.hill_rep(S) + self.beta_r - self.gamma * R
        return dS, dR

    def nullcline_S(self, R: np.ndarray, lam: float = 0.0) -> np.ndarray:
        """S as a function of R on the dS/dt = 0 nullcline."""
        return (self.alpha_s * (1.0 - lam) * self.hill_rep(R) + self.beta_s) / self.gamma

    def nullcline_R(self, S: np.ndarray, lam: float = 0.0) -> np.ndarray:
        """R as a function of S on the dR/dt = 0 nullcline."""
        return (self.alpha_r * (1.0 + lam) * self.hill_rep(S) + self.beta_r) / self.gamma

    def jacobian(self, state: Sequence[float], lam: float = 0.0) -> np.ndarray:
        S, R = state
        n, K = self.n, self.K

        def d_hill(x):
            x = max(x, 1e-12)
            return -n * (x / K) ** (n - 1) / K / (1.0 + (x / K) ** n) ** 2

        return np.array([
            [-self.gamma, self.alpha_s * (1.0 - lam) * d_hill(R)],
            [self.alpha_r * (1.0 + lam) * d_hill(S), -self.gamma],
        ])

    # ------------------------------------------------------------ analysis
    def fixed_points(self, lam: float = 0.0, r_max: float | None = None, n_grid: int = 4000) -> List[FixedPoint]:
        """All steady states for a given bias, with linear stability.

        The 2-D problem is reduced to a 1-D root search: on the S-nullcline
        S = S*(R); the remaining condition g(R) = R_nullcline(S*(R)) - R = 0 is
        bracketed on a fine grid and refined by bisection.
        """
        if r_max is None:
            r_max = (self.alpha_r * (1.0 + abs(lam)) + self.beta_r) / self.gamma * 1.2 + 0.5
        Rg = np.linspace(0.0, r_max, n_grid)
        g = self.nullcline_R(self.nullcline_S(Rg, lam), lam) - Rg
        out: List[FixedPoint] = []
        for i in range(len(Rg) - 1):
            if g[i] == 0.0 or g[i] * g[i + 1] < 0.0:
                lo, hi = Rg[i], Rg[i + 1]
                glo = g[i]
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    gm = self.nullcline_R(self.nullcline_S(mid, lam), lam) - mid
                    if glo * gm <= 0.0:
                        hi = mid
                    else:
                        lo, glo = mid, gm
                R0 = 0.5 * (lo + hi)
                S0 = float(self.nullcline_S(R0, lam))
                ev = np.linalg.eigvals(self.jacobian((S0, R0), lam))
                stable = bool(np.all(ev.real < 0.0))
                out.append(FixedPoint(S0, float(R0), stable, (float(ev[0].real), float(ev[1].real))))
        return out

    def is_bistable(self, lam: float = 0.0) -> bool:
        return sum(fp.stable for fp in self.fixed_points(lam)) == 2

    def bifurcation(self, lams: Sequence[float]) -> List[Tuple[float, FixedPoint]]:
        """Fixed points for every bias value in ``lams`` (for a bifurcation diagram)."""
        return [(float(l), fp) for l in lams for fp in self.fixed_points(float(l))]

    def bistable_range(self, lams: Sequence[float] | None = None) -> Tuple[float, float] | None:
        """(lam_min, lam_max) of the bistable interval sampled on ``lams``; None if never bistable."""
        if lams is None:
            lams = np.linspace(-0.95, 0.95, 381)
        bis = [float(l) for l in lams if self.is_bistable(float(l))]
        if not bis:
            return None
        return min(bis), max(bis)

    # ---------------------------------------------------------- simulation
    def simulate(self, t_end: float, dt: float = 0.01, lam: LamLike = 0.0,
                 state0: Sequence[float] = (2.0, 0.06)) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Integrate the ODEs with classical RK4.

        ``lam`` may be a constant or a callable lam(t).  Returns (t, S, R).
        """
        lam_fn = lam if callable(lam) else (lambda _t, _l=float(lam): _l)
        n_steps = int(round(t_end / dt))
        t = np.linspace(0.0, n_steps * dt, n_steps + 1)
        S = np.empty(n_steps + 1); R = np.empty(n_steps + 1)
        y = np.array(state0, dtype=float)
        S[0], R[0] = y
        for i in range(n_steps):
            ti = t[i]
            k1 = self.rhs(y, lam_fn(ti))
            k2 = self.rhs(y + 0.5 * dt * k1, lam_fn(ti + 0.5 * dt))
            k3 = self.rhs(y + 0.5 * dt * k2, lam_fn(ti + 0.5 * dt))
            k4 = self.rhs(y + dt * k3, lam_fn(ti + dt))
            y = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
            y = np.maximum(y, 0.0)
            S[i + 1], R[i + 1] = y
        return t, S, R

    def pulse(self, lam_pulse: float = 0.8, t_on: float = 20.0, t_off: float = 40.0,
              t_end: float = 80.0, dt: float = 0.01, state0: Sequence[float] = (2.0, 0.06)):
        """A transient catabolic pulse (lam = lam_pulse for t_on < t < t_off); returns (t, S, R, lam_t)."""
        lam_fn = lambda t: lam_pulse if t_on < t < t_off else 0.0  # noqa: E731
        t, S, R = self.simulate(t_end, dt, lam_fn, state0)
        lam_t = np.array([lam_fn(ti) for ti in t])
        return t, S, R, lam_t

    def switches_permanently(self, **kwargs) -> bool:
        """True if a pulse (see :meth:`pulse`) leaves the cell RUNX2-high after the pulse ends."""
        t, S, R, _ = self.pulse(**kwargs)
        return bool(R[-1] > S[-1])
