"""
A minimal reaction-diffusion model of the IHH-PTHrP feedback loop
(paper, Sections 2.3, 3.4, 5.3 and Figure 2).

On a one-dimensional axis x in [0, 1] from the periarticular end (x = 0) to the
chondro-osseous junction (x = 1):

    dP/dt = D_P d2P/dx2 + a_P s(x) [b + I^h / (K_I^h + I^h)] - k_P P
    dI/dt = D_I d2I/dx2 + a_I K_P^h / (K_P^h + P^h)             - k_I I

PTHrP (P) is produced near the periarticular end (Gaussian source s(x)) at a
rate enhanced by IHH (I); IHH is produced by cells in which PTHrP signaling is
low (prehypertrophic cells).  Hypertrophic commitment H(x) is taken as the same
decreasing Hill function of P.  The prehypertrophic boundary is where H crosses
0.5.  Zero-flux boundary conditions; explicit finite differences to steady state.

Dimensionless, illustrative parameters.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np


@dataclass(frozen=True)
class GrowthPlateSolution:
    x: np.ndarray
    P: np.ndarray
    I: np.ndarray
    H: np.ndarray
    converged: bool

    def boundary(self, threshold: float = 0.5) -> float:
        """Position of the prehypertrophic boundary (first x where H >= threshold)."""
        idx = np.argmax(self.H >= threshold)
        if self.H[idx] < threshold:  # never reaches threshold
            return float(self.x[-1])
        return float(self.x[idx])

    def proliferative_fraction(self, threshold: float = 0.5) -> float:
        """Fraction of the axis occupied by uncommitted (proliferative) cells."""
        return float(np.mean(self.H < threshold))


@dataclass
class GrowthPlateModel:
    D_P: float = 0.03
    D_I: float = 0.006
    a_P: float = 2.0
    k_P: float = 1.0
    a_I: float = 1.0
    k_I: float = 1.0
    K_P: float = 0.3
    K_I: float = 0.3
    basal_P: float = 0.25
    source_width: float = 0.08
    hill: float = 2.0
    n_grid: int = 201
    pthrp_production: bool = True

    def commitment(self, P: np.ndarray) -> np.ndarray:
        """Hypertrophic commitment as a decreasing Hill function of PTHrP."""
        h = self.hill
        return self.K_P ** h / (self.K_P ** h + np.maximum(P, 0.0) ** h)

    def steady_state(self, t_end: float = 40.0, tol: float = 1e-4) -> GrowthPlateSolution:
        """Integrate to steady state (stops early once max |dP/dt|, |dI/dt| < ``tol``).

        ``t_end`` is an upper bound; the slower Jansen-like variant (k_P x 0.35)
        needs t ~ 20-30 to meet ``tol``, so check ``converged`` on the result.
        """
        N = self.n_grid
        x = np.linspace(0.0, 1.0, N)
        dx = x[1] - x[0]
        dt = 0.2 * dx * dx / max(self.D_P, self.D_I)
        steps = int(t_end / dt)
        P = np.zeros(N); I = np.zeros(N)
        sP = np.exp(-x ** 2 / (2.0 * self.source_width ** 2))
        lapP = np.zeros(N); lapI = np.zeros(N)
        h = self.hill
        converged = False
        for step in range(steps):
            lapP[1:-1] = (P[2:] - 2 * P[1:-1] + P[:-2]) / dx ** 2
            lapI[1:-1] = (I[2:] - 2 * I[1:-1] + I[:-2]) / dx ** 2
            if self.pthrp_production:
                prodP = self.a_P * sP * (self.basal_P + I ** h / (self.K_I ** h + I ** h))
            else:
                prodP = 0.0
            prodI = self.a_I * self.K_P ** h / (self.K_P ** h + P ** h)
            dP = self.D_P * lapP + prodP - self.k_P * P
            dI = self.D_I * lapI + prodI - self.k_I * I
            P += dt * dP
            I += dt * dI
            P[0] = P[1]; P[-1] = P[-2]; I[0] = I[1]; I[-1] = I[-2]  # zero flux
            if step % 500 == 0 and step > 0:
                if max(np.abs(dP[1:-1]).max(), np.abs(dI[1:-1]).max()) < tol:
                    converged = True
                    break
        else:  # ran to t_end: check the final rate of change (interior nodes)
            converged = bool(max(np.abs(dP[1:-1]).max(), np.abs(dI[1:-1]).max()) < tol)
        H = self.commitment(P)
        return GrowthPlateSolution(x, P.copy(), I.copy(), H, converged)

    # ---------------------------------------------------------- scenarios
    def jansen_like(self, k_P_factor: float = 0.35) -> "GrowthPlateModel":
        """Longer effective PTHrP range, mimicking constitutive PTH1R signaling."""
        return replace(self, k_P=self.k_P * k_P_factor)

    def blomstrand_like(self) -> "GrowthPlateModel":
        """No PTHrP production/signaling, mimicking PTH1R or PTHrP loss."""
        return replace(self, pthrp_production=False)
