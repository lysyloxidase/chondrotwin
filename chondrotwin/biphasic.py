"""
A one-dimensional biphasic (poroelastic) cartilage model (paper, Section 5.5).

Confined compression of a cartilage layer of thickness h resting on an
impermeable base (z = 0) and loaded through a porous, free-draining platen
(z = h) with a constant stress sigma0 (the classical biphasic creep problem of
Mow et al., 1980).  With aggregate modulus H_A and permeability k, the solid
displacement u(z, t) obeys the consolidation equation

    du/dt = H_A k  d2u/dz2,     u(0, t) = 0,    H_A du/dz (h, t) = -sigma0

and the interstitial fluid pressure follows from 1-D equilibrium of the mixture,

    p(z, t) = sigma0 + H_A du/dz.

At t = 0+ the fluid carries the whole load (p = sigma0); as fluid exudes the
solid takes over and the equilibrium strain -sigma0/H_A is reached with the
characteristic time tau = h^2 / (H_A k).  Depth-resolved strain and pressure
are the cell-level mechanical stimuli used by :class:`chondrotwin.twin.MiniTwin`.

Explicit finite differences; SI units (Pa, m, s).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class BiphasicSolution:
    t: np.ndarray           # output times (s)
    z: np.ndarray           # depth coordinate (m), 0 = bone side, h = surface
    u: np.ndarray           # displacement (m), shape (n_t, n_z), negative = compression
    strain: np.ndarray      # du/dz, shape (n_t, n_z)
    pressure: np.ndarray    # fluid pressure (Pa), shape (n_t, n_z)
    sigma0: float
    tau: float

    @property
    def surface_displacement(self) -> np.ndarray:
        return self.u[:, -1]


@dataclass
class ConfinedCompression:
    H_A: float = 0.5e6      # aggregate modulus (Pa)
    k: float = 1.0e-15      # hydraulic permeability (m^4 / (N s))
    h: float = 1.5e-3       # thickness (m)
    n_z: int = 51

    @property
    def tau(self) -> float:
        """Characteristic consolidation time h^2 / (H_A k) in seconds."""
        return self.h ** 2 / (self.H_A * self.k)

    def equilibrium_strain(self, sigma0: float) -> float:
        return -sigma0 / self.H_A

    def solve(self, sigma0: float, t_end: float | None = None, n_out: int = 60,
              safety: float = 0.4) -> BiphasicSolution:
        """Creep response to a step stress ``sigma0`` (Pa, positive = compressive)."""
        if t_end is None:
            t_end = 5.0 * self.tau
        z = np.linspace(0.0, self.h, self.n_z)
        dz = z[1] - z[0]
        c = self.H_A * self.k
        dt = safety * dz * dz / c
        n_steps = int(np.ceil(t_end / dt))
        dt = t_end / n_steps
        out_idx = np.unique(np.round(np.linspace(0, n_steps, n_out)).astype(int))
        u = np.zeros(self.n_z)
        U = np.zeros((len(out_idx), self.n_z)); T = np.zeros(len(out_idx))
        j = 0
        if out_idx[0] == 0:
            U[0] = u; T[0] = 0.0; j = 1
        for step in range(1, n_steps + 1):
            lap = np.zeros_like(u)
            lap[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / dz ** 2
            u = u + dt * c * lap
            u[0] = 0.0                                     # fixed impermeable base
            u[-1] = u[-2] - sigma0 * dz / self.H_A         # traction boundary condition
            if j < len(out_idx) and step == out_idx[j]:
                U[j] = u; T[j] = step * dt; j += 1
        strain = np.gradient(U, dz, axis=1)
        pressure = sigma0 + self.H_A * strain
        pressure[:, -1] = 0.0                               # free-draining surface
        return BiphasicSolution(T, z, U, strain, np.maximum(pressure, 0.0), sigma0, self.tau)

    def solve_history(self, sigma_fn, t_end: float, dt_max: float = 0.02, n_out: int = 200,
                      safety: float = 0.4) -> BiphasicSolution:
        """Response to a time-varying applied stress ``sigma_fn(t)`` (Pa), e.g. a cyclic gait waveform.

        The time step is the smaller of the stability limit and ``dt_max`` (to resolve the loading).
        Returned ``pressure`` uses the stress at each output time; ``sigma0`` is the time-averaged stress.
        """
        z = np.linspace(0.0, self.h, self.n_z)
        dz = z[1] - z[0]
        c = self.H_A * self.k
        n_steps = int(np.ceil(t_end / min(safety * dz * dz / c, dt_max)))
        dt = t_end / n_steps
        out_idx = set(np.round(np.linspace(0, n_steps, n_out)).astype(int).tolist())
        u = np.zeros(self.n_z)
        U, T, S = [], [], []
        if 0 in out_idx:
            U.append(u.copy()); T.append(0.0); S.append(float(sigma_fn(0.0)))
        acc = 0.0
        for step in range(1, n_steps + 1):
            t = step * dt
            s = float(sigma_fn(t))
            acc += s
            u[1:-1] = u[1:-1] + dt * c * (u[2:] - 2 * u[1:-1] + u[:-2]) / dz ** 2
            u[0] = 0.0
            u[-1] = u[-2] - s * dz / self.H_A
            if step in out_idx:
                U.append(u.copy()); T.append(t); S.append(s)
        U = np.array(U); strain = np.gradient(U, dz, axis=1)
        pressure = np.array(S)[:, None] + self.H_A * strain
        pressure[:, -1] = 0.0
        return BiphasicSolution(np.array(T), z, U, strain, np.maximum(pressure, 0.0), acc / n_steps, self.tau)

    def strain_at(self, sigma0: float, depth_fraction: float = 0.5, t: float | None = None) -> float:
        """Compressive strain magnitude at a given depth fraction (0 = bone, 1 = surface)."""
        if t is None:
            return abs(self.equilibrium_strain(sigma0))
        sol = self.solve(sigma0, t_end=t, n_out=2)
        iz = int(round(depth_fraction * (self.n_z - 1)))
        return float(abs(sol.strain[-1, iz]))
