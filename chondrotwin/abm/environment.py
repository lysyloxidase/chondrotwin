"""The tissue environment shared by all cell agents.

Geometry
    A slab Lx x Ly x h (mm) with the articular surface at z = 0 and the
    tidemark/bone at z = h, divided into a regular grid of voxels.

Mechanics (quasi-static, 1-D equilibrium per column)
    Applied compressive stress sigma(x, y, t) is uniform through the depth of a
    column; the depth-dependent aggregate modulus H_A(z) (soft superficial
    zone, stiff deep zone) and the local matrix content M of each voxel give
    the local tissue strain  eps = sigma / (H_A(z) * M^q).  Matrix loss
    therefore raises local strain - the mechanical half of the vicious cycle.

Cytokine field (quasi-steady)
    An IL-1beta-like signal C obeys  D lap(C) - k C + P = 0  with production P
    from senescent (SASP) and hypertrophic/catabolic cells, a Dirichlet value
    C_syn at the articular surface (synovial fluid) and zero flux elsewhere.
    Cytokine kinetics (hours) are fast compared with cell-state kinetics (days),
    so the field is re-solved to steady state at every daily step with a
    pre-factorized sparse operator.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


@dataclass
class CartilageSlab:
    Lx: float = 1.0            # mm
    Ly: float = 1.0            # mm
    h: float = 2.0             # mm (thickness)
    nx: int = 10
    ny: int = 10
    nz: int = 20
    H_A_surface: float = 0.4e6  # Pa, aggregate modulus at the surface
    H_A_deep: float = 1.2e6     # Pa, at the tidemark
    H_A_power: float = 1.0
    matrix_power: float = 1.5   # H_eff = H_A(z) * M^q
    matrix_floor: float = 0.3   # M below this no longer softens the tissue (collagen network still carries load)
    strain_max: float = 0.6     # physical cap on compressive strain
    D_cytokine: float = 0.5     # mm^2/day
    k_cytokine: float = 5.0     # 1/day
    zone_bounds: Tuple[float, float] = (0.15, 0.6)   # depth fractions superficial|middle|deep

    def __post_init__(self):
        self.dx, self.dy, self.dz = self.Lx / self.nx, self.Ly / self.ny, self.h / self.nz
        self.voxel_volume = self.dx * self.dy * self.dz                     # mm^3
        self.zc = (np.arange(self.nz) + 0.5) * self.dz
        self.M_grid = np.ones((self.nx, self.ny, self.nz))
        self.C_grid = np.zeros((self.nx, self.ny, self.nz))
        self._lu = None

    # ------------------------------------------------------------------ geometry
    def voxel_index(self, x, y, z) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        ix = np.clip((x / self.dx).astype(int), 0, self.nx - 1)
        iy = np.clip((y / self.dy).astype(int), 0, self.ny - 1)
        iz = np.clip((z / self.dz).astype(int), 0, self.nz - 1)
        return ix, iy, iz

    def zone_of(self, z) -> np.ndarray:
        f = np.asarray(z) / self.h
        return np.where(f < self.zone_bounds[0], 0, np.where(f < self.zone_bounds[1], 1, 2)).astype(np.int8)

    # ------------------------------------------------------------------ mechanics
    def H_A(self, z) -> np.ndarray:
        f = np.clip(np.asarray(z, float) / self.h, 0, 1)
        return self.H_A_surface + (self.H_A_deep - self.H_A_surface) * f ** self.H_A_power

    def strain_grid(self, sigma_xy: np.ndarray) -> np.ndarray:
        """Tissue compressive strain per voxel for an applied stress map sigma_xy (nx x ny, Pa)."""
        H = self.H_A(self.zc)[None, None, :] * np.clip(self.M_grid, self.matrix_floor, 1.0) ** self.matrix_power
        return np.minimum(sigma_xy[:, :, None] / H, self.strain_max)

    def update_matrix_grid(self, ix, iy, iz, M_cells, alive, k_empty: float = 0.02) -> None:
        """Voxel matrix = mean territorial matrix of live cells in the voxel; empty voxels decay slowly."""
        idx = (ix[alive], iy[alive], iz[alive])
        s = np.zeros_like(self.M_grid); n = np.zeros_like(self.M_grid)
        np.add.at(s, idx, M_cells[alive]); np.add.at(n, idx, 1.0)
        has = n > 0
        self.M_grid[has] = s[has] / n[has]
        self.M_grid[~has] *= (1.0 - k_empty)

    # ------------------------------------------------------------------ cytokine
    def _build_operator(self):
        n = self.nx * self.ny * self.nz

        def lap1d(m, d):
            L = sparse.diags([np.ones(m - 1), -2 * np.ones(m), np.ones(m - 1)], [-1, 0, 1], format="lil")
            L[0, 0] = -1; L[m - 1, m - 1] = -1            # zero-flux boundaries
            return (L / d ** 2).tocsr()

        Ix, Iy, Iz = (sparse.identity(k, format="csr") for k in (self.nx, self.ny, self.nz))
        L = sparse.kron(sparse.kron(lap1d(self.nx, self.dx), Iy), Iz) + \
            sparse.kron(sparse.kron(Ix, lap1d(self.ny, self.dy)), Iz) + \
            sparse.kron(sparse.kron(Ix, Iy), lap1d(self.nz, self.dz))
        A = (-self.D_cytokine * L + self.k_cytokine * sparse.identity(n)).tolil()
        # Dirichlet at the articular surface (iz = 0): row -> identity
        surf = np.arange(n).reshape(self.nx, self.ny, self.nz)[:, :, 0].ravel()
        for r in surf:
            A.rows[r] = [r]; A.data[r] = [1.0]
        self._surface_rows = surf
        self._lu = splu(A.tocsc())

    def solve_cytokine(self, ix, iy, iz, production, C_synovial: float = 0.0, reference_density: float = 1.4e4) -> None:
        """Quasi-steady cytokine field from per-cell production rates.

        Production is expressed relative to a reference: a voxel in which every
        cell (at ``reference_density`` cells/mm^3) produces at rate 1 reaches
        C = 1 at steady state without diffusion.  ``C_synovial`` = 1 therefore
        means "as inflamed as a fully senescent tissue".
        """
        if self._lu is None:
            self._build_operator()
        P = np.zeros((self.nx, self.ny, self.nz))
        np.add.at(P, (ix, iy, iz), production)
        rhs = (P / self.voxel_volume * self.k_cytokine / reference_density).ravel()
        rhs[self._surface_rows] = C_synovial
        self.C_grid = np.maximum(self._lu.solve(rhs), 0.0).reshape(self.nx, self.ny, self.nz)
