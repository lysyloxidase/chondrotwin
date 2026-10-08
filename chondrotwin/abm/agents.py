"""Chondrocyte agents stored as a structure of arrays.

One agent is one chondrocyte.  All per-agent variables are NumPy arrays of a
fixed capacity so that every update rule is vectorized; ``alive`` marks the
slots in use.  Dead cells keep their slot (flagged ``dead``) so that cell loss
can be reported; proliferation fills free slots.

Per-agent state
---------------
x, y, z        position (mm); z = 0 at the articular surface, z = h at the bone
zone           0 superficial, 1 middle, 2 deep
S, R           SOX9 / RUNX2 activities of the toggle switch (continuous)
M              territorial matrix content (0-1): aggrecan/collagen around the cell
nrf2           antioxidant capacity (0-1), heterogeneous between cells
age            days since start (or since division)
senescent      SASP-producing, non-dividing, non-synthesizing
dead           removed from the tissue; ``death_cause`` 0 none, 1 ferroptosis, 2 apoptosis
cluster        id of the proliferative cluster the cell belongs to (0 = none)
damage         accumulated overload exposure (overload-days; drives the cumulative senescence hazard)
bias           cell-specific offset of the toggle-switch bias (stage-3 heterogeneity; 0 by default)
competent      can become senescent (stage-3 heterogeneity; True by default)
strain, lam, C diagnostic: last cell strain, toggle bias and local cytokine level
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict

import numpy as np

FATE_LABELS = ("articular", "hypertrophic", "senescent", "dead_ferroptosis", "dead_apoptosis")


@dataclass
class ChondrocyteAgents:
    capacity: int
    x: np.ndarray = field(init=False)
    y: np.ndarray = field(init=False)
    z: np.ndarray = field(init=False)
    zone: np.ndarray = field(init=False)
    S: np.ndarray = field(init=False)
    R: np.ndarray = field(init=False)
    M: np.ndarray = field(init=False)
    nrf2: np.ndarray = field(init=False)
    age: np.ndarray = field(init=False)
    senescent: np.ndarray = field(init=False)
    dead: np.ndarray = field(init=False)
    death_cause: np.ndarray = field(init=False)
    cluster: np.ndarray = field(init=False)
    damage: np.ndarray = field(init=False)
    bias: np.ndarray = field(init=False)
    competent: np.ndarray = field(init=False)
    alive: np.ndarray = field(init=False)
    strain: np.ndarray = field(init=False)
    lam: np.ndarray = field(init=False)
    C: np.ndarray = field(init=False)
    n_initial: int = 0

    def __post_init__(self):
        n = self.capacity
        f = lambda: np.zeros(n, dtype=np.float32)  # noqa: E731
        self.x, self.y, self.z = f(), f(), f()
        self.zone = np.zeros(n, dtype=np.int8)
        self.S, self.R, self.M, self.nrf2, self.age = f(), f(), f(), f(), f()
        self.senescent = np.zeros(n, dtype=bool)
        self.dead = np.zeros(n, dtype=bool)
        self.death_cause = np.zeros(n, dtype=np.int8)
        self.cluster = np.zeros(n, dtype=np.int32)
        self.damage = f()
        self.bias = f()
        self.competent = np.ones(n, dtype=bool)
        self.alive = np.zeros(n, dtype=bool)       # slot in use and not dead
        self.strain, self.lam, self.C = f(), f(), f()

    # ------------------------------------------------------------------ helpers
    @property
    def n_alive(self) -> int:
        return int(self.alive.sum())

    @property
    def n_used(self) -> int:
        return int((self.alive | self.dead).sum())

    def free_slots(self, k: int) -> np.ndarray:
        """Indices of up to ``k`` unused slots (never used, i.e. neither alive nor dead)."""
        free = np.where(~(self.alive | self.dead))[0]
        return free[:k]

    def fate(self) -> np.ndarray:
        """Categorical fate per slot: 0 articular, 1 hypertrophic, 2 senescent, 3 dead-ferroptosis, 4 dead-apoptosis, -1 unused."""
        out = np.full(self.capacity, -1, dtype=np.int8)
        used = self.alive | self.dead
        out[used & self.alive & ~self.senescent & (self.S >= self.R)] = 0
        out[used & self.alive & ~self.senescent & (self.R > self.S)] = 1
        out[used & self.alive & self.senescent] = 2
        out[used & self.dead & (self.death_cause == 1)] = 3
        out[used & self.dead & (self.death_cause == 2)] = 4
        return out

    def fractions(self) -> Dict[str, float]:
        """Fate fractions relative to the initial number of cells (so cell loss is visible)."""
        f = self.fate()
        n0 = max(self.n_initial, 1)
        return {lab: float(np.sum(f == i)) / n0 for i, lab in enumerate(FATE_LABELS)}

    def memory_bytes(self) -> int:
        return int(sum(getattr(self, a).nbytes for a in ("x", "y", "z", "zone", "S", "R", "M", "nrf2", "age", "senescent",
                                                          "dead", "death_cause", "cluster", "damage", "bias", "competent", "alive",
                                                          "strain", "lam", "C")))
