"""
A toy mechanistic "digital twin" loop (paper, Sections 5.6 and 7.1).

The loop closes the circle that existing models leave open: mechanics changes
the gene network, and the gene network changes the matrix that determines
mechanics.

    joint stress sigma0  --(biphasic tissue)-->  tissue strain = sigma0 / H_A
    tissue strain x PCM amplification          -->  cell strain
    cell strain --(loading.catabolic_bias)--> bias lam
    lam --(ToggleSwitch, one loading cycle)--> SOX9 (S), RUNX2 (R)
    S, R  --(matrix turnover)-->  H_A' = H_A (1 + k_syn S - k_deg R)

With the default parameters a stress of 0.05 MPa is physiological (cell strain
~0.15, i.e. m ~ 0.45, the anabolic optimum) and 0.15 MPa is injurious.
Under physiological stress the articular state maintains the matrix; under
overload RUNX2-driven degradation lowers H_A, which raises strain and feeds
back on the switch - an in-silico "vicious cycle".  Because the switch is
bistable, the outcome depends on history: a short overload episode can commit
the twin to progressive degradation even after the load returns to normal.

Illustrative parameters; cycles are abstract "loading epochs", not days.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Sequence, Union

import numpy as np

from .biphasic import ConfinedCompression
from .loading import catabolic_bias
from .toggle import ToggleSwitch

StressLike = Union[float, Sequence[float], Callable[[int], float]]


@dataclass(frozen=True)
class TwinRecord:
    cycle: int
    sigma0: float
    H_A: float
    tissue_strain: float
    cell_strain: float
    lam: float
    S: float
    R: float

    @property
    def fate(self) -> str:
        return "articular" if self.S > self.R else "hypertrophic"


@dataclass
class MiniTwin:
    toggle: ToggleSwitch = field(default_factory=ToggleSwitch)
    tissue: ConfinedCompression = field(default_factory=ConfinedCompression)
    pcm_amplification: float = 1.5     # cell strain / tissue strain (PCM strain transfer)
    strain_scale: float = 1.0 / 3.0    # cell strain that corresponds to m = 1 (injurious) in loading.py
    k_syn: float = 0.006               # matrix synthesis per unit SOX9 per cycle
    k_deg: float = 0.012               # matrix degradation per unit RUNX2 per cycle
    H_min: float = 0.05e6
    H_max: float = 0.80e6
    cycle_duration: float = 15.0       # toggle-switch time units per loading cycle
    dt: float = 0.02
    state: tuple = (2.0, 0.06)         # initial (S, R): articular

    def run(self, stress: StressLike, n_cycles: int = 60) -> List[TwinRecord]:
        """Run the loop for ``n_cycles``; ``stress`` is a constant (Pa), a sequence or a function of cycle."""
        if callable(stress):
            stress_fn = stress
        elif np.isscalar(stress):
            stress_fn = lambda c, s=float(stress): s  # noqa: E731
        else:
            seq = list(stress)
            stress_fn = lambda c, q=seq: q[min(c, len(q) - 1)]  # noqa: E731
        H_A = float(self.tissue.H_A)
        S, R = self.state
        records: List[TwinRecord] = []
        for c in range(n_cycles):
            sigma0 = float(stress_fn(c))
            tissue_strain = sigma0 / H_A
            cell_strain = self.pcm_amplification * tissue_strain
            lam = float(catabolic_bias(cell_strain / self.strain_scale))
            _, Ss, Rs = self.toggle.simulate(self.cycle_duration, self.dt, lam, (S, R))
            S, R = float(Ss[-1]), float(Rs[-1])
            S_mean, R_mean = float(Ss.mean()), float(Rs.mean())
            H_A = H_A * (1.0 + self.k_syn * S_mean - self.k_deg * R_mean)
            H_A = float(np.clip(H_A, self.H_min, self.H_max))
            records.append(TwinRecord(c, sigma0, H_A, tissue_strain, cell_strain, lam, S, R))
        return records

    @staticmethod
    def to_arrays(records: Sequence[TwinRecord]) -> dict:
        return {k: np.array([getattr(r, k) for r in records]) for k in
                ("cycle", "sigma0", "H_A", "tissue_strain", "cell_strain", "lam", "S", "R")}
