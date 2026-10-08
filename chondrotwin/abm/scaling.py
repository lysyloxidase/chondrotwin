"""How far can "one agent = one cell" be pushed?

Order-of-magnitude anatomy (adult human knee; values vary between studies and
with age, site and depth):

* cell density of adult human articular cartilage ~ 1e4 cells/mm^3 on average
  (9.6e3 mm^-3 full-depth in the medial femoral condyle, Hunziker et al. 2002;
  the model default 1.4e4 is at the upper end), with a strong depth gradient
  (superficial zone roughly 2-3x the mean)
* total knee cartilage volume (femoral + tibial + patellar) ~ 2e4 mm^3
  => of the order of 2-3e8 chondrocytes in one knee.

Memory (measured, see VALIDATION.md): the structure-of-arrays state takes
62 bytes per slot and slots are allocated for 1.5x the initial cell number
(~93 bytes per agent); the vectorized daily step allocates float64/int64
temporaries of the full capacity, ~260 bytes per agent at peak.  A whole knee
(2.8e8 agents) in one process therefore needs ~100 GB at peak - a large
workstation or a chunked / GPU implementation.  The time per step scales
linearly with the number of agents for a fixed voxel resolution (see
:func:`benchmark`); the extrapolation does not include the cost of the sparse
cytokine solve on a joint-sized voxel grid.  In
practice a *representative volume* (e.g. 1 x 1 x 2 mm, ~3e4 cells) is simulated
under the strain field of a tissue/joint-scale model, and whole-cartilage
behaviour is obtained by tiling such volumes over the joint surface with
different local loads (the multiscale coupling of the paper, Figure 7).
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Dict

import numpy as np

CELLS_PER_MM3_MEAN = 1.4e4
KNEE_CARTILAGE_VOLUME_MM3 = 2.0e4     # femoral + tibial + patellar, order of magnitude
BYTES_PER_AGENT = 93                  # persistent state: 62 B per slot x 1.5 capacity (ChondrocyteAgents)
PEAK_BYTES_PER_AGENT = 360            # state + ~260 B of temporaries during one vectorized step (measured)


def cells_for_volume(volume_mm3: float, density: float = CELLS_PER_MM3_MEAN) -> int:
    return int(volume_mm3 * density)


def memory_gb(n_agents: int, bytes_per_agent: int = BYTES_PER_AGENT) -> float:
    return n_agents * bytes_per_agent / 1e9


def tile_plan(volume_mm3: float = KNEE_CARTILAGE_VOLUME_MM3, rve_mm3: float = 2.0) -> Dict[str, float]:
    """Number of representative volumes (RVEs) needed to tile a tissue volume."""
    n_rve = int(np.ceil(volume_mm3 / rve_mm3))
    return {"tissue_volume_mm3": volume_mm3, "rve_mm3": rve_mm3, "n_rve": n_rve,
            "cells_per_rve": cells_for_volume(rve_mm3), "total_cells": cells_for_volume(volume_mm3)}


@dataclass
class Benchmark:
    n_agents: int
    seconds_per_step: float

    def extrapolate(self, n_agents: int, n_steps: int = 365) -> Dict[str, float]:
        per_step = self.seconds_per_step * n_agents / self.n_agents
        return {"n_agents": n_agents, "seconds_per_step": per_step, "hours_for_steps": per_step * n_steps / 3600,
                "memory_gb": memory_gb(n_agents), "peak_memory_gb": memory_gb(n_agents, PEAK_BYTES_PER_AGENT)}


def benchmark(n_steps: int = 5, Lx: float = 1.0, Ly: float = 1.0, h: float = 2.0, cells_per_mm3: float = CELLS_PER_MM3_MEAN) -> Benchmark:
    """Time the simulation of a slab for a few steps on this machine."""
    from .environment import CartilageSlab
    from .simulation import CartilageABM, Scenario
    sim = CartilageABM(slab=CartilageSlab(Lx=Lx, Ly=Ly, h=h), cells_per_mm3=cells_per_mm3)
    sc = Scenario(days=n_steps)
    t0 = time.perf_counter()
    sim.run(sc)
    dt = (time.perf_counter() - t0) / n_steps
    return Benchmark(sim.agents.n_initial, dt)


def scaling_table(bench: Benchmark) -> str:
    rows = [("1 mm^3 (RVE)", cells_for_volume(1.0)), ("1 x 1 x 2 mm slab", cells_for_volume(2.0)),
            ("1 cm^2 x 2 mm patch", cells_for_volume(200.0)), ("medial tibial plateau (~2.5e3 mm^3)", cells_for_volume(2.5e3)),
            ("whole knee (~2e4 mm^3)", cells_for_volume(KNEE_CARTILAGE_VOLUME_MM3))]
    lines = ["| tissue | agents | s / step (this machine) | hours / 365 steps | RAM state / peak |", "|---|---|---|---|---|"]
    for name, n in rows:
        e = bench.extrapolate(n)
        lines.append(f"| {name} | {n:.2e} | {e['seconds_per_step']:.3g} | {e['hours_for_steps']:.3g} | "
                     f"{e['memory_gb']:.3g} / {e['peak_memory_gb']:.3g} GB |")
    return "\n".join(lines)
