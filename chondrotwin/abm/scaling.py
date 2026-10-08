"""How far can "one agent = one cell" be pushed?

Order-of-magnitude anatomy (adult human knee; values vary between studies and
with age, site and depth):

* cell density of articular cartilage ~ 1e4 cells/mm^3 on average, with a
  strong depth gradient (roughly 2.4e4 superficial, 1.0e4 middle, 0.8e4 deep)
* total knee cartilage volume (femoral + tibial + patellar) ~ 2e4 mm^3
  => of the order of 1e8 chondrocytes in one knee.

With the structure-of-arrays representation used here (~60 bytes per agent)
1e8 agents need ~6 GB of RAM, which is feasible on a workstation; the time per
step scales linearly with the number of agents (see :func:`benchmark`).  In
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
BYTES_PER_AGENT = 60                  # float32/int8/bool fields of ChondrocyteAgents


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
                "memory_gb": memory_gb(n_agents)}


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
    lines = ["| tissue | agents | s / step (this machine) | hours / 365 steps | RAM |", "|---|---|---|---|---|"]
    for name, n in rows:
        e = bench.extrapolate(n)
        lines.append(f"| {name} | {n:.2e} | {e['seconds_per_step']:.3g} | {e['hours_for_steps']:.3g} | {e['memory_gb']:.3g} GB |")
    return "\n".join(lines)
