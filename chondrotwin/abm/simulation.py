"""The multi-agent cartilage model: one agent = one chondrocyte.

Scheduler (one step = one day)
    1. mechanics      applied stress map -> tissue strain per voxel (matrix feedback) -> cell strain
    2. signaling      cytokine field (quasi-steady) from SASP / catabolic cells
    3. fate dynamics  bias lam = catabolic_bias(cell strain) + inflammatory term;
                      stochastic toggle switch (SOX9/RUNX2) integrated for one day
    4. matrix         territorial matrix synthesis (SOX9) and degradation (RUNX2, cytokines)
    5. hazards        senescence, ferroptosis (overload x low Nrf2), apoptosis (overload x inflammation)
    6. proliferation  cells in degraded matrix divide and form clusters
    7. recording      fate fractions, matrix and strain by zone, cytokine, OARSI-like score

All rules are vectorized NumPy over the agent arrays; the environment grid
carries matrix content and the cytokine field.  Parameters are illustrative.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Union

import numpy as np

from ..loading import catabolic_bias
from ..toggle import ToggleSwitch
from .agents import FATE_LABELS, ChondrocyteAgents
from .environment import CartilageSlab

StressLike = Union[float, Callable[[int], float]]


@dataclass
class Scenario:
    name: str = "physiological"
    sigma0: StressLike = 0.05e6            # applied stress (Pa) or function of day; 0.05 MPa physiological, 0.12 MPa injurious
    focal_factor: float = 0.0              # extra stress in a focal region: sigma * (1 + f * gaussian)
    focal_width: float = 0.2               # mm
    C_synovial: float = 0.0                # cytokine level at the articular surface (inflamed joint)
    initial_senescent_fraction: float = 0.0
    initial_hypertrophic_fraction: float = 0.0
    days: int = 365

    def stress(self, day: int) -> float:
        return float(self.sigma0(day)) if callable(self.sigma0) else float(self.sigma0)


@dataclass
class Rates:
    """Illustrative daily rate constants."""
    pcm_amplification: float = 1.5
    strain_scale: float = 0.40             # cell strain corresponding to m = 1 (injurious) in loading.catabolic_bias
    inflammatory_gain: float = 0.6         # lam += gain * C / (K_C + C)
    K_C: float = 0.5
    toggle_noise: float = 0.10
    toggle_dt: float = 0.25
    toggle_substeps: int = 4
    k_syn: float = 0.008
    k_deg: float = 0.008
    cytokine_deg_gain: float = 2.0         # degradation x (1 + gain * C/(K_C+C))
    h_sen_base: float = 3e-5               # ageing hazard per day
    h_sen_overload: float = 0.012
    h_ferro: float = 0.006
    h_apop: float = 0.003
    overload_threshold: float = 0.8        # in m units
    overload_width: float = 0.08
    sasp_production: float = 1.0           # cytokine units/day per senescent cell
    catabolic_production: float = 0.25     # per hypertrophic cell (R > S)
    p_div: float = 0.003                   # division probability per day in degraded matrix
    M_div: float = 0.6
    cluster_max: int = 8
    k_empty: float = 0.02


@dataclass
class CartilageABM:
    slab: CartilageSlab = field(default_factory=CartilageSlab)
    cells_per_mm3: float = 1.4e4
    density_profile: tuple = (2.4e4, 1.0e4, 0.77e4)   # superficial, middle, deep (cells/mm^3), scaled to cells_per_mm3
    toggle: ToggleSwitch = field(default_factory=ToggleSwitch)
    rates: Rates = field(default_factory=Rates)
    max_growth: float = 1.5                # capacity = initial cells x max_growth
    seed: int = 0

    def __post_init__(self):
        self.rng = np.random.default_rng(self.seed)
        self.agents = self._seed_cells()
        self.history: List[Dict[str, float]] = []
        self._next_cluster = 1

    # ------------------------------------------------------------------ initialization
    def _seed_cells(self) -> ChondrocyteAgents:
        s = self.slab
        area = s.Lx * s.Ly
        bounds = (0.0, s.zone_bounds[0] * s.h, s.zone_bounds[1] * s.h, s.h)
        prof = np.array(self.density_profile, float)
        thick = np.diff(bounds)
        mean_density = np.sum(prof * thick) / s.h
        prof = prof * (self.cells_per_mm3 / mean_density)
        zs = []
        for k in range(3):
            n_k = int(round(prof[k] * area * thick[k]))
            zs.append(self.rng.uniform(bounds[k], bounds[k + 1], n_k))
        z = np.concatenate(zs)
        n = len(z)
        ag = ChondrocyteAgents(capacity=int(np.ceil(n * self.max_growth)))
        ag.n_initial = n
        idx = np.arange(n)
        ag.x[idx] = self.rng.uniform(0, s.Lx, n); ag.y[idx] = self.rng.uniform(0, s.Ly, n); ag.z[idx] = z
        ag.zone[idx] = s.zone_of(z)
        ag.S[idx] = 2.0 + 0.1 * self.rng.standard_normal(n); ag.R[idx] = 0.06 + 0.02 * self.rng.standard_normal(n)
        ag.M[idx] = np.clip(0.9 + 0.05 * self.rng.standard_normal(n), 0.5, 1.0)
        ag.nrf2[idx] = np.clip(self.rng.beta(4, 2, n), 0, 1)
        ag.alive[idx] = True
        return ag

    def seed_fates(self, senescent_fraction: float = 0.0, hypertrophic_fraction: float = 0.0,
                   fractions: Optional[Dict[str, float]] = None) -> None:
        """Set initial fate occupancy, e.g. from single-cell basin occupancy (chondrofate experiments)."""
        if fractions:
            senescent_fraction = fractions.get("senescent", senescent_fraction)
            hypertrophic_fraction = fractions.get("hypertrophic", hypertrophic_fraction)
        ag = self.agents
        idx = np.where(ag.alive)[0]
        pick = self.rng.permutation(idx)
        n_s = int(senescent_fraction * len(idx)); n_h = int(hypertrophic_fraction * len(idx))
        ag.senescent[pick[:n_s]] = True
        hyp = pick[n_s:n_s + n_h]
        ag.S[hyp], ag.R[hyp] = 0.06, 2.0

    # ------------------------------------------------------------------ one step
    def step(self, day: int, scenario: Scenario) -> Dict[str, float]:
        ag, s, r = self.agents, self.slab, self.rates
        alive = ag.alive
        ix, iy, iz = s.voxel_index(ag.x, ag.y, ag.z)
        # 1. mechanics -----------------------------------------------------------
        sigma = scenario.stress(day)
        xs = (np.arange(s.nx) + 0.5) * s.dx; ys = (np.arange(s.ny) + 0.5) * s.dy
        X, Y = np.meshgrid(xs, ys, indexing="ij")
        focal = np.exp(-((X - s.Lx / 2) ** 2 + (Y - s.Ly / 2) ** 2) / (2 * scenario.focal_width ** 2))
        sigma_xy = sigma * (1.0 + scenario.focal_factor * focal)
        s.update_matrix_grid(ix, iy, iz, ag.M, alive, r.k_empty)
        eps = s.strain_grid(sigma_xy)
        cell_strain = r.pcm_amplification * eps[ix, iy, iz]
        ag.strain[alive] = cell_strain[alive]
        m = cell_strain / r.strain_scale
        overload = 1.0 / (1.0 + np.exp(-(m - r.overload_threshold) / r.overload_width))
        # 2. cytokine field --------------------------------------------------------
        prod = np.zeros(ag.capacity, dtype=float)
        prod[alive & ag.senescent] = r.sasp_production
        prod[alive & ~ag.senescent & (ag.R > ag.S)] = r.catabolic_production
        s.solve_cytokine(ix[alive], iy[alive], iz[alive], prod[alive], scenario.C_synovial, reference_density=self.cells_per_mm3)
        C = s.C_grid[ix, iy, iz]
        ag.C[alive] = C[alive]
        infl = C / (r.K_C + C)
        # 3. fate dynamics ---------------------------------------------------------
        lam = np.clip(catabolic_bias(m) + r.inflammatory_gain * infl, -0.95, 0.95)
        ag.lam[alive] = lam[alive]
        dyn = alive & ~ag.senescent
        S, R = ag.S[dyn].astype(float), ag.R[dyn].astype(float); lam_d = lam[dyn]
        sq = r.toggle_noise * np.sqrt(r.toggle_dt)
        for _ in range(r.toggle_substeps):
            dS, dR = self.toggle.rhs_arrays(S, R, lam_d)
            S = np.maximum(S + r.toggle_dt * dS + sq * self.rng.standard_normal(len(S)), 0.0)
            R = np.maximum(R + r.toggle_dt * dR + sq * self.rng.standard_normal(len(R)), 0.0)
        ag.S[dyn], ag.R[dyn] = S, R
        # 4. matrix turnover -------------------------------------------------------
        syn = np.where(ag.senescent, 0.2, 1.0) * r.k_syn * ag.S * (1.0 - ag.M)
        deg = r.k_deg * ag.R * ag.M * (1.0 + r.cytokine_deg_gain * infl)
        ag.M[alive] = np.clip(ag.M[alive] + syn[alive] - deg[alive], 0.05, 1.0)
        # 5. hazards ---------------------------------------------------------------
        u = self.rng.uniform(size=ag.capacity)
        h_sen = r.h_sen_base + r.h_sen_overload * overload * (1.0 + 2.0 * infl)
        new_sen = alive & ~ag.senescent & (u < 1.0 - np.exp(-h_sen))
        ag.senescent[new_sen] = True
        u = self.rng.uniform(size=ag.capacity)
        h_ferro = r.h_ferro * overload * (1.0 - ag.nrf2)
        ferro = alive & (u < 1.0 - np.exp(-h_ferro))
        u = self.rng.uniform(size=ag.capacity)
        h_apop = r.h_apop * overload * (0.2 + infl)
        apop = alive & ~ferro & (u < 1.0 - np.exp(-h_apop))
        for mask, cause in ((ferro, 1), (apop, 2)):
            ag.dead[mask] = True; ag.alive[mask] = False; ag.death_cause[mask] = cause
        # 6. proliferation / clusters ----------------------------------------------
        cand = ag.alive & ~ag.senescent & (ag.M < r.M_div)
        u = self.rng.uniform(size=ag.capacity)
        div = np.where(cand & (u < r.p_div))[0]
        if len(div):
            free = ag.free_slots(len(div))
            div = div[:len(free)]
            for src, dst in zip(div, free):
                if ag.cluster[src] == 0:
                    ag.cluster[src] = self._next_cluster; self._next_cluster += 1
                members = np.sum(ag.cluster == ag.cluster[src])
                if members >= r.cluster_max:
                    continue
                ag.x[dst] = np.clip(ag.x[src] + self.rng.normal(0, 0.02), 0, s.Lx)
                ag.y[dst] = np.clip(ag.y[src] + self.rng.normal(0, 0.02), 0, s.Ly)
                ag.z[dst] = np.clip(ag.z[src] + self.rng.normal(0, 0.02), 0, s.h)
                ag.zone[dst] = s.zone_of(ag.z[dst])
                ag.S[dst], ag.R[dst], ag.M[dst], ag.nrf2[dst] = ag.S[src], ag.R[src], ag.M[src], ag.nrf2[src]
                ag.cluster[dst] = ag.cluster[src]; ag.alive[dst] = True; ag.age[dst] = 0.0
        ag.age[ag.alive] += 1.0
        # 7. record ----------------------------------------------------------------
        rec = {"day": day, "sigma_MPa": sigma / 1e6, **ag.fractions()}
        rec["n_alive"] = ag.n_alive
        for k, lab in enumerate(("superficial", "middle", "deep")):
            zm = ag.alive & (ag.zone == k)
            rec[f"matrix_{lab}"] = float(ag.M[zm].mean()) if zm.any() else float("nan")
            rec[f"strain_{lab}"] = float(ag.strain[zm].mean()) if zm.any() else float("nan")
        rec["cytokine_mean"] = float(s.C_grid.mean())
        rec["n_clusters"] = int(len(np.unique(ag.cluster[ag.alive & (ag.cluster > 0)])))
        rec["oarsi_like"] = self.oarsi_like()
        self.history.append(rec)
        return rec

    def oarsi_like(self) -> float:
        """0-6 degeneration score: surface matrix loss (0-3) + cell loss and clusters (0-3)."""
        ag = self.agents
        sup = ag.alive & (ag.zone == 0)
        m_sup = float(ag.M[sup].mean()) if sup.any() else 0.0
        cell_loss = 1.0 - ag.n_alive / max(ag.n_initial, 1)
        clusters = float(np.mean(ag.cluster[ag.alive] > 0)) if ag.n_alive else 1.0
        return float(3 * (1 - m_sup) + 1.5 * min(cell_loss / 0.3, 1) + 1.5 * min(clusters / 0.3, 1))

    def run(self, scenario: Scenario, record_every: int = 1, progress: bool = False) -> "SimulationResult":
        self.seed_fates(scenario.initial_senescent_fraction, scenario.initial_hypertrophic_fraction)
        for day in range(scenario.days):
            rec = self.step(day, scenario)
            if progress and day % 30 == 0:
                print(f"day {day:4d}  alive {rec['n_alive']:6d}  articular {rec['articular']:.2f}  hypertrophic {rec['hypertrophic']:.2f}  "
                      f"senescent {rec['senescent']:.2f}  M_sup {rec['matrix_superficial']:.2f}  OARSI-like {rec['oarsi_like']:.2f}")
        return SimulationResult(scenario.name, self.history, self.snapshot(), self.slab, self.agents.n_initial)

    def snapshot(self):
        import pandas as pd
        ag = self.agents
        used = ag.alive | ag.dead
        f = ag.fate()
        return pd.DataFrame({"x": ag.x[used], "y": ag.y[used], "z": ag.z[used], "zone": ag.zone[used], "S": ag.S[used],
                             "R": ag.R[used], "M": ag.M[used], "strain": ag.strain[used], "lam": ag.lam[used], "C": ag.C[used],
                             "senescent": ag.senescent[used], "dead": ag.dead[used], "cluster": ag.cluster[used],
                             "fate": [FATE_LABELS[k] if k >= 0 else "unused" for k in f[used]]})


@dataclass
class SimulationResult:
    name: str
    history: List[Dict[str, float]]
    cells: object            # pandas DataFrame snapshot at the end
    slab: CartilageSlab
    n_initial: int

    def frame(self):
        import pandas as pd
        return pd.DataFrame(self.history).set_index("day")

    def final(self) -> Dict[str, float]:
        return self.history[-1]
