"""
A small logical (Boolean) model of chondrocyte fate (paper, Sections 5.1 and 7).

The engine is generic: a :class:`BooleanNetwork` is a list of *internal* nodes
with update rules, plus *input* nodes that are clamped to a scenario.  With the
inputs fixed, all 2^N internal states are enumerated, the synchronous successor
of each state is computed, and the attractors (fixed points and cycles) are
collected together with their basin sizes.  For the illustrative chondrocyte
network (10 internal nodes, 1024 states) this takes well under a second.

The chondrocyte network encodes, in the simplest possible form, the wiring of
Figure 3 of the paper:

* mechanosensing: Piezo1 is activated by overload, or by physiological load in
  the presence of IL-1beta (inflammatory sensitization); TRPV4 by physiological
  load only;
* Ca2+ overload follows Piezo1; a stiff matrix drives nuclear YAP unless
  cytokines degrade it; NF-kB is driven by IL-1beta, Ca2+ overload or a
  senescent (SASP) state and is inhibited by YAP; HIF-2alpha follows NF-kB;
* the SOX9-RUNX2 toggle: SOX9 requires the absence of RUNX2, NF-kB, nuclear YAP
  and ferroptosis, and TGF-beta can override WNT; RUNX2 requires the absence of
  SOX9, PTHrP and ferroptosis, and a driver (WNT, HIF-2alpha, Ca2+ overload) or
  its own autoregulation;
* terminal exits: ferroptosis follows Ca2+ overload when Nrf2 antioxidant
  capacity is absent, senescence follows NF-kB plus Ca2+ overload; both are
  absorbing.

The attractors of this network under different inputs are the fate states of
Table 1 (articular, hypertrophic, dedifferentiated, senescent, ferroptotic).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import Callable, Dict, List, Mapping, Sequence, Tuple

State = Dict[str, int]
Rule = Callable[[State], bool]


@dataclass(frozen=True)
class Attractor:
    states: Tuple[Tuple[int, ...], ...]   # ordered states of the cycle (length 1 = fixed point)
    basin_size: int
    nodes: Tuple[str, ...]

    @property
    def is_fixed_point(self) -> bool:
        return len(self.states) == 1

    def as_dicts(self) -> List[State]:
        return [dict(zip(self.nodes, s)) for s in self.states]

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        kind = "fixed point" if self.is_fixed_point else f"cycle (period {len(self.states)})"
        active = ", ".join(n for n, v in zip(self.nodes, self.states[0]) if v) or "none"
        return f"{kind}; basin {self.basin_size}; active: {active}"


@dataclass
class BooleanNetwork:
    internal: Sequence[str]
    inputs: Sequence[str]
    rules: Mapping[str, Rule]
    default_inputs: Mapping[str, int] = field(default_factory=dict)
    absorbing: Sequence[str] = ()   # terminal flags (e.g. Senescent, Ferroptotic): once on, never off
    slow: Sequence[str] = ()        # fate-memory nodes allowed to vary in the default initial-condition set

    def __post_init__(self):
        missing = [n for n in self.internal if n not in self.rules]
        if missing:
            raise ValueError(f"rules missing for internal nodes: {missing}")

    # ---------------------------------------------------------------- core
    def step(self, state: State) -> State:
        """Synchronous update of the internal nodes (inputs are carried over)."""
        new = dict(state)
        for n in self.internal:
            new[n] = int(bool(self.rules[n](state)))
        return new

    def _inputs(self, scenario: Mapping[str, int] | None) -> Dict[str, int]:
        inp = {k: 0 for k in self.inputs}
        inp.update(self.default_inputs)
        if scenario:
            unknown = set(scenario) - set(self.inputs)
            if unknown:
                raise ValueError(f"unknown input nodes: {sorted(unknown)}")
            inp.update(scenario)
        return inp

    def attractors(self, scenario: Mapping[str, int] | None = None, initial: str = "slow") -> List[Attractor]:
        """Enumerate synchronous attractors for clamped inputs, with basin sizes.

        ``initial`` selects the set of initial states over which basins are counted:

        * ``"slow"`` (default): only the fate-memory nodes in ``self.slow`` vary
          (2^len(slow) initial states); all fast signaling nodes and absorbing
          flags start at 0.  This asks "where does a cell with a given prior
          fate go under these inputs?" and is the biologically meaningful count.
        * ``"viable"``: all internal states with absorbing flags off.
        * ``"all"``: all 2^N internal states (classical enumeration).

        Attractors with an empty basin for the chosen initial set are dropped.
        """
        if initial not in ("slow", "viable", "all"):
            raise ValueError("initial must be 'slow', 'viable' or 'all'")
        inp = self._inputs(scenario)
        nodes = tuple(self.internal)
        n = len(nodes)
        absorbing_idx = [nodes.index(a) for a in self.absorbing if a in nodes]
        slow_idx = [nodes.index(a) for a in self.slow if a in nodes]

        def counts(state: Tuple[int, ...]) -> bool:
            if initial == "all":
                return True
            if any(state[i] for i in absorbing_idx):
                return False
            if initial == "viable":
                return True
            return all(state[i] == 0 for i in range(n) if i not in slow_idx)
        # successor map over all 2^n internal states
        succ: Dict[Tuple[int, ...], Tuple[int, ...]] = {}
        for bits in product((0, 1), repeat=n):
            st = dict(zip(nodes, bits)); st.update(inp)
            nxt = self.step(st)
            succ[bits] = tuple(nxt[k] for k in nodes)
        # follow trajectories; memoize the attractor reached from each state
        reached: Dict[Tuple[int, ...], int] = {}
        attractors: List[Tuple[Tuple[int, ...], ...]] = []
        basins: List[int] = []
        for start in succ:
            path: List[Tuple[int, ...]] = []
            seen: Dict[Tuple[int, ...], int] = {}
            s = start
            while s not in reached and s not in seen:
                seen[s] = len(path); path.append(s); s = succ[s]
            if s in reached:
                aid = reached[s]
            else:  # closed a new cycle
                cyc = tuple(path[seen[s]:])
                # canonical rotation so the same cycle is recognised once
                k = min(range(len(cyc)), key=lambda i: cyc[i])
                cyc = cyc[k:] + cyc[:k]
                try:
                    aid = attractors.index(cyc)
                except ValueError:
                    attractors.append(cyc); basins.append(0); aid = len(attractors) - 1
            for p in path:
                reached[p] = aid
            basins[aid] += sum(1 for p in path if counts(p))
        return [Attractor(cyc, basins[i], nodes) for i, cyc in enumerate(attractors) if basins[i] > 0]

    def n_initial_states(self, initial: str = "slow") -> int:
        """Size of the initial-state set used for basin counting (see :meth:`attractors`)."""
        n_abs = sum(1 for a in self.absorbing if a in self.internal)
        if initial == "all":
            return 2 ** len(self.internal)
        if initial == "viable":
            return 2 ** (len(self.internal) - n_abs)
        return 2 ** sum(1 for a in self.slow if a in self.internal)

    def trajectory(self, state0: Mapping[str, int], scenario: Mapping[str, int] | None = None,
                   max_steps: int = 100) -> List[State]:
        """Synchronous trajectory from an initial internal state until it repeats."""
        inp = self._inputs(scenario)
        st = {k: 0 for k in self.internal}; st.update(state0); st.update(inp)
        out = [dict(st)]
        seen = {tuple(st[k] for k in self.internal)}
        for _ in range(max_steps):
            st = self.step(st)
            key = tuple(st[k] for k in self.internal)
            out.append(dict(st))
            if key in seen:
                break
            seen.add(key)
        return out


# ------------------------------------------------------------------ the chondrocyte network
INTERNAL = ("Piezo1", "TRPV4", "CaHigh", "YAP", "NFkB", "HIF2a", "SOX9", "RUNX2", "Senescent", "Ferroptotic")
INPUTS = ("LoadPhys", "LoadOver", "IL1b", "WNT", "TGFb", "Stiff", "Nrf2", "PTHrP")


def chondrocyte_network() -> BooleanNetwork:
    rules: Dict[str, Rule] = {
        "Piezo1": lambda s: s["LoadOver"] or (s["LoadPhys"] and s["IL1b"]),
        "TRPV4": lambda s: s["LoadPhys"] and not s["LoadOver"],
        "CaHigh": lambda s: s["Piezo1"],
        "YAP": lambda s: s["Stiff"] and not s["IL1b"],
        "NFkB": lambda s: (s["IL1b"] or s["CaHigh"] or s["Senescent"]) and not s["YAP"],
        "HIF2a": lambda s: s["NFkB"],
        "SOX9": lambda s: (not s["RUNX2"]) and (not s["NFkB"]) and (not s["YAP"]) and (not s["Ferroptotic"])
                          and (s["TGFb"] or not s["WNT"]),
        "RUNX2": lambda s: (not s["SOX9"]) and (not s["PTHrP"]) and (not s["Ferroptotic"])
                           and (s["WNT"] or s["HIF2a"] or s["CaHigh"] or s["RUNX2"]),
        "Senescent": lambda s: s["Senescent"] or (s["NFkB"] and s["CaHigh"] and not s["Ferroptotic"]),
        "Ferroptotic": lambda s: s["Ferroptotic"] or (s["CaHigh"] and not s["Nrf2"]),
    }
    return BooleanNetwork(INTERNAL, INPUTS, rules, default_inputs={"Nrf2": 1},
                          absorbing=("Senescent", "Ferroptotic"), slow=("SOX9", "RUNX2", "YAP"))


def classify_state(state: Mapping[str, int]) -> str:
    """Map a network state to a fate label of Table 1."""
    if state.get("Ferroptotic"):
        return "ferroptotic"
    if state.get("Senescent"):
        return "senescent"
    if state.get("SOX9") and not state.get("RUNX2"):
        return "articular (SOX9-high)"
    if state.get("RUNX2") and not state.get("SOX9"):
        return "hypertrophic (RUNX2-high)"
    if state.get("YAP") and not state.get("SOX9") and not state.get("RUNX2"):
        return "dedifferentiated (YAP-high)"
    if state.get("PTHrP") and not state.get("SOX9") and not state.get("RUNX2"):
        return "proliferative (PTHrP-held)"
    return "transit (saddle-like)"


SCENARIOS: Dict[str, Dict[str, int]] = {
    "basal": {},
    "physiological loading": {"LoadPhys": 1},
    "overload": {"LoadOver": 1},
    "overload, oxidative stress (Nrf2 off)": {"LoadOver": 1, "Nrf2": 0},
    "inflammation (IL-1beta)": {"IL1b": 1},
    "inflammation + physiological loading": {"IL1b": 1, "LoadPhys": 1},
    "stiff matrix (monolayer-like)": {"Stiff": 1},
    "WNT": {"WNT": 1},
    "WNT + TGF-beta": {"WNT": 1, "TGFb": 1},
    "growth plate, PTHrP present": {"PTHrP": 1, "WNT": 1},
}


def attractor_table(network: BooleanNetwork | None = None,
                    scenarios: Mapping[str, Mapping[str, int]] | None = None,
                    initial: str = "slow") -> List[Dict[str, object]]:
    """Attractors and their fate labels for every scenario (rows for a table).

    Note: a period-2 cycle alternating between SOX9/RUNX2 = (0, 0) and (1, 1) is
    the well-known artifact of synchronous updating of a mutual-inhibition motif;
    it is the discrete counterpart of the saddle of the ODE switch and is
    labelled "transit (saddle-like)".
    """
    net = network or chondrocyte_network()
    rows: List[Dict[str, object]] = []
    for name, inp in (scenarios or SCENARIOS).items():
        full_inp = net._inputs(inp)
        for a in net.attractors(inp, initial=initial):
            labels = sorted({classify_state({**d, **full_inp}) for d in a.as_dicts()})
            rows.append({
                "scenario": name,
                "type": "fixed point" if a.is_fixed_point else f"cycle (period {len(a.states)})",
                "fate": " / ".join(labels),
                "basin_fraction": a.basin_size / net.n_initial_states(initial),
                "active_nodes": ", ".join(n for n, v in zip(a.nodes, a.states[0]) if v) or "none",
            })
    return rows
