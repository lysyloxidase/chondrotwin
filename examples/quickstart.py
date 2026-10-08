"""Five-minute tour of chondrotwin (run from the repository root)."""
import numpy as np
from chondrotwin import (ToggleSwitch, GrowthPlateModel, chondrocyte_network, classify_state,
                         quasi_potential_langevin, ConfinedCompression, MiniTwin)
from chondrotwin.landscape import basin_depths

# 1. The SOX9-RUNX2 switch: fixed points, hysteresis window, memory of a pulse
sw = ToggleSwitch()
for fp in sw.fixed_points(0.0):
    print(f"lambda=0: {fp.label:28s} S={fp.S:.2f} R={fp.R:.2f}")
print("bistable window:", sw.bistable_range())
print("permanent switch after a lambda=0.8 pulse:", sw.switches_permanently(lam_pulse=0.8))

# 2. The growth plate: where is the prehypertrophic boundary?
gp = GrowthPlateModel()
print("boundary WT / Jansen-like / Blomstrand-like:",
      round(gp.steady_state().boundary(), 3), round(gp.jansen_like().steady_state().boundary(), 3),
      round(gp.blomstrand_like().steady_state().boundary(), 3))

# 3. The Boolean network: fate reached from the articular state under overload
net = chondrocyte_network()
print("overload ->", classify_state(net.trajectory({"SOX9": 1}, {"LoadOver": 1})[-1]))
print("overload, Nrf2 off ->", classify_state(net.trajectory({"SOX9": 1}, {"LoadOver": 1, "Nrf2": 0})[-1]))

# 4. A computable landscape: barrier heights from Langevin sampling
res = quasi_potential_langevin(lam=0.0, noise=0.3, n_walkers=1500, n_steps=3000, burn_in=1000)
print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in basin_depths(res["U"], res["S_edges"], res["R_edges"]).items()})

# 5. Tissue mechanics and the twin loop
tissue = ConfinedCompression()
print(f"consolidation time tau = {tissue.tau:.0f} s; equilibrium strain at 0.05 MPa = {tissue.equilibrium_strain(0.05e6):.3f}")
twin = MiniTwin()
a = MiniTwin.to_arrays(twin.run(lambda c: 0.15e6 if 10 <= c < 22 else 0.05e6, n_cycles=40))
print("H_A after an overload episode (MPa):", np.round(a["H_A"][[0, 9, 21, 39]] / 1e6, 3), "| final fate:",
      "hypertrophic" if a["R"][-1] > a["S"][-1] else "articular")
