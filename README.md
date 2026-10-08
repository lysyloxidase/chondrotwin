# chondrotwin

**Minimal computable models of chondrocyte fate** — code accompanying the review

> *The Chondrocyte Fate Landscape: A Mechanically Shaped, Computable Attractor Framework Integrating Signaling Networks, Mathematical Modeling, and Artificial Intelligence Toward a Digital Twin of the Chondrocyte.*
> International Journal of Molecular Sciences, Special Issue "Chondrocyte Fate: Signaling Networks in Development and Disease" (manuscript in preparation).

The paper argues that chondrocyte fate is best understood as a *mechanically shaped, computable attractor landscape* and proposes **CHONDRO-TWIN**, a mechanistic–AI hybrid digital twin of the chondrocyte. This repository makes the "computable" part concrete: every model below is small enough to read in an afternoon, runs in seconds on a laptop, depends only on NumPy (plus Matplotlib for figures), and is covered by tests. Two of the manuscript figures (Figure 2 and Figure 4) are produced directly from this code.

| Module | Model | Paper |
|---|---|---|
| `chondrotwin.toggle` | SOX9–RUNX2 bistable toggle switch (ODE): fixed points, stability, bifurcation diagram, hysteresis, memory of a transient insult | Sec. 3.1, Fig. 4, hypothesis H1 |
| `chondrotwin.growth_plate` | IHH–PTHrP reaction–diffusion loop (1-D PDE): prehypertrophic boundary; Jansen- and Blomstrand-like perturbations | Sec. 2.3, 3.4, 5.3, Fig. 2 |
| `chondrotwin.boolean` | Logical (Boolean) chondrocyte fate network with attractor enumeration; attractors = fate states of Table 1 | Sec. 5.1, Fig. 3, Table 1 |
| `chondrotwin.landscape` | Schematic landscape (Fig. 1a) and a *computed* quasi-potential U = −ln P from Langevin sampling of the switch, with basin depths and barrier heights | Sec. 2.1, 5.7, Fig. 1 |
| `chondrotwin.loading` | Inverted-U loading response, magnitude × frequency regime map, strain → network-bias mapping | Sec. 4.7, Fig. 6 |
| `chondrotwin.biphasic` | 1-D biphasic (poroelastic) confined-compression creep of cartilage: depth-resolved strain and fluid pressure | Sec. 5.5, Fig. 7 |
| `chondrotwin.twin` | Toy digital-twin loop: tissue mechanics → cell strain → switch bias → SOX9/RUNX2 → matrix turnover → mechanics | Sec. 5.6, 7.1, Fig. 9 |
| `chondrotwin.plotting` | Reproduces the computed manuscript figures and diagnostic figures for every module | — |

All parameters are illustrative (dimensionless or order-of-magnitude) and are **not** fitted to data. The purpose of the code is to make the qualitative logic of the framework explicit, testable, and reusable as a starting point for calibrated models.

## Installation

```bash
git clone https://github.com/<org>/chondrotwin.git
cd chondrotwin
pip install -e ".[dev]"        # numpy, matplotlib, pytest
pytest                          # ~30 tests, ~15 s
```

Python ≥ 3.9. The numerical modules need only NumPy; Matplotlib is imported lazily by `chondrotwin.plotting`.

## Quick start

```python
from chondrotwin import ToggleSwitch, GrowthPlateModel, chondrocyte_network, classify_state

sw = ToggleSwitch()                       # SOX9–RUNX2 switch, default parameters
for fp in sw.fixed_points(lam=0.0):       # lam: catabolic (+) / anabolic (−) bias
    print(fp.label, round(fp.S, 2), round(fp.R, 2))
print(sw.bistable_range())                # (-0.51, 0.51): the hysteresis window
print(sw.switches_permanently(lam_pulse=0.8))   # True: a transient insult commits the cell

gp = GrowthPlateModel()
print(gp.steady_state().boundary(), gp.jansen_like().steady_state().boundary())  # 0.205 -> 0.37

net = chondrocyte_network()
print(classify_state(net.trajectory({"SOX9": 1}, {"LoadOver": 1})[-1]))            # 'senescent'
print(classify_state(net.trajectory({"SOX9": 1}, {"LoadOver": 1, "Nrf2": 0})[-1]))  # 'ferroptotic'
```

See `examples/quickstart.py` for a five-minute tour of every module.

## Reproducing the figures

```bash
python scripts/make_figures.py figures            # Figure 2, Figure 4 + diagnostic figures
python scripts/manuscript_schematics.py figures/manuscript   # schematic Figures 1, 3, 5–10
python scripts/attractor_table.py > docs/attractors.md       # Boolean attractors per scenario
python scripts/sweep_toggle.py                    # hysteresis window vs cooperativity
```

`scripts/make_figures.py` writes:

* `figure2_growth_plate.png` — Figure 2 of the paper (zonation schematic + reaction–diffusion model + perturbations)
* `figure4_toggle_switch.png` — Figure 4 of the paper (phase plane, bifurcation diagram, memory of a pulse)
* `landscape_langevin.png` — schematic landscape next to Langevin-sampled quasi-potentials of the switch
* `loading_regime_map.png` — inverted-U response and magnitude × frequency map
* `boolean_attractors.png` — basin fractions of each fate state per input scenario
* `biphasic_creep.png` — confined-compression creep and interstitial pressure profiles
* `twin_loop.png` — the digital-twin loop under physiological loading, sustained overload, and an overload episode (hysteresis: the cell stays hypertrophic after the load returns to normal)

## What the models show

* **Bistability and hysteresis** (`toggle`). Mutual Hill-type repression of SOX9 and RUNX2 yields two attractors separated by a saddle for |λ| < 0.51 (defaults: α = 2, K = 0.5, n = 4, basal 0.05). A pulse λ = 0.8 switches the cell to the RUNX2-high state, which persists after the pulse ends (hypothesis H1). Without cooperativity (n = 1) there is no bistability.
* **Zonal patterning from one feedback loop** (`growth_plate`). A periarticular PTHrP source and IHH feedback produce a sharp prehypertrophic boundary; lengthening the PTHrP range (Jansen-like) moves it distally, removing PTHrP (Blomstrand-like) commits every cell.
* **Fate states as attractors** (`boolean`). Starting from the prior fate of a viable cell, the network reaches: articular or hypertrophic (basal; prior-fate-dependent), senescent (overload), ferroptotic (overload without Nrf2), hypertrophic (IL-1β or WNT), senescent (IL-1β + physiological load: Piezo1 sensitization), dedifferentiated (stiff matrix), proliferative (growth plate with PTHrP). See `docs/attractors.md`.
* **A computable landscape** (`landscape`). U = −ln P from Langevin sampling gives numerical basin depths and barrier heights that fall as the catabolic bias rises.
* **Closing the loop** (`twin`). Coupling the biphasic tissue to the switch and to matrix turnover produces the "vicious cycle" of overload → RUNX2 → degradation → higher strain, and shows that a short overload episode can commit the twin to progressive degradation even though the bias returns to anabolic values.

The equations and default parameters are documented in [`docs/models.md`](docs/models.md).

## Repository layout

```
chondrotwin/          package (numpy-only numerics; plotting optional)
tests/                pytest suite
scripts/              figure generation, attractor table, parameter sweep, manuscript schematics
examples/             quickstart tour
docs/                 model equations (models.md), attractor table (attractors.md)
figures/              generated figures (created by the scripts)
```

## Limitations

These are deliberately minimal models. The toggle switch has two variables; the growth-plate model is one-dimensional and steady-state; the Boolean network has ten internal nodes and uses synchronous updating (the period-2 cycle that appears for the SOX9/RUNX2 pair under some inputs is the well-known synchronous artifact of mutual inhibition and corresponds to the saddle of the ODE); the biphasic model is linear and confined; the twin loop treats "cycles" as abstract loading epochs. None of the parameters are calibrated. Extending any component toward the full CHONDRO-TWIN architecture (Section 7 of the paper) requires calibration against data, identifiability analysis, and uncertainty quantification as described there.

## Citing

If you use this code, please cite the paper (see `CITATION.cff`; the citation will be updated on publication).

## License

MIT — see `LICENSE`.
