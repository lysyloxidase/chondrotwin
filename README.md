# chondrotwin

**Minimal computable models of chondrocyte fate** — code accompanying the review

> *The Chondrocyte Fate Landscape: A Mechanically Shaped, Computable Attractor Framework Integrating Signaling Networks, Mathematical Modeling, and Artificial Intelligence Toward a Digital Twin of the Chondrocyte.*
> International Journal of Molecular Sciences, Special Issue "Chondrocyte Fate: Signaling Networks in Development and Disease" (manuscript in preparation).

The paper argues that chondrocyte fate is best understood as a *mechanically shaped, computable attractor landscape* and proposes **CHONDRO-TWIN**, a mechanistic–AI hybrid digital twin of the chondrocyte. This repository makes the "computable" part concrete: every model below is small enough to read in an afternoon, runs in seconds on a laptop, depends only on NumPy, SciPy and pandas (plus Matplotlib for figures), and is covered by tests and by a numerical verification report ([`VALIDATION.md`](VALIDATION.md)). The computed panels of two manuscript figures (Figures 2 and 4) can be regenerated from this code; no figures are stored in the repository.

| Module | Model | Paper |
|---|---|---|
| `chondrotwin.toggle` | SOX9–RUNX2 bistable toggle switch (ODE): fixed points, stability, bifurcation diagram, hysteresis, memory of a transient insult | Sec. 3.1, Fig. 4, hypothesis H1 |
| `chondrotwin.growth_plate` | IHH–PTHrP reaction–diffusion loop (1-D PDE): prehypertrophic boundary; Jansen- and Blomstrand-like perturbations | Sec. 2.3, 3.4, 5.3, Fig. 2 |
| `chondrotwin.boolean` | Logical (Boolean) chondrocyte fate network with attractor enumeration; attractors = fate states of Table 1 | Sec. 5.1, Fig. 3, Table 1 |
| `chondrotwin.landscape` | Schematic landscape (Fig. 1a) and a *computed* quasi-potential U = −ln P from Langevin sampling of the switch, with basin depths and barrier heights | Sec. 2.1, 5.7, Fig. 1 |
| `chondrotwin.loading` | Inverted-U loading response, magnitude × frequency regime map, strain → network-bias mapping | Sec. 4.7, Fig. 6 |
| `chondrotwin.biphasic` | 1-D biphasic (poroelastic) confined-compression creep of cartilage: depth-resolved strain and fluid pressure | Sec. 5.5, Fig. 7 |
| `chondrotwin.twin` | Toy digital-twin loop: tissue mechanics → cell strain → switch bias → SOX9/RUNX2 → matrix turnover → mechanics | Sec. 5.6, 7.1, Fig. 9 |
| `chondrotwin.abm` | **Multi-agent cartilage model**: one agent = one chondrocyte in a depth-zoned slab, with matrix-feedback mechanics, a cytokine field, per-cell stochastic toggle switches, senescence / ferroptosis / apoptosis hazards and proliferative clusters; scenarios, OARSI-like score, scaling helpers | Sec. 5.4, 5.6, 7.1, Fig. 7 |
| `chondrotwin.plotting` | Reproduces the computed manuscript figures and diagnostic figures for every module | — |

All parameters are illustrative (dimensionless or order-of-magnitude) and are **not** fitted to data. The purpose of the code is to make the qualitative logic of the framework explicit, testable, and reusable as a starting point for calibrated models.

## Installation

```bash
git clone https://github.com/lysyloxidase/chondrotwin.git
cd chondrotwin
pip install -e ".[dev]"        # numpy, scipy, pandas, matplotlib, pytest  (add ",data" for the knee pipeline)
pytest                          # 65 tests, ~35 s (2 skip without downloaded data)
```

Python ≥ 3.9. The ODE, PDE and Boolean modules need only NumPy; the agent-based model also uses SciPy (sparse cytokine solver) and pandas (outputs). Matplotlib is imported lazily by the plotting modules.

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

## Figures (generated locally, never committed)

```bash
python scripts/make_figures.py                     # -> output/figures: Figure 2 and 4 panels + diagnostic figures
python scripts/run_abm.py                          # -> output/abm: multi-agent scenarios, time courses, scaling table (~40 s)
python scripts/attractor_table.py > docs/attractors.md       # Boolean attractors per scenario
python scripts/sweep_toggle.py                     # hysteresis window vs cooperativity
python scripts/validate.py                         # numerical verification report (VALIDATION.md; ~8 min, --quick ~1 min)
```

The reference time courses of the multi-agent scenarios (`data/reference/abm/*.csv`) are kept in the repository and
re-checked against a fresh run by `scripts/validate.py`.

## Multi-agent model: can one agent be one cell?

Yes — and the package does exactly that. `CartilageABM` seeds a 1 × 1 × 2 mm slab with 28,000 agents at depth-dependent densities (superficial zone ≈ 2.4 × the middle zone; mean 1.4e4 cells mm⁻³, the upper end of reported adult human values — Hunziker et al. 2002 give 9.6e3 mm⁻³ full-depth for the medial femoral condyle, available via `CartilageABM(cells_per_mm3=9.6e3)`); each agent carries its own SOX9–RUNX2 switch, territorial matrix, antioxidant capacity and fate flags, and the environment carries the strain field (depth-dependent modulus softened by matrix loss) and a quasi-steady IL-1β-like cytokine field. A simulated year takes about 6–12 s on a laptop.

```python
from chondrotwin.abm import CartilageABM, Scenario
sim = CartilageABM()                                       # 1 x 1 x 2 mm, 28,000 chondrocytes
res = sim.run(Scenario("overload", sigma0=0.12e6, days=365))
print(res.final()["articular"], res.final()["oarsi_like"])  # fate fractions, matrix by zone, OARSI-like score
cells = res.cells                                          # per-cell DataFrame (position, S, R, M, strain, fate ...)
```

`python scripts/run_abm.py` runs six scenarios (physiological, sustained overload, overload episode, focal overload, inflamed joint, senescent seeding) and writes cross-sections, time courses, depth profiles and a scaling table. Tissue-level behaviour that is not written as a rule but follows from the model structure (attributions checked by ablation in [`VALIDATION.md`](VALIDATION.md), §8): degeneration starts at the surface — a direct consequence of the prescribed depth-dependent modulus (with a uniform modulus it disappears); a focal overload produces a focal lesion; a 60-day overload episode leaves progressive damage after the load returns to normal — driven mainly by the matrix → stiffness feedback (without it the superficial matrix partially recovers), with switch hysteresis preventing full recovery; and seeded senescent cells add a small, dose-dependent paracrine (SASP) contribution to matrix loss on top of their own reduced synthesis. Results are robust to the random seed (SD < 0.01 of every reported fraction over 5 seeds).

**Whole cartilage.** An adult knee has on the order of 2e4 mm³ of cartilage and ~1e4 cells mm⁻³ — roughly 2–3e8 chondrocytes. The agent state takes ~93 bytes per agent (62 B per slot, 1.5 × capacity for proliferation) and the vectorized daily step needs ~260 B per agent of temporaries; at ~20–35 ms per 28,000 agents per step (laptop, load-dependent), a whole knee (2.8e8 agents at the default density) at one agent per cell would need ~100 GB of RAM at peak and ~1 day per simulated year in this single-process NumPy implementation, before accounting for the cytokine solve on a joint-sized voxel grid. That calls for a chunked or GPU implementation (as in FLAME GPU or PhysiCell-class codes), and is rarely necessary. The intended use — and the paper's multiscale argument (Figure 7) — is a representative volume driven by the local strain field of a joint-scale finite element model, tiled over the joint surface; `chondrotwin.abm.tile_plan()` and `scaling_table()` give the numbers for any tissue volume, and `benchmark()` measures them on your machine.

## Knee cartilage simulation from open data

The agent-based model was also driven, calibrated and tested with open data (in vivo knee loads, 507 OAI cartilage
segmentations, single-cell and bulk transcriptomes of OA cartilage and synovium) in six stages, with pre-registered
independent tests. The full account is in [`docs/knee_study.md`](docs/knee_study.md); data sources, licenses and
citations in [`DATA.md`](DATA.md) (nothing is redistributed).

```bash
pip install -e ".[dev,data]"
python scripts/fetch_data.py                  # ~400 MB from the providers, checksummed
python scripts/knee_open_data.py              # stage 1; later stages: scripts/knee_stage{2..6}*.py (see docs/knee_study.md)
```

Main results (recomputed scorecard in [`docs/knee_stage6.md`](docs/knee_stage6.md)):

| test | initial model | final model |
|---|---|---|
| held-out patients GSE104782 OA6–10 (loss, lower is better) | 798 | 88 |
| GSE152805, damaged vs intact tibia (pre-registered) | 0 / 4 | 1 / 4 |
| GSE57218, 33 joint-matched pairs (pre-registered) | 0 / 4 | 2 / 4 |
| GSE117999, OA vs meniscus-tear knees (pre-registered) | 0 / 5 | 2 / 5 |
| GSE169077, OA vs normal (pre-registered, secondary) | 1 / 5 | 3 / 5 |

* **Supported:** medial-predominant cartilage loss follows from the measured in vivo medial load share; strain alone
  cannot flip the SOX9 / RUNX2 switch (identifiable, and robust across 48 analysis variants).
* **Not supported:** the model's senescence is early and near-complete in every dataset tested (also with a 124-gene
  SenMayo readout); it does not reproduce the OA-vs-normal-joint difference (SOX9 −2.2, MMP13 +3.6 log2); the synovial
  inflammation hypothesis could not be supported with open data.

## What the models show

* **Bistability and hysteresis** (`toggle`). Mutual Hill-type repression of SOX9 and RUNX2 yields two attractors separated by a saddle for |λ| < 0.51 (defaults: α = 2, K = 0.5, n = 4, basal 0.05). A pulse λ = 0.8 switches the cell to the RUNX2-high state, which persists after the pulse ends (hypothesis H1). Without cooperativity (n = 1) there is no bistability.
* **Zonal patterning from one feedback loop** (`growth_plate`). A periarticular PTHrP source and IHH feedback produce a sharp prehypertrophic boundary; lengthening the PTHrP range (Jansen-like) moves it distally, removing PTHrP (Blomstrand-like) commits every cell.
* **Fate states as attractors** (`boolean`). Starting from the prior fate of a viable cell, the network reaches: articular or hypertrophic (basal; prior-fate-dependent), senescent (overload), ferroptotic (overload without Nrf2), hypertrophic (IL-1β or WNT), senescent (IL-1β + physiological load: Piezo1 sensitization), dedifferentiated (stiff matrix), proliferative (growth plate with PTHrP: SOX9-high, RUNX2 held off by PTHrP). See `docs/attractors.md`.
* **A computable landscape** (`landscape`). U = −ln P from Langevin sampling gives numerical barrier heights; the escape barrier from the articular (SOX9-high) basin falls as the catabolic bias rises. At low noise walkers rarely cross between basins, so for λ ≠ 0 the *relative* depth of the two basins reflects the initial allocation of walkers rather than the stationary distribution (the function reports the crossing fraction; see `VALIDATION.md`, §4).
* **Tissue degeneration from cell rules** (`abm`). With one agent per cell, surface-first degeneration (from the depth-dependent modulus), focal lesions, and progressive damage after a transient overload (matrix → stiffness feedback plus switch hysteresis) arise from the per-cell switch, matrix feedback and hazards; senescent-cell seeding adds a small paracrine contribution (see `scripts/run_abm.py` and `VALIDATION.md`).
* **Closing the loop** (`twin`). Coupling the biphasic tissue to the switch and to matrix turnover produces the "vicious cycle" of overload → RUNX2 → degradation → higher strain, and shows that a short overload episode can commit the twin to progressive degradation even though the bias first returns to anabolic values (λ ≈ −0.14, inside the hysteresis window); continued matrix loss then raises the strain until the bias turns catabolic again (from cycle 26).

The equations and default parameters are documented in [`docs/models.md`](docs/models.md).

## Repository layout

```
chondrotwin/          package (numpy/scipy numerics; plotting optional); abm/ = multi-agent cartilage model; knee/ = open-data pipeline
tests/                pytest suite
scripts/              figures (to output/), attractor table, parameter sweep, verification (validate.py), data download and knee pipeline
examples/             quickstart tour
docs/                 model equations (models.md), attractor table (attractors.md), knee study (knee_study.md) and its stage reports
                      (knee_open_data.md, knee_stage2–6.md), pre-registrations and frozen parameters
data/reference/abm/   reference ABM time courses (CSV) checked by validate.py
VALIDATION.md         numerical verification: analytical solutions, convergence, reproducibility, seeds, ablations
DATA.md               open datasets of the knee pipeline: sources, sizes, licenses, citations
data/raw, data/derived downloaded and derived data - not tracked (output/ likewise)
```

## Limitations

These are deliberately minimal models. The toggle switch has two variables; the growth-plate model is one-dimensional and steady-state; the Boolean network has ten internal nodes and uses synchronous updating (the period-2 cycle that appears for the SOX9/RUNX2 pair under some inputs is the well-known synchronous artifact of mutual inhibition and corresponds to the saddle of the ODE); the biphasic model is linear and confined; the twin loop treats "cycles" as abstract loading epochs. Applied stresses in `twin` and `abm` (0.05 MPa "physiological", 0.12–0.15 MPa "injurious") are effective drained (solid-matrix) stresses chosen to give equilibrium tissue strains of ~10–30 %; they are not joint contact pressures, which are of the order of MPa but carried largely by interstitial fluid pressurization under physiological loading rates. The Langevin landscape is quasi-stationary within basins but not equilibrated between them at low noise. None of the parameters are calibrated. Extending any component toward the full CHONDRO-TWIN architecture (Section 7 of the paper) requires calibration against data, identifiability analysis, and uncertainty quantification as described there.

## Citing

If you use this code, please cite the accompanying article once it is published, and this repository.

## License

MIT — see `LICENSE`.
