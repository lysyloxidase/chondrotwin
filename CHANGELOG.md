# Changelog

## 0.9.0 — repository housekeeping
- No figures in the repository: the manuscript schematic script and the `figures/` folder were removed (the corrected
  schematic script is in the git history of the previous commit); scripts write figures to the ignored `output/`.
- Reference ABM time courses moved to `data/reference/abm/` (checked by `scripts/validate.py`).
- Author fields removed (`pyproject.toml`; `CITATION.cff` deleted); repository URL set.

## 0.9.0 (2026-10-08) — knee simulation from open data, stage 6 (uncertainty and robustness)
- `chondrotwin.knee.uq`, `scripts/knee_stage6_run.py`, `scripts/knee_stage6.py`: loss profile over the stage-2
  parameters with a misspecification-scaled 95 % region; prediction ranges inside it; a multiverse of 48 analysis
  choices; one scorecard of every test of stages 1-5, recomputed from the data.
- Findings: c_over_max ∈ [0, 0.2] (robustly below the saddle-node); k_cat not identified — the stage-2 claim that
  catabolism partly bypasses the switch is withdrawn; stage 2 beats stage 1 on held-out patients in 48/48 universes;
  near-complete senescence in degraded tissue holds over the whole parameter region.

## 0.8.0 (2026-10-08) — knee simulation from open data, stage 5
- ABM: optional joint-level cytokine background (`Scenario.C_joint`), inflammation-driven catabolism (`Rates.k_inf`) and a
  parameter for cytokine amplification of senescence (`Rates.sen_infl_gain`, default 2 as before); defaults bit-identical.
- Calibration on within-joint (GSE104782 OA1-5) plus between-joint (GSE114007 raw counts) contrasts
  (`calibrate.evaluate_joint`, `scripts/knee_stage5_calibrate.py`): not selected; the frozen model is the stage-2 model.
- Pre-registered premise test on synovium (`chondrotwin/knee/independent4.py`, GSE55235 + GSE55457): OA synovium not
  more inflamed than controls (valid RA positive control); post hoc caveat: the controls include joint-trauma patients.
- New data sources: hallmark inflammatory gene set, GSE55235/GSE55457; GSE178557 downloaded but unusable (no annotation).

## 0.7.0 (2026-10-08) — knee simulation from open data, stage 4
- Senescence re-measured with the SenMayo gene set (`chondrotwin.knee.senescence`): the stage 1-3 failure is confirmed
  in GSE104782, GSE152805 and GSE57218 with a readout that shares no gene with CDKN2A.
- ABM: optional inflammation-gated senescence (`Rates.sen_overload_weight`) and death of senescent cells
  (`Rates.senescent_death_rate`); calibration with a SenMayo target (`calibrate.PAIRS4`, `scripts/knee_stage4_calibrate.py`)
  did not select either — the frozen stage-4 model is the stage-2 model; the loss breakdown shows why (SASP is the only
  intra-cartilage inflammatory source, coupling senescence to catabolism).
- Pre-registered test (`chondrotwin/knee/independent3.py`, `docs/preregistration_stage4.md`): GSE117999 2/5 readouts
  (stage 1: 0/5), GSE169077 3/5; GSE114007 repaired with raw counts as a disclosed deviation (3/5).
- Cross-dataset synthesis: MMP13 is a joint-level, not lesion-level, OA feature — hypothesis for stage 5.
- New data sources SenMayo, GSE117999, GSE169077 (+ GPL96), GSE114007 raw counts.

## 0.6.0 (2026-10-08) — knee simulation from open data, stage 3
- ABM: optional cell-to-cell heterogeneity — per-cell switch offset (`Rates.switch_bias_sd`) and senescence competence
  (`Rates.senescence_competent_fraction`), drawn from a separate random stream (defaults bit-identical).
- Stage-3 calibration on GSE104782 OA1-5 (`scripts/knee_stage3_calibrate.py`): neither mechanism was selected; the
  frozen stage-3 parameters equal stage 2 (`docs/stage3_frozen_parameters.json`).
- Pre-registered bulk test (`chondrotwin/knee/independent2.py`, `docs/preregistration_stage3.md`): GSE57218 (33
  joint-matched pairs) 2/4 markers for the selected model (stage 1: 0/4); GSE114007 declared invalid post hoc (negative
  values in its "normalized counts"); deviations and diagnostics reported separately in `docs/knee_stage3.md`.
- New data sources GSE57218 (+ GPL6947 annotation) and GSE114007 in `chondrotwin.datasets` / `DATA.md`.

## 0.5.0 (2026-10-08) — knee simulation from open data, stage 2
- ABM: three optional mechanisms, each answering a stage-1 falsification; defaults reproduce the earlier model
  bit-for-bit: cumulative-damage senescence (`Rates.damage_D50`, per-cell `damage`), bounded mechanical bias of the
  switch (`Rates.c_over_max`), switch-independent catabolism (`Rates.k_cat`); `Rates.stage2(...)`.
- Calibration on GSE104782 patients OA1-5 with a patient split fixed in advance (`chondrotwin.knee.calibrate`,
  `scripts/knee_stage2_grid.py`, multi-seed selection `scripts/knee_stage2_select.py`); validation on OA6-10.
- Independent, pre-registered test on GSE152805 (damaged medial vs minimally damaged lateral tibia; protocol and
  its SHA-256 in `docs/preregistration_gse152805.md`, written before the data were opened).
- `scripts/knee_stage2.py` writes `docs/knee_stage2.md`. Results: held-out-patient loss 798 -> 88 (SOX9 partly
  retained, no switch by strain alone); senescence timing not improved; pre-registered test 1/4 (stage 1: 0/4), with
  the direction of change right for all markers; tissue prediction unchanged; the stage-2 model has a latent period
  before erosion.

## 0.4.0 (2026-10-08) — knee simulation from open data, stage 1
- `chondrotwin.datasets` + `scripts/fetch_data.py`: download of OrthoLoad knee loads, GEO GSE104782 and OAIZIB-CM
  segmentations from the providers, with pinned sizes and SHA-256; licenses and citations in `DATA.md`.
- `chondrotwin.knee`: in vivo medial / lateral force split (condylar distance calibrated to published peaks),
  tibial cartilage geometry of 507 OAI knees (central thickness, denuded area), mechanical bridge (time averaging
  verified with the new `ConfinedCompression.solve_history`), ABM dose-response, calibration on the train split and
  out-of-sample tests (lateral plate, held-out knees, cell states by degeneration stage).
- `scripts/knee_open_data.py` writes `docs/knee_open_data.md` (all numbers generated). Main findings: medial-
  predominant loss follows from the measured load share; the ABM's cell-state transitions are too early and too
  complete and its tissue response too threshold-like (see the report).
- Tests for the new modules (data-dependent ones are skipped when the data are not downloaded).

## 0.3.0 (2026-10-08) — verification release
- Numerical verification: `scripts/validate.py` and `VALIDATION.md` (analytical creep solution, time/grid convergence,
  Langevin sampling diagnostics, reproducibility of the published ABM time courses, 5-seed variability, ablations,
  measured memory per agent); new tests for the analytical biphasic solution, grid convergence, landscape barriers
  and the committed attractor table (43 tests).
- Boolean network: PTHrP now also permits SOX9 (PTHrP/PKA activates SOX9; Huang et al., PNAS 2001), so the
  growth-plate attractor "proliferative (PTHrP-held)" is SOX9-high as in the proliferative zone (was all-off).
  All other attractors, labels and basin fractions are unchanged.
- Growth plate: `steady_state` default horizon 12 -> 40 so the Jansen-like run actually converges (boundary unchanged, 0.37).
- Landscape: `quasi_potential_langevin` reports `crossed` (inter-basin crossing fraction) and accepts
  `start_fraction_sox9`; documented that basin depths are not equilibrated between basins at low noise.
- ABM scaling: corrected memory estimates (80 B per agent of state, ~340 B at peak, was 60 B) and whole-knee figures.
- Manuscript Figure 5: inflammatory Piezo1 sensitization corrected to IL-1α → p38 MAPK (HNF4, ATF2) → Piezo1↑
  (Lee et al., PNAS 2021; was "IL-1β → PI3K/AKT/mTOR").
- README / docs: dependency list, cell-density sources, meaning of the applied stresses (effective drained stresses,
  not joint contact pressures), ablation-based attribution of the ABM tissue-level behaviours, twin-loop bias after an
  overload episode.
- Figure images removed from the repository (they are published in the article; the scripts regenerate them);
  the ABM reference time courses (`figures/abm/*.csv`) are kept for the reproducibility check.

## 0.2.0 (2026-10-08)
- Multi-agent cartilage model (`chondrotwin.abm`): one agent = one chondrocyte; depth-zoned slab, matrix-feedback
  mechanics, quasi-steady cytokine field, stochastic toggle switch per cell, senescence / ferroptosis / apoptosis hazards,
  proliferative clusters, OARSI-like score; scenarios (physiological, overload, overload episode, focal, inflamed,
  senescent seeding); scaling helpers and benchmark for "whole cartilage" estimates; `scripts/run_abm.py`; tests.

## 0.1.0 (2026-10-08)
- Initial release accompanying the manuscript: toggle switch, growth-plate reaction-diffusion model,
  Boolean fate network, Langevin quasi-potential landscape, biphasic confined compression,
  digital-twin loop, figure scripts and tests.
