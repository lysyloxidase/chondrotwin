# Knee cartilage simulation from open data — the study, stage by stage

This document collects the knee study that was added to the package after the review's code: an agent-based cartilage model driven by open data, revised in stages, each stage calibrated on one part of the data and tested on another, with every independent test pre-registered (protocol and frozen parameters hashed before the test data were read). Each stage has a generated report in this folder; numbers below are copied from those reports.

## Method in brief

* **Data** — ten open datasets and two gene sets, downloaded from the providers and checksummed (`scripts/fetch_data.py`; sources, licenses and citations in [`DATA.md`](../DATA.md)); nothing is redistributed.
* **Calibration vs testing** — single-cell patients GSE104782 OA1–5 calibrate, OA6–10 validate (split fixed before any fitting); external datasets are pre-registered tests (`docs/preregistration_*.md`) or are labelled as previously seen.
* **Selection** — grid search, the best sets re-run with three seeds, selection by mean calibration loss only.
* **Precedence** — generated reports are never edited after the fact; where a later stage revises an earlier conclusion (e.g. the identifiability of `k_cat` in stage 6), the later stage takes precedence.
* **Defaults** — every added mechanism defaults to the published model; the published time courses reproduce bit for bit (`tests/test_abm.py`, `scripts/validate.py`).

## Stage 1 — open-data pipeline


The agent-based model is driven and tested with three open datasets (details, licenses and citations in
[`DATA.md`](../DATA.md); nothing is redistributed here):

* **loads** — in vivo tibial loads during walking in 8 subjects with instrumented knee replacements (OrthoLoad);
  the medial / lateral split uses a condylar distance calibrated to the published medial peaks (47 mm);
* **geometry** — medial and lateral tibial cartilage of 507 OAI knees with KL grade, sex, age and BMI (OAIZIB-CM),
  using the dataset's own train / test split;
* **cell states** — 1600 chondrocytes from 10 patients at five degeneration stages (GEO GSE104782).

```bash
pip install -e ".[dev,data]"
python scripts/fetch_data.py          # ~93 MB, checksummed
python scripts/knee_open_data.py      # writes docs/knee_open_data.md (~1 min; --recompute ~12 min)
```

One common load level is calibrated on the medial plate of the training knees; everything else is a prediction.
Results ([`docs/knee_open_data.md`](knee_open_data.md)):

* **Supported.** Medial-predominant cartilage loss follows from the measured medial force share (0.69 ± 0.11) without
  being fitted, and a valgus load distribution reverses it. MMP13+ and p16+ cells increase in the most degraded
  cartilage, as in the model.
* **Not supported.** In the model, senescence and the SOX9 → RUNX2 switch appear already at early degeneration and
  affect nearly all cells; in the data, p16 and RUNX2 do not rise before the last stage and > 90 % of cells keep SOX9.
  The tissue response is too threshold-like (no lateral loss predicted; 0.18 ± 0.13 mm observed), and early-OA
  thickening is absent from the model.

This is a first, cross-sectional calibration: it uses loads from replaced (not native) knees and no longitudinal
data.


## Stage 2

**Stage 2** ([`docs/knee_stage2.md`](knee_stage2.md)) revises the cell-state rules where stage 1 was falsified,
calibrates them on half of the single-cell patients (OA1–5) and tests them on the other half and on an independent
dataset (GSE152805, damaged medial vs intact lateral tibia) under a protocol registered before the data were opened
([`docs/preregistration_gse152805.md`](preregistration_gse152805.md)):

```bash
python scripts/knee_stage2_grid.py && python scripts/knee_stage2_grid.py --refine && python scripts/knee_stage2_select.py
python scripts/knee_stage2.py         # writes docs/knee_stage2.md
```

* held-out patients: loss 798 → 88 — strain alone no longer flips the switch (inflammation is needed) and SOX9 is
  partly retained in degraded tissue (41 % of intact vs 60 % in the data); the second revision (catabolism partly
  bypassing the switch) turned out not to be identifiable in stage 6;
* not fixed: senescence still appears too early and saturates (a cumulative-damage hazard was tested and not selected);
* independent pre-registered test: 1 / 4 markers within the per-patient range (stage 1: 0 / 4); the direction of change
  is right for all four, the magnitudes fail because the model has no RUNX2+ cells in intact cartilage;
* tissue: lateral loss still predicted ≈ 0; the revised model has a latent period before erosion, so longitudinal
  OAI data (NIH NDA agreement) are needed to fix its time scale.


## Stage 3

**Stage 3** ([`docs/knee_stage3.md`](knee_stage3.md)) tested cell-to-cell heterogeneity (a per-cell offset of the
switch, and limited senescence competence) and a new, pre-registered independent test on bulk transcriptomes
([`docs/preregistration_stage3.md`](preregistration_stage3.md); protocol and frozen parameters hashed before the
expression data were read):

```bash
python scripts/fetch_data.py gse57218 gse114007
python scripts/knee_stage3_calibrate.py && python scripts/knee_stage3.py     # writes docs/knee_stage3.md
```

* neither heterogeneity mechanism improved the calibration, so the frozen stage-3 model **is** the stage-2 model;
* primary test, GSE57218 (33 joint-matched preserved / OA-affected pairs): 2 / 4 markers (SOX9, RUNX2 consistent);
  fails MMP13 (lower in OA-affected tissue, model higher) and CDKN2A (flat but near the array background);
* the secondary test (GSE114007) was declared invalid after the fact: its "normalized counts" contain negative values;
* across all stages, senescence and end-stage catabolism are where the model's logic most likely fails.


## Stage 4

**Stage 4** ([`docs/knee_stage4.md`](knee_stage4.md); [`docs/preregistration_stage4.md`](preregistration_stage4.md))
asked whether the senescence failure is real and tested revised senescence logic:

```bash
python scripts/fetch_data.py senmayo gse117999 gse169077 gse114007
python scripts/knee_stage4_calibrate.py && python scripts/knee_stage4.py     # writes docs/knee_stage4.md
```

* the SenMayo senescence / SASP score (MSigDB version, 124 genes; no overlap with p16) confirms the failure in three datasets: senescence
  is late and moderate, the model's is early and near-complete;
* inflammation-gated senescence and death of senescent cells were not selected: in the model the senescent cells' SASP
  is the only inflammatory source inside cartilage, so less senescence also means less catabolism;
* pre-registered test (GSE117999, OA vs meniscus-tear knees): 2 / 5 readouts (stage 1: 0 / 5); secondary GSE169077 3 / 5;
  repaired GSE114007 (raw counts, a disclosed deviation) 3 / 5;
* across four bulk datasets MMP13 is ~12-fold higher in OA than in normal knees but not higher within OA joints —
  catabolism looks like a joint-level (synovial) rather than a lesion-level process, the hypothesis for stage 5.


## Stage 5

**Stage 5** ([`docs/knee_stage5.md`](knee_stage5.md); [`docs/preregistration_stage5.md`](preregistration_stage5.md))
tested joint-level (synovial) inflammation:

```bash
python scripts/fetch_data.py hallmark_inflammation gse55235
python scripts/knee_stage5_calibrate.py && python scripts/knee_stage5.py     # writes docs/knee_stage5.md
```

* joint-level inflammation and inflammation-driven catabolism improved the OA-vs-normal fit but worsened the
  within-joint cell fit; not selected — the model remains the stage-2 model;
* pre-registered premise test (synovium, GSE55235 + GSE55457): OA synovium is **not** more inflamed than control
  synovium (−0.09 [−0.24, 0.07]), with a valid RA positive control — but the controls include joint-trauma patients;
* the largest remaining gap is between OA and normal joints (SOX9 −2.2 / MMP13 +3.6 log2 in the data, −0.1 / +1.4
  in the model); further progress needs data that are not openly available (OAI longitudinal data, paired cohorts).


## Stage 6

**Stage 6** ([`docs/knee_stage6.md`](knee_stage6.md)) adds no hypothesis and no data; it quantifies how certain the
selected model is and collects every test in one recomputed scorecard:

```bash
python scripts/knee_stage6_run.py && python scripts/knee_stage6.py     # ~15 min; writes docs/knee_stage6.md
```

* identifiability: the strain-only switch bias is confined to c_over_max ∈ [0, 0.2], below the saddle-node (0.51) —
  "strain alone cannot flip the switch" is robust; k_cat is **not identified** (0 is inside the region), so the
  "catabolism bypasses the switch" part of stage 2 is withdrawn;
* robustness: the stage-2 model beats stage 1 on held-out patients in **48 / 48** analysis universes (simulation
  settings, stage bins, senescence readout, pseudo-count; loss ratio 1.8–10.3);
* near-complete senescence in degraded tissue holds over the whole parameter region — a structural failure.


## Overall

Only the stage-2 revision that strain alone cannot flip the SOX9 / RUNX2 switch survived; it is identifiable and robust (stage 6). The model's senescence (early, near-complete) and the gap between OA and normal joints are not fixed by any mechanism testable with these open data. The recomputed scorecard of all tests is in [`knee_stage6.md`](knee_stage6.md), section 4.
