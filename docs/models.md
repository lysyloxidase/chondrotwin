# Model equations and default parameters

All models are illustrative. Units are dimensionless unless stated otherwise.

## 1. SOX9–RUNX2 toggle switch (`chondrotwin.toggle`)

State: SOX9 activity S, RUNX2 activity R. Bias λ ∈ (−1, 1).

```
dS/dt = α_S (1 − λ) / (1 + (R/K)^n) + β_S − γ S
dR/dt = α_R (1 + λ) / (1 + (S/K)^n) + β_R − γ R
```

| Parameter | Default | Meaning |
|---|---|---|
| α_S, α_R | 2.0 | maximal production |
| K | 0.5 | half-maximal repression constant |
| n | 4 | Hill coefficient (cooperativity) |
| β_S, β_R | 0.05 | basal production |
| γ | 1.0 | first-order decay |

* λ > 0: catabolic bias (overload, inflammation, WNT, HIF-2α); λ < 0: anabolic bias (TGF-β, physiological loading, PTHrP).
* Fixed points are found by reducing to one dimension on the S-nullcline and bracketing the roots of `R_nullcline(S*(R)) − R`; stability from the Jacobian eigenvalues.
* Defaults give bistability for |λ| < 0.51 (saddle-node at λ* = ±0.5146 by bisection); the saddle-node bifurcations bound the hysteresis loop. A pulse λ = 0.8 for 20 time units switches a SOX9-high cell permanently to RUNX2-high (`ToggleSwitch.pulse`); the smallest 20-unit pulse that does so is λ ≈ 0.53.
* `scripts/sweep_toggle.py` shows that the hysteresis window widens with n and shrinks with K.

## 2. IHH–PTHrP reaction–diffusion loop (`chondrotwin.growth_plate`)

Axis x ∈ [0, 1] from the periarticular end (0) to the chondro-osseous junction (1). PTHrP P(x, t), IHH I(x, t).

```
∂P/∂t = D_P ∂²P/∂x² + a_P s(x) [ b + I^h / (K_I^h + I^h) ] − k_P P
∂I/∂t = D_I ∂²I/∂x² + a_I K_P^h / (K_P^h + P^h)              − k_I I
H(x)  = K_P^h / (K_P^h + P^h)          (hypertrophic commitment)
s(x)  = exp(−x² / (2 w²))              (periarticular PTHrP source)
```

| Parameter | Default |
|---|---|
| D_P, D_I | 0.03, 0.006 |
| a_P, k_P | 2.0, 1.0 |
| a_I, k_I | 1.0, 1.0 |
| K_P, K_I | 0.3, 0.3 |
| b (basal PTHrP production) | 0.25 |
| w (source width) | 0.08 |
| h (Hill coefficient) | 2 |

Zero-flux boundaries; explicit finite differences (201 nodes) integrated to steady state (t ≤ 40, stopping when max |∂P/∂t|, |∂I/∂t| < 10⁻⁴; `converged` reports whether this was met). The prehypertrophic boundary is the first x where H ≥ 0.5 (≈ 0.205 for the defaults; 0.204–0.205 by linear interpolation on 201 and 401 nodes). `jansen_like()` multiplies k_P by 0.35 (longer PTHrP range → boundary ≈ 0.37); `blomstrand_like()` removes PTHrP production (H ≈ 1 everywhere).

## 3. Boolean chondrocyte fate network (`chondrotwin.boolean`)

Internal nodes: Piezo1, TRPV4, CaHigh, YAP, NFkB, HIF2a, SOX9, RUNX2, Senescent, Ferroptotic.
Inputs (clamped per scenario): LoadPhys, LoadOver, IL1b, WNT, TGFb, Stiff, Nrf2 (default 1), PTHrP.

```
Piezo1      = LoadOver ∨ (LoadPhys ∧ IL1b)                     # inflammatory sensitization
TRPV4       = LoadPhys ∧ ¬LoadOver
CaHigh      = Piezo1
YAP         = Stiff ∧ ¬IL1b                                     # cytokines degrade YAP/TAZ
NFkB        = (IL1b ∨ CaHigh ∨ Senescent) ∧ ¬YAP                # YAP ⊣ NF-κB
HIF2a       = NFkB
SOX9        = ¬RUNX2 ∧ ¬NFkB ∧ ¬YAP ∧ ¬Ferroptotic ∧ (TGFb ∨ PTHrP ∨ ¬WNT)   # PTHrP/PKA activates SOX9
RUNX2       = ¬SOX9 ∧ ¬PTHrP ∧ ¬Ferroptotic ∧ (WNT ∨ HIF2a ∨ CaHigh ∨ RUNX2)
Senescent   = Senescent ∨ (NFkB ∧ CaHigh ∧ ¬Ferroptotic)        # absorbing
Ferroptotic = Ferroptotic ∨ (CaHigh ∧ ¬Nrf2)                    # absorbing
```

Synchronous updating. `attractors(scenario, initial="slow")` counts basins over the 2³ = 8 initial states of the fate-memory nodes (SOX9, RUNX2, YAP) with all fast nodes and absorbing flags off; `initial="viable"` (256 states) and `initial="all"` (1024 states) are also available. `classify_state` maps a state to the fate labels of Table 1 (articular, hypertrophic, dedifferentiated, proliferative/PTHrP-held, senescent, ferroptotic, transit); "proliferative (PTHrP-held)" is SOX9 on, RUNX2 off with PTHrP present, as in the proliferative zone of the growth plate (PTHrP signaling phosphorylates and activates SOX9 through PKA; Huang et al., PNAS 2001).

Sources for the wiring: Piezo1/2 in injurious loading (Lee et al., PNAS 2014) and inflammatory sensitization of Piezo1 by IL-1α via p38 MAPK, HNF4 and ATF2/CREBP1 (Lee et al., PNAS 2021); TRPV4 in physiological loading (O'Conor et al., PNAS 2014); reciprocal YAP/TAZ–NF-κB inhibition and cytokine-induced YAP/TAZ degradation (Deng et al., Nat Commun 2018); HIF-2α as a catabolic/hypertrophic regulator (Saito et al., Nat Med 2010; Yang et al., Nat Med 2010); overload-induced, GPX4-regulated ferroptosis via Piezo1-mediated Ca²⁺ influx (Wang et al., J Adv Res 2022).

## 4. Quasi-potential landscape (`chondrotwin.landscape`)

Langevin version of the toggle switch (Euler–Maruyama, additive noise σ = 0.25, dt = 0.01):

```
dS = f_S(S, R; λ) dt + σ dW_S ,   dR = f_R(S, R; λ) dt + σ dW_R ,   S, R ≥ 0
```

3000 walkers (half started near each attractor; `start_fraction_sox9` changes the split) are integrated for 6000 steps; positions after a 1500-step burn-in are histogrammed on a 60 × 60 grid of [0, 3]². The quasi-potential is U = −ln P, shifted so that min U = 0 (Wang et al., PNAS 2011).

*Stationarity.* At low noise few walkers cross between basins during the run (`crossed` in the output: 0 % at λ = 0, σ = 0.18; 1.5 % at σ = 0.25), so P is quasi-stationary within each basin but the split between basins is set by the initial allocation. At λ = 0 this split is exact by symmetry. For λ ≠ 0 the escape barrier from the SOX9-high basin is insensitive to the split (λ = 0.35: 3.2–3.5 for 50/50, 90/10 and 10/90 starts at σ = 0.25), whereas the relative depth of the two basins and the barrier from the RUNX2-high side are not (VALIDATION.md, §4). Barrier heights in units of −ln P scale roughly as 1/σ². `basin_depths` reports the minimum of U in each half-plane and the lowest U on the diagonal S = R (a proxy for the ridge); if the diagonal was never visited the ridge is a lower bound and is flagged.

The schematic landscape of Figure 1a is a sum of Gaussian wells at the positions of the fate states (`DEFAULT_BASINS`), U = 2.4 − Σ depth·exp(−r²/2w²).

## 5. Loading regimes (`chondrotwin.loading`)

Net response to loading magnitude m (dimensionless; 1 ≈ injurious strain):

```
R(m) = a_anab exp(−((m − m_opt)/w_opt)²) − a_over σ((m − m_over)/w_over) − a_disuse exp(−(m/w_disuse)²)
```

with m_opt = 0.45, w_opt = 0.2, m_over = 0.82, w_over = 0.05, w_disuse = 0.12, a = (1.0, 0.85, 0.4). `regime_map` gates the anabolic term by a log-normal frequency window centred at 1 Hz and adds a static-compression penalty. `catabolic_bias(strain)` maps a cell-level strain magnitude to the switch bias λ ∈ (−0.95, 0.95) with the same shape (c = 0.5, 0.9, 0.3).

## 6. Biphasic confined compression (`chondrotwin.biphasic`)

Cartilage layer of thickness h on an impermeable base (z = 0), loaded through a free-draining platen (z = h) by a step stress σ₀:

```
∂u/∂t = H_A k ∂²u/∂z²,   u(0, t) = 0,   H_A ∂u/∂z (h, t) = −σ₀
p(z, t) = σ₀ + H_A ∂u/∂z,   p(h, t) = 0
```

Defaults: H_A = 0.5 MPa, k = 10⁻¹⁵ m⁴ N⁻¹ s⁻¹, h = 1.5 mm, 51 nodes → τ = h²/(H_A k) = 4500 s; equilibrium strain −σ₀/H_A. Explicit finite differences with dt = 0.4 dz²/(H_A k). Verified against the analytical creep solution (Mow et al., 1980): surface displacement within 1 % and pore pressure within 4 % of σ₀ at 51 nodes, first-order convergence in dz (VALIDATION.md, §3).

## 7. Digital-twin loop (`chondrotwin.twin`)

Per loading cycle c with applied stress σ₀(c):

```
ε_tissue = σ₀ / H_A
ε_cell   = A_PCM · ε_tissue                      (A_PCM = 1.5)
λ        = catabolic_bias(ε_cell / ε_scale)      (ε_scale = 1/3)
(S, R)   ← integrate the toggle switch for T = 15 time units at bias λ
H_A      ← H_A · (1 + k_syn ⟨S⟩ − k_deg ⟨R⟩),  clipped to [0.05, 0.8] MPa   (k_syn = 0.006, k_deg = 0.012)
```

With the defaults, σ₀ = 0.05 MPa is physiological (ε_cell ≈ 0.15, m ≈ 0.45, the anabolic optimum) and σ₀ = 0.15 MPa is injurious. These are effective drained (solid-matrix) stresses chosen to give equilibrium tissue strains of 10 % and 30 %; they are not joint contact pressures (of the order of MPa, carried largely by interstitial fluid pressurization at physiological loading rates). An overload episode of 12 cycles switches the cell to RUNX2-high; after the load returns to 0.05 MPa the bias first becomes anabolic again (λ ≈ −0.14) but lies inside the hysteresis window, so the cell stays hypertrophic and the matrix keeps degrading; the softer matrix raises the strain until the bias turns catabolic again from cycle 26 — the in-silico vicious cycle of Figure `twin_loop.png`.

## 8. Multi-agent cartilage model (`chondrotwin.abm`) — one agent = one chondrocyte

**Agents.** Each chondrocyte is one agent with position (x, y, z; z = 0 at the articular surface), zone (superficial / middle / deep), toggle-switch state (S, R), territorial matrix content M ∈ [0.05, 1], antioxidant capacity nrf2 ~ Beta(4, 2), age, flags (senescent, dead with cause, cluster id). The population is a structure of arrays (float32/int8/bool, 62 bytes per slot; slots for 1.5 × the initial cells, ≈ 93 bytes per agent) so every rule is a vectorized NumPy operation.

**Environment.** A slab Lx × Ly × h (default 1 × 1 × 2 mm) on a voxel grid (10 × 10 × 20). Cells are seeded at depth-dependent densities with relative values 2.4 : 1.0 : 0.77 (superficial / middle / deep), scaled to a mean of 1.4e4 mm⁻³ → 28,000 agents per slab. The mean is at the upper end of reported adult human values (Hunziker et al., Osteoarthritis Cartilage 2002: 9.6e3 mm⁻³ full-depth, medial femoral condyle); set `cells_per_mm3` to change it.

**Mechanics (per day).** Applied stress σ(x, y) (optionally with a Gaussian focal factor) is uniform through each column; the aggregate modulus increases with depth, H_A(z) = 0.4 → 1.2 MPa, and is softened by matrix loss, H_eff = H_A(z) · max(M_voxel, 0.3)^1.5, where M_voxel is the mean territorial matrix of the live cells in the voxel (empty voxels decay at 2 %/day). Tissue strain ε = min(σ / H_eff, 0.6); cell strain = 1.5 ε (pericellular amplification). Matrix loss therefore raises local strain.

**Cytokine field (quasi-steady).** D ∇²C − k C + P = 0 with D = 0.5 mm² day⁻¹, k = 5 day⁻¹, Dirichlet C = C_syn at the surface (synovial fluid), zero flux elsewhere; P from senescent cells (rate 1) and hypertrophic cells (0.25), normalized so that a fully senescent voxel at the reference density gives C = 1. Solved once per day with a pre-factorized sparse operator (cytokine kinetics are fast relative to cell-state kinetics).

**Fate dynamics.** λ_i = catabolic_bias(cell strain / 0.40) + 0.6 · C/(0.5 + C), clipped to (−0.95, 0.95); the toggle switch of Section 1 is integrated for one day (4 Euler–Maruyama substeps, Δt = 0.25, noise 0.10) per cell with its own λ_i. Senescent cells are frozen.

**Matrix turnover.** dM = 0.008 · S · (1 − M) · (0.2 if senescent) − 0.008 · R · M · (1 + 2 · C/(0.5 + C)) per day.

**Hazards (per day).** With overload o = sigmoid((m − 0.8)/0.08), m = cell strain / 0.40: senescence 3e-5 + 0.012 · o · (1 + 2 C/(0.5 + C)); ferroptosis 0.006 · o · (1 − nrf2); apoptosis 0.003 · o · (0.2 + C/(0.5 + C)). Dead cells keep their slot (cell loss is reported).

**Stage-2 options** (`Rates`; defaults give the rules above exactly, see `docs/knee_stage2.md`). Senescence hazard
× D / (D + `damage_D50`), D = accumulated overload-days per cell; strain bias `catabolic_bias(m, c_over = c_over_max)`;
degradation 0.008 · (R + `k_cat` · o) · M · (1 + 2 C/(0.5 + C)). Calibrated on GSE104782 patients OA1–5:
`damage_D50` = 0, `c_over_max` = 0.15, `k_cat` = 0.5 (`Rates.stage2(0, 0.15, 0.5)`). Stage 6 (`docs/knee_stage6.md`)
showed that only `c_over_max` is identified (95 % region [0, 0.2]); `k_cat` is not (the region includes 0).

**Stage-3 options** (not selected by calibration; defaults 0 / 1 keep the model unchanged): `switch_bias_sd` — a fixed
per-cell offset δ ~ N(0, sd) added to λ; `senescence_competent_fraction` — only this fraction of cells can become
senescent. Both are drawn from a separate random stream, so the main simulation stream is unaffected.

**Stage-4 options** (not selected by calibration): `sen_overload_weight` w — the overload-driven senescence hazard is
multiplied by w + (1 − w) · C/(0.5 + C); `senescent_death_rate` — daily apoptosis hazard of senescent cells (an extra
random draw only when > 0).

**Stage-5 options** (not selected): `Scenario.C_joint` — joint-level cytokine background added to C throughout the
tissue; `Rates.k_inf` — catabolic drive += k_inf · C/(0.5 + C); `Rates.sen_infl_gain` — the factor 2 in the senescence
hazard's (1 + 2 C/(0.5 + C)), now a parameter.

**Proliferation.** Live, non-senescent cells with M < 0.6 divide with probability 0.003 per day (cluster size ≤ 8); daughters are placed within ≈ 20 µm and inherit the parent state → chondrocyte clusters in degraded matrix.

**Outputs.** Fate fractions (articular, hypertrophic, senescent, dead by ferroptosis / apoptosis, relative to the initial cell number), matrix and strain by zone, mean cytokine, number of clusters, and an OARSI-like score 0–6 (surface matrix loss + cell loss + clusters). Scenarios: physiological (σ = 0.05 MPa), overload (0.12 MPa), overload episode, focal overload, inflamed joint (C_syn = 1), senescent seeding.

**Mechanical stimuli.** As in Section 7, the applied stresses (0.05 MPa physiological, 0.12 MPa injurious) are effective drained stresses giving equilibrium strains of ~10–30 %, not joint contact pressures.

**Robustness and attribution.** Over 5 random seeds the SD of every reported fraction and zonal matrix content is < 0.01. Ablations (VALIDATION.md, §8): surface-first degeneration disappears with a uniform modulus (it follows from H_A(z)); after an overload episode the superficial matrix partially recovers without the matrix → stiffness feedback (it drives the progression), and still progresses with a monostable switch (n = 1); seeded senescent cells do not induce extra senescence, and their effect on the matrix is small (10 % seeding: −0.014 superficial, −0.021 deep), of which the paracrine SASP part is ~1/3 (superficial) to ~2/3 (deep, further from the cytokine-clearing surface); the rest is their own reduced synthesis.

**Scaling.** Time per step is linear in the number of agents for a fixed voxel resolution (≈ 20–35 ms per 28,000 agents per day on a laptop core, depending on load). The agent state needs ≈ 93 B per agent and one vectorized step ≈ 260 B per agent of temporaries. A whole human knee (~2e4 mm³ of cartilage, ~2.8e8 cells at the default density) would need ~100 GB at peak and ~1 day per simulated year with this single-process NumPy implementation (excluding the cytokine solve on a joint-sized grid) — a job for a chunked or GPU implementation; the intended use is a representative volume driven by a joint-scale strain field (Figure 7 of the paper); `tile_plan()` and `scaling_table()` give the numbers.
