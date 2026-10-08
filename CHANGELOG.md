# Changelog

## 0.2.0 (2026-10-08)
- Multi-agent cartilage model (`chondrotwin.abm`): one agent = one chondrocyte; depth-zoned slab, matrix-feedback
  mechanics, quasi-steady cytokine field, stochastic toggle switch per cell, senescence / ferroptosis / apoptosis hazards,
  proliferative clusters, OARSI-like score; scenarios (physiological, overload, overload episode, focal, inflamed,
  senescent seeding); scaling helpers and benchmark for "whole cartilage" estimates; `scripts/run_abm.py`; tests.

## 0.1.0 (2026-10-08)
- Initial release accompanying the manuscript: toggle switch, growth-plate reaction-diffusion model,
  Boolean fate network, Langevin quasi-potential landscape, biphasic confined compression,
  digital-twin loop, figure scripts and tests.
