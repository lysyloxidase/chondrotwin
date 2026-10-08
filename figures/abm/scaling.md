Benchmark on this machine: 28000 agents, 31.3 ms per step

| tissue | agents | s / step (this machine) | hours / 365 steps | RAM state / peak |
|---|---|---|---|---|
| 1 mm^3 (RVE) | 1.40e+04 | 0.0156 | 0.00159 | 0.00112 / 0.00476 GB |
| 1 x 1 x 2 mm slab | 2.80e+04 | 0.0313 | 0.00317 | 0.00224 / 0.00952 GB |
| 1 cm^2 x 2 mm patch | 2.80e+06 | 3.13 | 0.317 | 0.224 / 0.952 GB |
| medial tibial plateau (~2.5e3 mm^3) | 3.50e+07 | 39.1 | 3.96 | 2.8 / 11.9 GB |
| whole knee (~2e4 mm^3) | 2.80e+08 | 313 | 31.7 | 22.4 / 95.2 GB |

Tiling plan for a whole knee: {'tissue_volume_mm3': 20000.0, 'rve_mm3': 2.0, 'n_rve': 10000, 'cells_per_rve': 28000, 'total_cells': 280000000}
