Benchmark on this machine: 28000 agents, 19.4 ms per step

| tissue | agents | s / step (this machine) | hours / 365 steps | RAM |
|---|---|---|---|---|
| 1 mm^3 (RVE) | 1.40e+04 | 0.00968 | 0.000982 | 0.00084 GB |
| 1 x 1 x 2 mm slab | 2.80e+04 | 0.0194 | 0.00196 | 0.00168 GB |
| 1 cm^2 x 2 mm patch | 2.80e+06 | 1.94 | 0.196 | 0.168 GB |
| medial tibial plateau (~2.5e3 mm^3) | 3.50e+07 | 24.2 | 2.45 | 2.1 GB |
| whole knee (~2e4 mm^3) | 2.80e+08 | 194 | 19.6 | 16.8 GB |

Tiling plan for a whole knee: {'tissue_volume_mm3': 20000.0, 'rve_mm3': 2.0, 'n_rve': 10000, 'cells_per_rve': 28000, 'total_cells': 280000000}
