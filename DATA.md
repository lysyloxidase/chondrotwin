# Open data used by the knee pipeline

The knee simulation (`chondrotwin.knee`, `scripts/knee_open_data.py`, `scripts/knee_stage2*.py`, `scripts/knee_stage3*.py`, `scripts/knee_stage4*.py`, `scripts/knee_stage5*.py`) uses ten open datasets and two gene sets (~400 MB). **None of them is
redistributed in this repository.** `python scripts/fetch_data.py` downloads them from the original providers into
`data/raw/` (or `$CHONDROTWIN_DATA`) and checks the size and SHA-256 of every file against the values pinned in
`chondrotwin/datasets.py`. Derived tables are written to `data/derived/`; both folders are ignored by git.

| Source | Files (size) | What is used | License / terms | Cite |
|---|---|---|---|---|
| **OrthoLoad** — standardized knee loads | `AverageSubject_AllLoadLevels_Data_Diagrams.zip` (4.7 MB), `K1L…K9L_Walking.xlsx` (8 × ~0.8 MB) | 6-component tibial loads during level walking in 8 subjects with instrumented knee replacements → medial / lateral force split | non-commercial use; cite the article **and** the OrthoLoad database | Bergmann G, et al. Standardized loads acting in knee implants. *PLoS ONE* 2014;9:e86035; OrthoLoad, orthoload.com |
| **GEO GSE104782** — scRNA-seq of OA cartilage | `GSE104782_allcells_UMI_count.txt.gz` (10.2 MB), cell table `.xlsx` (0.1 MB), series matrix (3 kB) | 1600 chondrocytes, 10 patients × 5 degeneration stages; published clusters; marker transcripts | NCBI GEO public data | Ji Q, et al. Single-cell RNA-seq analysis reveals the progression of human osteoarthritis. *Ann Rheum Dis* 2019;78:100–110 |
| **GEO GSE152805** — scRNA-seq of OA cartilage, independent test | 6 cartilage samples (3 patients × damaged medial / minimally damaged lateral tibia): `GSM4626766…771_*.matrix.mtx.gz` (31–48 MB each), `barcodes.tsv.gz`, `genes.tsv.gz` (≈ 228 MB in total) | pre-registered stage-2 test (`docs/preregistration_gse152805.md`) | NCBI GEO public data | Chou CH, et al. Synovial cell cross-talk with cartilage plays a major role in the pathogenesis of osteoarthritis. *Sci Rep* 2020;10:10868 |
| **GEO GSE57218** — RAAK, paired preserved / OA-affected cartilage | `GSE57218_series_matrix.txt.gz` (18.7 MB), platform annotation `GPL6947.annot.gz` (7.0 MB) | pre-registered stage-3 bulk test, 33 joint-matched pairs | NCBI GEO public data | Ramos YFM, et al. Genes involved in the osteoarthritis process identified through genome wide expression analysis in articular cartilage; the RAAK study. *PLoS ONE* 2014;9:e103056 |
| **GEO GSE114007** — normal vs OA cartilage, RNA-seq | `GSE114007_normal_normalized.counts.txt.gz` (1.5 MB), `GSE114007_OA_normalized.counts.txt.gz` (1.7 MB), two series matrices (4 kB) | secondary stage-3 test — **declared invalid**: the "normalized counts" are batch-corrected values with negative entries; repaired in stage 4 with `GSE114007_raw_counts.xlsx` (4.1 MB), reported as a deviation | NCBI GEO public data | Fisch KM, et al. *Osteoarthritis Cartilage* 2018;26:1531–1538 |
| **GEO GSE117999** — OA vs meniscus-tear knee cartilage (Agilent) | `GSE117999_normalized_matrix.xlsx` (19.8 MB), series matrix (7.1 MB) | pre-registered stage-4 test (primary) | NCBI GEO public data | Rai MF, et al. *Osteoarthritis Cartilage* 2019;27:945–955 |
| **GEO GSE169077** — OA vs normal knee cartilage pools (Affymetrix U133A) | series matrix (1.0 MB), `GPL96.annot.gz` (4.5 MB) | pre-registered stage-4 test (secondary) | NCBI GEO public data | GEO GSE169077 (Attur M, Abramson SB; no primary publication linked) |
| **SenMayo** gene set (MSigDB `SAUL_SEN_MAYO`) | `SAUL_SEN_MAYO.json` (1.3 kB) | 124-gene senescence / SASP readout (stage 4) | CC BY 4.0 (article); MSigDB terms | Saul D, et al. *Nat Commun* 2022;13:4827 |
| **GEO GSE55235 + GSE55457** — synovium, healthy / OA / RA (Affymetrix U133A) | two series matrices (1.4 + 1.7 MB) | pre-registered stage-5 premise test (is OA synovium inflamed?); the control group includes joint-trauma patients | NCBI GEO public data | Woetzel D, et al. *Arthritis Res Ther* 2014;16:R84 |
| **MSigDB HALLMARK_INFLAMMATORY_RESPONSE** | `HALLMARK_INFLAMMATORY_RESPONSE.json` (1.9 kB) | inflammation score (stage 5) | CC BY 4.0 (article); MSigDB terms | Liberzon A, et al. *Cell Syst* 2015;1:417–425 |
| **GEO GSE178557** — OA vs fracture-control cartilage (Agilent) | 8 sample files (17.6 MB) + series matrix (1.1 MB) | **downloaded, not used**: no gene annotation short of a 1.7 GB platform file | NCBI GEO public data | GEO GSE178557 |
| **OAIZIB-CM** (Zenodo 10.5281/zenodo.14934086) | `labelsTr.zip` (57.5 MB), `labelsTs.zip` (14.1 MB), `info.zip` (52 kB) — segmentation masks only, no MR images | medial / lateral tibial cartilage of 507 OAI knees; KL grade, sex, age, BMI; the dataset's train / test split | **CC BY-NC 4.0**; citing CartiMorph is compulsory | Yao Y, et al. CartiMorph. *Med Image Anal* 2024;91:103035; Ambellan F, et al. *Med Image Anal* 2019;52:109–118; Osteoarthritis Initiative |

Also used, from the published text only (no data file): the mean medial force peaks of the OrthoLoad cohort
(Kutzner I, et al. Knee adduction moment and medial contact force — facts about their correlation during gait.
*PLoS ONE* 2013;8:e81036) to calibrate the condylar distance.

Because two of the sources are non-commercial, results derived from them inherit that restriction; the code in this
repository remains MIT-licensed.

Not used (and why): the Grand Challenge knee-load data (SimTK login required), the OAI images and longitudinal
quantitative MRI (NIH NDA data use agreement required) — the latter is the main target for stage 2.
