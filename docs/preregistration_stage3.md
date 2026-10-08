# Pre-registration: stage-3 independent test on bulk transcriptomes (GSE57218, GSE114007)

The protocol (datasets, preprocessing, statistics, model counterparts, criteria C1/C2) is the module docstring and
code of `chondrotwin/knee/independent2.py`.

* Written and saved: **2026-10-08T13:37:06Z** (UTC).
* SHA-256 of `chondrotwin/knee/independent2.py` at registration: `ec648c11f77b71619797540e08e445d546620e95eb8d143d943ac3c867158f1d`
* What had been seen before registration (to write a parser): the sample metadata lines of the GSE57218 series
  matrix (titles, accessions, sex, age, disease state, processing description) and its table header; the GPL6947
  probe annotation for SOX9, RUNX2, MMP13, CDKN2A; the header line and the first five gene symbols of the two
  GSE114007 count files and the sample metadata of its series matrices. **No expression value was read.**
* Order of work: registration → implementation and calibration of the stage-3 model on GSE104782 patients OA1–5
  → parameters frozen in `docs/stage3_frozen_parameters.json` (with its own SHA-256 below, added at freezing) →
  only then are the expression values read. The result is reported whatever it is.

Limitation: local registration without an external time stamp; the git history of the release is the only
independent record.

## Frozen parameters

* Frozen at **2026-10-08T13:45:00Z** (UTC), before any expression value of GSE57218 / GSE114007 was read.
* SHA-256 of `docs/stage3_frozen_parameters.json`: `f3064df284136bc0417bf964e27536a8a79bb723cbbdf476ba56548d09ad321b`
