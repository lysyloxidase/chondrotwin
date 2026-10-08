# Pre-registration: stage-4 independent test (GSE117999 primary, GSE169077 secondary; SenMayo readout)

The protocol is the module docstring and code of `chondrotwin/knee/independent3.py`.

* Written and saved: **2026-10-08T14:20:11Z** (UTC).
* SHA-256 of `chondrotwin/knee/independent3.py` at registration: `790b58e3131d93f2cf8684011700e6ea9ec60e09ee9c175b93bced790a4dd4a8`
* Seen before registration (to write a parser): the SenMayo gene list; the column names, sheet name and identifier
  columns of `GSE117999_normalized_matrix.xlsx`; the sample metadata lines and table headers of both series
  matrices; the GPL96 annotation rows of the four markers. **No expression value was read.**
* Order: registration → exploratory re-measurement of senescence on the previously seen datasets (GSE104782,
  GSE152805, GSE57218) → stage-4 model revisions calibrated on GSE104782 OA1–5 → parameters frozen in
  `docs/stage4_frozen_parameters.json` (hash appended below at freezing) → only then the expression values of
  GSE117999 / GSE169077 are read. Reported whatever the result.

Limitation: local registration without an external time stamp.

## Frozen parameters

* Frozen at **2026-10-08T14:27:32Z** (UTC), before any expression value of GSE117999 / GSE169077 was read.
* SHA-256 of `docs/stage4_frozen_parameters.json`: `0205c1e81fd0b5e89a3819ba080e858934b4ca1de36edb56ad215adfbc81dd4d`
