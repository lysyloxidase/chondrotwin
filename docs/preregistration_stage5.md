# Pre-registration: stage-5 premise test (synovial inflammation in OA; GSE55235 + GSE55457)

The protocol is the module docstring and code of `chondrotwin/knee/independent4.py`.

* Written and saved: **2026-10-08T14:44:05Z** (UTC).
* SHA-256 of `chondrotwin/knee/independent4.py` at registration: `d9fb31bc6fca15ef8674a81b024085877797f4953f4f76115d74979a99df2735`
* Seen before registration: the hallmark gene list, the sample metadata lines and table headers of both series
  matrices. **No expression value was read.**
* Scope: this tests the *premise* of the stage-5 model (OA joints carry a joint-level inflammatory input), not its
  numbers. The cartilage model itself is calibrated and checked only on previously seen data in stage 5 (no fresh
  cartilage dataset could be used: GSE178557 was downloaded but has no usable gene annotation short of a 1.7 GB file).
* Order: registration → stage-5 model calibrated on GSE104782 OA1–5 + GSE114007 raw counts → parameters frozen in
  `docs/stage5_frozen_parameters.json` (hash appended below) → only then the synovial expression values are read.

Limitation: local registration without an external time stamp.

## Frozen parameters

* Frozen at **2026-10-08T14:50:11Z** (UTC), before any expression value of GSE55235 / GSE55457 was read.
* SHA-256 of `docs/stage5_frozen_parameters.json`: `e418804508775819f216afbd479fc2f12939a07738ec17679fa89ac087a50900`
