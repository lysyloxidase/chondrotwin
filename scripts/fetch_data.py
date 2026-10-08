#!/usr/bin/env python
"""Download the open datasets used by the knee pipeline into data/raw (or $CHONDROTWIN_DATA).

Usage:  python scripts/fetch_data.py [orthoload] [gse104782] [oaizib_cm]

Licenses and required citations are listed in chondrotwin/datasets.py and DATA.md.
"""
import sys

from chondrotwin.datasets import SOURCES, fetch

if __name__ == "__main__":
    names = sys.argv[1:] or list(SOURCES)
    for n in names:
        print(f"[{n}] {SOURCES[n].license}")
    fetch(names)
