#!/usr/bin/env python
"""Regenerate the computed manuscript figures and the diagnostic figures.

Usage:  python scripts/make_figures.py [OUTDIR]     (default: figures/)
"""
import sys
from chondrotwin.plotting import make_all_figures

if __name__ == "__main__":
    outdir = sys.argv[1] if len(sys.argv) > 1 else "figures"
    for name, path in make_all_figures(outdir).items():
        print(f"{name:24s} -> {path}")
