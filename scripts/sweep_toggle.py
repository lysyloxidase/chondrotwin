#!/usr/bin/env python
"""Parameter sweep: how the bistable (hysteresis) window of the SOX9-RUNX2 switch depends on cooperativity n and K.

Usage:  python scripts/sweep_toggle.py
"""
import numpy as np
from chondrotwin import ToggleSwitch

if __name__ == "__main__":
    print("n     K     bistable range of lambda")
    for n in (1.5, 2.0, 3.0, 4.0, 6.0):
        for K in (0.3, 0.5, 0.8):
            rng = ToggleSwitch(n=n, K=K).bistable_range(np.linspace(-0.95, 0.95, 191))
            print(f"{n:<5} {K:<5} {'none' if rng is None else f'[{rng[0]:+.2f}, {rng[1]:+.2f}]'}")
