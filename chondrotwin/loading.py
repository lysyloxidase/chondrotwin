"""
Loading regimes (paper, Section 4.7 and Figure 6).

Phenomenological functions describing the inverted-U relationship between
loading magnitude and the net anabolic/catabolic response of chondrocytes, the
magnitude x frequency "loading-regime map", and the mapping from cell-level
strain to the catabolic bias ``lam`` that drives the toggle switch in the
digital-twin loop.  Shapes are illustrative; all magnitudes are dimensionless
(roughly, 1.0 ~ injurious strain).
"""
from __future__ import annotations

import numpy as np


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def inverted_u(m, m_opt: float = 0.45, w_opt: float = 0.2, m_over: float = 0.82, w_over: float = 0.05,
               disuse_w: float = 0.12, a_anab: float = 1.0, a_over: float = 0.85, a_disuse: float = 0.4):
    """Net anabolic (+) / catabolic (-) response versus loading magnitude ``m``."""
    m = np.asarray(m, dtype=float)
    anabolic = a_anab * np.exp(-((m - m_opt) / w_opt) ** 2)
    overload = a_over * _sigmoid((m - m_over) / w_over)
    disuse = a_disuse * np.exp(-(m / disuse_w) ** 2)
    return anabolic - overload - disuse


def regime_map(M, F, f_opt: float = 1.0, f_width: float = 0.8, static_w: float = 0.15, a_static: float = 0.35, **kw):
    """Net response over a magnitude (M) x frequency (F, Hz) grid.

    Dynamic loading near ``f_opt`` is required for the anabolic response; very
    low frequencies (static compression) are catabolic above the disuse range.
    """
    M = np.asarray(M, dtype=float); F = np.asarray(F, dtype=float)
    freq_window = np.exp(-((np.log(F + 0.05) - np.log(f_opt)) / f_width) ** 2)
    base = inverted_u(M, **kw)
    anabolic_part = kw.get("a_anab", 1.0) * np.exp(-((M - kw.get("m_opt", 0.45)) / kw.get("w_opt", 0.2)) ** 2)
    # replace the anabolic term by its frequency-gated version
    Z = base - anabolic_part + anabolic_part * freq_window
    Z -= a_static * np.exp(-(F / static_w) ** 2) * (M > 0.15)
    return Z


def catabolic_bias(strain, m_opt: float = 0.45, w_opt: float = 0.2, m_over: float = 0.82, w_over: float = 0.05,
                   disuse_w: float = 0.12, c_anab: float = 0.5, c_over: float = 0.9, c_disuse: float = 0.3,
                   clip: float = 0.95):
    """Map cell-level strain magnitude to the toggle-switch bias lam in (-clip, clip).

    Physiological strains give a negative (anabolic) bias, overload a positive
    (catabolic) bias, disuse a mild positive bias.
    """
    m = np.abs(np.asarray(strain, dtype=float))
    lam = (c_over * _sigmoid((m - m_over) / w_over)
           - c_anab * np.exp(-((m - m_opt) / w_opt) ** 2)
           + c_disuse * np.exp(-(m / disuse_w) ** 2))
    return np.clip(lam, -clip, clip)
