"""Small helpers for the generated Markdown reports of the knee pipeline."""
from __future__ import annotations

import numpy as np
import pandas as pd


def table(df: pd.DataFrame, floatfmt: str = "{:.3f}") -> str:
    """Markdown table; itertuples keeps each column's dtype, so integer columns print as integers."""
    cols = [df.index.name or ""] + [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for row in df.itertuples():
        idx, vals = row[0], row[1:]
        cells = [str(idx)] + [floatfmt.format(v) if isinstance(v, (float, np.floating)) else str(v) for v in vals]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)
