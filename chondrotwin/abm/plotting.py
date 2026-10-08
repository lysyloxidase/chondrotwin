"""Figures for the agent-based cartilage model."""
from __future__ import annotations

from typing import Sequence

import numpy as np

from .agents import FATE_LABELS
from .simulation import SimulationResult

FATE_COLORS = {"articular": "#1d4ed8", "hypertrophic": "#f59e0b", "senescent": "#7c3aed",
               "dead_ferroptosis": "#dc2626", "dead_apoptosis": "#6b7280"}


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7,
                         "ytick.labelsize": 7, "legend.fontsize": 7})
    return plt


def figure_cross_sections(results: Sequence[SimulationResult], path: str | None = None, y_slice: float = 0.15):
    """x-z cross-section of every scenario (cells within +-y_slice mm of the mid-plane), colored by fate, plus matrix."""
    plt = _plt()
    n = len(results)
    fig, axes = plt.subplots(2, n, figsize=(3.2 * n, 6.0), squeeze=False)
    for j, res in enumerate(results):
        c = res.cells; s = res.slab
        m = np.abs(c["y"] - s.Ly / 2) < y_slice
        ax = axes[0, j]
        for lab in FATE_LABELS:
            mm = m & (c["fate"] == lab)
            ax.scatter(c.loc[mm, "x"], c.loc[mm, "z"], s=5, color=FATE_COLORS[lab], label=lab if j == 0 else None, alpha=0.8, linewidths=0)
        ax.set_ylim(s.h, 0); ax.set_xlim(0, s.Lx); ax.set_title(f"{res.name}: fates (day {len(res.history)})", loc="left")
        ax.set_xlabel("x (mm)"); ax.set_ylabel("depth (mm), 0 = surface")
        for zb in (s.zone_bounds[0] * s.h, s.zone_bounds[1] * s.h):
            ax.axhline(zb, color="gray", lw=0.5, ls=":")
        ax2 = axes[1, j]
        alive = m & ~c["dead"]
        sc = ax2.scatter(c.loc[alive, "x"], c.loc[alive, "z"], s=5, c=c.loc[alive, "M"], cmap="YlGn", vmin=0, vmax=1, linewidths=0)
        ax2.set_ylim(s.h, 0); ax2.set_xlim(0, s.Lx); ax2.set_title("territorial matrix M", loc="left"); ax2.set_xlabel("x (mm)")
        if j == n - 1:
            fig.colorbar(sc, ax=ax2, fraction=0.05).set_label("M")
    axes[0, 0].legend(loc="lower right", fontsize=6, markerscale=2)
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig


def figure_time_courses(results: Sequence[SimulationResult], path: str | None = None):
    plt = _plt()
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 5.6))
    pal = plt.cm.tab10(np.linspace(0, 1, max(len(results), 3)))
    for res, col in zip(results, pal):
        f = res.frame()
        axes[0, 0].plot(f.index, f["articular"], color=col, lw=1.5, label=res.name)
        axes[0, 1].plot(f.index, f["hypertrophic"], color=col, lw=1.5)
        axes[0, 2].plot(f.index, f["senescent"], color=col, lw=1.5)
        axes[1, 0].plot(f.index, f["dead_ferroptosis"] + f["dead_apoptosis"], color=col, lw=1.5)
        axes[1, 1].plot(f.index, f["matrix_superficial"], color=col, lw=1.5, label=f"{res.name} superficial")
        axes[1, 1].plot(f.index, f["matrix_deep"], color=col, lw=1.0, ls="--")
        axes[1, 2].plot(f.index, f["oarsi_like"], color=col, lw=1.5)
    titles = ["articular fraction", "hypertrophic fraction", "senescent fraction", "dead fraction (ferroptosis + apoptosis)",
              "matrix content (solid: superficial, dashed: deep)", "OARSI-like score (0-6)"]
    for ax, t in zip(axes.ravel(), titles):
        ax.set_title(t, loc="left"); ax.set_xlabel("day")
    axes[0, 0].legend(fontsize=6); fig.tight_layout()
    if path:
        fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig


def figure_depth_profiles(results: Sequence[SimulationResult], path: str | None = None, n_bins: int = 20):
    """Final depth profiles of matrix, cell strain, cytokine and fate composition."""
    plt = _plt()
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.2))
    pal = plt.cm.tab10(np.linspace(0, 1, max(len(results), 3)))
    for res, col in zip(results, pal):
        c = res.cells; s = res.slab
        bins = np.linspace(0, s.h, n_bins + 1); zc = 0.5 * (bins[:-1] + bins[1:])
        idx = np.clip(np.digitize(c["z"], bins) - 1, 0, n_bins - 1)
        alive = ~c["dead"].values
        def prof(v, mask):
            out = np.full(n_bins, np.nan)
            for b in range(n_bins):
                mm = mask & (idx == b)
                if mm.any(): out[b] = float(np.mean(v[mm]))
            return out
        axes[0].plot(prof(c["M"].values, alive), zc, color=col, lw=1.5, label=res.name)
        axes[1].plot(prof(c["strain"].values, alive), zc, color=col, lw=1.5)
        axes[2].plot(prof(c["C"].values, alive), zc, color=col, lw=1.5)
        axes[3].plot(prof((c["dead"].values | c["senescent"].values).astype(float), np.ones(len(c), bool)), zc, color=col, lw=1.5)
    for ax, t in zip(axes, ["matrix M", "cell strain", "cytokine C", "fraction senescent or dead"]):
        ax.set_ylim(results[0].slab.h, 0); ax.set_title(t, loc="left"); ax.set_xlabel(t)
    axes[0].set_ylabel("depth (mm)"); axes[0].legend(fontsize=6); fig.tight_layout()
    if path:
        fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig
