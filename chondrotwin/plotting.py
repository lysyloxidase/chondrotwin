"""
Figure generation from the package models.

``make_all_figures(outdir)`` reproduces the two *computed* figures of the paper
(Figure 2: growth plate; Figure 4: toggle switch) and produces additional
figures that document the other modules (Langevin landscape, loading-regime
map, Boolean attractor map, biphasic creep, digital-twin loop).

Matplotlib is imported lazily so that the numerical modules stay dependency-light.
"""
from __future__ import annotations

import os
from typing import Dict, List

import numpy as np

from .biphasic import ConfinedCompression
from .boolean import SCENARIOS, attractor_table, chondrocyte_network
from .growth_plate import GrowthPlateModel
from .landscape import basin_depths, quasi_potential_langevin, schematic_landscape, DEFAULT_BASINS
from .loading import inverted_u, regime_map
from .toggle import ToggleSwitch
from .twin import MiniTwin


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7})
    return plt


# ---------------------------------------------------------------- Figure 4
def figure_toggle(model: ToggleSwitch | None = None, lam_plane: float = 0.3, path: str | None = None):
    plt = _plt()
    model = model or ToggleSwitch()
    fig = plt.figure(figsize=(8.2, 3.2))
    gs = fig.add_gridspec(1, 3, wspace=0.38)
    # (a) phase plane
    ax = fig.add_subplot(gs[0])
    S, R = np.meshgrid(np.linspace(0, 3, 30), np.linspace(0, 3, 30))
    dS = model.alpha_s * (1 - lam_plane) * model.hill_rep(R) + model.beta_s - model.gamma * S
    dR = model.alpha_r * (1 + lam_plane) * model.hill_rep(S) + model.beta_r - model.gamma * R
    ax.streamplot(S, R, dS, dR, color="#9ca3af", density=0.9, linewidth=0.5, arrowsize=0.6)
    r = np.linspace(0, 3, 400)
    ax.plot(model.nullcline_S(r, lam_plane), r, color="#1d4ed8", lw=1.6, label="dSOX9/dt = 0")
    ax.plot(r, model.nullcline_R(r, lam_plane), color="#dc2626", lw=1.6, label="dRUNX2/dt = 0")
    for fp in model.fixed_points(lam_plane):
        ax.plot(fp.S, fp.R, "o", ms=6, mfc="k" if fp.stable else "white", mec="k")
    ax.set_xlabel("SOX9 activity"); ax.set_ylabel("RUNX2 activity")
    ax.set_title(f"(a) phase plane (λ = {lam_plane})", loc="left"); ax.legend(fontsize=6, loc="upper right")
    ax.set_xlim(0, 3); ax.set_ylim(0, 3)
    # (b) bifurcation
    ax2 = fig.add_subplot(gs[1])
    lams = np.linspace(-0.95, 0.95, 300)
    pts = model.bifurcation(lams)
    st = [(l, fp.R) for l, fp in pts if fp.stable]; un = [(l, fp.R) for l, fp in pts if not fp.stable]
    if st:
        ax2.plot(*zip(*st), ".", color="#111827", ms=2, label="stable")
    if un:
        ax2.plot(*zip(*un), ".", color="#9ca3af", ms=1.5, label="unstable (saddle)")
    rng = model.bistable_range(lams)
    if rng:
        ax2.axvspan(rng[0], rng[1], color="#fde68a", alpha=0.5)
        ax2.text(sum(rng) / 2, 3.75, "bistable\n(hysteresis)", ha="center", fontsize=6)
    ax2.set_xlabel("λ: catabolic bias →   (← anabolic bias)"); ax2.set_ylabel("steady-state RUNX2")
    ax2.set_title("(b) bifurcation diagram", loc="left"); ax2.legend(fontsize=6, loc="center left"); ax2.set_ylim(0, 4.2)
    # (c) pulse
    ax3 = fig.add_subplot(gs[2])
    t, S_t, R_t, lam_t = model.pulse()
    ax3.fill_between(t, 0, 4.0, where=lam_t > 0, color="#fde68a", alpha=0.5, label="transient overload (λ = 0.8)")
    ax3.plot(t, S_t, color="#1d4ed8", lw=1.6, label="SOX9"); ax3.plot(t, R_t, color="#dc2626", lw=1.6, label="RUNX2")
    ax3.set_xlabel("time (a.u.)"); ax3.set_ylabel("activity"); ax3.set_title("(c) memory of a transient insult", loc="left")
    ax3.legend(fontsize=6, loc="center right"); ax3.set_ylim(0, 4.0)
    if path:
        fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig


# ---------------------------------------------------------------- Figure 2
def figure_growth_plate(model: GrowthPlateModel | None = None, path: str | None = None):
    plt = _plt()
    from matplotlib.patches import Circle, Ellipse, Rectangle
    model = model or GrowthPlateModel()
    fig = plt.figure(figsize=(8.2, 4.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.1, 1.1], wspace=0.35)
    ax = fig.add_subplot(gs[0])
    zones = [("Resting zone\n(SOX9; PTHrP-\nresponsive)", 0.80, 1.0, "#dbeafe"),
             ("Proliferative zone\n(SOX9, COL2A1;\ncolumns, IHH target)", 0.44, 0.80, "#bfdbfe"),
             ("Prehypertrophic\n(IHH, PTH1R)", 0.31, 0.44, "#fde68a"),
             ("Hypertrophic\n(RUNX2, COL10A1,\nMMP13, VEGF)", 0.08, 0.31, "#fdba74"),
             ("Bone / marrow\n(transdifferentiation;\nosteoclasts)", 0.0, 0.08, "#d1d5db")]
    for name, y0, y1, c in zones:
        ax.add_patch(Rectangle((0, y0), 1, y1 - y0, facecolor=c, edgecolor="k", lw=0.5))
        ax.text(1.06, (y0 + y1) / 2, name, va="center", fontsize=5.8)
    rng = np.random.default_rng(1)
    for _ in range(16):
        ax.add_patch(Circle((rng.uniform(0.1, 0.9), rng.uniform(0.84, 0.97)), 0.02, fc="white", ec="k", lw=0.4))
    for col in np.linspace(0.15, 0.85, 5):
        for y in np.linspace(0.47, 0.77, 7):
            ax.add_patch(Ellipse((col, y), 0.11, 0.03, fc="white", ec="k", lw=0.4))
        for y in (0.345, 0.405):
            ax.add_patch(Circle((col, y), 0.025, fc="white", ec="k", lw=0.4))
    for col in np.linspace(0.16, 0.84, 4):
        for y in (0.135, 0.245):
            ax.add_patch(Circle((col, y), 0.048, fc="white", ec="k", lw=0.4))
    ax.annotate("", xy=(-0.09, 0.33), xytext=(-0.09, 0.97), arrowprops=dict(arrowstyle="-|>", color="#1d4ed8", lw=1.6))
    ax.text(-0.16, 0.66, "PTHrP gradient", rotation=90, va="center", ha="center", color="#1d4ed8", fontsize=6.5)
    ax.annotate("", xy=(-0.36, 0.97), xytext=(-0.36, 0.33), arrowprops=dict(arrowstyle="-|>", color="#b45309", lw=1.6))
    ax.text(-0.44, 0.66, "IHH (feedback)", rotation=90, va="center", ha="center", color="#b45309", fontsize=6.5)
    ax.set_xlim(-0.55, 2.0); ax.set_ylim(-0.02, 1.02); ax.axis("off"); ax.set_title("(a) growth-plate zonation", loc="left")
    sol = model.steady_state()
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(sol.x, sol.P / sol.P.max(), color="#1d4ed8", lw=1.6, label="PTHrP (normalized)")
    ax2.plot(sol.x, sol.I / sol.I.max(), color="#b45309", lw=1.6, label="IHH (normalized)")
    ax2.plot(sol.x, sol.H, color="#dc2626", lw=1.6, ls="--", label="hypertrophic\ncommitment")
    xb = sol.boundary()
    ax2.axvline(xb, color="gray", lw=0.8, ls=":"); ax2.text(xb + 0.02, 0.05, "prehypertrophic\nboundary", fontsize=6, color="gray")
    ax2.set_xlabel("distance from periarticular end (a.u.)"); ax2.set_ylabel("level (a.u.)")
    ax2.set_title("(b) minimal reaction–diffusion model", loc="left")
    ax2.legend(loc="center right", bbox_to_anchor=(1.0, 0.33), fontsize=5.8)
    ax3 = fig.add_subplot(gs[2])
    for lab, m, c in [("wild type", model, "#111827"),
                      ("PTH1R gain of function\n(Jansen-like: longer PTHrP range)", model.jansen_like(), "#1d4ed8"),
                      ("PTH1R / PTHrP loss\n(Blomstrand-like)", model.blomstrand_like(), "#dc2626")]:
        s = m.steady_state(); ax3.plot(s.x, s.H, lw=1.6, color=c, label=lab)
    ax3.set_xlabel("distance from periarticular end (a.u.)"); ax3.set_ylabel("hypertrophic commitment")
    ax3.set_title("(c) perturbing the loop shifts the boundary", loc="left"); ax3.legend(fontsize=5.6, loc="center right")
    if path:
        fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig


# ---------------------------------------------------------------- extras
def figure_landscape(path: str | None = None, lam_values=(0.0, 0.35), noise: float = 0.18, seed: int = 0):
    """Schematic landscape (Figure 1a) next to computed Langevin quasi-potentials of the toggle switch."""
    plt = _plt()
    fig, axes = plt.subplots(1, 1 + len(lam_values), figsize=(3.6 * (1 + len(lam_values)), 3.4))
    X, Y, U = schematic_landscape()
    ax = axes[0]
    ax.contourf(X, Y, U, levels=26, cmap="viridis_r")
    for name, cx, cy, *_ in DEFAULT_BASINS:
        ax.text(cx, cy, name.replace(" / ", "\n"), ha="center", va="center", fontsize=4.8, color="white", fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title("(a) schematic landscape (Figure 1a)", loc="left")
    for ax, lam in zip(axes[1:], lam_values):
        res = quasi_potential_langevin(lam=lam, noise=noise, seed=seed)
        Sc = 0.5 * (res["S_edges"][:-1] + res["S_edges"][1:]); Rc = 0.5 * (res["R_edges"][:-1] + res["R_edges"][1:])
        Um = np.ma.masked_invalid(res["U"].T)
        cs = ax.contourf(Sc, Rc, Um, levels=20, cmap="viridis_r")
        fig.colorbar(cs, ax=ax, fraction=0.045, pad=0.02).set_label("U = −ln P", fontsize=7)
        d = basin_depths(res["U"], res["S_edges"], res["R_edges"])
        sym = "≥" if d["ridge_is_lower_bound"] else "≈"
        ax.set_title(f"(λ = {lam}) barrier SOX9→RUNX2 {sym} {d['barrier_from_SOX9']:.1f}", loc="left", fontsize=7.5)
        ax.set_xlabel("SOX9"); ax.set_ylabel("RUNX2")
    fig.suptitle("Computable landscape: quasi-potential of the toggle switch from Langevin sampling", fontsize=8)
    if path:
        fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig


def figure_loading(path: str | None = None):
    plt = _plt()
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.2), gridspec_kw=dict(width_ratios=[1, 1.25], wspace=0.3))
    m = np.linspace(0, 1.1, 400); r = inverted_u(m)
    ax.plot(m, r, color="#111827", lw=1.8); ax.axhline(0, color="gray", lw=0.6)
    ax.fill_between(m, r, 0, where=r > 0, color="#bbf7d0", alpha=0.7); ax.fill_between(m, r, 0, where=r < 0, color="#fecaca", alpha=0.7)
    ax.set_xlabel("loading magnitude (a.u.)"); ax.set_ylabel("net anabolic (+) / catabolic (−)"); ax.set_title("(a) inverted-U response", loc="left")
    M, F = np.meshgrid(np.linspace(0, 1.1, 300), np.linspace(0, 3, 300))
    Z = regime_map(M, F)
    cs = ax2.contourf(M, F, Z, levels=20, cmap="RdYlGn", vmin=-1, vmax=1); ax2.contour(M, F, Z, levels=[0], colors="k", linewidths=0.8)
    fig.colorbar(cs, ax=ax2, fraction=0.04, pad=0.02).set_label("net response", fontsize=7)
    ax2.set_xlabel("loading magnitude (a.u.)"); ax2.set_ylabel("frequency (Hz)"); ax2.set_title("(b) loading-regime map", loc="left")
    if path:
        fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig


def figure_boolean(path: str | None = None):
    """Attractor map: fate of the largest-basin attractor in every scenario."""
    plt = _plt()
    rows = attractor_table()
    scen = list(SCENARIOS)
    fates = ["articular (SOX9-high)", "hypertrophic (RUNX2-high)", "dedifferentiated (YAP-high)", "proliferative (PTHrP-held)",
             "senescent", "ferroptotic", "transit (saddle-like)"]
    mat = np.zeros((len(scen), len(fates)))
    for r in rows:
        i = scen.index(r["scenario"])
        for f in str(r["fate"]).split(" / "):
            j = fates.index(f) if f in fates else len(fates) - 1
            mat[i, j] += r["basin_fraction"]
    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    im = ax.imshow(mat, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(fates))); ax.set_xticklabels([f.split(" (")[0] for f in fates], rotation=25, ha="right")
    ax.set_yticks(range(len(scen))); ax.set_yticklabels(scen)
    for i in range(len(scen)):
        for j in range(len(fates)):
            if mat[i, j] > 0:
                ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=6.5, color="white" if mat[i, j] > 0.5 else "black")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02).set_label("basin fraction", fontsize=7)
    ax.set_title("Boolean chondrocyte network: attractor basins per scenario (initial states: prior fate of a viable cell)", loc="left")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig


def figure_biphasic(path: str | None = None, sigma0: float = 0.1e6):
    plt = _plt()
    tissue = ConfinedCompression()
    sol = tissue.solve(sigma0)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.2))
    ax.plot(sol.t / sol.tau, -sol.surface_displacement / tissue.h, color="#111827", lw=1.6)
    ax.axhline(sigma0 / tissue.H_A, color="#dc2626", ls="--", lw=1, label="equilibrium strain σ₀/H_A")
    ax.set_xlabel("time / τ"); ax.set_ylabel("surface compressive strain"); ax.set_title("(a) biphasic creep (confined compression)", loc="left"); ax.legend(fontsize=6)
    for frac in (0.02, 0.1, 0.3, 1.0, 3.0):
        i = int(np.argmin(np.abs(sol.t / sol.tau - frac)))
        ax2.plot(sol.pressure[i] / sigma0, sol.z / tissue.h, lw=1.4, label=f"t = {frac:g} τ")
    ax2.set_xlabel("fluid pressure / σ₀"); ax2.set_ylabel("depth (0 = bone, 1 = surface)"); ax2.set_title("(b) interstitial fluid pressure", loc="left"); ax2.legend(fontsize=6)
    if path:
        fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig


def figure_twin(path: str | None = None):
    plt = _plt()
    twin = MiniTwin()
    phys = twin.run(0.05e6, n_cycles=60)
    over = twin.run(0.15e6, n_cycles=60)
    episode = twin.run(lambda c: 0.15e6 if 10 <= c < 22 else 0.05e6, n_cycles=60)
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 3.0))
    for recs, lab, c in [(phys, "physiological (0.05 MPa)", "#1d4ed8"), (over, "sustained overload (0.15 MPa)", "#dc2626"),
                         (episode, "overload episode (cycles 10–22)", "#b45309")]:
        a = MiniTwin.to_arrays(recs)
        axes[0].plot(a["cycle"], a["H_A"] / 1e6, lw=1.6, color=c, label=lab)
        axes[1].plot(a["cycle"], a["lam"], lw=1.6, color=c)
        axes[2].plot(a["cycle"], a["R"] - a["S"], lw=1.6, color=c)
    axes[0].set_ylabel("aggregate modulus H_A (MPa)"); axes[0].set_title("(a) matrix", loc="left"); axes[0].legend(fontsize=5.8)
    axes[1].set_ylabel("bias λ"); axes[1].set_title("(b) mechanical → network input", loc="left"); axes[1].axhline(0, color="gray", lw=0.6)
    axes[2].set_ylabel("RUNX2 − SOX9"); axes[2].set_title("(c) fate (>0: hypertrophic)", loc="left"); axes[2].axhline(0, color="gray", lw=0.6)
    for ax in axes:
        ax.set_xlabel("loading cycle")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    return fig


def make_all_figures(outdir: str = "figures") -> Dict[str, str]:
    os.makedirs(outdir, exist_ok=True)
    paths = {
        "figure2_growth_plate": os.path.join(outdir, "figure2_growth_plate.png"),
        "figure4_toggle_switch": os.path.join(outdir, "figure4_toggle_switch.png"),
        "landscape_langevin": os.path.join(outdir, "landscape_langevin.png"),
        "loading_regime_map": os.path.join(outdir, "loading_regime_map.png"),
        "boolean_attractors": os.path.join(outdir, "boolean_attractors.png"),
        "biphasic_creep": os.path.join(outdir, "biphasic_creep.png"),
        "twin_loop": os.path.join(outdir, "twin_loop.png"),
    }
    figure_growth_plate(path=paths["figure2_growth_plate"])
    figure_toggle(path=paths["figure4_toggle_switch"])
    figure_landscape(path=paths["landscape_langevin"])
    figure_loading(path=paths["loading_regime_map"])
    figure_boolean(path=paths["boolean_attractors"])
    figure_biphasic(path=paths["biphasic_creep"])
    figure_twin(path=paths["twin_loop"])
    return paths
