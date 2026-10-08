"""Schematic (non-computed) figures of the manuscript: Figures 1, 3, 5-10, drawn with matplotlib.

Figures 2 and 4 are also reproduced here for completeness, but their canonical,
package-based implementation is ``chondrotwin.plotting``.

Usage:  python scripts/manuscript_schematics.py [OUTDIR]   (default: figures/manuscript)
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle, Rectangle, Polygon
from PIL import Image

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join('figures', 'manuscript')
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 9,
                     'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7})
dims = {}


def save(fig, n):
    p = f'{OUT}/fig{n}.png'
    fig.savefig(p, dpi=220, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    dims[str(n)] = Image.open(p).size


def box(ax, x, y, w, h, text, fc='#f3f4f6', ec='#374151', fs=7, lw=0.8, bold=False, tc='black', pad=0.02, style='round,pad=0.02,rounding_size=0.08'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs, color=tc,
            fontweight='bold' if bold else 'normal', wrap=True)


def arrow(ax, a, b, color='#111827', lw=1.0, style='-|>', ls='-', shrinkA=2, shrinkB=2, ms=9, rad=0.0):
    ax.annotate('', xy=b, xytext=a,
                arrowprops=dict(arrowstyle=style, color=color, lw=lw, ls=ls, shrinkA=shrinkA, shrinkB=shrinkB,
                                mutation_scale=ms, connectionstyle=f'arc3,rad={rad}'))


# ---------------------------------------------------------------- Figure 1: attractor landscape
def fig1():
    fig = plt.figure(figsize=(8.2, 4.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1], wspace=0.38)
    ax = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1])
    X, Y = np.meshgrid(np.linspace(0, 10, 500), np.linspace(0, 10, 500))
    basins = [('Chondro-\nprogenitor', 1.3, 7.9, 1.4, 1.0), ('Proliferative', 3.7, 8.5, 1.5, 0.9),
              ('Prehypertrophic /\nhypertrophic', 6.3, 7.7, 1.8, 1.0), ('Osteoblast', 8.8, 8.4, 1.5, 0.8),
              ('Permanent /\narticular', 3.0, 5.0, 2.2, 1.25), ('Cluster', 1.0, 3.5, 0.9, 0.7),
              ('Dedifferentiated /\nfibroblastic', 6.1, 4.7, 1.2, 0.9), ('Senescent', 4.3, 1.8, 1.9, 1.0),
              ('Ferroptotic', 7.1, 1.5, 1.7, 0.85), ('Apoptotic', 8.8, 3.5, 1.6, 0.85)]
    U = np.full_like(X, 2.4)
    for name, cx, cy, d, w in basins:
        U -= d * np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / (2 * w ** 2))
    cs = ax.contourf(X, Y, U, levels=26, cmap='viridis_r')
    ax.contour(X, Y, U, levels=13, colors='white', linewidths=0.3, alpha=0.5)
    for name, cx, cy, d, w in basins:
        ax.text(cx, cy, name, ha='center', va='center', fontsize=6.2, color='white', fontweight='bold')
    dev = [((1.3, 7.9), (3.7, 8.5)), ((3.7, 8.5), (6.3, 7.7)), ((6.3, 7.7), (8.8, 8.4)), ((1.3, 7.9), (3.0, 5.0))]
    for a, b in dev:
        arrow(ax, a, b, color='white', lw=1.3, shrinkA=16, shrinkB=16)
    path = [((3.0, 5.0), (6.1, 4.7)), ((3.0, 5.0), (6.3, 7.7)), ((3.0, 5.0), (4.3, 1.8)), ((6.1, 4.7), (4.3, 1.8)),
            ((4.3, 1.8), (7.1, 1.5)), ((6.1, 4.7), (8.8, 3.5)), ((3.0, 5.0), (1.0, 3.5)), ((6.3, 7.7), (6.1, 4.7))]
    for a, b in path:
        arrow(ax, a, b, color='#ff7b7b', lw=1.2, ls='--', shrinkA=16, shrinkB=16)
    ax.annotate('', xy=(9.55, 6.4), xytext=(9.55, 0.5),
                arrowprops=dict(arrowstyle='-|>', lw=2.2, color='white', mutation_scale=14))
    ax.text(9.85, 3.4, 'mechanical load, inflammation, ageing tilt the surface', rotation=90, ha='center', va='center', fontsize=6, color='white')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel('differentiation axis: SOX9-high  →  RUNX2-high')
    ax.set_ylabel('terminal / pathological  ←→  developmental')
    ax.set_title('(a) chondrocyte fate as a quasi-potential landscape', loc='left')
    cb = fig.colorbar(cs, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label('U (a.u.)', fontsize=7); cb.ax.tick_params(labelsize=6)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], color='white', lw=1.3, label='developmental transition'),
                       Line2D([0], [0], color='#ff7b7b', lw=1.2, ls='--', label='pathological transition')],
              loc='lower left', fontsize=6, framealpha=0.85)
    # panel b
    s = np.linspace(-1.7, 1.7, 400)
    U0 = s ** 4 - 2 * s ** 2
    U1 = s ** 4 - 2 * s ** 2 - 1.0 * s
    ax2.plot(s, U0, color='#1d4ed8', lw=1.6, label='physiological loading')
    ax2.plot(s, U1, color='#dc2626', lw=1.6, label='overload / inflammation')
    ax2.plot([-1, 1], [U0[np.argmin(np.abs(s + 1))], U0[np.argmin(np.abs(s - 1))]], 'o', color='#1d4ed8', ms=5)
    sm = s[np.argmin(U1)]
    ax2.plot([sm], [U1.min()], 'o', color='#dc2626', ms=5)
    ax2.annotate('', xy=(0, 0), xytext=(0, -1), arrowprops=dict(arrowstyle='<->', color='k', lw=1))
    ax2.text(0.08, -0.5, 'ΔU\n(barrier)', fontsize=6.5, va='center')
    ax2.text(-1, -1.35, 'SOX9-high\n(articular)', ha='center', fontsize=6.5, color='#1d4ed8')
    ax2.text(1, -1.35, 'RUNX2-high\n(hypertrophic)', ha='center', fontsize=6.5, color='#1d4ed8')
    ax2.set_xlabel('position along the permanent–hypertrophic axis')
    ax2.set_ylabel('U (a.u.)')
    ax2.set_ylim(-2.6, 2.2)
    ax2.legend(loc='upper center', fontsize=6)
    ax2.set_title('(b) load tilts the basins and lowers one barrier', loc='left')
    save(fig, 1)


# ---------------------------------------------------------------- Figure 2: growth plate + IHH–PTHrP model
def simulate(kP=1.0, aP=2.0, noP=False, T=40.0):  # T=40: the Jansen-like case (kP=0.35) is not yet steady at T=12
    N = 201
    x = np.linspace(0, 1, N); dx = x[1] - x[0]
    DP, DI = 0.03, 0.006
    KI, KP = 0.3, 0.3
    dt = 0.2 * dx * dx / DP
    steps = int(T / dt)
    P = np.zeros(N); I = np.zeros(N)
    sP = np.exp(-x ** 2 / (2 * 0.08 ** 2))
    lapP = np.zeros(N); lapI = np.zeros(N)
    for _ in range(steps):
        lapP[1:-1] = (P[2:] - 2 * P[1:-1] + P[:-2]) / dx ** 2
        lapI[1:-1] = (I[2:] - 2 * I[1:-1] + I[:-2]) / dx ** 2
        prodP = 0.0 if noP else aP * sP * (0.25 + I ** 2 / (KI ** 2 + I ** 2))
        P += dt * (DP * lapP + prodP - kP * P)
        I += dt * (DI * lapI + 1.0 * KP ** 2 / (KP ** 2 + P ** 2) - 1.0 * I)
        P[0] = P[1]; P[-1] = P[-2]; I[0] = I[1]; I[-1] = I[-2]
    H = KP ** 2 / (KP ** 2 + P ** 2)
    return x, P, I, H


def fig2():
    fig = plt.figure(figsize=(8.2, 4.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.1, 1.1], wspace=0.35)
    ax = fig.add_subplot(gs[0])
    zones = [('Resting zone\n(SOX9; PTHrP-\nresponsive)', 0.80, 1.0, '#dbeafe'),
             ('Proliferative zone\n(SOX9, COL2A1;\ncolumns, IHH target)', 0.44, 0.80, '#bfdbfe'),
             ('Prehypertrophic\n(IHH, PTH1R)', 0.31, 0.44, '#fde68a'),
             ('Hypertrophic\n(RUNX2, COL10A1,\nMMP13, VEGF)', 0.08, 0.31, '#fdba74'),
             ('Bone / marrow\n(transdifferentiation;\nosteoclasts)', 0.0, 0.08, '#d1d5db')]
    for name, y0, y1, c in zones:
        ax.add_patch(Rectangle((0, y0), 1, y1 - y0, facecolor=c, edgecolor='k', lw=0.5))
        ax.text(1.06, (y0 + y1) / 2, name, va='center', fontsize=5.8)
    rng = np.random.default_rng(1)
    for _ in range(16):
        ax.add_patch(Circle((rng.uniform(0.1, 0.9), rng.uniform(0.84, 0.97)), 0.02, fc='white', ec='k', lw=0.4))
    for col in np.linspace(0.15, 0.85, 5):
        for y in np.linspace(0.47, 0.77, 7):
            ax.add_patch(Ellipse((col, y), 0.11, 0.03, fc='white', ec='k', lw=0.4))
    for col in np.linspace(0.15, 0.85, 5):
        for y in [0.345, 0.405]:
            ax.add_patch(Circle((col, y), 0.025, fc='white', ec='k', lw=0.4))
    for col in np.linspace(0.16, 0.84, 4):
        for y in [0.135, 0.245]:
            ax.add_patch(Circle((col, y), 0.048, fc='white', ec='k', lw=0.4))
    ax.annotate('', xy=(-0.09, 0.33), xytext=(-0.09, 0.97), arrowprops=dict(arrowstyle='-|>', color='#1d4ed8', lw=1.6))
    ax.text(-0.16, 0.66, 'PTHrP gradient', rotation=90, va='center', ha='center', color='#1d4ed8', fontsize=6.5)
    ax.annotate('', xy=(-0.36, 0.97), xytext=(-0.36, 0.33), arrowprops=dict(arrowstyle='-|>', color='#b45309', lw=1.6))
    ax.text(-0.44, 0.66, 'IHH (feedback)', rotation=90, va='center', ha='center', color='#b45309', fontsize=6.5)
    ax.set_xlim(-0.55, 2.0); ax.set_ylim(-0.02, 1.02); ax.axis('off')
    ax.set_title('(a) growth-plate zonation', loc='left')
    # simulations
    x, P, I, H = simulate()
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(x, P / P.max(), color='#1d4ed8', lw=1.6, label='PTHrP (normalized)')
    ax2.plot(x, I / I.max(), color='#b45309', lw=1.6, label='IHH (normalized)')
    ax2.plot(x, H, color='#dc2626', lw=1.6, ls='--', label='hypertrophic\ncommitment')
    xb = x[np.argmin(np.abs(H - 0.5))]
    ax2.axvline(xb, color='gray', lw=0.8, ls=':')
    ax2.text(xb + 0.02, 0.05, 'prehypertrophic\nboundary', fontsize=6, color='gray')
    ax2.set_xlabel('distance from periarticular end (a.u.)'); ax2.set_ylabel('level (a.u.)')
    ax2.set_title('(b) minimal reaction–diffusion model', loc='left'); ax2.legend(loc='center right', bbox_to_anchor=(1.0, 0.33), fontsize=5.8)
    ax3 = fig.add_subplot(gs[2])
    for lab, kw, c in [('wild type', dict(), '#111827'), ('PTH1R gain of function\n(Jansen-like: longer PTHrP range)', dict(kP=0.35), '#1d4ed8'),
                       ('PTH1R / PTHrP loss\n(Blomstrand-like)', dict(noP=True), '#dc2626')]:
        xx, _, _, HH = simulate(**kw)
        ax3.plot(xx, HH, lw=1.6, color=c, label=lab)
    ax3.set_xlabel('distance from periarticular end (a.u.)'); ax3.set_ylabel('hypertrophic commitment')
    ax3.set_title('(c) perturbing the loop shifts the boundary', loc='left'); ax3.legend(fontsize=5.6, loc='center right')
    save(fig, 2)


# ---------------------------------------------------------------- Figure 3: signaling network map
def fig3():
    fig, ax = plt.subplots(figsize=(8.4, 6.0))
    ax.set_xlim(0, 13.6); ax.set_ylim(0, 9.2); ax.axis('off')
    N = {}
    def node(name, x, y, w=1.5, h=0.5, fc='#f3f4f6', fs=6.4, bold=False):
        N[name] = (x, y, w, h)
        box(ax, x - w / 2, y - h / 2, w, h, name, fc=fc, fs=fs, bold=bold)
    # inputs
    for nm, x, w in [('Mechanical load', 1.6, 2.0), ('IL-1β / TNF-α', 4.4, 2.0), ('Hypoxia', 7.2, 1.6), ('Oxidative stress / ageing', 10.3, 2.6)]:
        node(nm, x, 8.6, w=w, fc='#fee2e2', bold=True)
    # sensors
    for nm, x, w in [('Piezo1 / Piezo2', 0.9, 1.5), ('TRPV4', 2.4, 1.1), ('Primary cilium', 3.8, 1.4), ('Integrin–FAK–Src', 5.4, 1.65)]:
        node(nm, x, 7.3, w=w, fc='#fef3c7')
    # transducers
    for nm, x, w in [('Ca²⁺', 1.3, 1.0), ('YAP / TAZ', 3.3, 1.2), ('NF-κB', 5.0, 1.0), ('PI3K/AKT–mTOR', 6.75, 1.6), ('HIF-1α', 8.2, 0.9), ('HIF-2α', 9.3, 0.9), ('cGAS–STING', 10.7, 1.3)]:
        node(nm, x, 6.0, w=w, fc='#e0f2fe')
    node('Nrf2', 12.2, 6.0, w=1.0, fc='#e0f2fe')
    # growth factor pathways
    for nm, x in [('WNT / β-catenin', 1.3), ('TGF-β / SMAD2/3', 3.1), ('BMP / SMAD1/5/8', 4.9), ('FGFR3–MAPK', 6.6), ('IHH', 8.1), ('PTHrP / PTH1R', 9.6)]:
        node(nm, x, 4.6, w=1.55, fc='#dcfce7')
    # core switch
    node('SOX9 (+SOX5/6)', 3.6, 3.0, w=1.9, h=0.6, fc='#bfdbfe', fs=7, bold=True)
    node('RUNX2 (+MEF2C)', 7.6, 3.0, w=1.9, h=0.6, fc='#fed7aa', fs=7, bold=True)
    # outputs
    for nm, x, fc in [('Permanent / articular', 1.4, '#bfdbfe'), ('Hypertrophic', 4.0, '#fed7aa'), ('Dedifferentiated', 6.2, '#ddd6fe'),
                      ('Senescent', 8.4, '#e9d5ff'), ('Apoptotic / ferroptotic', 10.6, '#fecaca')]:
        node(nm, x, 1.1, w=1.9, h=0.55, fc=fc, fs=6.6, bold=True)

    def E(a, b, kind='act', rad=0.0, color=None):
        xa, ya, wa, ha = N[a]; xb, yb, wb, hb = N[b]
        # start/end at box edges (vertical direction)
        sa = (xa, ya - ha / 2) if yb < ya else (xa, ya + ha / 2)
        sb = (xb, yb + hb / 2) if yb < ya else (xb, yb - hb / 2)
        if kind == 'act':
            arrow(ax, sa, sb, color=color or '#1f2937', lw=0.9, style='-|>', ms=8, rad=rad)
        else:
            arrow(ax, sa, sb, color=color or '#dc2626', lw=0.9, style='-[,widthB=0.3,lengthB=0.01', ms=8, rad=rad)
    E('Mechanical load', 'Piezo1 / Piezo2'); E('Mechanical load', 'TRPV4'); E('Mechanical load', 'Primary cilium'); E('Mechanical load', 'Integrin–FAK–Src')
    E('Piezo1 / Piezo2', 'Ca²⁺'); E('TRPV4', 'Ca²⁺'); E('Primary cilium', 'Ca²⁺'); E('Integrin–FAK–Src', 'YAP / TAZ')
    E('IL-1β / TNF-α', 'NF-κB'); E('IL-1β / TNF-α', 'PI3K/AKT–mTOR'); E('IL-1β / TNF-α', 'YAP / TAZ', 'inh', rad=0.15)
    E('PI3K/AKT–mTOR', 'Piezo1 / Piezo2', rad=-0.35, color='#6b7280')
    E('Hypoxia', 'HIF-1α'); E('Hypoxia', 'HIF-2α'); E('Oxidative stress / ageing', 'cGAS–STING'); E('Oxidative stress / ageing', 'Nrf2')
    E('cGAS–STING', 'NF-κB', rad=0.25); E('YAP / TAZ', 'NF-κB', 'inh', rad=0.1)
    E('Ca²⁺', 'SOX9 (+SOX5/6)', color='#2563eb'); E('Ca²⁺', 'Apoptotic / ferroptotic', rad=0.42, color='#dc2626')
    E('YAP / TAZ', 'Dedifferentiated', rad=0.25); E('NF-κB', 'SOX9 (+SOX5/6)', 'inh'); E('NF-κB', 'Senescent', rad=-0.2)
    E('HIF-1α', 'SOX9 (+SOX5/6)', rad=0.2); E('HIF-2α', 'RUNX2 (+MEF2C)'); E('cGAS–STING', 'Senescent', rad=0.1); E('cGAS–STING', 'Apoptotic / ferroptotic')
    E('Nrf2', 'Apoptotic / ferroptotic', 'inh')
    E('WNT / β-catenin', 'RUNX2 (+MEF2C)', rad=0.1); E('WNT / β-catenin', 'SOX9 (+SOX5/6)', 'inh')
    E('TGF-β / SMAD2/3', 'SOX9 (+SOX5/6)'); E('TGF-β / SMAD2/3', 'RUNX2 (+MEF2C)', 'inh', rad=0.1)
    E('BMP / SMAD1/5/8', 'SOX9 (+SOX5/6)'); E('BMP / SMAD1/5/8', 'RUNX2 (+MEF2C)')
    E('FGFR3–MAPK', 'RUNX2 (+MEF2C)'); E('IHH', 'PTHrP / PTH1R', rad=-0.3); E('IHH', 'RUNX2 (+MEF2C)'); E('PTHrP / PTH1R', 'RUNX2 (+MEF2C)', 'inh')
    E('PI3K/AKT–mTOR', 'SOX9 (+SOX5/6)', rad=0.15, color='#6b7280')
    # toggle
    arrow(ax, (4.55, 3.1), (6.65, 3.1), color='#dc2626', lw=1.4, style='-[,widthB=0.35,lengthB=0.01', ms=9)
    arrow(ax, (6.65, 2.9), (4.55, 2.9), color='#dc2626', lw=1.4, style='-[,widthB=0.35,lengthB=0.01', ms=9)
    ax.text(5.6, 3.35, 'mutual inhibition\n(bistable switch)', ha='center', fontsize=6, color='#dc2626')
    E('SOX9 (+SOX5/6)', 'Permanent / articular'); E('RUNX2 (+MEF2C)', 'Hypertrophic', rad=0.1)
    # legend
    arrow(ax, (0.2, 0.25), (1.0, 0.25), color='#1f2937', lw=0.9, style='-|>', ms=8); ax.text(1.1, 0.25, 'activation', va='center', fontsize=6.5)
    arrow(ax, (2.3, 0.25), (3.1, 0.25), color='#dc2626', lw=0.9, style='-[,widthB=0.3,lengthB=0.01', ms=8); ax.text(3.2, 0.25, 'inhibition', va='center', fontsize=6.5)
    arrow(ax, (4.4, 0.25), (5.2, 0.25), color='#6b7280', lw=0.9, style='-|>', ms=8); ax.text(5.3, 0.25, 'context-dependent / indirect', va='center', fontsize=6.5)
    for y, lab in [(8.6, 'inputs'), (7.3, 'mechanosensors'), (6.0, 'transducers'), (4.6, 'growth-factor pathways'), (3.0, 'core switch'), (1.1, 'fate outputs')]:
        ax.text(13.5, y, lab, ha='right', va='center', fontsize=6, color='#6b7280', style='italic')
    save(fig, 3)


# ---------------------------------------------------------------- Figure 4: SOX9–RUNX2 toggle switch
aS = aR = 2.0; K = 0.5; n = 4; bS = bR = 0.05


def f1(S, R, lam=0.0):
    return aS * (1 - lam) / (1 + (R / K) ** n) + bS - S


def f2(S, R, lam):
    return aR * (1 + lam) / (1 + (S / K) ** n) + bR - R


def fixed_points(lam):
    Rg = np.linspace(0, 4.5, 4000)
    Sg = aS * (1 - lam) / (1 + (Rg / K) ** n) + bS
    g = aR * (1 + lam) / (1 + (Sg / K) ** n) + bR - Rg
    out = []
    for i in range(len(Rg) - 1):
        if g[i] == 0 or g[i] * g[i + 1] < 0:
            lo, hi = Rg[i], Rg[i + 1]
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                Sm = aS * (1 - lam) / (1 + (mid / K) ** n) + bS
                gm = aR * (1 + lam) / (1 + (Sm / K) ** n) + bR - mid
                if g[i] * gm <= 0: hi = mid
                else: lo = mid
            R0 = 0.5 * (lo + hi); S0 = aS * (1 - lam) / (1 + (R0 / K) ** n) + bS
            d1R = -aS * (1 - lam) * n * (R0 / K) ** (n - 1) / K / (1 + (R0 / K) ** n) ** 2
            d2S = -aR * (1 + lam) * n * (S0 / K) ** (n - 1) / K / (1 + (S0 / K) ** n) ** 2
            det = 1 - d1R * d2S
            out.append((S0, R0, det > 0))
    return out


def fig4():
    fig = plt.figure(figsize=(8.2, 3.2))
    gs = fig.add_gridspec(1, 3, wspace=0.38)
    ax = fig.add_subplot(gs[0]); lam = 0.3
    S, R = np.meshgrid(np.linspace(0, 3.0, 30), np.linspace(0, 3.0, 30))
    ax.streamplot(S, R, f1(S, R, lam), f2(S, R, lam), color='#9ca3af', density=0.9, linewidth=0.5, arrowsize=0.6)
    r = np.linspace(0, 3.0, 400)
    ax.plot(aS * (1 - lam) / (1 + (r / K) ** n) + bS, r, color='#1d4ed8', lw=1.6, label='dSOX9/dt = 0')
    s = np.linspace(0, 3.0, 400)
    ax.plot(s, aR * (1 + lam) / (1 + (s / K) ** n) + bR, color='#dc2626', lw=1.6, label='dRUNX2/dt = 0')
    for S0, R0, stable in fixed_points(lam):
        ax.plot(S0, R0, 'o', ms=6, mfc='k' if stable else 'white', mec='k')
    ax.set_xlabel('SOX9 activity'); ax.set_ylabel('RUNX2 activity'); ax.set_title('(a) phase plane (λ = 0.3)', loc='left')
    ax.legend(fontsize=6, loc='upper right'); ax.set_xlim(0, 3.0); ax.set_ylim(0, 3.0)
    ax.text(1.75, 0.65, 'articular\nattractor', fontsize=6, ha='center'); ax.text(1.0, 2.25, 'hypertrophic\nattractor', fontsize=6, ha='center')
    # bifurcation
    ax2 = fig.add_subplot(gs[1])
    lams = np.linspace(-0.95, 0.95, 400)
    st_x, st_y, un_x, un_y = [], [], [], []
    for l in lams:
        for S0, R0, stable in fixed_points(l):
            (st_x if stable else un_x).append(l); (st_y if stable else un_y).append(R0)
    ax2.plot(st_x, st_y, '.', color='#111827', ms=2, label='stable')
    ax2.plot(un_x, un_y, '.', color='#9ca3af', ms=1.5, label='unstable (saddle)')
    # find bistable range
    bist = [l for l in lams if sum(1 for *_, s_ in fixed_points(l) if s_) == 2]
    if bist:
        ax2.axvspan(min(bist), max(bist), color='#fde68a', alpha=0.5)
        ax2.text((min(bist) + max(bist)) / 2, 3.75, 'bistable\n(hysteresis)', ha='center', fontsize=6)
        ax2.annotate('', xy=(max(bist), 0.3), xytext=(max(bist) - 0.6, 0.3), arrowprops=dict(arrowstyle='-|>', color='#dc2626', lw=1))
        ax2.annotate('', xy=(min(bist), 2.3), xytext=(min(bist) + 0.6, 2.3), arrowprops=dict(arrowstyle='-|>', color='#1d4ed8', lw=1))
    ax2.set_xlabel('λ: catabolic bias (overload, inflammation)  →\n←  anabolic bias (TGF-β, physiological load)')
    ax2.set_ylabel('steady-state RUNX2'); ax2.set_title('(b) bifurcation diagram', loc='left'); ax2.legend(fontsize=6, loc='center left'); ax2.set_ylim(0, 4.2)
    # time course
    ax3 = fig.add_subplot(gs[2])
    T = 80; dt = 0.01; t = np.arange(0, T, dt)
    lam_t = np.where((t > 20) & (t < 40), 0.8, 0.0)
    S0, R0 = 2.0, 0.06; Ss, Rs = [], []
    for l in lam_t:
        S0 += dt * f1(S0, R0, l); R0 += dt * f2(S0, R0, l); Ss.append(S0); Rs.append(R0)
    ax3.fill_between(t, 0, 4.0, where=lam_t > 0, color='#fde68a', alpha=0.5, label='transient overload (λ = 0.8)')
    ax3.plot(t, Ss, color='#1d4ed8', lw=1.6, label='SOX9'); ax3.plot(t, Rs, color='#dc2626', lw=1.6, label='RUNX2')
    ax3.set_xlabel('time (a.u.)'); ax3.set_ylabel('activity'); ax3.set_title('(c) memory of a transient insult', loc='left')
    ax3.legend(fontsize=6, loc='center right'); ax3.set_ylim(0, 4.0)
    save(fig, 4)


# ---------------------------------------------------------------- Figure 5: mechanotransduction map
def fig5():
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.set_xlim(0, 12.2); ax.set_ylim(0, 8); ax.axis('off')
    # loading cascade strip
    casc = ['joint load', 'tissue strain &\nhydrostatic pressure', 'interstitial fluid flow', 'osmotic / ionic shifts', 'membrane tension &\ncell deformation']
    for i, c in enumerate(casc):
        box(ax, 0.3 + i * 2.35, 7.1, 2.1, 0.7, c, fc='#fee2e2', fs=6.5)
        if i < len(casc) - 1:
            arrow(ax, (0.3 + i * 2.35 + 2.1, 7.45), (0.3 + (i + 1) * 2.35, 7.45), lw=1.0, ms=8)
    ax.text(0.3, 6.75, 'the loading cascade: one joint-level load becomes several cell-level signals', fontsize=6.5, style='italic', color='#6b7280')
    cx, cy = 4.4, 3.4
    ax.add_patch(Ellipse((cx, cy), 7.4, 5.6, fc='#fdf2e9', ec='#c2410c', lw=0.8, ls='--'))
    ax.text(cx - 3.4, cy + 2.55, 'pericellular matrix (perlecan, collagen VI, aggrecan)', fontsize=6.3, color='#c2410c')
    ax.add_patch(Ellipse((cx, cy), 5.0, 3.6, fc='#e8f4fd', ec='#1e3a8a', lw=1.0))
    ax.add_patch(Ellipse((cx, cy - 0.1), 2.0, 1.4, fc='#cfd8ea', ec='#1e3a8a', lw=0.8))
    xs = np.linspace(cx - 0.8, cx + 0.8, 200)
    ax.plot(xs, cy - 0.1 + 0.25 * np.sin(14 * (xs - cx)), color='#1e3a8a', lw=0.7)
    ax.text(cx, cy - 0.1 - 0.5, 'nucleus: chromatin\nmechanoresponse', ha='center', va='top', fontsize=5.8, color='#1e3a8a')
    a, b = 2.5, 1.8
    def mem(theta):
        th = np.deg2rad(theta); return (cx + a * np.cos(th), cy + b * np.sin(th))
    # channels and sensors
    for theta, lab, fc in [(100, 'Piezo1 / Piezo2\n(high strain)', '#fde68a'), (145, 'TRPV4\n(physiological, osmotic)', '#fde68a'),
                           (300, 'Cx43 hemichannels /\nP2 receptors (ATP)', '#e0f2fe')]:
        x, y = mem(theta)
        ax.add_patch(Rectangle((x - 0.18, y - 0.14), 0.36, 0.28, fc=fc, ec='k', lw=0.6))
        off = (0, 0.55) if theta < 180 else (0, -0.55)
        ax.text(x + off[0], y + off[1], lab, ha='center', va='center', fontsize=5.8)
    # integrins
    x, y = mem(220)
    for dx in (-0.25, 0, 0.25):
        ax.plot([x + dx, x + dx], [y - 0.2, y + 0.2], color='k', lw=1.2)
        ax.plot([x + dx - 0.08, x + dx, x + dx + 0.08], [y - 0.28, y - 0.2, y - 0.28], color='k', lw=1.0)
    ax.text(x - 0.1, y - 0.75, 'integrin α5β1–FAK–Src\n(matrix composition, stiffness)', ha='center', fontsize=5.8)
    # stress fibers to nucleus / LINC
    nx_, ny_ = cx - 0.95, cy - 0.3
    for dx in (-0.25, 0, 0.25):
        ax.plot([x + dx, nx_], [y + 0.2, ny_], color='#4b5563', lw=0.8)
    ax.add_patch(Rectangle((nx_ - 0.12, ny_ - 0.12), 0.24, 0.24, fc='#9ca3af', ec='k', lw=0.5))
    ax.text(1.05, 3.1, 'actin stress\nfibers → LINC /\nlamin A/C', fontsize=5.6, ha='center', va='center')
    # cilium
    x, y = mem(20)
    ax.plot([x, x + 1.4], [y, y + 0.9], color='#065f46', lw=2.2); ax.add_patch(Circle((x, y), 0.12, fc='#065f46', ec='k', lw=0.4))
    ax.text(x + 0.3, y + 1.15, 'primary cilium\n(IFT88, polycystins,\nciliary TRPV4, Hedgehog)', fontsize=5.8, va='bottom', ha='center', color='#065f46')
    # intracellular signaling labels inside cell
    ax.text(cx + 0.3, cy + 1.2, 'Ca²⁺ transients', fontsize=6, color='#b45309', ha='center')
    ax.text(cx + 1.6, cy + 0.4, 'YAP / TAZ', fontsize=6, color='#7c3aed', ha='center')
    ax.text(cx - 0.5, cy - 1.25, 'RhoA / ROCK', fontsize=6, color='#4b5563', ha='center')
    ax.text(cx + 1.5, cy - 1.0, 'PI3K/AKT/mTOR\nNF-κB', fontsize=6, color='#b91c1c', ha='center')
    # right-hand outcome boxes
    outs = [('physiological window\nTRPV4 → oscillatory Ca²⁺ → anabolic genes\n(ACAN, COL2A1); SOX9 maintained', '#dcfce7'),
            ('overload\nPiezo1/2 → Ca²⁺ overload → mitochondrial\nstress, cGAS–STING, GPX4↓ → ferroptosis,\napoptosis, senescence', '#fecaca'),
            ('stiff / degraded matrix\nintegrin–FAK → RhoA → nuclear YAP/TAZ →\ndedifferentiation; YAP ⊣ NF-κB', '#ede9fe'),
            ('inflammation\nIL-1α → p38 MAPK (HNF4, ATF2) → Piezo1↑;\nTAK1 → YAP/TAZ degradation (feed-forward)', '#fee2e2')]
    y0 = 6.1
    for t, fc in outs:
        box(ax, 8.75, y0 - 1.25, 3.15, 1.15, t, fc=fc, fs=5.5)
        y0 -= 1.45
    arrow(ax, (cx + 2.5, cy + 0.9), (8.75, 5.5), lw=0.8, ms=7, color='#6b7280')
    arrow(ax, (cx + 2.5, cy + 0.4), (8.75, 4.05), lw=0.8, ms=7, color='#6b7280')
    arrow(ax, (cx + 2.5, cy - 0.2), (8.75, 2.6), lw=0.8, ms=7, color='#6b7280')
    arrow(ax, (cx + 2.5, cy - 0.8), (8.75, 1.15), lw=0.8, ms=7, color='#6b7280')
    ax.text(0.3, 0.25, 'Mechanosensors operate on different time scales: ion channels (ms–s) → cytoskeleton/YAP (min–h) → chromatin and matrix remodeling (h–d).',
            fontsize=6, style='italic', color='#6b7280')
    save(fig, 5)


# ---------------------------------------------------------------- Figure 6: loading-regime map
def response(m):
    return 1.0 * np.exp(-((m - 0.45) / 0.2) ** 2) - 0.85 / (1 + np.exp(-(m - 0.82) / 0.05)) - 0.4 * np.exp(-(m / 0.12) ** 2)


def fig6():
    fig = plt.figure(figsize=(8.2, 3.3))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.3)
    ax = fig.add_subplot(gs[0])
    m = np.linspace(0, 1.1, 400); r = response(m)
    ax.plot(m, r, color='#111827', lw=1.8)
    ax.axhline(0, color='gray', lw=0.6)
    ax.fill_between(m, r, 0, where=r > 0, color='#bbf7d0', alpha=0.7); ax.fill_between(m, r, 0, where=r < 0, color='#fecaca', alpha=0.7)
    ax.axvspan(0, 0.15, color='#e5e7eb', alpha=0.5); ax.axvspan(0.78, 1.1, color='#e5e7eb', alpha=0.5)
    ax.text(0.075, 0.75, 'disuse /\nimmobilization', ha='center', fontsize=6); ax.text(0.45, 1.12, 'physiological\nwindow (TRPV4)', ha='center', fontsize=6)
    ax.text(0.94, 0.75, 'overload\n(Piezo1/2)', ha='center', fontsize=6)
    ax.set_xlabel('loading magnitude (strain, a.u.)'); ax.set_ylabel('net anabolic (+) / catabolic (−) response')
    ax.set_title('(a) inverted-U dose–response', loc='left'); ax.set_ylim(-1.1, 1.4)
    ax2 = fig.add_subplot(gs[1])
    M, F = np.meshgrid(np.linspace(0, 1.1, 300), np.linspace(0, 3, 300))
    fw = np.exp(-((np.log(F + 0.05) - np.log(1.0)) / 0.8) ** 2)
    Z = 1.0 * np.exp(-((M - 0.45) / 0.2) ** 2) * fw - 0.85 / (1 + np.exp(-(M - 0.82) / 0.05)) - 0.4 * np.exp(-(M / 0.12) ** 2) - 0.35 * np.exp(-(F / 0.15) ** 2) * (M > 0.15)
    cs = ax2.contourf(M, F, Z, levels=20, cmap='RdYlGn', vmin=-1, vmax=1)
    ax2.contour(M, F, Z, levels=[0], colors='k', linewidths=0.8)
    ax2.text(0.42, 1.0, 'anabolic\nwindow', ha='center', fontsize=6.5, fontweight='bold')
    ax2.text(0.95, 1.5, 'overload:\ncatabolic', ha='center', fontsize=6, color='white', fontweight='bold')
    ax2.text(0.07, 1.5, 'disuse', ha='center', fontsize=6, color='white', fontweight='bold', rotation=90)
    ax2.text(0.5, 0.12, 'static compression: catabolic', ha='center', fontsize=5.8)
    # patient-specific shifted window (aged / early OA)
    from matplotlib.patches import Ellipse as El
    ax2.add_patch(El((0.36, 0.95), 0.34, 1.4, fc='none', ec='#1d4ed8', lw=1.2, ls='--'))
    ax2.text(0.36, 2.0, 'narrower window:\naged / early-OA\n(stiff PCM, sensitized Piezo1)', ha='center', fontsize=5.6, color='#1d4ed8')
    ax2.set_xlabel('loading magnitude (strain, a.u.)'); ax2.set_ylabel('loading frequency (Hz)')
    ax2.set_title('(b) loading-regime map (magnitude × frequency)', loc='left')
    cb = fig.colorbar(cs, ax=ax2, fraction=0.04, pad=0.02); cb.set_label('net response', fontsize=7); cb.ax.tick_params(labelsize=6)
    save(fig, 6)


# ---------------------------------------------------------------- Figure 7: multiscale hierarchy
def fig7():
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.set_xlim(0, 12.9); ax.set_ylim(0, 10); ax.axis('off')
    levels = [('Whole joint', 'poroviscoelastic FE on patient / atlas geometry (OAI); gait-derived loads\n→ stress, strain, fluid pressure fields', '#fee2e2'),
              ('Tissue', 'biphasic, fibril-reinforced cartilage; collagen / aggrecan / fixed charge\n→ depth-dependent strain, flow, osmotic environment', '#ffedd5'),
              ('Chondron / PCM', 'micro-FE of cell + pericellular matrix (perlecan, collagen VI)\n→ cell strain, membrane tension, local osmolarity', '#fef3c7'),
              ('Cell / gene network', 'Boolean–ODE hybrid GRN (SOX9–RUNX2, WNT, TGF-β, YAP/TAZ, NF-κB …)\n→ attractor occupancy, matrix synthesis / degradation rates', '#dcfce7'),
              ('Fate landscape', 'quasi-potential surface, basin depths, barrier heights, transition rates\n→ fate probabilities, progression risk', '#dbeafe')]
    anchors = ['MRI morphometry, gait analysis, joint-space width', 'mechanical testing, histology (OARSI), T2 / dGEMRIC',
               'AFM of PCM, chondron imaging, osmotic challenge', 'scRNA-seq, reporter assays, Ca²⁺ imaging, explant loading',
               'trajectory inference, lineage tracing, longitudinal cohorts']
    y = 8.6
    for (name, desc, fc), anc in zip(levels, anchors):
        box(ax, 1.6, y - 0.75, 6.5, 1.5, '', fc=fc)
        ax.text(1.75, y + 0.35, name, fontsize=7.5, fontweight='bold', va='center')
        ax.text(1.75, y - 0.25, desc, fontsize=5.7, va='center')
        ax.text(8.4, y, 'validation anchors:\n' + anc, fontsize=5.7, va='center', color='#374151')
        y -= 1.85
    arrow(ax, (0.9, 9.2), (0.9, 1.0), color='#b91c1c', lw=2.0, ms=14)
    ax.text(0.45, 5.1, 'load transmission (top-down):\njoint load → cell-level stimuli', rotation=90, ha='center', va='center', fontsize=6.3, color='#b91c1c')
    arrow(ax, (12.0, 1.0), (12.0, 9.2), color='#1d4ed8', lw=2.0, ms=14)
    ax.text(12.45, 5.1, 'emergence (bottom-up):\nmatrix turnover → stiffness, thickness, geometry', rotation=90, ha='center', va='center', fontsize=6.3, color='#1d4ed8')
    ax.text(0.3, 0.25, 'Each level has its own formalism (Table 6) and its own validation anchor; CHONDRO-TWIN couples them in both directions.', fontsize=6.2, style='italic', color='#6b7280')
    save(fig, 7)


# ---------------------------------------------------------------- Figure 8: AI/ML landscape
def fig8():
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.set_xlim(0, 12); ax.set_ylim(0, 9); ax.axis('off')
    quads = [(0.3, 4.85, 'IMAGING  (learns structure & progression)', '#dbeafe',
              '• U-Net / V-Net, transformers, SAM-2 (SAMRI-2),\n   KneeXNet-2.5D, KneeSeg-U\n• cartilage Dice ≈ 0.86–0.90; OA diagnosis, TKR and\n   progression prediction (AUC ≈ 0.7–0.97)\n• data: OAI / MOST radiographs & 3T MRI; T2 maps, radiomics\n• failure modes: domain shift, weak symptom coupling,\n   limited interpretability'),
             (6.1, 4.85, 'OMICS  (learns cell states & networks)', '#dcfce7',
              '• scRNA-seq / snRNA-seq / spatial atlases;\n   RNA velocity, pseudotime\n• GRN inference (SCENIC, CellOracle); scVI latent spaces\n• foundation models: Geneformer, scGPT, scFoundation,\n   UCE; perturbation models (GEARS)\n• caution: perturbation prediction not yet better than\n   linear baselines'),
             (0.3, 0.4, 'HYBRID / PHYSICS-INFORMED  (learns dynamics)', '#fef3c7',
              '• PINNs for biphasic / poroelastic inverse problems\n• neural ODEs; FE surrogates and mesh-based graph\n   networks (real-time joint mechanics)\n• NN load prediction improving FE-based OA\n   classification (OAI)\n• Bayesian calibration, global sensitivity analysis,\n   uncertainty quantification'),
             (6.1, 0.4, 'REASONING & STRUCTURE  (assists)', '#ede9fe',
              '• LLM agents: literature curation, hypothesis\n   generation, provenance tracking\n• AlphaFold 2/3 and cryo-EM: Piezo1/2, TRPV4 gating\n   and drug-binding sites\n• explainable AI for clinical traceability\n• human oversight mandatory; immature for\n   autonomous use')]
    for x, y, title, fc, body in quads:
        box(ax, x, y, 5.6, 3.75, '', fc=fc)
        ax.text(x + 0.2, y + 3.4, title, fontsize=7.3, fontweight='bold', va='center')
        ax.text(x + 0.2, y + 1.6, body, fontsize=6.0, va='center')
    box(ax, 4.9, 4.28, 2.2, 0.44, 'CHONDRO-TWIN components', fc='white', fs=6, bold=True)
    save(fig, 8)


# ---------------------------------------------------------------- Figure 9: CHONDRO-TWIN workflow
def fig9():
    fig, ax = plt.subplots(figsize=(8.4, 4.9))
    ax.set_xlim(0, 12); ax.set_ylim(0, 9); ax.axis('off')
    # inputs
    box(ax, 0.3, 1.3, 2.3, 6.6, '', fc='#f3f4f6')
    ax.text(1.45, 7.55, 'INPUTS', ha='center', fontsize=7.5, fontweight='bold')
    for i, t in enumerate(['qMRI / radiographs\n(geometry, T2, thickness)', 'gait & activity\n(loading history)', 'scRNA-seq / biopsy\n(state priors)',
                           'mechanical tests, AFM\n(tissue / PCM properties)', 'demographics, biomarkers,\ninflammatory status']):
        box(ax, 0.45, 6.4 - i * 1.15, 2.0, 0.95, t, fc='white', fs=5.8)
    # core
    box(ax, 3.3, 1.3, 5.3, 6.6, '', fc='#eff6ff')
    ax.text(5.95, 7.55, 'DIGITAL TWIN', ha='center', fontsize=7.5, fontweight='bold')
    box(ax, 3.5, 5.1, 2.4, 2.1, 'MECHANISTIC CORE\ngene-regulatory network\n(Boolean–ODE hybrid;\nattractors = fate states)', fc='#dcfce7', fs=6)
    box(ax, 6.2, 5.1, 2.2, 2.1, 'biphasic poro-\nviscoelastic FE\n(joint → chondron;\nload → cell stimuli)', fc='#fee2e2', fs=6)
    arrow(ax, (5.9, 6.4), (6.2, 6.4), lw=1.0, ms=8); arrow(ax, (6.2, 5.9), (5.9, 5.9), lw=1.0, ms=8)
    ax.text(6.05, 6.75, 'stimuli', fontsize=5.2, ha='center'); ax.text(6.05, 5.55, 'matrix\nturnover', fontsize=5.2, ha='center')
    box(ax, 3.5, 3.2, 4.9, 1.6, 'LEARNED LAYER\nGRN priors (single-cell foundation models, GNNs)\nimaging-derived phenotypes • FE surrogates & PINNs\ncalibrated uncertainty on every component', fc='#fef3c7', fs=5.8)
    box(ax, 3.5, 1.5, 4.9, 1.4, 'BAYESIAN DATA ASSIMILATION & UQ\nidentifiability • global sensitivity • credible intervals\ncontinuous updating from follow-up data', fc='#ede9fe', fs=5.8)
    arrow(ax, (5.95, 5.1), (5.95, 4.8), lw=1.0, ms=8); arrow(ax, (5.95, 3.2), (5.95, 2.9), lw=1.0, ms=8)
    # outputs
    box(ax, 9.3, 1.3, 2.5, 6.6, '', fc='#f3f4f6')
    ax.text(10.55, 7.55, 'OUTPUTS', ha='center', fontsize=7.5, fontweight='bold')
    for i, t in enumerate(['patient-specific fate\nlandscape (basins, barriers)', 'progression risk\n(with credible interval)', 'in silico trial /\nvirtual cohort enrichment',
                           'loading prescription\n(magnitude, frequency)', 'drug / combination\nrecommendation & timing']):
        box(ax, 9.45, 6.4 - i * 1.15, 2.2, 0.95, t, fc='white', fs=5.8)
    for y in (6.9, 5.3, 3.9, 2.2):
        arrow(ax, (2.6, y), (3.3, y), lw=1.0, ms=8)
    for y in (6.4, 4.5, 2.3):
        arrow(ax, (8.6, y), (9.3, y), lw=1.0, ms=8)
    # feedback
    arrow(ax, (10.55, 1.3), (1.45, 1.3), lw=1.4, ms=10, color='#b91c1c', rad=0.0)
    ax.text(6.0, 0.95, 'feedback: follow-up imaging, outcomes, adverse events → recalibration (ASME V&V 40 / FDA credibility loop)', ha='center', fontsize=6.2, color='#b91c1c')
    ax.text(0.3, 0.35, 'Design principles: fate states as explicit attractors • mechanics as a modeled input, not a covariate • every learned component carries calibrated uncertainty.',
            fontsize=6.0, style='italic', color='#6b7280')
    save(fig, 9)


# ---------------------------------------------------------------- Figure 10: roadmap
def fig10():
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.set_xlim(0, 12.6); ax.set_ylim(0, 8.6); ax.axis('off')
    cols = [('NEAR TERM (1–2 years)', 3.5), ('MID TERM (3–5 years)', 6.6), ('LONG TERM (5–10 years)', 9.7)]
    for t, x in cols:
        box(ax, x - 1.45, 7.6, 2.9, 0.6, t, fc='#111827', tc='white', fs=6.8, bold=True)
    lanes = [('Data standards', '#dbeafe', ['FAIR cartilage atlases;\nSBML/CellML/SED-ML\nmodel deposition', 'minimum reporting\nstandards for chondrocyte\nmechanobiology data', 'federated multi-cohort\ndata commons with\nprovenance']),
             ('Benchmarks', '#dcfce7', ['OAI/MOST segmentation &\nprogression leaderboards;\nstandardized OA scRNA atlas', 'held-out validation sets for\nGRN and FE predictions;\nexplant loading benchmark', 'prospective registry\nlinking twins to outcomes']),
             ('Credibility', '#fef3c7', ['identifiability & UQ as\npublication requirements;\nV&V 40 context of use', 'external validation of\nhybrid twins; credibility\nevidence packages', 'regulatory qualification\nof twin-based endpoints\nand enrichment']),
             ('Clinical translation', '#fee2e2', ['in silico re-analysis of\nFORWARD / lorecivivint\ncohorts; H1–H5 tests', 'in silico trials & virtual\ncohorts for DMOAD design;\nmodel-guided rehabilitation', 'patient-specific twins\nin care pathways;\nmechanotherapeutic timing'])]
    y = 6.7
    for name, fc, items in lanes:
        ax.add_patch(Rectangle((0.2, y - 0.8), 12.2, 1.6, fc=fc, ec='none', alpha=0.45))
        ax.text(1.1, y, name.replace(' ', '\n') if name == 'Clinical translation' else name, fontsize=6.8, fontweight='bold', va='center', ha='center')
        for (t, x), it in zip(cols, items):
            box(ax, x - 1.4, y - 0.65, 2.8, 1.3, it, fc='white', fs=5.6)
        y -= 1.7
    # gates
    for x in (5.05, 8.15):
        ax.plot([x, x], [0.5, 7.4], color='#b91c1c', lw=1.2, ls='--')
        ax.add_patch(Polygon([[x, 7.55], [x + 0.22, 7.3], [x, 7.05], [x - 0.22, 7.3]], fc='#b91c1c', ec='k', lw=0.5))
        ax.text(x, 0.3, 'external validation gate', ha='center', fontsize=6, color='#b91c1c')
    save(fig, 10)


for fn in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10):
    fn()
    print(fn.__name__, 'ok')
json.dump(dims, open(f'{OUT}/dims.json', 'w'))
print(dims)
