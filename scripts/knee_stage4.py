#!/usr/bin/env python
"""Stage 4: is the senescence failure real? Re-measurement with SenMayo, revised senescence logic, calibration on
GSE104782 OA1-5 (SenMayo target), pre-registered test on GSE117999 / GSE169077, and the repaired (deviation) GSE114007
analysis with raw counts.

Usage:  python scripts/knee_stage4_calibrate.py      # once; freezes the parameters
        python scripts/knee_stage4.py [--out docs/knee_stage4.md]
"""
import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from chondrotwin import datasets
from chondrotwin.abm import Rates
from chondrotwin.knee import calibrate as C
from chondrotwin.knee import independent3 as I3
from chondrotwin.knee import loads as L
from chondrotwin.knee import senescence as S
from chondrotwin.knee.analysis import mixture_test, observed_losses
from chondrotwin.knee.erosion import dose_response
from chondrotwin.knee.geometry import tibial_geometry
from chondrotwin.knee.independent2 import model_log2 as bulk_model_log2
from chondrotwin.knee.report import table

ROOT = Path(__file__).resolve().parents[1]
DERIVED = datasets.data_dir().parent / "derived"
H_MED, H_LAT, YEAR = 2.2, 2.6, 5
KEYS = ["damage_D50", "c_over_max", "k_cat", "switch_bias_sd", "senescence_competent_fraction",
        "sen_overload_weight", "senescent_death_rate", "baseline_runx2"]
STAGE1 = dict(damage_D50=0.0, c_over_max=0.9, k_cat=0.0, switch_bias_sd=0.0, senescence_competent_fraction=1.0,
              sen_overload_weight=1.0, senescent_death_rate=0.0, baseline_runx2=0.0)
STAGE2 = dict(STAGE1, c_over_max=0.15, k_cat=0.5)
out = []


def md(s=""):
    print(s)
    out.append(s)


def pooled(tables, p):
    m = np.ones(len(tables), bool)
    for k, v in p.items():
        m &= np.isclose(tables[k].astype(float), v)
    d = tables[m]
    cols = ["runx2_high", "senescent", "clustered", "catabolic", "mean_S", "mean_R", "mean_cat"]
    w = d["n"].fillna(0.0)
    res = d[cols].mul(w, axis=0).groupby(d["stage"]).sum().div(w.groupby(d["stage"]).sum(), axis=0)
    res["n"] = w.groupby(d["stage"]).sum().astype(int)
    return res.reindex(range(5))


def registration():
    pre = (ROOT / "docs" / "preregistration_stage4.md").read_text(encoding="utf-8")
    reg = re.search(r"independent3\.py` at registration: `([0-9a-f]{64})`", pre).group(1)
    frz = re.search(r"stage4_frozen_parameters\.json`: `([0-9a-f]{64})`", pre).group(1)
    h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()   # noqa: E731
    return reg == h(ROOT / "chondrotwin" / "knee" / "independent3.py"), frz == h(ROOT / "docs" / "stage4_frozen_parameters.json")


def gse114007_raw(senmayo, n_boot=2000, seed=0):
    """Deviation from the stage-3 protocol: raw counts instead of the (batch-corrected) 'normalized counts'."""
    x = pd.ExcelFile(datasets.require("gse114007", "GSE114007_raw_counts.xlsx"))
    nrm = x.parse("Normal").set_index("symbol").groupby(level=0).sum()
    oa = x.parse("OA").set_index("symbol").groupby(level=0).sum()
    common = nrm.index.intersection(oa.index)
    cpm_n = nrm.loc[common] / nrm.loc[common].sum() * 1e6
    cpm_o = oa.loc[common] / oa.loc[common].sum() * 1e6
    rng = np.random.default_rng(seed)
    rows = {}
    for g in ("SOX9", "RUNX2", "MMP13", "CDKN2A"):
        a, b = cpm_o.loc[g].to_numpy(float), cpm_n.loc[g].to_numpy(float)
        bt = [np.log2(rng.choice(a, len(a)).mean() / rng.choice(b, len(b)).mean()) for _ in range(n_boot)]
        rows[g] = dict(diff=float(np.log2(a.mean() / b.mean())), ci_low=float(np.percentile(bt, 2.5)),
                       ci_high=float(np.percentile(bt, 97.5)), units="log2")
    lg = np.log2(pd.concat([cpm_o, cpm_n], axis=1) + 1)
    genes = [g for g in senmayo if g in lg.index and lg.loc[g].std() > 0]
    z = lg.loc[genes]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1, ddof=1), axis=0)
    sc = z.mean(axis=0)
    d, lo, hi = I3.welch(sc[cpm_o.columns], sc[cpm_n.columns])
    rows["SenMayo"] = dict(diff=d, ci_low=lo, ci_high=hi, units="z")
    res = pd.DataFrame(rows).T
    res.attrs.update(n_oa=cpm_o.shape[1], n_normal=cpm_n.shape[1], min_value=float(min(nrm.min().min(), oa.min().min())))
    return res


def c1_table(readout: pd.DataFrame, models: dict) -> pd.DataFrame:
    rows = []
    for g, r in readout.iterrows():
        row = {"readout": g, "data": float(r["diff"]), "95 % low": float(r["ci_low"]), "95 % high": float(r["ci_high"]),
               "units": r["units"]}
        for name, m in models.items():
            row[f"{name} model log2"] = m[g]
            row[f"{name} C1"] = I3.criterion_c1(float(r["diff"]), float(r["ci_low"]), float(r["ci_high"]), m[g])
        rows.append(row)
    return pd.DataFrame(rows).set_index("readout")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "docs" / "knee_stage4.md"))
    a = ap.parse_args()
    t0 = time.perf_counter()
    proto_ok, frozen_ok = registration()
    P4 = json.loads((ROOT / "docs" / "stage4_frozen_parameters.json").read_text())["stage4_parameters"]
    sel = pd.read_csv(DERIVED / "stage4_selection.csv")
    tabs = pd.read_csv(DERIVED / "stage4_selection_tables.csv")
    grid = pd.read_csv(DERIVED / "stage4_grid.csv")
    summ = sel.groupby(KEYS)[["loss_cal", "loss_val"]].agg(["mean", "std"]).sort_values(("loss_cal", "mean"))
    same = all(np.isclose(P4[k], STAGE2[k]) for k in KEYS)
    PARAMS = {"stage 1": STAGE1, "stage 2 (= 3)": STAGE2} if same else {"stage 1": STAGE1, "stage 2 (= 3)": STAGE2, "stage 4": P4}
    st = {n: pooled(tabs, p) for n, p in PARAMS.items()}
    sel_name = "stage 2 (= 3)" if same else "stage 4"

    md("# Knee cartilage simulation from open data — stage 4 results")
    md()
    md("_Generated by `scripts/knee_stage4.py`; do not edit by hand. Pre-registration: `docs/preregistration_stage4.md`._")
    md()
    md(f"Registration check: protocol unchanged — **{proto_ok}**; frozen parameters unchanged — **{frozen_ok}**.")
    md()

    # ---------------------------------------------------------------- 1. re-measurement
    md("## 1. Is the senescence failure real? Re-measurement with SenMayo (previously seen data; exploratory)")
    md()
    senmayo = I3.senmayo_genes()
    ji = S.ji_scores()
    md(f"SenMayo (Saul et al., Nat Commun 2022): {len(senmayo)} genes, mostly SASP; it contains neither CDKN2A nor CDKN1A, "
       f"so it is an independent readout of the one used in stages 1–3. A cell is senescence-high if its score exceeds "
       f"the 90th percentile of the reference cells.")
    md()
    md("GSE104782, senescence-high fraction relative to stage 0 (95 % patient bootstrap):")
    md()
    md(table(S.stage_ratios(ji), "{:.2f}"))
    md()
    md("GSE152805, damaged medial vs minimally damaged lateral (reference: lateral cells of the same patient):")
    md()
    md(table(S.gse152805_senmayo(), "{:.3f}"))
    md()
    (m57, lo57, hi57, n57), _ = S.gse57218_senmayo()
    md(f"GSE57218, paired OA − preserved SenMayo score: {m57:+.3f} z [{lo57:+.3f}, {hi57:+.3f}] ({n57} pairs).")
    md()
    md("All three agree with the single-marker picture: senescence does not rise before the last stage in GSE104782, "
       "rises 2–5-fold in damaged vs intact tissue in GSE152805 and not detectably in GSE57218 — the model's early, "
       "near-complete senescence is a real failure, not a readout artefact.")
    md()

    # ---------------------------------------------------------------- 2. calibration
    md("## 2. Revised senescence logic and calibration (GSE104782 OA1–5, SenMayo target)")
    md()
    md("* `sen_overload_weight` w — overload-driven senescence is multiplied by w + (1 − w) · inflammation: with w < 1 "
       "overload alone induces less senescence and inflammatory co-signals are needed (cf. stage 2 for the switch);")
    md("* `senescent_death_rate` — senescent cells die (apoptosis) at this daily hazard.")
    md()
    md(f"Grid of {len(grid)} sets (with `c_over_max` 0.15 / 0.3, `k_cat` 0.5); 5 best re-run with seeds 2–3 together with "
       f"the stage-1 and stage-2 models; all scored with the same SenMayo-based loss:")
    md()
    rows = []
    for n, p in PARAMS.items():
        r = summ.loc[tuple(float(p[k]) for k in KEYS)]
        rows.append({"model": n, "loss OA1–5": f"{r[('loss_cal', 'mean')]:.0f} ± {r[('loss_cal', 'std')]:.0f}",
                     "loss OA6–10": f"{r[('loss_val', 'mean')]:.0f} ± {r[('loss_val', 'std')]:.0f}"})
    md(table(pd.DataFrame(rows).set_index("model")))
    md()
    md("Selected (frozen) parameters: " + ", ".join(f"{k} = {P4[k]:g}" for k in KEYS) + ".")
    md()
    if same:
        md("**Neither senescence revision was selected; the stage-4 model is the stage-2 model.**")
        md()
        cal = C.data_targets(C.ji_cell_table(), C.CAL_PATIENTS, pairs=C.PAIRS4)
        rows = []
        for name, p in (("stage 2", STAGE2), ("gated, w = 0.3", dict(STAGE2, sen_overload_weight=0.3))):
            t = pooled(tabs, p)
            err = (C.model_log_ratios(t, C.PAIRS4) - cal["log_ratio"]) ** 2 / cal["var"]
            rows.append({"model": name, **{f"loss {k}": float(v) for k, v in err.sum().items()},
                         "senescent, stage 0": t.loc[0, "senescent"], "senescent, stage 4": t.loc[4, "senescent"],
                         "RUNX2-high, stage 4": t.loc[4, "runx2_high"]})
        md("Why (calibration patients, loss by readout):")
        md()
        md(table(pd.DataFrame(rows).set_index("model"), "{:.3g}"))
        md()
        md("Gating lowers senescence in intact tissue but not in degraded tissue, where cytokines are high anyway, so the "
           "damaged/intact ratio gets *worse*; and with less senescence there is less cytokine, so fewer cells switch "
           "to RUNX2-high and the catabolic (MMP13) fit deteriorates. In this architecture the senescent cells' SASP is "
           "the only inflammatory source inside cartilage, which couples senescence to catabolism; the data show "
           "catabolism without massive senescence, pointing to an inflammatory source outside cartilage (synovium).")
        md()
    cells = C.ji_cell_table()
    val = C.data_targets(cells, C.VAL_PATIENTS, pairs=C.PAIRS4)
    rows = []
    for gname in ("SENMAYO", "RUNX2", "MMP13", "SOX9"):
        for s in (1, 2, 3, 4):
            rows.append({"readout / stage": f"{gname} {s}", "data OA6–10": val["log_ratio"].loc[s, gname],
                         "± sd": float(np.sqrt(val["var"].loc[s, gname])),
                         **{n: C.model_log_ratios(t, C.PAIRS4).loc[s, gname] for n, t in st.items()}})
    md("Validation patients, log-ratio stage s / stage 0:")
    md()
    md(table(pd.DataFrame(rows).set_index("readout / stage"), "{:.2f}"))
    md()
    md("Live cells by local stage (selected model):")
    md()
    md(table(st[sel_name][["runx2_high", "senescent", "mean_S", "mean_R", "n"]]))
    md()

    # ---------------------------------------------------------------- 3. independent test
    md("## 3. Pre-registered independent test — OA vs non-OA knee cartilage")
    md()
    models = {}
    for n, t in st.items():
        m = bulk_model_log2(t)
        m["SenMayo"] = m["CDKN2A"]
        models[n] = m
    res = {}
    for name, loader in (("GSE117999 (primary)", I3.gse117999), ("GSE169077 (secondary)", I3.gse169077)):
        expr, groups = loader()
        r = I3.readouts(expr, groups, senmayo)
        tb = c1_table(r, models)
        res[name] = tb
        md(f"**{name}**: {r.attrs['n_oa']} OA vs {r.attrs['n_non']} non-OA samples; SenMayo genes present "
           f"{r.attrs['senmayo_present']} / {r.attrs['senmayo_total']}.")
        md()
        md(table(tb, "{:.2f}"))
        md()
    for name, tb in res.items():
        md(f"* {name}, C1: " + ", ".join(f"{n} {int(tb[f'{n} C1'].sum())} / {len(tb)}" for n in models))
    md()
    md("GSE117999's non-OA group are knees with a meniscus tear (cartilage from the non-weight-bearing notch); "
       "GSE169077 compares pools of five donors each.")
    md()

    # ---------------------------------------------------------------- 4. tissue
    md("## 4. Tissue test (OAIZIB-CM)")
    md()
    if same:
        md("Unchanged: the selected model is the stage-2 model (see `docs/knee_stage2.md`).")
        md()
    else:
        cache = DERIVED / "abm_dose_response_stage4.csv"
        if cache.exists():
            dose4 = pd.read_csv(cache)
        else:
            r4 = Rates(**{k: P4[k] for k in KEYS if k != "baseline_runx2"})
            sig = np.round(np.arange(0.06, 0.3001, 0.01), 3)
            dose4 = pd.concat([dose_response(sig, h, years=YEAR, seeds=(0, 1), rates=r4,
                                             initial_hypertrophic_fraction=P4["baseline_runx2"]) for h in (H_MED, H_LAT)],
                              ignore_index=True)
            dose4.to_csv(cache, index=False)
        geo = tibial_geometry()
        loads = L.load_all()
        mr = L.summary(loads, L.calibrate_condylar_distance(loads)).medial_ratio.values
        rows = []
        for name, dose in (("stage 2", pd.read_csv(DERIVED / "abm_dose_response_stage2.csv")), ("stage 4", dose4)):
            for yr in (2, YEAR):
                r = mixture_test(geo, mr, dose, H_MED, H_LAT, yr)
                ok = r.get("calibrated", False)
                rows.append({"model": name, "ABM years": yr, "load scale (MPa)": r["load_scale_MPa"] if ok else float("nan"),
                             "predicted lateral loss (mm)": r["predicted_mean_lateral_loss_mm"] if ok else float("nan"),
                             "note": "" if ok else "no erosion at any stress yet"})
        md(table(pd.DataFrame(rows).set_index("model")))
        md()
        obs = observed_losses(geo).set_index(["split", "compartment"])
        md(f"Observed lateral loss: train {obs.loc[('train', 'lateral'), 'loss_mm']:.3f} ± {obs.loc[('train', 'lateral'), 'se_mm']:.3f}, "
           f"test {obs.loc[('test', 'lateral'), 'loss_mm']:.3f} ± {obs.loc[('test', 'lateral'), 'se_mm']:.3f} mm.")
        md()

    # ---------------------------------------------------------------- 5. GSE114007 repair (deviation)
    md("## 5. Repair of the invalid stage-3 secondary test (GSE114007) — deviation from the stage-3 protocol")
    md()
    raw = gse114007_raw(senmayo)
    md(f"Raw counts (`GSE114007_raw_counts.xlsx`; smallest value {raw.attrs['min_value']:.0f}, i.e. no negative entries) "
       f"instead of the registered 'normalized counts'; CPM over shared genes; {raw.attrs['n_oa']} OA vs "
       f"{raw.attrs['n_normal']} normal. The directions of the four markers had been glimpsed in stage 3 (normalized "
       f"values), so this is a repaired secondary analysis, not an independent test.")
    md()
    md(table(c1_table(raw, models), "{:.2f}"))
    md()

    # ---------------------------------------------------------------- cross-dataset view (post hoc synthesis)
    from chondrotwin.knee import independent2 as I2
    pairs = I2.gse57218_pairs()
    g57 = {g: I2.t_interval(pairs[g].to_numpy())[:3] for g in ("SOX9", "RUNX2", "MMP13", "CDKN2A")}
    g57["SenMayo"] = (m57, lo57, hi57)
    views = {"GSE57218: OA-affected vs preserved, same joint": {k: v for k, v in g57.items()},
             "GSE117999: OA vs meniscus-tear knees": {g: tuple(res["GSE117999 (primary)"].loc[g, ["data", "95 % low", "95 % high"]]) for g in res["GSE117999 (primary)"].index},
             "GSE169077: OA vs normal (pools)": {g: tuple(res["GSE169077 (secondary)"].loc[g, ["data", "95 % low", "95 % high"]]) for g in res["GSE169077 (secondary)"].index},
             "GSE114007: OA vs normal (raw counts)": {g: (float(raw.loc[g, "diff"]), float(raw.loc[g, "ci_low"]), float(raw.loc[g, "ci_high"])) for g in raw.index}}

    def fmt(t):
        d, lo, hi = t
        star = " *" if (lo > 0 or hi < 0) else ""
        return f"{d:+.2f}{star}"
    cross = pd.DataFrame({k: {g: fmt(v[g]) for g in ("SOX9", "RUNX2", "MMP13", "CDKN2A", "SenMayo")} for k, v in views.items()}).T
    cross["model (stage 2) log2"] = ""
    m2 = models["stage 2 (= 3)"]
    cross.loc["model: damaged vs intact", :] = [f"{m2[g]:+.2f}" for g in ("SOX9", "RUNX2", "MMP13", "CDKN2A", "SenMayo")] + [""]
    cross = cross.drop(columns="model (stage 2) log2")
    md("## 6. Cross-dataset view (post hoc synthesis)")
    md()
    md("Differences (log2; SenMayo in z) with * where the 95 % interval excludes 0:")
    md()
    md(table(cross))
    md()
    md("MMP13 is strongly higher only when OA knees are compared with *normal* knees; within OA joints it is not higher "
       "(slightly lower in affected than in preserved cartilage) and against meniscus-tear knees it does not differ. Catabolism behaves as a property of the whole OA "
       "joint, not of the local lesion — consistent with a joint-level (synovial) inflammatory source and inconsistent "
       "with the model, which localises catabolism to damaged tissue via the SASP of senescent cells. No senescence "
       "readout shows the large increase the model predicts.")
    md()

    # ---------------------------------------------------------------- 7. conclusions
    md("## 7. Conclusions")
    md()
    prim = res["GSE117999 (primary)"]
    md("* **The senescence failure is real.** A 124-gene SASP/senescence score that shares no gene with CDKN2A gives "
       "the same picture as p16 in three datasets: late, moderate senescence in damaged tissue.")
    if same:
        md("* **The two revisions of the senescence logic were not selected** by the calibration; the model is unchanged "
           "from stage 2.")
    else:
        md("* **Selected revision:** " + ", ".join(f"{k} = {P4[k]:g}" for k in ("sen_overload_weight", "senescent_death_rate")) +
           f"; validation loss {summ.loc[tuple(float(STAGE2[k]) for k in KEYS)][('loss_val', 'mean')]:.0f} (stage 2) → "
           f"{summ.loc[tuple(float(P4[k]) for k in KEYS)][('loss_val', 'mean')]:.0f} (stage 4) under the SenMayo loss.")
    md("* **Mechanism:** senescence and catabolism are coupled through the SASP, the model's only intra-cartilage "
       "inflammatory source; reducing senescence removes the catabolic co-signal. Across four datasets catabolism is a "
       "joint-level feature (OA vs normal knees), not a lesion-level one (within OA joints). The next hypothesis to test "
       "is a joint-level (synovial) inflammatory input that drives catabolism throughout the joint without senescence.")
    md("* **Pre-registered test (GSE117999):** " + ", ".join(f"{n} {int(prim[f'{n} C1'].sum())} / {len(prim)}" for n in models) +
       " readouts consistent (C1).")
    md()
    md(f"_Runtime {time.perf_counter() - t0:.0f} s._")
    Path(a.out).write_text("\n".join(out) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    main()
