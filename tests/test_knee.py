import numpy as np
import pandas as pd
import pytest

from chondrotwin import ConfinedCompression, datasets
from chondrotwin.knee import erosion
from chondrotwin.knee.bridge import stress_ratio_lateral_to_medial
from chondrotwin.knee.geometry import plate_metrics
from chondrotwin.knee.loads import GaitLoads

HAVE_ORTHOLOAD = all(datasets.verify(f) for f in datasets.FILES if f.source == "orthoload")
HAVE_GSE = all(datasets.verify(f) for f in datasets.FILES if f.source == "gse104782")


def test_plate_metrics_on_synthetic_plate():
    # R x A x S volume, 1 x 1 x 0.5 mm voxels; a 20 x 20 plate, 4 voxels (2 mm) thick, with a 4 x 4 hole
    lab = np.zeros((30, 30, 20), dtype=np.uint8)
    lab[5:25, 5:25, 8:12] = 4
    lab[13:17, 13:17, :] = 0
    m = plate_metrics(lab, (1.0, 1.0, 0.5), ("R", "A", "S"), 4)
    assert m["area_mm2"] == 400 - 16
    assert m["denuded_mm2"] == 16
    assert np.isclose(m["thickness_mean_mm"], 2.0)
    assert m["central_thickness_mm"] < 2.0           # the hole is in the central region and counts as zero


def test_medial_lateral_split():
    t = np.linspace(0, 1, 11); fz = np.full(11, 2000.0)
    centred = GaitLoads("x", 1000.0, 0.0, t, fz, np.zeros(11))
    f_med, f_lat = centred.split(0.047)
    assert np.allclose(f_med, 1000.0) and np.allclose(f_lat, 1000.0)
    medial_shift = GaitLoads("x", 1000.0, 5.0, t, fz, np.full(11, -20.0))    # My < 0: force shifted medially
    assert medial_shift.medial_ratio(0.047) > 0.5
    extreme = GaitLoads("x", 1000.0, 5.0, t, fz, np.full(11, -500.0))
    assert np.allclose(extreme.split(0.047)[0], fz)                         # clipped to the axial force


def test_stress_ratio():
    assert np.isclose(stress_ratio_lateral_to_medial(0.5, 1000, 1000), 1.0)
    assert stress_ratio_lateral_to_medial(0.7, 1100, 1000) < 0.5


def test_invert_dose_response():
    s = np.array([0.06, 0.07, 0.08, 0.09]); e = np.array([0.0, 0.0, 0.2, 0.6])
    assert np.isclose(erosion.invert(s, e, 0.4), 0.085)
    assert np.isnan(erosion.invert(s, e, 1.0))


def test_cyclic_loading_equals_time_averaged_stress():
    """The bridge assumption: cartilage consolidation integrates the gait cycle (tau >> 1 s)."""
    tis = ConfinedCompression(h=0.5e-3)                     # tau = 500 s, keeps the test fast
    wave = lambda t: 1e5 * (1 + np.sin(2 * np.pi * t))      # noqa: E731  1 Hz, mean 0.1 MPa
    cyc = tis.solve_history(wave, t_end=3 * tis.tau, dt_max=0.025, n_out=20)
    const = tis.solve(1e5, t_end=3 * tis.tau, n_out=20)
    assert np.isclose(cyc.sigma0, 1e5, rtol=1e-3)
    assert abs(cyc.u[-1, -1] - const.u[-1, -1]) < 0.02 * abs(const.u[-1, -1])


@pytest.mark.skipif(not HAVE_ORTHOLOAD, reason="OrthoLoad data not downloaded (python scripts/fetch_data.py orthoload)")
def test_orthoload_condylar_distance_is_plausible():
    from chondrotwin.knee import loads as L
    ld = L.load_all()
    l = L.calibrate_condylar_distance(ld)
    assert 0.040 <= l <= 0.055
    s = L.summary(ld, l)
    assert 0.6 < s.medial_ratio.mean() < 0.75 and s.loc["K6L", "medial_ratio"] == s.medial_ratio.min()


@pytest.mark.skipif(not HAVE_GSE, reason="GSE104782 not downloaded (python scripts/fetch_data.py gse104782)")
def test_ji2019_design():
    from chondrotwin.knee.cells import ji2019_cells
    c = ji2019_cells()
    assert len(c) == 1600 and c.patient.nunique() == 10
    assert (c.groupby(["patient", "stage"]).size() == 32).all()


def test_uq_threshold_scales_with_misspecification():
    from chondrotwin.knee import uq
    assert np.isclose(uq.chi2_threshold(16, 14.0), 5.991, atol=1e-3)        # well specified: plain chi2_2
    assert np.isclose(uq.chi2_threshold(16, 140.0), 5.991 * 10, atol=1e-2)  # 10x over-dispersed


def test_uq_stage_binning():
    from chondrotwin.knee import uq
    cells = pd.DataFrame({"M": [0.95, 0.7, 0.5, 0.3, 0.1], "runx2_high": [0, 0, 1, 1, 1], "senescent": [0, 1, 1, 1, 1],
                          "clustered": 0, "catabolic": [0, 0, 1, 1, 1], "S": [2, 2, 1, 1, 0.1], "R": [0.1, 0.1, 1, 2, 2], "cat": 0.0})
    st = uq.stage_table(cells, uq.BINS["0.2/0.4/0.6/0.8 (default)"])
    assert list(st["n"]) == [1, 1, 1, 1, 1] and st.loc[4, "mean_S"] == 0.1
