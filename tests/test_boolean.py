from chondrotwin import chondrocyte_network, classify_state, SCENARIOS
from chondrotwin.boolean import attractor_table


def fates(scenario):
    net = chondrocyte_network()
    return {classify_state(d) for a in net.attractors(scenario) for d in a.as_dicts()}


def test_basal_has_articular_and_hypertrophic_attractors():
    f = fates({})
    assert "articular (SOX9-high)" in f and "hypertrophic (RUNX2-high)" in f


def test_physiological_loading_keeps_articular_state():
    net = chondrocyte_network()
    traj = net.trajectory({"SOX9": 1}, {"LoadPhys": 1})
    assert classify_state(traj[-1]) == "articular (SOX9-high)"
    assert traj[-1]["TRPV4"] == 1 and traj[-1]["Piezo1"] == 0


def test_overload_leads_to_senescence_with_antioxidant_capacity():
    net = chondrocyte_network()
    traj = net.trajectory({"SOX9": 1}, {"LoadOver": 1})
    assert classify_state(traj[-1]) == "senescent"


def test_overload_without_nrf2_leads_to_ferroptosis():
    net = chondrocyte_network()
    traj = net.trajectory({"SOX9": 1}, {"LoadOver": 1, "Nrf2": 0})
    assert classify_state(traj[-1]) == "ferroptotic"


def test_stiff_matrix_dedifferentiates():
    net = chondrocyte_network()
    traj = net.trajectory({"SOX9": 1}, {"Stiff": 1})
    assert classify_state(traj[-1]) == "dedifferentiated (YAP-high)"


def test_wnt_drives_hypertrophy_and_tgfb_rescues():
    net = chondrocyte_network()
    assert classify_state(net.trajectory({"SOX9": 1}, {"WNT": 1})[-1]) == "hypertrophic (RUNX2-high)"
    assert classify_state(net.trajectory({"SOX9": 1}, {"WNT": 1, "TGFb": 1})[-1]) == "articular (SOX9-high)"


def test_inflammation_sensitizes_piezo1_to_physiological_load():
    net = chondrocyte_network()
    traj = net.trajectory({"SOX9": 1}, {"IL1b": 1, "LoadPhys": 1})
    assert traj[-1]["Piezo1"] == 1
    assert classify_state(traj[-1]) == "senescent"


def test_basin_fractions_sum_to_one():
    net = chondrocyte_network()
    for scenario in SCENARIOS.values():
        assert sum(a.basin_size for a in net.attractors(scenario)) == net.n_initial_states("slow") == 8
        assert sum(a.basin_size for a in net.attractors(scenario, initial="viable")) == 256
        assert sum(a.basin_size for a in net.attractors(scenario, initial="all")) == 2 ** 10


def test_growth_plate_scenario_holds_cells_proliferative():
    net = chondrocyte_network()
    traj = net.trajectory({"SOX9": 1}, {"PTHrP": 1, "WNT": 1})
    assert classify_state(traj[-1]) == "proliferative (PTHrP-held)"


def test_basal_viable_cells_cannot_reach_terminal_states():
    f = fates({})
    assert "senescent" not in f and "ferroptotic" not in f


def test_attractor_table_has_rows():
    rows = attractor_table()
    assert len(rows) >= len(SCENARIOS)
