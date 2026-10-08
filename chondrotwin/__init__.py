"""
chondrotwin
===========

Minimal, fully reproducible computational models accompanying the review
"The Chondrocyte Fate Landscape: A Mechanically Shaped, Computable Attractor
Framework Integrating Signaling Networks, Mathematical Modeling, and Artificial
Intelligence Toward a Digital Twin of the Chondrocyte" (IJMS, Special Issue
"Chondrocyte Fate: Signaling Networks in Development and Disease").

The package implements, in plain NumPy, the building blocks of the CHONDRO-TWIN
framework at the smallest scale at which each idea becomes computable:

* :class:`ToggleSwitch` - the SOX9-RUNX2 bistable switch (ODE; Figure 4 of the paper)
* :class:`GrowthPlateModel` - the IHH-PTHrP reaction-diffusion loop (PDE; Figure 2)
* :class:`BooleanNetwork` / :func:`chondrocyte_network` - a logical fate network whose
  attractors are the fate states of Table 1
* :func:`quasi_potential_langevin` - a computable Waddington landscape (U = -ln P)
* :class:`ConfinedCompression` - a 1-D biphasic (poroelastic) cartilage model
* :class:`MiniTwin` - a toy mechanistic "digital twin" loop coupling tissue mechanics,
  the gene switch, and matrix turnover
* :mod:`chondrotwin.abm` - a multi-agent cartilage model in which one agent is one
  chondrocyte (vectorized; ~3e4 cells per mm^2 x 2 mm slab run in seconds)

All models are illustrative: parameters are dimensionless or order-of-magnitude
and are not fitted to data. See ``docs/models.md`` for the equations.
"""

__version__ = "0.2.0"

from .toggle import ToggleSwitch, FixedPoint
from .growth_plate import GrowthPlateModel, GrowthPlateSolution
from .boolean import BooleanNetwork, Attractor, chondrocyte_network, classify_state, SCENARIOS
from .landscape import schematic_landscape, cross_section, quasi_potential_langevin, DEFAULT_BASINS
from .loading import inverted_u, regime_map, catabolic_bias
from .biphasic import ConfinedCompression, BiphasicSolution
from .twin import MiniTwin, TwinRecord
from .abm import CartilageABM, CartilageSlab, Scenario

__all__ = [
    "ToggleSwitch", "FixedPoint",
    "GrowthPlateModel", "GrowthPlateSolution",
    "BooleanNetwork", "Attractor", "chondrocyte_network", "classify_state", "SCENARIOS",
    "schematic_landscape", "cross_section", "quasi_potential_langevin", "DEFAULT_BASINS",
    "inverted_u", "regime_map", "catabolic_bias",
    "ConfinedCompression", "BiphasicSolution",
    "MiniTwin", "TwinRecord",
    "CartilageABM", "CartilageSlab", "Scenario",
    "__version__",
]
