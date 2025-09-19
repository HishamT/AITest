"""Light-weight vortex panel solver for airfoil analysis."""
from .airfoil import AirfoilGeometry, load_dat
from .panel import Panel, build_panels
from .solver import SimulationParameters, SimulationResult, VortexPanelSolver

__all__ = [
    "AirfoilGeometry",
    "SimulationParameters",
    "SimulationResult",
    "VortexPanelSolver",
    "build_panels",
    "load_dat",
]
