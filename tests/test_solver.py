import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from airfoil_cfd.solver import SimulationParameters, VortexPanelSolver


def test_lift_increases_with_angle():
    solver = VortexPanelSolver.from_dat_file("data/naca2412.dat")
    params_zero = SimulationParameters(
        freestream_speed=40.0,
        angle_of_attack_deg=0.0,
        chord=1.0,
    )
    params_high = SimulationParameters(
        freestream_speed=40.0,
        angle_of_attack_deg=5.0,
        chord=1.0,
    )
    result_zero = solver.solve(params_zero)
    result_high = solver.solve(params_high)
    assert result_high.cl_circulation > result_zero.cl_circulation
    assert result_high.cl_pressure > result_zero.cl_pressure
    assert abs(result_high.cd_pressure) < 0.1
