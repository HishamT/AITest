"""Command-line interface for running the vortex panel solver."""
from __future__ import annotations

import argparse
from pathlib import Path

from airfoil_cfd.solver import SimulationParameters, VortexPanelSolver


def positive_float(value: str) -> float:
    num = float(value)
    if num <= 0:
        raise argparse.ArgumentTypeError("Value must be positive.")
    return num


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Solve potential flow over an airfoil using a vortex panel method.",
    )
    parser.add_argument("dat_file", type=Path, help="Path to the DAT coordinate file")
    parser.add_argument(
        "--speed",
        type=positive_float,
        default=50.0,
        help="Freestream speed magnitude (m/s).",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=5.0,
        help="Angle of attack in degrees (positive nose-up).",
    )
    parser.add_argument(
        "--rho",
        type=positive_float,
        default=1.225,
        help="Freestream density (kg/m^3).",
    )
    parser.add_argument(
        "--mu",
        type=positive_float,
        default=1.7894e-5,
        help="Dynamic viscosity (Pa·s).",
    )
    parser.add_argument(
        "--chord",
        type=positive_float,
        default=None,
        help="Reference chord length (m). Defaults to the geometric chord.",
    )
    parser.add_argument(
        "--reynolds",
        type=positive_float,
        default=None,
        help="Override the Reynolds number used in the skin-friction model.",
    )
    parser.add_argument(
        "--skin-model",
        choices=["laminar", "turbulent"],
        default="turbulent",
        help="Skin-friction model used to estimate viscous drag.",
    )
    parser.add_argument(
        "--quadrature",
        type=int,
        default=16,
        help="Order of Gauss-Legendre quadrature used for panel integrations.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    solver = VortexPanelSolver.from_dat_file(args.dat_file, quadrature_order=args.quadrature)
    params = SimulationParameters(
        freestream_speed=args.speed,
        angle_of_attack_deg=args.alpha,
        rho=args.rho,
        mu=args.mu,
        chord=args.chord,
        reynolds_number=args.reynolds,
        skin_friction_model=args.skin_model,
    )
    result = solver.solve(params)

    print("Results")
    print("========")
    print(f"Reynolds number: {result.reynolds_number:,.1f}")
    print(f"Circulation (per unit span): {result.circulation:.4f} m^2/s")
    print(f"Lift coefficient (circulation): {result.cl_circulation:.4f}")
    print(f"Lift coefficient (pressure):    {result.cl_pressure:.4f}")
    print(f"Pressure drag coefficient:      {result.cd_pressure:.5f}")
    print(f"Profile drag coefficient:       {result.cd_profile:.5f}")
    print(f"Total drag coefficient:         {result.cd_total:.5f}")
    print(f"Moment coefficient (c/4):       {result.cm_quarter_chord:.5f}")

    cp_min = min(result.pressure_coefficient)
    cp_max = max(result.pressure_coefficient)
    print(f"Cp min/max: {cp_min:.3f} / {cp_max:.3f}")


if __name__ == "__main__":
    main()
