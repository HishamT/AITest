"""A simple vortex panel solver tailored for airfoil DAT files."""
from __future__ import annotations

from dataclasses import dataclass
from math import atan2, cos, log, log10, pi, sin, sqrt
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from .airfoil import AirfoilGeometry, load_dat
from .panel import Panel, build_panels


@dataclass
class SimulationParameters:
    freestream_speed: float
    angle_of_attack_deg: float
    rho: float = 1.225
    mu: float = 1.7894e-5
    chord: float | None = None
    reynolds_number: float | None = None
    skin_friction_model: str = "turbulent"


@dataclass
class SimulationResult:
    panels: List[Panel]
    gamma: List[float]
    tangential_velocity: List[float]
    pressure_coefficient: List[float]
    circulation: float
    cl_circulation: float
    cl_pressure: float
    cd_pressure: float
    cd_profile: float
    cd_total: float
    cm_quarter_chord: float
    reynolds_number: float


class VortexPanelSolver:
    """Solves potential flow about an airfoil using a vortex panel method."""

    def __init__(self, panels: Iterable[Panel], quadrature_order: int = 16) -> None:
        self.panels = list(panels)
        if len(self.panels) < 4:
            raise ValueError("At least four panels are required for the solver.")
        self._precompute_geometry()
        self._compute_influence_coefficients()

    @classmethod
    def from_geometry(
        cls, geometry: AirfoilGeometry, quadrature_order: int = 16
    ) -> "VortexPanelSolver":
        return cls(build_panels(geometry), quadrature_order)

    @classmethod
    def from_dat_file(
        cls, path: str | Path, quadrature_order: int = 16
    ) -> "VortexPanelSolver":
        geometry = load_dat(path)
        return cls.from_geometry(geometry, quadrature_order)

    def _precompute_geometry(self) -> None:
        self._lengths = [panel.length for panel in self.panels]
        self._control_points = [panel.control_point for panel in self.panels]
        self._tangents = [panel.tangent for panel in self.panels]
        self._normals = [panel.normal for panel in self.panels]
        self._betas = [panel.beta for panel in self.panels]
        self._upper_te_index = 0
        self._lower_te_index = len(self.panels) - 1
        xs = [pt[0] for pt in self._control_points]
        self._chord = max(xs) - min(xs)

    def _compute_influence_coefficients(self) -> None:
        n = len(self.panels)
        self._source_normal = [[0.0 for _ in range(n)] for _ in range(n)]
        self._source_tangential = [[0.0 for _ in range(n)] for _ in range(n)]
        self._vortex_normal = [[0.0 for _ in range(n)] for _ in range(n)]
        self._vortex_tangential = [[0.0 for _ in range(n)] for _ in range(n)]

        for i in range(n):
            xi, yi = self._control_points[i]
            beta_i = self._betas[i]
            for j in range(n):
                if i == j:
                    self._source_normal[i][j] = 0.5
                    self._source_tangential[i][j] = 0.0
                    self._vortex_normal[i][j] = 0.0
                    self._vortex_tangential[i][j] = -0.5
                    continue
                panel_j = self.panels[j]
                beta_j = self._betas[j]
                dx = xi - panel_j.xa
                dy = yi - panel_j.ya
                cos_beta = cos(beta_j)
                sin_beta = sin(beta_j)
                A = -(dx * cos_beta + dy * sin_beta)
                B = dx * dx + dy * dy
                C = sin(beta_i - beta_j)
                D = cos(beta_i - beta_j)
                E_sq = B - A * A
                if E_sq <= 0.0:
                    E = 0.0
                else:
                    E = sqrt(E_sq)
                Sj = panel_j.length
                log_argument = (Sj * Sj + 2.0 * Sj * A + B) / B if B != 0.0 else 1.0
                if log_argument <= 0.0:
                    log_argument = 1e-12
                term1 = log(log_argument)
                term2 = 0.0
                denominator = B + Sj * A
                term2 = atan2(E * Sj, denominator)
                common = 0.5 / pi
                cn = common * (C * term1 + D * term2)
                ct = common * (D * term1 - C * term2)
                self._source_normal[i][j] = cn
                self._source_tangential[i][j] = ct
                self._vortex_normal[i][j] = -ct
                self._vortex_tangential[i][j] = cn

    @staticmethod
    def _vector_dot(a: Tuple[float, float], b: Tuple[float, float]) -> float:
        return a[0] * b[0] + a[1] * b[1]

    def solve(self, params: SimulationParameters) -> SimulationResult:
        v_inf = params.freestream_speed
        if v_inf <= 0.0:
            raise ValueError("Freestream speed must be positive.")
        alpha = params.angle_of_attack_deg * pi / 180.0
        u_inf = (v_inf * cos(alpha), v_inf * sin(alpha))
        n = len(self.panels)

        system_matrix = [row[:] for row in self._vortex_normal]
        rhs = [-self._vector_dot(normal, u_inf) for normal in self._normals]

        kutta_row = [0.0 for _ in range(n)]
        kutta_row[self._upper_te_index] = 1.0
        kutta_row[self._lower_te_index] = 1.0
        system_matrix[-1] = kutta_row
        rhs[-1] = 0.0

        gamma = self._solve_linear_system(system_matrix, rhs)

        tangential_velocity: List[float] = []
        for i in range(n):
            induced = 0.0
            for j in range(n):
                induced += self._vortex_tangential[i][j] * gamma[j]
            tangential_velocity.append(induced + self._vector_dot(self._tangents[i], u_inf))

        cp = [1.0 - (vt / v_inf) ** 2 for vt in tangential_velocity]

        chord = params.chord if params.chord is not None else self._chord
        chord = float(chord)
        circulation = sum(g * l for g, l in zip(gamma, self._lengths))
        cl_circ = 2.0 * circulation / (v_inf * chord)

        cx = 0.0
        cy = 0.0
        for cp_i, normal, length in zip(cp, self._normals, self._lengths):
            cx += -cp_i * normal[0] * length / chord
            cy += -cp_i * normal[1] * length / chord
        cl_pressure = cy * cos(alpha) - cx * sin(alpha)
        cd_pressure = cx * cos(alpha) + cy * sin(alpha)

        xs = [pt[0] for pt in self._control_points]
        x_min = min(xs)
        x_ref = x_min + 0.25 * chord
        y_ref = 0.0
        moment = 0.0
        for cp_i, normal, length, point in zip(
            cp, self._normals, self._lengths, self._control_points
        ):
            arm_x = point[0] - x_ref
            arm_y = point[1] - y_ref
            moment += -cp_i * (arm_x * normal[1] - arm_y * normal[0]) * length
        cm_quarter_chord = moment / (chord**2)

        reynolds = (
            params.reynolds_number
            if params.reynolds_number is not None
            else params.rho * v_inf * chord / params.mu
        )
        if reynolds <= 0.0:
            raise ValueError("Reynolds number must be positive.")

        skin_model = params.skin_friction_model.lower()
        if skin_model == "laminar":
            cf = 1.328 / sqrt(reynolds)
        elif skin_model == "turbulent":
            cf = 0.455 / (log10(reynolds) ** 2.58)
        else:
            raise ValueError(
                "Skin friction model must be either 'laminar' or 'turbulent'."
            )
        cd_profile = cf
        cd_total = cd_pressure + cd_profile

        return SimulationResult(
            panels=self.panels,
            gamma=gamma,
            tangential_velocity=tangential_velocity,
            pressure_coefficient=cp,
            circulation=circulation,
            cl_circulation=cl_circ,
            cl_pressure=cl_pressure,
            cd_pressure=cd_pressure,
            cd_profile=cd_profile,
            cd_total=cd_total,
            cm_quarter_chord=cm_quarter_chord,
            reynolds_number=float(reynolds),
        )

    @staticmethod
    def _solve_linear_system(
        matrix: List[List[float]], rhs: Sequence[float]
    ) -> List[float]:
        n = len(matrix)
        aug = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
        for col in range(n):
            pivot = col
            max_value = abs(aug[pivot][col])
            for row in range(col + 1, n):
                value = abs(aug[row][col])
                if value > max_value:
                    max_value = value
                    pivot = row
            if max_value == 0.0:
                raise ValueError("Singular linear system encountered.")
            if pivot != col:
                aug[col], aug[pivot] = aug[pivot], aug[col]
            pivot_value = aug[col][col]
            for j in range(col, n + 1):
                aug[col][j] /= pivot_value
            for row in range(n):
                if row == col:
                    continue
                factor = aug[row][col]
                for j in range(col, n + 1):
                    aug[row][j] -= factor * aug[col][j]
        return [aug[i][n] for i in range(n)]
