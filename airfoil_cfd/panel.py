"""Panel geometry utilities for the vortex panel solver."""
from __future__ import annotations

from dataclasses import dataclass
from math import atan2
from typing import Iterable, List, Tuple

from .airfoil import AirfoilGeometry


@dataclass
class Panel:
    xa: float
    ya: float
    xb: float
    yb: float

    def __post_init__(self) -> None:
        dx = self.xb - self.xa
        dy = self.yb - self.ya
        self.length = (dx**2 + dy**2) ** 0.5
        if self.length == 0.0:
            raise ValueError("Degenerate panel with zero length encountered.")
        tx = dx / self.length
        ty = dy / self.length
        self.tangent = (tx, ty)
        self.beta = atan2(ty, tx)
        # The outward normal is obtained by rotating the tangent vector
        # clockwise (negative 90 degrees). The input coordinate ordering
        # therefore determines the orientation automatically.
        self.normal = (ty, -tx)
        self.control_point = (
            0.5 * (self.xa + self.xb),
            0.5 * (self.ya + self.yb),
        )

    def vector_to(self, point: Tuple[float, float]) -> Tuple[float, float]:
        return (point[0] - self.control_point[0], point[1] - self.control_point[1])


def build_panels(geometry: AirfoilGeometry) -> List[Panel]:
    points = geometry.points
    panels: List[Panel] = []
    for start, end in zip(points[:-1], points[1:]):
        panels.append(Panel(start[0], start[1], end[0], end[1]))
    return panels


def panel_centroid_array(panels: Iterable[Panel]) -> List[Tuple[float, float]]:
    return [panel.control_point for panel in panels]
