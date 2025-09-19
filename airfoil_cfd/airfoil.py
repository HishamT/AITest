"""Utilities for loading and preprocessing airfoil geometry files."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


@dataclass
class AirfoilGeometry:
    """Represents the discretised surface of an airfoil."""

    x: List[float]
    y: List[float]

    @property
    def points(self) -> List[Tuple[float, float]]:
        return list(zip(self.x, self.y))

    @property
    def chord(self) -> float:
        return max(self.x) - min(self.x)

    def closed(self) -> "AirfoilGeometry":
        if self.points[0] != self.points[-1]:
            x = list(self.x)
            y = list(self.y)
            x.append(self.x[0])
            y.append(self.y[0])
            return AirfoilGeometry(x, y)
        return self


def _read_coordinate_lines(lines: Iterable[str]) -> Tuple[List[float], List[float]]:
    xs: List[float] = []
    ys: List[float] = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        numeric: List[float] = []
        for part in parts:
            try:
                numeric.append(float(part))
            except ValueError:
                continue
        if len(numeric) < 2:
            continue
        xs.append(numeric[0])
        ys.append(numeric[1])
    if len(xs) < 4:
        raise ValueError("Airfoil file must contain at least four coordinate pairs.")
    return xs, ys


def load_dat(path: Path | str) -> AirfoilGeometry:
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"Airfoil file not found: {file_path}")

    with file_path.open("r", encoding="utf8") as handle:
        xs, ys = _read_coordinate_lines(handle)

    geometry = AirfoilGeometry(xs, ys).closed()
    return geometry
