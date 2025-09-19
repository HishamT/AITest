import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from airfoil_cfd.airfoil import load_dat


def test_load_dat_closes_geometry():
    geometry = load_dat("data/naca2412.dat")
    assert geometry.points[0] == geometry.points[-1]
    assert len(geometry.x) == len(geometry.y)
    assert geometry.chord > 0.9
    assert geometry.chord < 1.1
