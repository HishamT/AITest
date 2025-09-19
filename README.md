# Airfoil CFD Solver

This project implements a two-dimensional vortex panel solver that can ingest
airfoil coordinate files in DAT format. It provides utilities for loading
geometry, running a potential-flow simulation, and estimating aerodynamic
coefficients such as lift, drag, and moments.

## Features

- DAT reader tolerant to headers and comments.
- Constant-strength vortex panel method with Kutta condition enforcement.
- CLI for running analyses against supplied airfoil files.
- Basic estimation of viscous (profile) drag using skin-friction correlations.
- Pure-Python implementation without third-party dependencies.

## Usage

1. Place your airfoil coordinate file (``*.dat``) inside the repository.
2. Run the solver from the command line:

```bash
python scripts/run_solver.py data/naca2412.dat --speed 50 --alpha 5
```

The script prints Reynolds number, circulation, lift and drag coefficients, and
basic surface pressure information.

## Testing

Install no additional dependencies are required. Execute the tests with:

```bash
pytest
```
