# Dish Antenna Sidelobe Minimization

This project explores methods for minimizing the Side-Lobe Levels (SLL) of a dish antenna by placing and weighting auxiliary dipoles. The approach is based on linear algebra and uses the LSQR algorithm to solve an underdetermined system of equations that models the electromagnetic field contributions from these dipoles. The goal is to cancel or suppress sidelobes while preserving the main beam characteristics.

The MATLAB code simulates the full radiation pattern, evaluates dipole placements, and applies numerical solvers to optimize dipole currents.

📚 Full MATLAB documentation is available here:  
👉 [https://gurenig.github.io/PO_matlab_docs/](https://gurenig.github.io/PO_matlab_docs/)

---

## Features
- Far-field calculation of a parabolic dish antenna
- Dipole array superposition with custom weights
- SLL minimization using LSQR
- Beamwidth and peak analysis
- 2D and 3D radiation pattern visualization

## Requirements
- MATLAB R2023 or later (tested)
- No proprietary toolboxes required (uses built-in functions)

## License
MIT
