# Quantum Fisher Information Analysis

This repository provides Julia code to analyze quantum Fisher information (QFI) for quantum states, with a focus on mixed states arising from depolarization and their corresponding metrics like the Bures metric. It includes utilities for eigen-decomposition and QFI calculation using various methods.

## Files

### `qfish_eigendecomp.jl`

This script performs:
- Eigen-decomposition of quantum states
- Analytical and numerical computation of the Quantum Fisher Information (QFI)
- Analysis of mixed quantum states composed of a pure entangled component and a maximally mixed component (depolarization model)

Functions include:
- `qfi_eigendecomp(rho, L)`: Computes QFI using the eigen-decomposition method.
- `depolarized_state(psi, p)`: Constructs a depolarized state from a pure state `psi` and mixing parameter `p`.
- `generate_pure_state(...)`: (If present) Utility for generating pure quantum states, such as GHZ or product states.

### `qfish_bures.jl`

This script includes:
- Computation of the Bures metric and fidelity between quantum states
- Derivation of QFI from the Bures metric
- Perturbative and numerical analysis of parameterized quantum states

Key features:
- `bures_distance(rho, sigma)`: Calculates the Bures distance between two quantum states.
- `fidelity(rho, sigma)`: Computes the Uhlmann fidelity.
- `qfi_bures(rho_func, theta)`: Computes the QFI via second-order expansion of the Bures metric.

## Requirements

- Julia 1.6 or later
- Dependencies:
  - `LinearAlgebra`
  - `SparseArrays`
  - `QuantumOptics` (optional, if any part of your code uses it – otherwise, remove mention)

Install dependencies using:
```julia
using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("SparseArrays")
```
