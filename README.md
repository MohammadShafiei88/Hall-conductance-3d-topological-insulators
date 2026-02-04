# 3D Topological Insulator Hall Conductance Calculator

This repository contains a **C program** to compute the **Hall conductance of 3D topological insulators** using the **Green's Function formalism**. The program constructs the system Hamiltonian, calculates Green’s functions including lead self-energies, and computes transmission and Hall conductance for a cubic lattice system.

---

## Features

- Simulates a **3D topological insulator** on a cubic lattice.
- Computes **Green's functions** for the device region with left and right leads.
- Implements **lead self-energy contributions**.
- Constructs **gamma matrices** and calculates **transmission probabilities**.
- Outputs **energy-resolved Hall conductance**.
- Fully written in **C** for high performance.

---

## System Description

- **Lattice dimensions:** `Nx × Ny × Nz`
- **Basis:** 4-component per site
- **Hamiltonian terms:**
  - Onsite energies (`MM0`)
  - Nearest-neighbor hopping in x, y, z directions (`Tx`, `Ty`, `Tz`)
  - Zeeman splitting at top and bottom surfaces (`delta_top`, `delta_bottom`)
  - Spin-orbit and SIA effects (`v_sia`)
- **Leads:** Left and right semi-infinite leads via self-energies
- **Transmission:** Computed as `Tr(Γ_left * G * Γ_right * G†)`

---

## Requirements

- **Compiler:** GCC or any C compiler with `<complex.h>` support  
- **BLAS/LAPACK:** Required for matrix inversion (`cgetrf_`, `cgetri_`)  
- **OS:** Linux, macOS, Windows (with LAPACK linkage)


