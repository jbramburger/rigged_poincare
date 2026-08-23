# Structure, Dynamical Organization, and Transport in Chaotic Dynamics through Koopman Spectral Analysis

This repository contains MATLAB scripts for the data-driven analysis of nonlinear dynamical systems using Koopman operator techniques, with particular emphasis on **rigged Dynamic Mode Decomposition (riggedDMD)**. The examples use spectral measures and smoothed generalized eigenfunctions to uncover coherent regions, symbolic dynamics, unstable periodic-orbit organization, and macroscopic transport in chaotic systems.

The examples range from a one-dimensional chaotic map to Poincaré maps obtained from low- and high-dimensional continuous-time systems. Although the governing equations are known, the Koopman analyses are performed directly from trajectory or Poincaré-section data and do not require closed-form expressions for the corresponding return maps.

---

## Overview

The repository currently features five primary scripts:

- **Logistic Map (`logistic.m`)**
  - Application of EDMD, mpEDMD, and riggedDMD to a canonical one-dimensional chaotic map.
  - Computation of kneading sequences and visualization of eigenfunctions.

- **Rössler Poincaré Section (`rossler.m`)**
  - Construction of a Poincaré return map from trajectory data.
  - Data-driven spectral analysis using riggedDMD.
  - Identification of fixed points, preimages, and geometric structure.

- **Kuramoto–Sivashinsky Poincaré Section (`kuramoto_sivashinsky.m`)**
  - Analysis of a 32-mode Galerkin truncation of the Kuramoto–Sivashinsky equation.
  - Projection of high-dimensional data onto principal components.
  - Application of riggedDMD to uncover low-dimensional structure.
  - Integration with symbolic dynamics and unstable periodic-orbit analysis.

- **Forced Duffing Oscillator: Regime I (`duffing_regime_I.m`)**
  - Analysis of the parameter regime $\beta=1$, $\delta=0.3$, $\omega=1.2$, and $\gamma=0.5$.
  - Identification of a dominant spectral component near $2\pi/3$.
  - Construction of a three-region cyclic partition of the Poincaré section.
  - Comparison of the generalized eigenfunction with unstable periodic orbits of map periods one, two, and three.

- **Forced Duffing Oscillator: Regime II (`duffing_regime_II.m`)**
  - Analysis of the parameter regime $\beta=0.25$, $\delta=0.1$, $\omega=2$, and $\gamma=2.5$.
  - Identification of related spectral maxima near $2\pi/7$, $4\pi/7$, and $6\pi/7$, with the $6\pi/7$ component dominant.
  - Construction of a seven-region transport model with the progression $R_j\mapsto R_{j+3}\pmod 7$.
  - Robustness tests across smoothing parameters and analysis of transport errors near the low-modulus skeleton.

---

## Repository Structure

- **`logistic.m`** — Logistic map analysis.
- **`rossler.m`** — Rössler Poincaré-section analysis.
- **`kuramoto_sivashinsky.m`** — Kuramoto–Sivashinsky Poincaré-section analysis.
- **`duffing_regime_I.m`** — Forced Duffing oscillator, Regime I.
- **`duffing_regime_II.m`** — Forced Duffing oscillator, Regime II.

- **`ks_upos/`** — Unstable periodic-orbit data for the Kuramoto–Sivashinsky example.
  - `per1sym.mat`
  - `per1asym.mat`
  - `per2asym.mat`
  - `per3asym.mat`
  - `per4asym.mat`

- **`ks_utils/`** — Auxiliary scripts and data for the Kuramoto–Sivashinsky analysis.
  - `ks_adjacency_cycles.m` — Compute symbolic cycles predicted by the adjacency matrix.
  - `ks_find_per_3.m` — Identify period-three orbits and their preimages using generalized-eigenfunction minima.
  - `ks_find_cycles.m` — Search the Poincaré-section data for candidate cycles using return errors.
  - `ks_plot_upos.m` — Plot continuous-time unstable periodic orbits and their projections onto the Poincaré section.
  - `ks_psec_0298_x2=0.mat` — Poincaré-section data gathered from a long trajectory.

- **`duffing_utils/`** — Auxiliary scripts and Poincaré data for the forced Duffing oscillator.
  - `duffing_upo_search.m` — Locate low-period unstable periodic orbits in Regime I using close returns and shooting.
  - `duffing_upo_plot.m` — Overlay the periodic-orbit intersections on the Regime I phase partition and log-modulus field.
  - `duffing_psec_I.mat` — Poincaré-section data for Regime I.
  - `duffing_psec_II.mat` — Poincaré-section data for Regime II.

- **`duffing_results/`** — Generated Duffing results and figures. This directory is created automatically when the Duffing scripts are run.

---

## Dependencies

This repository relies on the **riggedDMD** implementation developed by Colbrook et al.:

<https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition>

To use the scripts in this repository:

1. Clone or download the riggedDMD repository.
2. Copy its `main_routines` folder into the root of this repository, or add that folder to your MATLAB path:

   ```matlab
   addpath('path_to_riggedDMD/main_routines')
   ```

3. Open MATLAB in the root directory of this repository.

The **Optimization Toolbox** is optional for `duffing_upo_search.m`. The script uses `fsolve` when it is available and otherwise falls back to `fminsearch`. The Kuramoto–Sivashinsky example may require the **Statistics and Machine Learning Toolbox** for principal component analysis.

---

## Running the Examples

The five primary scripts can be run independently from the repository root:

```matlab
logistic
rossler
kuramoto_sivashinsky
duffing_regime_I
duffing_regime_II
```

The Duffing scripts automatically load their corresponding Poincaré data:

- `duffing_regime_I.m` loads `duffing_utils/duffing_psec_I.mat`.
- `duffing_regime_II.m` loads `duffing_utils/duffing_psec_II.mat`.

To reproduce the Duffing unstable periodic-orbit overlays, run:

```matlab
duffing_regime_I
run('duffing_utils/duffing_upo_search.m')
run('duffing_utils/duffing_upo_plot.m')
```

The first script computes the Regime I generalized eigenfunction, the second locates and saves the unstable periodic orbits, and the third generates the phase and modulus overlays. All generated Duffing files are written to `duffing_results/`.
