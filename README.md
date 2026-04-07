# Koopman Operator Methods for Poincare Maps of Chaotic Dynamics

This repository contains MATLAB scripts for the data-driven analysis of nonlinear dynamical systems using Koopman operator techniques, with a particular focus on rigged Dynamic Mode Decomposition (riggedDMD).

The codes are designed to study spectral properties, generalized eigenfunctions, and symbolic dynamics arising in chaotic systems, with applications ranging from low-dimensional maps to high-dimensional PDEs.

---

## Overview

The repository currently features three primary scripts:

- **Logistic Map**
  - Application of EDMD, mpEDMD, and riggedDMD to a canonical one-dimensional chaotic map.
  - Includes computation of kneading sequences and visualization of eigenfunctions.

- **Rössler Poincaré Section**
  - Construction of a Poincaré return map from trajectory data.
  - Data-driven spectral analysis using riggedDMD.
  - Identification of fixed points, preimages, and geometric structure.

- **Kuramoto–Sivashinsky (KS) Poincaré Section**
  - Analysis of a 32-mode Galerkin truncation of the KS equation.
  - Projection of high-dimensional data onto principal components.
  - Application of riggedDMD to uncover low-dimensional structure.
  - Integration with symbolic dynamics and periodic orbit analysis.

---

## Repository Structure

- logistic.m        (Logistic map analysis)
- rossler.m        (Rössler Poincaré section analysis)
- kuramoto_sivashinsky.m             (KS Poincaré section analysis)

- ks_upos/                   (Unstable periodic orbit data for KS)
  - per1sym.mat
  - per1asym.mat
  - per2asym.mat
  - per3asym.mat
  - per4asym.math

- ks_utils/                  (Auxiliary scripts for KS analysis)
  - symbolic dynamics tools
  - transition matrix construction
  - periodic orbit detection
  - visualization utilities

## Dependencies

This repository relies on the **riggedDMD** implementation developed by Colbrook et al.:

https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition

To use the scripts in this repository:

1. Clone or download the riggedDMD repository.
2. Ensure that the `main_routines` folder from that repository is on your MATLAB path:
   ```matlab
   addpath('path_to_riggedDMD/main_routines')
