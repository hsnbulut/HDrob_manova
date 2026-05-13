This repository contains the R code, simulation outputs, summary tables, and graphical materials associated with the manuscript:

**A Robust High-Dimensional MANOVA Test Based on Weighted MRCD Estimation**

## Overview

The manuscript proposes a weighted minimum regularized covariance determinant (MRCD)-based robust Wilks' Lambda test for one-way high-dimensional MANOVA. The proposed method is designed for settings where the number of variables is comparable to, or larger than, the sample size and where the data may contain outlying observations.

The method combines MRCD-based robust location and scatter estimation with a robust distance-based reweighting step and uses permutation calibration to obtain p-values. The performance of the proposed test is evaluated through Monte Carlo simulations in terms of Type-I error control, power, and robustness under structured contamination.

## Repository Structure

```text
HDrob_manova/
├── code/       # R scripts and functions used in the simulation study
├── results/    # Raw and processed simulation results, tables, and figures
└── README.md
