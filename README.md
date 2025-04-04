# STAT 689 Final Project - Fall 2024

This project implements and compares multiple multiple-testing correction methods using simulated gene expression data. It provides a performance-based evaluation framework and visualization tools to assess FDR control and statistical power.

## Author

Chaoan Li  
Email: [chaoan@tamu.edu](mailto:chaoan@tamu.edu)

## Overview

The script simulates gene expression testing under both the null and alternative hypotheses, applies multiple FDR-controlling methods, and compares their effectiveness using various performance metrics and visualizations. A custom scoring function, `RejectScore`, is used to balance false discoveries and statistical power.

## Features

- Simulates p-values for 10,000 genes
- Incorporates both traditional and flexible FDR-controlling procedures:
  - Benjamini-Hochberg (BH)
  - Barber-Candès (BC)
  - Adaptive Flexible BC (FBC)
  - Adaptive Flexible BH (FBH)
  - Storey's method (ST)
- Evaluates:
  - True/False Positive/Negative counts
  - FDR
  - Statistical power
  - Rejection rate
  - Custom `RejectScore`
- Plots flexible weight functions used in FBC/FBH
- Visualizes threshold behavior of Storey’s method with varying lambda
- Faceted p-value distribution plots across methods

## Dependencies

The script requires the following R packages:

- `nloptr` – for numerical optimization
- `tidyverse` – for data manipulation and visualization

These packages will be installed automatically if not already present.

## How to Run

1. Open the script in R or RStudio.
2. Run the entire file. It will:
   - Generate simulated data
   - Apply each FDR method
   - Output performance metrics
   - Display several ggplot2-based visualizations

## Custom Functions

### `RejectScore(false_positive, true_positive)`
A logistic-type function scoring the trade-off between false and true discoveries.

### `Flexfunction(function_type, x)`
Defines transformations for p-values based on various functional forms:
- `"log"`
- `"sqrt"`
- `"linear"`
- `"square"`
- `"exp"`

### `threshold(alpha, method, ...)`
Determines the optimal rejection threshold for each method using constrained optimization.

### `calculate_performance(results, T)`
Calculates FDR, power, rejection rate, and other metrics based on decisions at threshold `T`.

## Outputs

- Summary tables of performance metrics
- Threshold plots for Storey’s method vs lambda
- Faceted visualizations of p-values colored by discovery category (TP, FP, TN, FN)

## Notes

- The custom `RejectScore` is useful for methods evaluation when balancing Type I and Type II errors.
- Thresholds are computed using `nloptr` to satisfy FDR constraints.
- Visualizations highlight how methods differ in terms of rejected hypotheses and their accuracy.

## License

This code is free to use, modify, and distribute with proper attribution.

---
For questions or feedback, please contact [chaoan@tamu.edu](mailto:chaoan@tamu.edu).
