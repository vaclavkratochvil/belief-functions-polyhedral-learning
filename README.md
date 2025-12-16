Learning Belief Functions from Data via Polyhedral Methods ?
Computational Supplement
================

# Overview

This repository contains supplementary materials accompanying the paper

**M. Daniel, R. Jiroušek, V. Kratochvíl**  
*On Learning Belief Functions from Data Using Confidence Bounds*

The goal of this repository is to ensure **computational transparency
and reproducibility** of the numerical examples and optimization
procedures presented in the paper.

------------------------------------------------------------------------

# Content

The repository includes:

- **R Markdown (`.Rmd`) file** [Kybernetika
  manuscript](https://github.com/vaclavkratochvil/belief-functions-polyhedral-learning/blob/main/Kybernetika.md)
  implementing:
  - computation of Jeffreys confidence interval lower bounds,
  - construction of pseudo-belief functions from data,
  - upper approximation of pseudo-belief functions by belief functions,
  - linear programming formulations for selecting representative belief
    functions.
- Compiled **HTML/PDF outputs** corresponding to the supplementary
  material.

------------------------------------------------------------------------

# Implementation

All computations are implemented in **R** using standard libraries.

Main technologies: - **R (version ≥ 4.2)** - **R Markdown** - Base R
statistical functions (`qbeta`) - Linear programming solvers
(e.g. `lpSolve`)

All computations are deterministic; no stochastic components are used.

------------------------------------------------------------------------

# How to Run

1.  Open the project in **RStudio**.
2.  Open the main supplementary `.Rmd` file.
3.  Click **Knit** (HTML or PDF), or run:

Alternatively, individual code blocks can be copied and executed
directly in R. Users may provide their own input data or organize the
code into custom test cases depending on their experimental needs.

``` r
rmarkdown::render("Kybernetika.Rmd")
```

The compiled output will be generated automatically.

------------------------------------------------------------------------

# Reproducibility

All numerical results reported in the paper can be reproduced by running
the provided R Markdown files without any manual intervention.

------------------------------------------------------------------------

# License and Usage

The code is provided for **research and academic use**. If you use or
adapt this material, please cite the associated paper.
