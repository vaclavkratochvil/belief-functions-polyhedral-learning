From Confidence Bounds to Belief Functions: Computational Supplement
================

# Overview

This repository contains supplementary materials accompanying the paper

**M. Daniel, R. Jiroušek, V. Kratochvíl**  
*From Confidence Bounds to Belief Functions*, submitted to Kybernetika
Journal, December 2025. <https://www.kybernetika.cz/>

The goal of this repository is to ensure **computational transparency
and reproducibility** of the numerical examples and optimization
procedures presented in the paper.

------------------------------------------------------------------------

# Content

The repository includes:

- **R Markdown (`.Rmd`) file**
  [analysis.Rmd](https://github.com/vaclavkratochvil/belief-functions-polyhedral-learning/blob/main/analysis.Rmd)
  (compiled into user frilendly
  [analysis.md](https://github.com/vaclavkratochvil/belief-functions-polyhedral-learning/blob/main/analysis.md)
  implementing:
  - computation of Jeffreys confidence interval lower bounds,
  - construction of pseudo-belief functions from data,
  - upper approximation of pseudo-belief functions by belief functions,
  - linear programming formulations for selecting representative belief
    functions.
- Compiled **HTML outputs** corresponding to the supplementary material.

------------------------------------------------------------------------

# Implementation

All computations are implemented in **R** using standard libraries.

Main technologies: - **R (version ≥ 4.2)** - **R Markdown** - Base R
statistical functions (`qbeta`) - Linear programming solver
(`lpSolve`) - Tools supporting polyhedral geometry computations (`rcdd`)

All computations are deterministic; no stochastic components are used.

------------------------------------------------------------------------

# How to Run

You can use one of the following two options.

## Local execution using RStudio

1.  Open **RStudio**.
2.  Download and open the main supplementary `analysis.Rmd` file.
3.  Click **Knit** (HTML or PDF), or run the document using
    `rmarkdown::render()`.

``` r
rmarkdown::render("analysis.Rmd")
```

The compiled output will be generated automatically.

Alternatively, individual code chunks can be executed interactively.
This allows users to: - run selected parts of the analysis step by
step, - modify input values or confidence levels, - organize the code
into custom test scenarios.

Each code chunk is self-contained and can be evaluated independently in
RStudio using the *Run Current Chunk* or *Run All Chunks* options. Note,
however, that some chunks rely on helper functions defined earlier in
the document and may therefore require prior evaluation of those chunks.

## Binder

An interactive **RStudio environment** for reproducing the results is
also available via Binder:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/vaclavkratochvil/belief-functions-polyhedral-learning/main?urlpath=rstudio&filepath=analysis.Rmd)

The Binder link opens RStudio with the main manuscript file
(`analysis.Rmd`) preloaded. The results can be reproduced by clicking
*Knit* or by running individual code chunks interactively.

------------------------------------------------------------------------

# Reproducibility

All numerical results reported in the paper can be reproduced by running
the provided R Markdown files without any manual intervention.

------------------------------------------------------------------------

# License and Usage

The code is provided for **research and academic use**. If you use or
adapt this material, please cite the associated paper.
