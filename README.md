From Confidence Bounds to Belief Functions: Computational Supplement
================

# Overview

The supplement is a self-contained, step-by-step R analysis accompanying
*From Confidence Bounds to Belief Functions* (M. Daniel, R. Jiroušek and
V. Kratochvíl).

- **analysis.Rmd**: explanation and executable R code in one document.
- **analysis.R**: the same analysis as an ordinary R script, with the
  explanation retained as comments.
- **analysis.md**: rendered version for reading.

The analysis needs no additional scripts or input data files.
Worked-example inputs and recorded benchmark summaries are included in
the document.

# Run the analysis

Install the packages used by the worked examples:

``` r
install.packages(c("lpSolve", "rcdd", "knitr"))
source("analysis.R")
```

Alternatively, open analysis.Rmd in RStudio and Knit it. Rendering also
requires rmarkdown and Pandoc. By default, the worked examples are
calculated and the included benchmark summaries are displayed.

To repeat the time and memory measurements:

``` r
install.packages(c("highs", "Matrix", "callr", "ps"))
options(bf.run_benchmarks = TRUE)
source("analysis.R")
options(bf.run_benchmarks = FALSE)
```

The same option works when knitting the Rmd. All calculation and summary
functions are defined in the analysis itself. For a quick trial, after
loading the script use run_instance(5) or benchmark_frames(5:6,
repetitions = 2).

# What is measured

The supplement covers Jeffreys lower bounds, UAP, the polyhedral
formulation, the three LP criteria, and vertex enumeration for the small
examples. The resource study uses five repetitions per small-example
task, 20 random inputs per frame size 5–10, and one exploratory input
for 11–15 and 20.

Each measurement runs in a fresh R process. Peak process memory includes
native solver allocations. Matrix-size, time and memory guards stop
overly large calculations; a guarded run is not a proof of
infeasibility. Seeds and RNG settings fix the inputs. Times, memory and
the choice between alternative optimal solutions can depend on the
environment and solver version.

# Maintenance

These utilities are optional and are not dependencies of the analysis:

- run_r_experiments.R –check: consistency and guard tests.
- run_r_experiments.R –extract: regenerate analysis.R from the Rmd.
- run_r_experiments.R –render: regenerate analysis.md.
- run_r_experiments.R –run: rerun the benchmarks and save raw results in
  r_experiment_results after the complete run.
- check_self_contained.R: test execution from an empty working
  directory.
- update_article_tables.R: regenerate the red manuscript TeX tables from
  the saved raw measurements using the analysis’s summary functions.

# Usage

If you use or adapt this material, please cite the associated paper.
