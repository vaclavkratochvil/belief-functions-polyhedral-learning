From Confidence Bounds to Belief Functions: Computational Supplement
================
Václav Kratochvíl, Radim Jiroušek, Milan Daniel

# Computational supplement

This repository accompanies the paper *From Confidence Bounds to Belief
Functions*, submitted to **Kybernetika**. Its main file is
[`analysis.Rmd`](analysis.md), a self-contained and executable account
of the worked examples, polyhedral calculations, and scaling experiment
reported in the paper. The explanations and the R code are kept together
so that a reader can follow the construction step by step and reproduce
individual results in RStudio.

# Contents of `analysis.Rmd`

## Statistical construction of the lower bounds

The first part constructs a pseudo-belief function from observations. An
observation assigned to a subset $B\subseteq\Omega$ supports every event
$A\supseteq B$. These counts are aggregated and converted to lower
bounds by Jeffreys intervals. A small example with precise and ambiguous
observations shows the complete calculation.

## Upper Approximation Procedure

The notebook implements the Upper Approximation Procedure (UAP),
including explicit input, successful output, and failure detection. It
reproduces the ternary example and the quaternary calculations
corresponding to Tables 1 and 4–6 of the paper. The examples include the
original bounds, modified bounds, a change of the confidence level, and
Shafer discounting.

## Polyhedral representation and vertices

The dominance conditions are written as linear inequalities for the
basic probability assignment. This gives an H-representation of the
feasible polytope

$$
\mathcal M_g=\{m:\operatorname{Bel}_m(A)\geq g(A),\ m\geq0,\
m(\emptyset)=0,\ \sum_A m(A)=1\}.
$$

The notebook constructs this representation and uses `rcdd`/`cddlib` to
enumerate its vertices for the ternary and quaternary examples. It
reproduces the vertex counts used in the article.

## Linear-programming upper approximations

Three linear criteria select representative belief functions from the
polytope:

- **L1**, minimizing the total upward correction of the belief values;
- **HD**, maximizing Dubois–Prade nonspecificity;
- **CW**, using inverse-cardinality weights.

The corresponding LP solutions are compared with UAP in Tables 7–9.
Table 10 then reports the number of polytope vertices attaining the
optimum of each criterion. These computations use `lpSolve` for the
small, transparent worked examples.

## Scaling experiment

Calculations for $|\Omega|=3,4$ take negligible time, so the last part
illustrates how far the LP formulation can be taken on a standard
computer. All scaling calculations use **HiGHS** and sparse matrices.

For each $|\Omega|=5,\ldots,10$, the experiment generates 20 data sets.
For each $|\Omega|=11,12,13$, it generates five. Every data set contains
200 observations distributed among all singletons and $k$ randomly
selected nonsingleton subsets, where $k$ is sampled separately for every
data set. Every selected subset first receives at least one observation;
the remaining observations are allocated multinomially.

All calculations through $|\Omega|=12$ finished. None of the five
calculations for $|\Omega|=13$ finished within the 60-second limit per
LP; larger frames were therefore not attempted in this study.

For every data set, the code constructs Jeffreys lower bounds, prepares
the sparse dominance matrix, and solves the L1, HD, and CW programs.
Preparation time and the three optimization times are measured
separately with `system.time()`. Their sum is also retained in
`results.csv`. The memory column is a simple R-managed-memory indicator
obtained with `gc(reset = TRUE)`; it does not include memory allocated
internally by HiGHS. The final chunk groups the measurements by
$|\Omega|$ and reports arithmetic mean and maximum.

# Running the supplement

Install the required packages once:

``` r
install.packages(c(
  "lpSolve", "rcdd", "knitr", "rmarkdown", "highs", "Matrix"
))
```

Open `analysis.Rmd` in RStudio and run the chunks in order, or knit the
complete document. The worked examples and Tables 1 and 4–10 are
evaluated directly.

The longer scaling experiment is disabled by default. To run it, change

``` r
run_scaling_experiment <- FALSE
```

to `TRUE` in the `run-scaling-experiment` chunk, run that chunk, and
then run the `scaling-table` chunk. Results are saved after every
attempted data set in `results.csv` in the repository root. A previous
file is overwritten. The summary is saved as `summary.csv`, and the TeX
version of the table as `table_resources_frames.tex`.

The fixed seeds make the generated inputs reproducible. Computation
times, memory measurements, and the particular solution selected on an
optimal face can still depend on the computer and solver version.

# Files

- `analysis.Rmd`: main computational supplement;
- `analysis.md`: rendered GitHub version of the main supplement;
- `README.Rmd`: source of this repository description;
- `README.md`: rendered repository description.

The code is provided for research and academic reproducibility. When
using or adapting it, please cite the accompanying article.
