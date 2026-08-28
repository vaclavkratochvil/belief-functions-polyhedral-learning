From Confidence Bounds to Belief Functions: Experimental Supplement
================
Václav Kratochvíl, Radim Jiroušek, Milan Daniel
August 28, 2026

- [Introduction](#introduction)
- [Jeffreys Lower Bounds](#jeffreys-lower-bounds)
  - [Example: Ambiguous observations](#example-ambiguous-observations)
- [Upper Approximation Procedure
  (UAP)](#upper-approximation-procedure-uap)
  - [Table 1 (ternary, original $g$)](#table-1-ternary-original-g)
  - [Table 1 (modified $\tilde g$)](#table-1-modified-tilde-g)
  - [UAP Reproduction of Tables 4, 5,
    6](#uap-reproduction-of-tables-4-5-6)
    - [Table 4 — Jeffreys lower bounds with
      $\alpha = 0.05$](#table-4--jeffreys-lower-bounds-with-alpha--005)
    - [Table 5 — Jeffreys lower bounds with
      $\alpha = 0.02$](#table-5--jeffreys-lower-bounds-with-alpha--002)
    - [Table 6 — Discounted pseudo-belief function
      ($\delta = 0.02$)](#table-6--discounted-pseudo-belief-function-delta--002)
- [Polyhedral Representation of the Upper
  Approximation](#polyhedral-representation-of-the-upper-approximation)
  - [Constructing the Polytope in R](#constructing-the-polytope-in-r)
  - [Vertex Enumeration of
    $\mathcal{M}_g$](#vertex-enumeration-of-mathcalm_g)
- [Linear Programming
  Approximations](#linear-programming-approximations)
  - [Call the LP solver](#call-the-lp-solver)
  - [LP Reproduction of Tables 7, 8,
    9](#lp-reproduction-of-tables-7-8-9)
    - [Table 7](#table-7)
    - [Table 8](#table-8)
    - [Table 9](#table-9)
    - [Table 10](#table-10)
- [Larger frames and computational
  cost](#larger-frames-and-computational-cost)
  - [1. Generate a random input](#1-generate-a-random-input)
  - [2. Keep only the necessary
    constraints](#2-keep-only-the-necessary-constraints)
  - [3. Solve and check all three
    criteria](#3-solve-and-check-all-three-criteria)
  - [4. Measure time and memory](#4-measure-time-and-memory)
  - [5. Repeat the experiment](#5-repeat-the-experiment)
  - [6. Display the results](#6-display-the-results)
    - [Small examples](#small-examples)
    - [Twenty random inputs per frame
      size](#twenty-random-inputs-per-frame-size)
    - [Exploratory larger frames](#exploratory-larger-frames)

# Introduction

This notebook provides a compact computational companion to the paper  
**“From Confidence Bounds to Belief Functions”, submitted to
Kybernetika.**

Its purpose is to **reproduce and document all worked examples from the
paper**:

- construction of **Jeffreys-based pseudo-belief functions** $g$ from
  count data,
- application of the **Upper Approximation Procedure (UAP)** to these
  functions (Examples 1–3, Tables 1, 4–6),
- and the **polyhedral LP approach** used to select dominating belief
  functions (Examples 4–5, Tables 7–9).

In particular, the notebook:

1.  Shows how Jeffreys’ intervals turn simple frequency tables into a
    normalized, monotone set function $g$ on $\mathcal{P}(\Omega)$.
2.  Applies the Upper Approximation Procedure to $g$ and its
    modifications, exactly reproducing Tables 1, 4, 5, and 6.
3.  Constructs the dominance polytope $\mathcal{M}_g$ and uses several
    linear objectives (L1, HD, CW) to obtain LP-based upper
    approximations, reproducing Tables 7–9.
4.  For the **quaternary Example 5**, compares the LP solutions with the
    outcome of UAP in a single summary table.

The notebook therefore serves as an executable “log file” for all
numerical tables in the paper and as a reference implementation of the
UAP and LP-based upper approximation methods.

------------------------------------------------------------------------

``` r
# Required packages
library(lpSolve)
library(knitr)

# ===========================
# Auxiliary functions
# ===========================

# Generate all subsets
powerset <- function(Omega) {
  n <- length(Omega)
  unlist(lapply(0:n, function(i) combn(Omega, i, simplify = FALSE)), recursive = FALSE)
}

subset_name <- function(S) {
  if (length(S)==0) return("{}")
  paste0("{", paste(S, collapse=","), "}")
}

mobius_transform <- function(subsets, g_values) {
  k <- length(subsets)
  m_g <- numeric(k)
  
  for (i in seq_len(k)) {
    A <- subsets[[i]]
    
    idx <- sapply(subsets, function(B)
      all(B %in% A))
    
    signs <- sapply(subsets[idx], function(B)
      (-1)^(length(A) - length(B)))
    
    m_g[i] <- sum(signs * g_values[idx])
  }
  
  m_g
}
```

# Jeffreys Lower Bounds

To construct the statistical lower approximation $g(A)$ on the Boolean
algebra $\mathcal{P}(\Omega)$, we follow the Jeffreys–based procedure
described in Section 2 of the paper.

Assume that observations are recorded on subsets $B \subseteq \Omega$,
possibly ambiguous (e.g. $\{w_1,w_2\}$). Let $n_B$ denote the number of
observations associated with $B$, and let

$$
N = \sum_{B \subseteq \Omega} n_B
$$

be the total number of observed cases.

An observation assigned to $B$ supports the event $A$ only when
$B \subseteq A$. The effective number of “successes’’ for event $A$ is
therefore

$$
s(A) = \sum_{B \subseteq A} n_B .
$$

Under Jeffreys’ prior $\mathrm{Beta}(1/2,1/2)$, the posterior
distribution of the binomial parameter $p(A)$ is

$$
p(A) \mid \mathrm{data} \sim
\mathrm{Beta}\!\left(s(A)+\tfrac12,\; N-s(A)+\tfrac12 \right),
$$

and the corresponding lower $(1-\alpha)$-credible bound is

$$
g(A) = L(A) =
\mathrm{qbeta}\!\left(\tfrac{\alpha}{2},\, s(A)+\tfrac12,\,
N-s(A)+\tfrac12 \right).
$$

Following the article, we additionally impose the normalization
condition $g(\Omega) = 1 ,$ which ensures that the resulting set
function behaves as a *pseudo-belief function* compatible with the
dominance relation used later in the Upper Approximation Procedure.

All computations of $g(A)$ in this notebook use the function
`jeffreys_pseudo_belief_function()` defined below, which implements the
mapping $\{n_B\} \mapsto g(\cdot)$ exactly as described above.

For convenience, we include the empty set $\emptyset$ in the internal
representation (with $g(\emptyset)=0$, but we omit it from the printed
tables.

``` r
# ------------------------------------------------------------
# Jeffreys confidence interval for the binomial model
# Returns the lower and upper (1 - alpha) credible bounds.
# ------------------------------------------------------------
jeffreys_interval <- function(successes, total, alpha = 0.05,
                              a0 = 0.5, b0 = 0.5) {
  
  lower <- qbeta(alpha / 2, successes + a0, total - successes + b0)
  upper <- qbeta(1 - alpha / 2, successes + a0, total - successes + b0)
  
  list(lower = lower, upper = upper)
}



# ------------------------------------------------------------
# Normalize a named vector of counts:
# Ensures all subsets exist; missing subsets are assigned zero.
# ------------------------------------------------------------
normalize_counts <- function(counts_named, subsets) {
  
  subset_names <- sapply(subsets, subset_name)
  out <- setNames(rep(0, length(subsets)), subset_names)
  
  if (is.null(names(counts_named)))
    stop("counts_named must be a named vector, e.g. c('{w1}'=5, '{w1,w2}'=2).")
  
  # fill matching names
  out[names(counts_named)] <- as.numeric(counts_named)
  
  out
}



# ------------------------------------------------------------
# Construct Jeffreys-based pseudo-belief function g(A):
#   g(A) = Jeffreys lower bound for f(A)
#   where "successes" for A is the total number of observations
#   assigned to subsets B \subseteq A.
#
# The function automatically enforces g(Omega) = 1, consistently with the paper
# ------------------------------------------------------------
jeffreys_pseudo_belief_function <- function(Omega, counts_named, alpha = 0.05) {
  
  # all subsets in cardinality order
  subsets <- powerset(Omega)
  subset_names <- sapply(subsets, subset_name)
  
  # full count vector for all subsets
  counts_full <- normalize_counts(counts_named, subsets)
  N <- sum(counts_full)
  
  # compute "successes" s(A) = sum_{B \subseteq A} n_B
  successes <- sapply(subsets, function(A) {
    idx <- sapply(subsets, function(B) all(B %in% A))
    sum(counts_full[idx])
  })
  
  # allocate vector for g(A)
  g_values <- numeric(length(subsets))
  
  for (i in seq_along(subsets)) {
    A <- subsets[[i]]
    
    # by convention g(emptyset) = 0
    if (length(A) == 0) {
      g_values[i] <- 0
    } else {
      ji <- jeffreys_interval(successes[i], N, alpha = alpha)
      g_values[i] <- ji$lower
    }
  }
  
  # normalize: g(Omega) must equal 1
  g_values[length(g_values)] <- 1
  
  list(
    g            = g_values,
    subsets      = subsets,
    subset_names = subset_names,
    counts       = counts_full,
    N            = N
  )
}
```

## Example: Ambiguous observations

This example illustrates how the Jeffreys-based construction processes
**partially ambiguous measurements**.

We consider $30$ observations in total:

- $25$ observations are fully specified (assigned to the singletons
  $\{w_1\}, \{w_2\}, \{w_3\}$),
- $3$ observations cannot distinguish between states $w_1$ and $w_2$ and
  are therefore recorded as the set $\{w_1,w_2\}$,
- $2$ observations provide no discriminatory information at all,
  represented by the full set $\{w_1,w_2,w_3\}$.

In this setting, an ambiguous observation recorded as $B$ contributes
evidence to every event $A$ such that $B \subseteq A$.  
Consequently, the resulting Jeffreys lower bounds $g(A)$ naturally
incorporate both precise and imprecise information without forcing an
artificial resolution of uncertainty.

``` r
## Frame of discernment
Omega <- c("w1", "w2", "w3")

## Example: 30 observations (25 precise, 3 ambiguous {w1,w2}, 2 fully ambiguous)
counts <- c(
  "{w1}"        = 10,
  "{w2}"        = 8,
  "{w3}"        = 7,
  "{w1,w2}"     = 3,
  "{w1,w2,w3}"  = 2
)

## Compute Jeffreys pseudo-belief function g(A)
g_result <- jeffreys_pseudo_belief_function(Omega, counts, alpha = 0.05)

## Create output table (drop empty set)
df <- data.frame(
  subset = g_result$subset_names,
  g       = round(g_result$g, 3),
  row.names = NULL
)

df <- df[df$subset != "{}", ]   # hide the empty set

df
```

    ##       subset     g
    ## 2       {w1} 0.186
    ## 3       {w2} 0.135
    ## 4       {w3} 0.111
    ## 5    {w1,w2} 0.523
    ## 6    {w1,w3} 0.390
    ## 7    {w2,w3} 0.328
    ## 8 {w1,w2,w3} 1.000

# Upper Approximation Procedure (UAP)

To obtain a belief function $f$ that dominates the statistically derived
pseudo-belief function $g$, the paper employs the **Upper Approximation
Procedure** (UAP) introduced in Section 3.  
The method constructs a mass assignment $m$ and the corresponding belief
function $f$ by processing the subsets of $\Omega$ in increasing order
of cardinality.

For every set $A$, the value of $f(A)$ is determined so that

$$
f(A) \ge g(A)
\quad\text{and}\quad
f(A) \ge \sum_{B \subset A} m(B),
$$

where the second inequality ensures the internal consistency of the
belief function. The increment

$$
m(A) = f(A) - \sum_{B \subset A} m(B)
$$

defines the basic probability assigned to $A$.

If the procedure reaches the top element $\Omega$ and produces
$f(\Omega) > 1$, the construction fails; otherwise the output is a valid
belief function dominating $g$.  
This behaviour matches the examples in the article, including the
failure for Table 4 and the successful runs for Table 1 (modified),
Table 5, and Table 6.

All computations below use the implementation `upper_approximation()`,
which follows the pseudocode of the article exactly.The procedure
processes subsets in increasing cardinality, ensuring the Möbius
recursion is respected.

``` r
# ------------------------------------------------------------
# Upper Approximation Procedure (UAP)
#
# Input:
#   Omega     ... vector of frame elements
#   g_values  ... numeric vector of g(A) in cardinality order,
#                 or named by subset_name().
#
# Output:
#   list(success, subsets, f, m)
#
# The procedure assumes that g(Omega) = 1 (enforced here).
# ------------------------------------------------------------
upper_approximation <- function(Omega, g_values) {
  
  # 1) Generate subsets in canonical order (cardinality-first)
  subsets <- powerset(Omega)
  subsets <- subsets[order(sapply(subsets, length))]
  subset_names <- sapply(subsets, subset_name)
  
  # 2) Reorder / name g_values
  if (!is.null(names(g_values))) {
    g_values <- g_values[subset_names]
  } else {
    names(g_values) <- subset_names
  }
  
  # 3) Enforce normalization
  full_idx <- which(subset_names == subset_name(Omega))
  g_values[full_idx] <- 1
  
  # 4) Compute m_g BEFORE UAP
  m_g <- mobius_transform(subsets, g_values)
  
  k <- length(subsets)
  f <- numeric(k)
  m_f <- numeric(k)
  
  # 5) Main UAP loop
  for (i in seq_len(k)) {
    A <- subsets[[i]]
    
    idx <- sapply(subsets, function(B)
      length(B) < length(A) && all(B %in% A))
    
    sum_prev <- sum(m_f[idx])
    
    f[i]   <- max(g_values[i], sum_prev)
    m_f[i] <- f[i] - sum_prev
  }
  
  # 6) Failure check
  if (f[k] > 1 + 1e-12) {
    return(list(
      success = FALSE,
      message = "UAP failed: f(Omega) > 1",
      subsets = subset_names,
      g = g_values,
      m_g = m_g,
      f = f,
      m_f = m_f
    ))
  }
  
  list(
    success = TRUE,
    subsets = subset_names,
    g = g_values,
    m_g = m_g,
    f = f,
    m_f = m_f
  )
}



# ============================================
# Helper: Print UAP table exactly like in paper
# ============================================

uap_table <- function(result) {
  df <- data.frame(
    A        = result$subsets,
    g        = round(result$g, 6),
    m_g      = round(result$m_g, 6),
    sum_prev = round(result$f - result$m_f, 6),
    f        = round(result$f, 6),
    m_f      = round(result$m_f, 6),
    row.names = NULL
  )
  
  # remove empty set
  df <- df[-1, ]
  
  df
}
```

## Table 1 (ternary, original $g$)

``` r
Omega3 <- c("w1","w2","w3")

g_ternary <- c(
  "{}"=0,
  "{w1}"=0.1, "{w2}"=0.1, "{w3}"=0.2,
  "{w1,w2}"=0.5, "{w1,w3}"=0.4, "{w2,w3}"=0.7,
  "{w1,w2,w3}"=1.0
)

res <- upper_approximation(Omega3, g_ternary)
res$success  # expected: FALSE (as in Table 1)
```

    ## [1] FALSE

``` r
table <- uap_table(res)

# print table
knitr::kable(
  table,
  caption = "Table 1 — Upper Approximation for the original pseudo-belief function g on a ternary frame.",
  align = "lrrrr",
  row.names = FALSE
)
```

| A          |   g |  m_g | sum_prev |   f | m_f |
|:-----------|----:|-----:|---------:|----:|:----|
| {w1}       | 0.1 |  0.1 |      0.0 | 0.1 | 0.1 |
| {w2}       | 0.1 |  0.1 |      0.0 | 0.1 | 0.1 |
| {w3}       | 0.2 |  0.2 |      0.0 | 0.2 | 0.2 |
| {w1,w2}    | 0.5 |  0.3 |      0.2 | 0.5 | 0.3 |
| {w1,w3}    | 0.4 |  0.1 |      0.3 | 0.4 | 0.1 |
| {w2,w3}    | 0.7 |  0.4 |      0.3 | 0.7 | 0.4 |
| {w1,w2,w3} | 1.0 | -0.2 |      1.2 | 1.2 | 0.0 |

Table 1 \<U+2014\> Upper Approximation for the original pseudo-belief
function g on a ternary frame.

For the original pseudo-belief function $g$, the UAP fails because the
lower bounds on several subsets jointly force the final value
$f(\Omega)$ to exceed 1.

## Table 1 (modified $\tilde g$)

``` r
g_ternary_tilde <- c(
  "{}"=0,
  "{w1}"=0.2, "{w2}"=0.2, "{w3}"=0.3,
  "{w1,w2}"=0.5, "{w1,w3}"=0.4, "{w2,w3}"=0.7,
  "{w1,w2,w3}"=1.0
)

res_tilde <- upper_approximation(Omega3, g_ternary_tilde)
res_tilde$success  # expected: TRUE
```

    ## [1] TRUE

``` r
table <- uap_table(res_tilde)

# print table
knitr::kable(
  table,
  caption = "Table 1 (modified) — UAP succeeds for the adjusted pseudo-belief.",
  align = "lrrrr",
  row.names = FALSE
)
```

| A          |   g |  m_g | sum_prev |   f | m_f |
|:-----------|----:|-----:|---------:|----:|:----|
| {w1}       | 0.2 |  0.2 |      0.0 | 0.2 | 0.2 |
| {w2}       | 0.2 |  0.2 |      0.0 | 0.2 | 0.2 |
| {w3}       | 0.3 |  0.3 |      0.0 | 0.3 | 0.3 |
| {w1,w2}    | 0.5 |  0.1 |      0.4 | 0.5 | 0.1 |
| {w1,w3}    | 0.4 | -0.1 |      0.5 | 0.5 | 0.0 |
| {w2,w3}    | 0.7 |  0.2 |      0.5 | 0.7 | 0.2 |
| {w1,w2,w3} | 1.0 |  0.1 |      1.0 | 1.0 | 0.0 |

Table 1 (modified) \<U+2014\> UAP succeeds for the adjusted
pseudo-belief.

The modified function $\tilde{g}$, which only increases singleton
values, removes this conflict: the accumulated mass from proper subsets
stays below 1, allowing the UAP to construct a valid dominating belief
function. The example illustrates that the feasibility of UAP is
strongly dependent on the internal coherence of the statistical lower
bounds.

## UAP Reproduction of Tables 4, 5, 6

In the quaternary case $|\Omega| = 4$, the data originate from a simple
frequency table (Table 2 in the article).  
Applying the Jeffreys construction yields a pseudo-belief function $g$
whose properties depend sensitively on the confidence level $\alpha$.  
This section reproduces the behaviour of the Upper Approximation
Procedure under three different conditions corresponding to Tables 4, 5,
and 6.

``` r
Omega4 <- c("w1","w2","w3","w4")

counts_raw <- c(
  "{w1}"=14, "{w2}"=17, "{w3}"=15, "{w4}"=6
)
```

### Table 4 — Jeffreys lower bounds with $\alpha = 0.05$

For the standard $95\%$ Jeffreys interval, many lower bounds become
relatively tight, especially on two– and three–element subsets.  
When fed into the UAP, the cumulative mass assigned to proper subsets
grows too large, and the procedure attempts to set

$$
f(\Omega) > 1,
$$

which violates the normalization constraint.  
Consequently, **UAP fails** for this choice of $\alpha$, matching the
result reported in Table 4 of the article.

This failure illustrates that the lower statistical bounds may be too
optimistic for a belief function to dominate them.

``` r
# Compute pseudo-belief function g(A)
g4_res  <- jeffreys_pseudo_belief_function(Omega4, counts_raw, alpha = 0.05)
g4_vals <- g4_res$g
names(g4_vals) <- g4_res$subset_names

# Apply UAP
res4 <- upper_approximation(Omega4, g4_vals)

# Expected: FALSE (matches article Table 4)
res4$success
```

    ## [1] FALSE

``` r
# Construct table
tab4 <- uap_table(res4)

# print table
knitr::kable(
  tab4,
  caption = "Table 4 — UAP failure for Jeffreys lower bounds with alpha = 0.05.",
  align = "lrrrr",
  row.names = FALSE
)
```

| A             |        g |       m_g | sum_prev |        f | m_f      |
|:--------------|---------:|----------:|---------:|---------:|:---------|
| {w1}          | 0.163456 |  0.163456 | 0.000000 | 0.163456 | 0.163456 |
| {w2}          | 0.211449 |  0.211449 | 0.000000 | 0.211449 | 0.211449 |
| {w3}          | 0.179203 |  0.179203 | 0.000000 | 0.179203 | 0.179203 |
| {w4}          | 0.049625 |  0.049625 | 0.000000 | 0.049625 | 0.049625 |
| {w1,w2}       | 0.460575 |  0.085670 | 0.374905 | 0.460575 | 0.085670 |
| {w1,w3}       | 0.422637 |  0.079978 | 0.342659 | 0.422637 | 0.079978 |
| {w1,w4}       | 0.261527 |  0.048447 | 0.213081 | 0.261527 | 0.048447 |
| {w2,w3}       | 0.479846 |  0.089194 | 0.390652 | 0.479846 | 0.089194 |
| {w2,w4}       | 0.313479 |  0.052405 | 0.261074 | 0.313479 | 0.052405 |
| {w3,w4}       | 0.278643 |  0.049815 | 0.228828 | 0.278643 | 0.049815 |
| {w1,w2,w3}    | 0.777559 | -0.031391 | 0.808950 | 0.808950 | 0.000000 |
| {w1,w2,w4}    | 0.579468 | -0.031584 | 0.611051 | 0.611051 | 0.000000 |
| {w1,w3,w4}    | 0.538932 | -0.031591 | 0.570524 | 0.570524 | 0.000000 |
| {w2,w3,w4}    | 0.600113 | -0.031578 | 0.631691 | 0.631691 | 0.000000 |
| {w1,w2,w3,w4} | 1.000000 |  0.116902 | 1.009241 | 1.009241 | 0.000000 |

Table 4 \<U+2014\> UAP failure for Jeffreys lower bounds with alpha =
0.05.

### Table 5 — Jeffreys lower bounds with $\alpha = 0.02$

Reducing the confidence level to $\alpha = 0.02$ produces wider Jeffreys
intervals and therefore **smaller lower bounds** $g(A)$.  
These relaxed bounds reduce the dominance pressure on the UAP recursion,
allowing it to satisfy all inequalities while keeping

$$
f(\Omega) = 1.
$$

The UAP therefore **succeeds**, giving the belief function listed in
Table 5. This demonstrates how the statistical confidence level directly
influences the feasibility of constructing a dominating belief function.

``` r
## Table 5 — Jeffreys lower bounds (alpha = 0.02)

# Compute pseudo-belief function g(A)
g5_res  <- jeffreys_pseudo_belief_function(Omega4, counts_raw, alpha = 0.02)
g5_vals <- g5_res$g
names(g5_vals) <- g5_res$subset_names

# Apply UAP
res5 <- upper_approximation(Omega4, g5_vals)

# Expected: TRUE (matches article Table 5)
res5$success
```

    ## [1] TRUE

``` r
# Construct table
tab5 <- uap_table(res5)

# Print table

# Knit-friendly table
knitr::kable(
  tab5,
  caption = "Table 5 — UAP succeeds for Jeffreys lower bounds with alpha = 0.02.",
  align = "lrrrr",
  row.names = FALSE
)
```

| A             |        g |       m_g | sum_prev |        f | m_f      |
|:--------------|---------:|----------:|---------:|---------:|:---------|
| {w1}          | 0.146595 |  0.146595 | 0.000000 | 0.146595 | 0.146595 |
| {w2}          | 0.192394 |  0.192394 | 0.000000 | 0.192394 | 0.192394 |
| {w3}          | 0.161570 |  0.161570 | 0.000000 | 0.161570 | 0.161570 |
| {w4}          | 0.040873 |  0.040873 | 0.000000 | 0.040873 | 0.040873 |
| {w1,w2}       | 0.435508 |  0.096519 | 0.338989 | 0.435508 | 0.096519 |
| {w1,w3}       | 0.398030 |  0.089866 | 0.308165 | 0.398030 | 0.089866 |
| {w1,w4}       | 0.240625 |  0.053157 | 0.187469 | 0.240625 | 0.053157 |
| {w2,w3}       | 0.454600 |  0.100637 | 0.353964 | 0.454600 | 0.100637 |
| {w2,w4}       | 0.291045 |  0.057777 | 0.233267 | 0.291045 | 0.057777 |
| {w3,w4}       | 0.257197 |  0.054754 | 0.202443 | 0.257197 | 0.054754 |
| {w1,w2,w3}    | 0.754524 | -0.033056 | 0.787580 | 0.787580 | 0.000000 |
| {w1,w2,w4}    | 0.553884 | -0.033431 | 0.587315 | 0.587315 | 0.000000 |
| {w1,w3,w4}    | 0.513368 | -0.033446 | 0.546814 | 0.546814 | 0.000000 |
| {w2,w3,w4}    | 0.574583 | -0.033421 | 0.608005 | 0.608005 | 0.000000 |
| {w1,w2,w3,w4} | 1.000000 |  0.139214 | 0.994141 | 1.000000 | 0.005859 |

Table 5 \<U+2014\> UAP succeeds for Jeffreys lower bounds with alpha =
0.02.

### Table 6 — Discounted pseudo-belief function ($\delta = 0.02$)

Another way to restore feasibility is **Shafer’s discounting**.  
Applying a discount factor $\delta = 0.02$ multiplies every proper
subset by $(1-\delta)$ while keeping $g(\Omega)=1$.  
This uniform reduction weakens the dominance requirements across all
levels of the lattice:

$$
\hat{g}(A) = (1-\delta)\, g(A), \qquad A \neq \Omega.
$$

The discounted function $\hat{g}$ becomes compatible with a valid belief
function, and the UAP again **succeeds**, reproducing Table 6 of the
article.

These three examples highlight that the UAP is highly sensitive to the
shape of the pseudo-belief function, and even minimal adjustments—either
statistical ($\alpha$) or structural ($\delta$)—may determine whether a
feasible upper approximation exists.

``` r
## Table 6 — Discounting delta = 0.02 applied to Jeffreys g from Table 4

delta <- 0.02

# Start from g of Table 4 (alpha = 0.05)
g4_vals <- g4_res$g
names(g4_vals) <- g4_res$subset_names

# Apply Shafer discounting to all proper subsets
g6_vals <- g4_vals
is_top <- names(g6_vals) == "{w1,w2,w3,w4}"

g6_vals[!is_top] <- (1 - delta) * g6_vals[!is_top]
g6_vals[is_top]  <- 1   # enforce g(Omega) = 1

# Run UAP
res6 <- upper_approximation(Omega4, g6_vals)

# Expected: TRUE
res6$success
```

    ## [1] TRUE

``` r
# Construct table
tab6 <- uap_table(res6)

# print table
knitr::kable(
  tab6,
  caption = "Table 6 — UAP succeeds after Shafer discounting with delta = 0.02.",
  align = "lrrrr",
  row.names = FALSE
)
```

| A             |        g |       m_g | sum_prev |        f | m_f      |
|:--------------|---------:|----------:|---------:|---------:|:---------|
| {w1}          | 0.160187 |  0.160187 | 0.000000 | 0.160187 | 0.160187 |
| {w2}          | 0.207220 |  0.207220 | 0.000000 | 0.207220 | 0.207220 |
| {w3}          | 0.175619 |  0.175619 | 0.000000 | 0.175619 | 0.175619 |
| {w4}          | 0.048632 |  0.048632 | 0.000000 | 0.048632 | 0.048632 |
| {w1,w2}       | 0.451363 |  0.083957 | 0.367407 | 0.451363 | 0.083957 |
| {w1,w3}       | 0.414184 |  0.078378 | 0.335806 | 0.414184 | 0.078378 |
| {w1,w4}       | 0.256297 |  0.047478 | 0.208819 | 0.256297 | 0.047478 |
| {w2,w3}       | 0.470249 |  0.087410 | 0.382839 | 0.470249 | 0.087410 |
| {w2,w4}       | 0.307209 |  0.051357 | 0.255852 | 0.307209 | 0.051357 |
| {w3,w4}       | 0.273071 |  0.048819 | 0.224251 | 0.273071 | 0.048819 |
| {w1,w2,w3}    | 0.762008 | -0.030763 | 0.792771 | 0.792771 | 0.000000 |
| {w1,w2,w4}    | 0.567878 | -0.030952 | 0.598830 | 0.598830 | 0.000000 |
| {w1,w3,w4}    | 0.528154 | -0.030959 | 0.559113 | 0.559113 | 0.000000 |
| {w2,w3,w4}    | 0.588110 | -0.030947 | 0.619057 | 0.619057 | 0.000000 |
| {w1,w2,w3,w4} | 1.000000 |  0.134564 | 0.989056 | 1.000000 | 0.010944 |

Table 6 \<U+2014\> UAP succeeds after Shafer discounting with delta =
0.02.

# Polyhedral Representation of the Upper Approximation

The Upper Approximation of pseudo-belief function $g$ defined in the
article is the set of all basic probability assignments $m$ whose
induced belief function $f$ dominates the pseudo-belief function $g$.  
For every subset $A \subseteq \Omega$,

$$
f(A) = \sum_{B \subseteq A} m(B)
\qquad\text{and}\qquad
f(A) \ge g(A).
$$

Together with the standard belief-function constraints

$$
m(\emptyset)=0,\qquad m(A)\ge 0,\qquad \sum_{A} m(A)=1,
$$

the feasible region becomes a **polytope** in
$\mathbb{R}^{2^{|\Omega|}}$, described entirely by linear inequalities:

$$
\mathcal{M}_g = \bigl\{ m \,\mid\, 
\mathbf{M}_{\mathrm{bel}}\, m \ge g,\;
m \ge 0,\;
m(\emptyset)=0,\;
\sum_A m(A)=1
\bigr\}.
$$

Here, $\mathbf{M}_{\mathrm{bel}}$ is the *incidence matrix* of the
subset lattice: its entry is 1 exactly when $B \subseteq A$.  
This matrix linearly maps a BPA $m$ to the corresponding belief function
$f$.

The polyhedral formulation serves two purposes:

1.  it provides a global geometric interpretation of the dominance
    condition $f \ge g$;

2.  it allows the use of dedicated tools (such as **cddlib**) to
    enumerate the **vertices** of $\mathcal{M}_g$ and analyse its
    structure.

In particular, the ternary example in the article yields $\;\mathbf{28}$
extreme points,  
while the quaternary example expands to $\;\mathbf{13\,889}$ vertices,  
illustrating the rapid combinatorial growth of the polytope as the
dimension increases.

## Constructing the Polytope in R

The following function constructs this polyhedron directly in terms of
its linear inequality and equality constraints. This representation is
*unified*: it is suitable both for **vertex enumeration** (`rcdd`) and
for **linear programming** (LP-based approximations).

``` r
### Constructing the polytope M_g using a unified LP representation
build_lp_constraints <- function(Omega, g_values) {
  
  subsets       <- powerset(Omega)
  subset_names  <- sapply(subsets, subset_name)
  n_vars        <- length(subsets)
  
  stopifnot(length(g_values) == n_vars)
  
  # -------------------------------
  # Dominance constraints: M_bel m >= g
  # -------------------------------
  M_dom <- matrix(0, nrow=n_vars, ncol=n_vars)
  for (i in seq_along(subsets)) {
    A <- subsets[[i]]
    M_dom[i, sapply(subsets, function(B) all(B %in% A))] <- 1
  }
  b_dom <- g_values
  
  # -------------------------------
  # Nonnegativity: m >= 0
  # -------------------------------
  M_nonneg <- diag(n_vars)
  b_nonneg <- rep(0, n_vars)
  
  # Stack all inequalities
  M_ineq <- rbind(M_dom, M_nonneg)
  b_ineq <- c(b_dom, b_nonneg)
  
  # -------------------------------
  # Equalities:
  #   sum m(A) = 1
  #   m(emptyset)   = 0
  # -------------------------------
  M_eq <- matrix(0, nrow=2, ncol=n_vars)
  M_eq[1, ] <- 1         # normalization
  M_eq[2, 1] <- 1        # empty set mass = 0
  b_eq <- c(1, 0)
  
  colnames(M_ineq) <- subset_names
  colnames(M_eq)   <- subset_names
  
  list(
    subsets = subsets,
    subset_names = subset_names,
    M_ineq = M_ineq,
    b_ineq = b_ineq,
    M_eq   = M_eq,
    b_eq   = b_eq
  )
}

# Example build
Omega <- c("w1","w2","w3")
g_example <- c(0, 0.1,0.2,0.3, 0.5,0.6,0.8, 1.0)

lp_data <- build_lp_constraints(Omega, g_example)
lp_data[3:6]
```

    ## $M_ineq
    ##       {} {w1} {w2} {w3} {w1,w2} {w1,w3} {w2,w3} {w1,w2,w3}
    ##  [1,]  1    0    0    0       0       0       0          0
    ##  [2,]  1    1    0    0       0       0       0          0
    ##  [3,]  1    0    1    0       0       0       0          0
    ##  [4,]  1    0    0    1       0       0       0          0
    ##  [5,]  1    1    1    0       1       0       0          0
    ##  [6,]  1    1    0    1       0       1       0          0
    ##  [7,]  1    0    1    1       0       0       1          0
    ##  [8,]  1    1    1    1       1       1       1          1
    ##  [9,]  1    0    0    0       0       0       0          0
    ## [10,]  0    1    0    0       0       0       0          0
    ## [11,]  0    0    1    0       0       0       0          0
    ## [12,]  0    0    0    1       0       0       0          0
    ## [13,]  0    0    0    0       1       0       0          0
    ## [14,]  0    0    0    0       0       1       0          0
    ## [15,]  0    0    0    0       0       0       1          0
    ## [16,]  0    0    0    0       0       0       0          1
    ## 
    ## $b_ineq
    ##  [1] 0.0 0.1 0.2 0.3 0.5 0.6 0.8 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ## 
    ## $M_eq
    ##      {} {w1} {w2} {w3} {w1,w2} {w1,w3} {w2,w3} {w1,w2,w3}
    ## [1,]  1    1    1    1       1       1       1          1
    ## [2,]  1    0    0    0       0       0       0          0
    ## 
    ## $b_eq
    ## [1] 1 0

## Vertex Enumeration of $\mathcal{M}_g$

The unified LP description of $\mathcal{M}_g$ can be converted into the
$H$-representation required by `rcdd`.  
Each vertex of this polytope corresponds to an *extreme* belief function
that still dominates the pseudo-belief function $g$.

The number of vertices grows rapidly, as illustrated on polytopes
related to above mentioned examples on $|\Omega| = 3,4$. We also list
first 10 vertices for each polytope.

This illustrates the exponential geometric complexity of the dominance
constraints.

``` r
library(rcdd)

# Convert LP constraints to rcdd H-representation
lp_to_H <- function(lp) {
  H_ineq <- cbind(0, -lp$b_ineq, lp$M_ineq)  # inequalities: 0 + Mx - b >= 0
  H_eq   <- cbind(1, -lp$b_eq,   lp$M_eq)    # equalities:  type=1
  d2q(rbind(H_ineq, H_eq))
}

# Extract vertices returned by rcdd
extract_vertices <- function(Vout) {
  out <- Vout$output
  V   <- out[out[,1]=="0", , drop=FALSE]   # rows marked as vertices
  if (nrow(V) == 0) return(NULL)
  q2d(V[,-c(1,2), drop=FALSE])             # drop type and constant
}

# Wrapper for computing vertices directly from g-values
enumerate_vertices <- function(Omega, g_values) {
  lp <- build_lp_constraints(Omega, g_values)
  H  <- lp_to_H(lp)
  V  <- scdd(H)
  verts <- extract_vertices(V)
  colnames(verts) <- names(g_values)
  list(
    subsets = lp$subset_names,
    vertices = verts,
    n_vertices = if (is.null(verts)) 0L else nrow(verts)
  )
}

# --- Ternary example ---
res <- enumerate_vertices(Omega3, g_ternary)
res$n_vertices  # expected: 28
```

    ## [1] 28

``` r
res$vertices[1:10, ]
```

    ##       {} {w1} {w2} {w3}      {w1,w2} {w1,w3} {w2,w3} {w1,w2,w3}
    ##  [1,]  0  0.1  0.4  0.3 5.551115e-17     0.0     0.0        0.2
    ##  [2,]  0  0.1  0.2  0.5 2.000000e-01     0.0     0.0        0.0
    ##  [3,]  0  0.2  0.5  0.2 1.000000e-01     0.0     0.0        0.0
    ##  [4,]  0  0.1  0.4  0.3 2.000000e-01     0.0     0.0        0.0
    ##  [5,]  0  0.2  0.2  0.2 1.000000e-01     0.0     0.3        0.0
    ##  [6,]  0  0.1  0.2  0.3 2.000000e-01     0.0     0.2        0.0
    ##  [7,]  0  0.1  0.3  0.2 1.000000e-01     0.1     0.2        0.0
    ##  [8,]  0  0.1  0.5  0.2 1.000000e-01     0.1     0.0        0.0
    ##  [9,]  0  0.1  0.4  0.2 0.000000e+00     0.1     0.1        0.1
    ## [10,]  0  0.1  0.4  0.2 0.000000e+00     0.2     0.1        0.0

``` r
# --- Quaternary example ---
res4 <- enumerate_vertices(Omega4, g5_vals)
res4$n_vertices  # expected: 13889
```

    ## [1] 13889

``` r
res4$vertices[1:10, ]
```

    ##       {}      {w1}      {w2}      {w3}       {w4}    {w1,w2} {w1,w3} {w1,w4}
    ##  [1,]  0 0.1997516 0.2501711 0.2415141 0.04087347 0.06308764       0       0
    ##  [2,]  0 0.1997516 0.1923940 0.2992912 0.04087347 0.06308764       0       0
    ##  [3,]  0 0.1465951 0.1923940 0.3190163 0.11837567 0.09651909       0       0
    ##  [4,]  0 0.1465951 0.1923940 0.3190163 0.24547552 0.09651909       0       0
    ##  [5,]  0 0.1465951 0.1923940 0.3190163 0.11837567 0.09651909       0       0
    ##  [6,]  0 0.1465951 0.1923940 0.3190163 0.09402999 0.09651909       0       0
    ##  [7,]  0 0.1465951 0.1923940 0.3190163 0.09402999 0.09651909       0       0
    ##  [8,]  0 0.1465951 0.1923940 0.2727424 0.09402999 0.14279295       0       0
    ##  [9,]  0 0.1465951 0.1923940 0.2727424 0.09402999 0.09651909       0       0
    ## [10,]  0 0.1465951 0.1923940 0.2727424 0.09402999 0.09651909       0       0
    ##       {w2,w3}    {w2,w4}   {w3,w4} {w1,w2,w3} {w1,w2,w4} {w1,w3,w4} {w2,w3,w4}
    ##  [1,]       0 0.00000000 0.0000000 0.00000000          0 0.03122834  0.1733737
    ##  [2,]       0 0.05777713 0.0000000 0.00000000          0 0.14682492  0.0000000
    ##  [3,]       0 0.00000000 0.0000000 0.00000000          0 0.00000000  0.1270999
    ##  [4,]       0 0.00000000 0.0000000 0.00000000          0 0.00000000  0.0000000
    ##  [5,]       0 0.00000000 0.1270999 0.00000000          0 0.00000000  0.0000000
    ##  [6,]       0 0.02434568 0.0000000 0.00000000          0 0.00000000  0.1270999
    ##  [7,]       0 0.15144554 0.0000000 0.00000000          0 0.00000000  0.0000000
    ##  [8,]       0 0.15144554 0.0000000 0.00000000          0 0.00000000  0.0000000
    ##  [9,]       0 0.02434568 0.0000000 0.04627386          0 0.00000000  0.1270999
    ## [10,]       0 0.15144554 0.0000000 0.04627386          0 0.00000000  0.0000000
    ##       {w1,w2,w3,w4}
    ##  [1,]             0
    ##  [2,]             0
    ##  [3,]             0
    ##  [4,]             0
    ##  [5,]             0
    ##  [6,]             0
    ##  [7,]             0
    ##  [8,]             0
    ##  [9,]             0
    ## [10,]             0

# Linear Programming Approximations

The polytope $\mathcal{M}_g$ contains all belief functions $f$
dominating the pseudo-belief function $g$.  
To select a *single representative* belief function, we minimise a
linear objective $$
\min_{m \in \mathcal{M}_g} c^\top m,
$$ where the vector $c$ encodes a preference for certain focal sets.

We now implement three criteria used in the article:

- **L1** (minimal upward correction)  
- **HD** (negative Dubois–Prade entropy; LP minimisation = entropy
  maximisation)
- **CW** (cardinality-weighted)

These objectives encode different preferences. L1 minimizes the total
increase of the belief values above the input bounds and is therefore a
data-fidelity criterion. HD maximizes expected log-cardinality and
directly favours nonspecificity. CW also favours larger focal elements,
but uses inverse cardinality rather than logarithmic weights.
Consequently, the three LPs may select different optimal faces, although
ties or coincident solutions are also possible.

The following unified solver builds the LP constraints and computes the
optimal $m$ and induced belief function $f$ for any objective.

``` r
# A list of subsets is used in the examples. For larger frames it is enough
# to pass their cardinalities, without constructing the subsets themselves.
subset_sizes <- function(subsets) {
  if (is.numeric(subsets)) subsets else lengths(subsets)
}

objective_L1 <- function(subsets, Omega) {
  2^(length(Omega) - subset_sizes(subsets))
}

objective_HD <- function(subsets, Omega) {
  -log2(pmax(1, subset_sizes(subsets)))
}

objective_CW <- function(subsets, Omega) {
  sizes <- subset_sizes(subsets)
  ifelse(sizes == 0, 0, 1 / pmax(1, sizes))
}
```

## Call the LP solver

``` r
# ================================================================
# solve_lp_for_objective()
# ---------------------------------------------------------------
# Solves the linear program defining the upper approximation f >= g
# for a chosen linear objective function in the belief space.
#
# The LP is formulated over the basic probability assignment vector m,
# subject to the standard belief-function constraints:
#
#   (1)  M_bel * m   >= g    (dominance: f(A) >= g(A))
#   (2)  m(A)        >= 0    (non-negativity)
#   (3)  m(emptyset)        =  0    (empty set has zero mass)
#   (4)  sum_A m(A)  =  1    (normalization)
#
# The objective is always of the form:
#          minimize   sum_A  w(A) * m(A)
# where w(A) is supplied by obj_fun (L1, HD, CW, …).
#
# Input:
#   Omega     ... frame of discernment
#   g_values  ... pseudo-belief function g(A)
#   obj_fun   ... objective-generating function, e.g. objective_L1
#
# Output:
#   A list with:
#       subset ... subset names in lexicographic order
#       f      ... resulting belief function f(A)
#       m      ... resulting BPA m(A)
#
# ================================================================
solve_lp_for_objective <- function(Omega, g_values, obj_fun) {
  
  # Construct all LP constraint matrices for f >= g
  lp_data <- build_lp_constraints(Omega, g_values)
  
  subsets      <- lp_data$subsets
  subset_names <- lp_data$subset_names
  nA <- length(subsets)
  
  # ------------------------------------------------------------
  # Objective vector: w(A) for each subset A
  # (determined by objective_L1, objective_HD, objective_CW, ...)
  # ------------------------------------------------------------
  obj <- obj_fun(subsets, Omega)
  
  # ------------------------------------------------------------
  # Full LP system:
  #   [M_ineq]  m   >= b_ineq
  #   [M_eq  ]  m   == b_eq
  # ------------------------------------------------------------
  M_full <- rbind(lp_data$M_ineq, lp_data$M_eq)
  dir_full <- c(rep(">=", nrow(lp_data$M_ineq)), rep("=", nrow(lp_data$M_eq)))
  rhs_full <- c(lp_data$b_ineq, lp_data$b_eq)
  
  # ------------------------------------------------------------
  # Solve LP using lpSolve
  # ------------------------------------------------------------
  sol <- lp(
    direction    = "min",
    objective.in = obj,
    const.mat    = M_full,
    const.dir    = dir_full,
    const.rhs    = rhs_full
  )
  
  # Extract BPA m(A)
  m <- sol$solution
  names(m) <- subset_names
  
  # Reconstruct belief f(A) = sum_{B \subseteq A} m(B)
  f <- sapply(subsets, function(A) {
    idx <- sapply(subsets, function(B) all(B %in% A))
    sum(m[idx])
  })
  
  list(
    subset = subset_names,
    f = f,
    m = m
  )
}
```

## LP Reproduction of Tables 7, 8, 9

The following code block reproduces Tables 7, 8, and 9 from the paper
by:

1.  constructing the polytope $\mathcal{M}_g$,
2.  solving the LP for each objective,
3.  computing the corresponding belief and mass assignments,
4.  presenting all results in a unified table.

As in the article, we omit the empty set from printed tables.

``` r
# ================================================================
# compute_table()
# ---------------------------------------------------------------
# Computes all mass and belief values for several approximation
# methods (L1, CW, HD, UAP) and returns a unified table.
#
# Input:
#   Omega     - vector of states ("w1","w2",...)
#   g_values  - pseudo-belief function g(A), already ordered and normalized
#
# Output:
#   A list with:
#     $table   - data frame containing g(A), f(A), and m(A) for all methods
#     $subsets - list of subsets of Omega (used by criteria function)
#
# This function does NOT compute criterion values; that is delegated
# to compute_criteria(), keeping responsibilities cleanly separated.
# ================================================================
compute_table <- function(Omega, g_values) {
  
  lp_data <- build_lp_constraints(Omega, g_values)
  subsets <- lp_data$subsets
  subset_names <- lp_data$subset_names
  idx <- 2:length(subset_names)   # remove empty set
  
  # ---- Solve all LP variants ----
  sol_L1  <- solve_lp_for_objective(Omega, g_values, objective_L1)
  sol_HD  <- solve_lp_for_objective(Omega, g_values, objective_HD)
  sol_CW  <- solve_lp_for_objective(Omega, g_values, objective_CW)
  uap     <- upper_approximation(Omega, g_values)
  
  # ---- Build main table (body) ----
  # Each row corresponds to a subset A in lexicographic order.
  # Columns contain:
  #   g(A)       - pseudo-belief
  #   f_method   - belief produced by each approximation method
  #   m_method   - corresponding mass assignments
  table_body <- data.frame(
    subset = subset_names,
    g      = g_values,
    # m_g    = uap$m_g,
    f_L1   = sol_L1$f,
    f_HD   = sol_HD$f,
    f_CW   = sol_CW$f,
    f_UAP  = uap$f,
    
    m_L1   = sol_L1$m,
    m_HD   = sol_HD$m,
    m_CW   = sol_CW$m,
    m_UAP  = uap$m_f,
    
    row.names = NULL
  )
  
  list(
    table_body = table_body,
    subsets = subsets
  )
}

# ================================================================
# evaluate_lp_criteria()
# ---------------------------------------------------------------
# Computes values of each criterion (L1, CW, HD) for each mass
# vector in the table produced by compute_table_explicit().
#
# Input:
#   table_body - data frame containing m_* columns
#   subsets    - list of subsets (same order as masses)
#   Omega      - frame of discernment
#
# Output:
#   A matrix with rows = criteria (L1, CW, HD)
#                 columns = mass vectors (m_L1, m_CW, m_HD, m_rSP, m_UAP ...)
#
# The formula is always:
#     criterion_value = sum_A  w(A) * m(A)
# where w(A) is given by objective_* functions.
# ================================================================
evaluate_lp_criteria <- function(table_body, subsets, Omega) {
  
  # Find all mass columns automatically (m_L1, m_CW, m_HD, m_rSP, m_UAP, ...)
  m_cols <- grep("^m_", names(table_body), value = TRUE)
  
  # Available criteria (objective functions defined earlier)
  crit_map <- list(
    L1 = objective_L1,
    CW = objective_CW,
    HD = function(subsets, Omega) {-objective_HD(subsets, Omega)}
  )
  
  # Initialize (criteria × mass vectors)
  crit_matrix <- matrix(
    NA,
    nrow = length(crit_map),
    ncol = length(m_cols),
    dimnames = list(names(crit_map), m_cols)
  )
  
  # Calculate weighted sum for each criterion–method pair
  for (crit_name in names(crit_map)) {
    w <- crit_map[[crit_name]](subsets, Omega)  # objective weights
    for (m_col in m_cols) {
      m_vec <- table_body[[m_col]]
      crit_matrix[crit_name, m_col] <- sum(w * m_vec)
    }
  }
  
  crit_matrix
}
```

### Table 7

``` r
table <- compute_table(g_values = g_ternary, Omega = Omega3)
knitr::kable(
  table$table_body[-1,],
  digits = 3,
  caption = "Upper approximations (L1, CW, HD, UAP) for the ternary example \\(g\\)"
)
```

|     | subset     |   g | f_L1 | f_HD | f_CW | f_UAP | m_L1 | m_HD | m_CW | m_UAP |
|:----|:-----------|----:|-----:|-----:|-----:|------:|-----:|-----:|-----:|------:|
| 2   | {w1}       | 0.1 |  0.2 |  0.1 |  0.2 |   0.1 |  0.2 |  0.1 |  0.2 |   0.1 |
| 3   | {w2}       | 0.1 |  0.2 |  0.3 |  0.2 |   0.1 |  0.2 |  0.3 |  0.2 |   0.1 |
| 4   | {w3}       | 0.2 |  0.2 |  0.2 |  0.2 |   0.2 |  0.2 |  0.2 |  0.2 |   0.2 |
| 5   | {w1,w2}    | 0.5 |  0.5 |  0.5 |  0.5 |   0.5 |  0.1 |  0.1 |  0.1 |   0.3 |
| 6   | {w1,w3}    | 0.4 |  0.4 |  0.4 |  0.4 |   0.4 |  0.0 |  0.1 |  0.0 |   0.1 |
| 7   | {w2,w3}    | 0.7 |  0.7 |  0.7 |  0.7 |   0.7 |  0.3 |  0.2 |  0.3 |   0.4 |
| 8   | {w1,w2,w3} | 1.0 |  1.0 |  1.0 |  1.0 |   1.2 |  0.0 |  0.0 |  0.0 |   0.0 |

Upper approximations (L1, CW, HD, UAP) for the ternary example $g$

``` r
evaluate_lp_criteria(table_body = table$table_body, subsets = table$subsets, Omega = Omega3)
```

    ##    m_L1 m_HD m_CW m_UAP
    ## L1  3.2  3.2  3.2   3.2
    ## CW  0.8  0.8  0.8   0.8
    ## HD  0.4  0.4  0.4   0.8

### Table 8

``` r
table <- compute_table(g_values = g_ternary_tilde, Omega = Omega3)
knitr::kable(
  table$table_body[-1,],
  digits = 3,
  caption = "Upper approximations (L1, CW, HD, UAP) for the ternary example \\(\\tilde{g}\\)"
)
```

|     | subset     |   g | f_L1 | f_HD | f_CW | f_UAP | m_L1 | m_HD | m_CW | m_UAP |
|:----|:-----------|----:|-----:|-----:|-----:|------:|-----:|-----:|-----:|------:|
| 2   | {w1}       | 0.2 |  0.2 |  0.2 |  0.2 |   0.2 |  0.2 |  0.2 |  0.2 |   0.2 |
| 3   | {w2}       | 0.2 |  0.2 |  0.2 |  0.2 |   0.2 |  0.2 |  0.2 |  0.2 |   0.2 |
| 4   | {w3}       | 0.3 |  0.3 |  0.3 |  0.3 |   0.3 |  0.3 |  0.3 |  0.3 |   0.3 |
| 5   | {w1,w2}    | 0.5 |  0.5 |  0.5 |  0.5 |   0.5 |  0.1 |  0.1 |  0.1 |   0.1 |
| 6   | {w1,w3}    | 0.4 |  0.5 |  0.5 |  0.5 |   0.5 |  0.0 |  0.0 |  0.0 |   0.0 |
| 7   | {w2,w3}    | 0.7 |  0.7 |  0.7 |  0.7 |   0.7 |  0.2 |  0.2 |  0.2 |   0.2 |
| 8   | {w1,w2,w3} | 1.0 |  1.0 |  1.0 |  1.0 |   1.0 |  0.0 |  0.0 |  0.0 |   0.0 |

Upper approximations (L1, CW, HD, UAP) for the ternary example
$\tilde{g}$

``` r
evaluate_lp_criteria(table_body = table$table_body, subsets = table$subsets, Omega = Omega3)
```

    ##    m_L1 m_HD m_CW m_UAP
    ## L1 3.40 3.40 3.40  3.40
    ## CW 0.85 0.85 0.85  0.85
    ## HD 0.30 0.30 0.30  0.30

### Table 9

``` r
table <- compute_table(g_values = g5_vals, Omega = Omega4)
knitr::kable(
  table$table_body[-1,],
  digits = 3,
  caption = "Upper approximations (L1, CW, HD, UAP) for the quaternary example"
)
```

|     | subset        |     g |  f_L1 |  f_HD |  f_CW | f_UAP |  m_L1 |  m_HD |  m_CW | m_UAP |
|:----|:--------------|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| 2   | {w1}          | 0.147 | 0.158 | 0.158 | 0.147 | 0.147 | 0.158 | 0.158 | 0.147 | 0.147 |
| 3   | {w2}          | 0.192 | 0.203 | 0.203 | 0.192 | 0.192 | 0.203 | 0.203 | 0.192 | 0.192 |
| 4   | {w3}          | 0.162 | 0.173 | 0.173 | 0.162 | 0.162 | 0.173 | 0.173 | 0.162 | 0.162 |
| 5   | {w4}          | 0.041 | 0.052 | 0.052 | 0.041 | 0.041 | 0.052 | 0.052 | 0.041 | 0.041 |
| 6   | {w1,w2}       | 0.436 | 0.436 | 0.436 | 0.436 | 0.436 | 0.074 | 0.074 | 0.097 | 0.097 |
| 7   | {w1,w3}       | 0.398 | 0.398 | 0.398 | 0.398 | 0.398 | 0.068 | 0.068 | 0.090 | 0.090 |
| 8   | {w1,w4}       | 0.241 | 0.241 | 0.241 | 0.241 | 0.241 | 0.031 | 0.031 | 0.053 | 0.053 |
| 9   | {w2,w3}       | 0.455 | 0.455 | 0.455 | 0.455 | 0.455 | 0.079 | 0.079 | 0.101 | 0.101 |
| 10  | {w2,w4}       | 0.291 | 0.291 | 0.291 | 0.291 | 0.291 | 0.035 | 0.035 | 0.058 | 0.058 |
| 11  | {w3,w4}       | 0.257 | 0.257 | 0.257 | 0.257 | 0.257 | 0.032 | 0.032 | 0.055 | 0.055 |
| 12  | {w1,w2,w3}    | 0.755 | 0.755 | 0.755 | 0.788 | 0.788 | 0.000 | 0.000 | 0.000 | 0.000 |
| 13  | {w1,w2,w4}    | 0.554 | 0.554 | 0.554 | 0.587 | 0.587 | 0.000 | 0.000 | 0.000 | 0.000 |
| 14  | {w1,w3,w4}    | 0.513 | 0.513 | 0.513 | 0.547 | 0.547 | 0.000 | 0.000 | 0.000 | 0.000 |
| 15  | {w2,w3,w4}    | 0.575 | 0.575 | 0.575 | 0.608 | 0.608 | 0.000 | 0.000 | 0.000 | 0.000 |
| 16  | {w1,w2,w3,w4} | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 0.095 | 0.095 | 0.006 | 0.006 |

Upper approximations (L1, CW, HD, UAP) for the quaternary example

``` r
evaluate_lp_criteria(table_body = table$table_body, subsets = table$subsets, Omega = Omega4)
```

    ##         m_L1      m_HD      m_CW     m_UAP
    ## L1 6.0592484 6.0592484 6.1481516 6.1481516
    ## CW 0.7692513 0.7692513 0.7692513 0.7692513
    ## HD 0.5088785 0.5088785 0.4644269 0.4644269

### Table 10

Number of vertrices with the same value

``` r
# ================================================================
# compute_vertex_minima()
# ---------------------------------------------------------------
# Computes, for a pseudo-belief function g, how many extreme points
# of the polytope M_g minimize each LP criterion (L1, HD, CW).
#
# The function:
#   1. constructs the polytope M_g via enumerate_vertices(),
#   2. enumerates all extreme points (vertices),
#   3. evaluates each LP objective over all vertices,
#   4. counts how many vertices achieve the minimum value.
#
# Input:
#   g_values ... numeric vector of g(A) in the same subset order
#   Omega    ... frame of discernment, e.g. c("w1","w2","w3")
#
# Output:
#   data.frame with one row and columns:
#       all_vertices ... total number of extreme points
#       L1           ... count of vertices minimizing L1 objective
#       HD           ... count of vertices minimizing HD objective
#       CW           ... count of vertices minimizing CW objective
#
# The result corresponds to Table 10 in the article.
# ================================================================
compute_vertex_minima <- function(g_values, Omega) {
  
  # enumerate all vertices of M_g
  verts <- enumerate_vertices(Omega, g_values)$vertices
  subsets <- powerset(Omega)
  
  # number of vertices
  n_verts <- nrow(verts)
  
  # criterion → objective function
  crit_map <- list(
    L1 = objective_L1,
    HD = objective_HD,
    CW = objective_CW
  )
  
  # initialize result vector
  out <- c(
    all_vertices = n_verts,
    L1 = NA,
    HD = NA,
    CW = NA
  )
  
  # compute for each criterion
  for (crit in names(crit_map)) {
    w <- crit_map[[crit]](subsets, Omega)
    vals <- verts %*% w
    min_val <- min(vals)
    out[crit] <- sum(abs(vals - min_val) < 1e-14)
  }
  
  as.data.frame(t(out))
}


## ============================================================
## Construct Table 10 — counts of vertices minimizing each criterion
## ============================================================

row1 <- compute_vertex_minima(g_ternary,       Omega3)
row2 <- compute_vertex_minima(g_ternary_tilde, Omega3)
row3 <- compute_vertex_minima(g5_vals,         Omega4)

# Add descriptive labels (first column)
Table10 <- rbind(
  "Ternary original g"= row1,
  "Ternary modified g_tilde" = row2,
  "Quaternary example"= row3
)

# Convert row names into a proper column
Table10 <- cbind(
  Example = rownames(Table10),
  Table10,
  row.names = NULL
)

# Print table
knitr::kable(
  Table10,
  caption = "Table 10 — Counts of extreme points and numbers achieving minimum value under each objective.",
  align = "lcrrr",
  digits = 0
)
```

| Example                  | all_vertices |  L1 |  HD |  CW |
|:-------------------------|:------------:|----:|----:|----:|
| Ternary original g       |      28      |   3 |   3 |   3 |
| Ternary modified g_tilde |      27      |   4 |   4 |   4 |
| Quaternary example       |    13889     |   1 |   1 |  14 |

Table 10 \<U+2014\> Counts of extreme points and numbers achieving
minimum value under each objective.

# Larger frames and computational cost

We now repeat the same calculation on randomly generated data. The steps
are unchanged: obtain Jeffreys lower bounds, construct a dominating
belief function, and compare L1, HD and CW. We reuse the functions
already introduced: the interval calculation, the three objective
functions, and the small-example routines for LP/UAP and vertex
enumeration.

For educational clarity, the small-frame examples use a straightforward
dense-matrix formulation with the R package `lpSolve`. For larger
frames, we use HiGHS, an open-source optimization solver accessed
through the R package `highs`, which supports large sparse constraint
matrices. The data provide a lower bound for each subset of the frame,
so the number of potential dominance constraints grows exponentially
with the frame size. Their coefficients are nonzero only when one subset
is contained in another, making sparse storage particularly useful here.

We therefore remove redundant constraints, represent subsets by integer
bit masks, and store only the nonzero matrix entries. This changes the
implementation, not the optimization problem or the three criteria.

## 1. Generate a random input

For a frame of size n, choose all n singletons and n random nonsingleton
focal sets. Assign normalized Gamma(1,1) masses and draw 200
observations. The counts determine the Jeffreys bounds at alpha = 0.05;
as before, the empty-set and full-frame bounds are fixed to 0 and 1.

The following helper computes all subset sums without building a matrix.
The same helper will later reconstruct belief values from a mass vector.
The random seed fixes the generated input.

``` r
bf_cardinality <- function(n) {
  card <- 0L
  for (bit in seq_len(n)) card <- c(card, card + 1L)
  card
}

bf_zeta <- function(x, n) {
  for (bit in 0:(n - 1L)) {
    step <- 2^bit
    blocks <- matrix(x, nrow = 2L * step)
    blocks[(step + 1L):(2L * step), ] <-
      blocks[(step + 1L):(2L * step), , drop = FALSE] +
      blocks[seq_len(step), , drop = FALSE]
    x <- as.vector(blocks)
  }
  x
}

bf_generate <- function(n, seed, N = 200L, alpha = 0.05) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  card <- bf_cardinality(n)
  selected <- c(as.integer(2^(0:(n - 1L))),
                sample(which(card >= 2L) - 1L, n, replace = FALSE))
  probability <- rgamma(2L * n, shape = 1)
  probability <- probability / sum(probability)
  counts <- integer(2^n)
  counts[selected + 1L] <- as.vector(rmultinom(1L, N, probability))
  successes <- bf_zeta(counts, n)
  lookup <- jeffreys_interval(0:N, N, alpha)$lower
  g <- lookup[successes + 1L]
  g[c(1L, length(g))] <- c(0, 1)
  list(card = card, counts = counts, successes = successes, g = g,
       selected = selected, probability = probability)
}
```

## 2. Keep only the necessary constraints

If a nonempty B is contained in A and both have the same lower bound,
the constraint for A follows from that for B and monotonicity. It is
enough to check immediate predecessors. The empty set is excluded from
this test: a nonempty zero-count singleton has a positive Jeffreys lower
bound.

The objective weights below come from the same L1, HD and CW functions
used in the small examples. Only subset cardinalities are needed.

``` r
bf_keep_rows <- function(successes, n) {
  masks <- 0:(length(successes) - 1L)
  keep <- rep(TRUE, length(successes))
  keep[c(1L, length(keep))] <- FALSE
  for (bit in 0:(n - 1L)) {
    b <- as.integer(2^bit)
    a <- masks[bitwAnd(masks, b) != 0L]
    predecessor <- bitwXor(a, b)
    redundant <- predecessor != 0L &
      successes[a + 1L] == successes[predecessor + 1L]
    keep[a[redundant] + 1L] <- FALSE
  }
  which(keep) - 1L
}

bf_sparse_matrix <- function(rows, card, nvars) {
  loadNamespace("Matrix")
  sizes <- 2^card[rows + 1L]
  offsets <- c(0, cumsum(sizes))
  indices <- integer(tail(offsets, 1L))
  for (i in seq_along(rows)) {
    bits <- as.integer(2^(which(as.logical(intToBits(rows[i]))) - 1L))
    submasks <- 0L
    for (b in bits) submasks <- c(submasks, submasks + b)
    indices[seq.int(offsets[i] + 1L, offsets[i + 1L])] <- submasks
  }
  row_matrix <- methods::new(
    "dgRMatrix", Dim = as.integer(c(length(rows), nvars)),
    p = as.integer(offsets), j = indices, x = rep(1, length(indices)))
  methods::as(row_matrix, "CsparseMatrix")
}

bf_objectives <- function(card, n) {
  criteria <- list(L1 = objective_L1, HD = objective_HD, CW = objective_CW)
  lapply(criteria, function(criterion) criterion(card, seq_len(n)))
}

bf_violation <- function(m, g, n) {
  if (any(!is.finite(m))) return(Inf)
  max(0, -min(m), abs(m[1L]), abs(sum(m) - 1),
      max(g - bf_zeta(m, n)))
}
```

## 3. Solve and check all three criteria

The next function builds the sparse LP once and calls HiGHS for each
objective. Each solution is checked against **all** lower bounds,
including the rows removed above, as well as nonnegativity,
normalization and zero empty-set mass. We retain times, focal counts and
the largest constraint violation. A temporary checkpoint allows partial
results to survive if a large calculation has to be stopped.

``` r
bf_scaling_worker <- function(job, checkpoint) {
  started <- as.numeric(Sys.time())
  n <- job$n
  input <- bf_generate(n, job$seed)
  rows <- bf_keep_rows(input$successes, n)
  nnz <- sum(2^input$card[rows + 1L])
  result <- list(
    Omega_size = n, repetition = job$repetition, seed = job$seed,
    variables = 2^n, dominance_constraints = length(rows),
    dominance_nonzeros = nnz,
    estimated_matrix_MiB = (12 * nnz + 4 * (2^n + 1)) / 2^20,
    attempted = FALSE, successful = FALSE, outcome = "preparing",
    seconds_preparation = NA_real_, seconds_three_LPs = NA_real_,
    seconds_total = NA_real_, completed_criteria = 0L,
    seconds_L1 = NA_real_, seconds_HD = NA_real_, seconds_CW = NA_real_,
    focal_L1 = NA_integer_, focal_HD = NA_integer_, focal_CW = NA_integer_,
    objective_L1 = NA_real_, objective_HD = NA_real_, objective_CW = NA_real_,
    n_distinct_solutions = NA_integer_, max_violation = NA_real_, message = "")
  saveRDS(result, checkpoint)
  if (nnz > job$max_nonzeros) {
    result$outcome <- "matrix_limit"
    result$message <- "Dominance nnz exceeds pre-allocation safety limit"
    result$seconds_preparation <- as.numeric(Sys.time()) - started
    result$seconds_total <- result$seconds_preparation
    return(result)
  }
  result$attempted <- TRUE
  saveRDS(result, checkpoint)
  tryCatch({
    dominance <- bf_sparse_matrix(rows, input$card, 2^n)
    A <- rbind(dominance, Matrix::Matrix(1, nrow = 1L, ncol = 2^n, sparse = TRUE))
    upper <- rep(Inf, 2^n)
    upper[1L] <- 0
    lhs <- c(input$g[rows + 1L], 1)
    rhs <- c(rep(Inf, length(rows)), 1)
    objectives <- bf_objectives(input$card, n)
    result$seconds_preparation <- as.numeric(Sys.time()) - started
    result$outcome <- "solving"
    saveRDS(result, checkpoint)
    solutions <- list()
    for (criterion in names(objectives)) {
      tic <- as.numeric(Sys.time())
      fit <- highs::highs_solve(
        L = objectives[[criterion]], lower = 0, upper = upper,
        A = A, lhs = lhs, rhs = rhs,
        control = highs::highs_control(
          threads = 1L, time_limit = job$lp_time_limit,
          primal_feasibility_tolerance = 1e-9,
          dual_feasibility_tolerance = 1e-9,
          ipm_optimality_tolerance = 1e-9))
      result[[paste0("seconds_", criterion)]] <- as.numeric(Sys.time()) - tic
      if (!identical(fit$status_message, "Optimal"))
        stop(paste(criterion, fit$status_message))
      m <- fit$primal_solution
      solutions[[criterion]] <- m
      result[[paste0("focal_", criterion)]] <- sum(m > 1e-8)
      result[[paste0("objective_", criterion)]] <- sum(objectives[[criterion]] * m)
      violation <- bf_violation(m, input$g, n)
      result$max_violation <- max(c(result$max_violation, violation), na.rm = TRUE)
      result$completed_criteria <- length(solutions)
      saveRDS(result, checkpoint)
      if (violation > 1e-7) stop("Independent feasibility validation failed")
    }
    result$seconds_three_LPs <- sum(unlist(result[paste0("seconds_", names(objectives))]))
    result$n_distinct_solutions <- nrow(unique(round(t(do.call(cbind, solutions)), 8)))
    result$successful <- TRUE
    result$outcome <- "solved"
  }, error = function(e) {
    result$outcome <<- if (grepl("bad_alloc|cannot allocate|memory",
                                conditionMessage(e), ignore.case = TRUE))
      "allocation_failure" else "solver_failure"
    result$message <<- conditionMessage(e)
  })
  result$seconds_total <- as.numeric(Sys.time()) - started
  result
}
```

## 4. Measure time and memory

A measurement uses a fresh R process, so memory from earlier
calculations is not carried over. Solver time covers the three LP calls;
total computation time also includes preparation and validation. Process
startup and package loading are excluded from these computation times.

Peak RSS measures resident process memory, including R, loaded packages
and native solver allocations. It is not the size of the R objects
alone. On Windows we also record peak committed memory and use the
operating system’s peak counters. Elsewhere RSS is sampled and is a
lower bound on the true peak. All memory values are in MiB (2^20 bytes).

For safety, a child is stopped above 1,024 MiB process memory or below
128 MiB available system memory. Polling every 20 ms can allow a short
allocation burst to overshoot the threshold. The matrix is not allocated
above 30 million nonzeros. These are resource guards, not tests of
mathematical feasibility.

<details>
<summary>
Implementation of the measurement helper
</summary>

This helper is independent of the belief-function calculation. It loads
only the functions needed by a child process and returns the result with
memory measurements. It can be read separately from the experiment
itself.

``` r
bf_memory <- function(handle = ps::ps_handle()) {
  info <- ps::ps_memory_info(handle)
  peak <- if ("peak_wset" %in% names(info)) info[["peak_wset"]] else info[["rss"]]
  committed <- if ("peak_pagefile" %in% names(info)) info[["peak_pagefile"]] else NA_real_
  current <- if ("mem_private" %in% names(info)) info[["mem_private"]] else info[["rss"]]
  c(rss = info[["rss"]], peak_rss = peak, peak_commit = committed, guard = current) / 2^20
}

bf_child <- function(job, functions, checkpoint) {
  e <- list2env(functions, parent = globalenv())
  for (name in names(functions))
    if (is.function(e[[name]])) environment(e[[name]]) <- e
  if (job$kind == "small") {
    library(lpSolve)
    library(rcdd)
  } else {
    loadNamespace("highs")
    loadNamespace("Matrix")
  }
  gc()
  baseline <- e$bf_memory()
  if (job$kind == "small") {
    tic <- as.numeric(Sys.time())
    value <- if (job$task == "LP_UAP")
      e$compute_table(job$Omega, job$g) else
      e$compute_vertex_minima(job$g, job$Omega)
    result <- list(example = job$example, task = job$task,
                   repetition = job$repetition, successful = TRUE,
                   outcome = "solved",
                   seconds_total = as.numeric(Sys.time()) - tic)
    if (job$task == "vertices")
      stopifnot(value$all_vertices == job$expected_vertices)
  } else result <- e$bf_scaling_worker(job, checkpoint)
  result$baseline_RSS_MiB <- unname(baseline["rss"])
  final <- e$bf_memory()
  result$peak_RSS_MiB <- unname(final["peak_rss"])
  result$peak_commit_MiB <- unname(final["peak_commit"])
  result
}

bf_measure <- function(job, functions, memory_limit_MiB = 1024,
                       min_available_MiB = 128, poll_seconds = 0.02) {
  # Do not serialize the caller's whole notebook environment into the child.
  # bf_child reconstructs a dedicated environment containing only this bundle.
  functions <- lapply(functions, function(f) {
    if (is.function(f)) environment(f) <- baseenv()
    f
  })
  child_entry <- bf_child
  environment(child_entry) <- baseenv()
  checkpoint <- tempfile(fileext = ".rds")
  logs <- c(tempfile(), tempfile())
  on.exit(unlink(c(checkpoint, logs)), add = TRUE)
  child <- callr::r_bg(
    child_entry, args = list(job, functions, checkpoint), libpath = .libPaths(),
    user_profile = FALSE, system_profile = FALSE, supervise = TRUE,
    stdout = logs[1], stderr = logs[2])
  on.exit(if (child$is_alive()) child$kill(), add = TRUE)
  handle <- ps::ps_handle(child$get_pid())
  peak_rss <- 0
  peak_commit <- NA_real_
  reason <- ""
  tic <- as.numeric(Sys.time())
  while (child$is_alive()) {
    mem <- tryCatch(bf_memory(handle), error = function(e) NULL)
    if (!is.null(mem)) {
      peak_rss <- max(peak_rss, mem["peak_rss"])
      if (is.finite(mem["peak_commit"]))
        peak_commit <- max(c(peak_commit, mem["peak_commit"]), na.rm = TRUE)
      available <- ps::ps_system_memory()[["avail"]] / 2^20
      if (mem["guard"] > memory_limit_MiB || mem["rss"] > memory_limit_MiB)
        reason <- "process_memory_limit"
      else if (available < min_available_MiB)
        reason <- "system_memory_limit"
    }
    if (as.numeric(Sys.time()) - tic > job$wall_limit)
      reason <- "wall_time_limit"
    if (nzchar(reason)) {
      child$kill()
      break
    }
    Sys.sleep(poll_seconds)
  }
  result <- tryCatch(child$get_result(), error = function(e) {
    partial <- if (file.exists(checkpoint))
      tryCatch(readRDS(checkpoint), error = function(e) list()) else list()
    partial$successful <- FALSE
    partial$outcome <- if (nzchar(reason)) reason else "process_failure"
    partial$message <- if (nzchar(reason)) reason else conditionMessage(e)
    partial
  })
  result$peak_RSS_MiB <- max(c(result$peak_RSS_MiB, peak_rss), na.rm = TRUE)
  result$peak_commit_MiB <- if (all(is.na(c(result$peak_commit_MiB, peak_commit))))
    NA_real_ else max(c(result$peak_commit_MiB, peak_commit), na.rm = TRUE)
  if (is.null(result$baseline_RSS_MiB)) result$baseline_RSS_MiB <- NA_real_
  result$wall_seconds <- as.numeric(Sys.time()) - tic
  result$memory_limit_MiB <- memory_limit_MiB
  result$min_available_MiB <- min_available_MiB
  result$peak_method <- if (.Platform$OS.type == "windows")
    "OS lifetime peak working set" else "sampled RSS lower bound"
  result
}

combine_results <- function(records) {
  columns <- unique(unlist(lapply(records, names)))
  result <- do.call(rbind, lapply(records, function(x) {
    x[setdiff(columns, names(x))] <- NA
    as.data.frame(x[columns], stringsAsFactors = FALSE)
  }))
  result
}
```

</details>

With that helper in place, a complete experiment has a short interface.
The function bundle contains only the definitions already given above.

``` r
experiment_functions <- function() {
  mget(c("powerset", "subset_name", "subset_sizes", "mobius_transform",
    "upper_approximation", "build_lp_constraints", "lp_to_H",
    "extract_vertices", "enumerate_vertices", "objective_L1", "objective_HD",
    "objective_CW", "solve_lp_for_objective", "compute_table",
    "compute_vertex_minima", "jeffreys_interval", "bf_cardinality", "bf_zeta",
    "bf_generate", "bf_keep_rows", "bf_sparse_matrix", "bf_objectives",
    "bf_violation", "bf_scaling_worker", "bf_memory"),
    envir = environment(experiment_functions))
}

run_instance <- function(n, repetition = 1L,
                         seed = 20260826L + 1000L * n + repetition,
                         time_limit = 60) {
  job <- list(kind = "scaling", n = n, repetition = repetition, seed = seed,
              max_nonzeros = 3e7, lp_time_limit = time_limit,
              wall_limit = 3 * time_limit + 120)
  result <- bf_measure(job, experiment_functions())
  # Preserve identifiers even if a process stops before its first checkpoint.
  result$Omega_size <- n
  result$repetition <- repetition
  result$seed <- seed
  result
}
```

For example, the following runs one instance on a five-element frame.
The full set of returned fields is available in result.

``` r
result <- run_instance(n = 5)
result[c("successful", "seconds_three_LPs", "peak_RSS_MiB",
         "focal_L1", "focal_HD", "focal_CW")]
```

## 5. Repeat the experiment

For n = 5,…,10, we use 20 independently seeded inputs per size. For n =
11,…,15 and 20, we use one input per size and a 180-second limit per LP.
These larger cases are exploratory probes, not a universal dimension
limit. Success or failure also depends on the generated counts.

The small examples reuse compute_table() and compute_vertex_minima(),
with five repetitions per task. The same process-level measurement is
used for all experiments.

``` r
benchmark_frames <- function(sizes, repetitions = 20L, exploratory = FALSE) {
  jobs <- expand.grid(repetition = seq_len(repetitions), n = sizes)
  records <- lapply(seq_len(nrow(jobs)), function(i) {
    n <- jobs$n[i]
    r <- jobs$repetition[i]
    seed <- 20260826L + if (exploratory) n else 1000L * n + r
    run_instance(n, r, seed, time_limit = if (exploratory) 180 else 60)
  })
  combine_results(records)
}

benchmark_small <- function(repetitions = 5L) {
  inputs <- list(
    ternary = list(Omega = Omega3, g = g_ternary, nv = row1$all_vertices),
    modified = list(Omega = Omega3, g = g_ternary_tilde, nv = row2$all_vertices),
    quaternary = list(Omega = Omega4, g = g5_vals, nv = row3$all_vertices))
  jobs <- expand.grid(repetition = seq_len(repetitions),
                      task = c("LP_UAP", "vertices"), example = names(inputs),
                      stringsAsFactors = FALSE)
  records <- lapply(seq_len(nrow(jobs)), function(i) {
    job <- as.list(jobs[i, ])
    input <- inputs[[job$example]]
    job <- c(job, list(kind = "small", Omega = input$Omega, g = input$g,
                      expected_vertices = input$nv, wall_limit = 120))
    result <- bf_measure(job, experiment_functions())
    result$example <- job$example
    result$task <- job$task
    result$repetition <- job$repetition
    result
  })
  combine_results(records)
}
```

The following summaries use median \[maximum\] for times and per-process
memory peaks. Failed or guarded runs are not included in
successful-solve time and focal-count summaries. Their observed memory
remains reported. A missing three-LP time means that the full
calculation did not finish; the memory observed before a stop is not the
memory required to finish.

``` r
bf_pair <- function(x, digits = 3L) {
  if (!length(x) || all(is.na(x))) return("--")
  sprintf(paste0("%.", digits, "f [%.", digits, "f]"),
          median(x, na.rm = TRUE), max(x, na.rm = TRUE))
}
summarise_small <- function(small_r) {
  small_table_r <- do.call(rbind, lapply(split(small_r,
      interaction(small_r$example, small_r$task, drop = TRUE)), function(d) {
    data.frame(Example = d$example[1], Task = d$task[1],
      Success = paste0(sum(d$successful), "/", nrow(d)),
      "Time (s)" = bf_pair(d$seconds_total[d$successful]),
      "Peak RSS (MiB)" = bf_pair(d$peak_RSS_MiB, 1),
      "Peak commit (MiB)" = bf_pair(d$peak_commit_MiB, 1),
      check.names = FALSE)
  }))
  small_table_r <- small_table_r[order(
    match(small_table_r$Example, c("ternary", "modified", "quaternary")),
    match(small_table_r$Task, c("LP_UAP", "vertices"))), ]
  small_table_r$Example <- c(ternary = "Original ternary",
    modified = "Modified ternary", quaternary = "Quaternary")[small_table_r$Example]
  small_table_r
}

summarise_scaling <- function(scaling_r) {
  scaling_table_r <- do.call(rbind, lapply(split(scaling_r,
      factor(scaling_r$Omega_size, levels = sort(unique(scaling_r$Omega_size)))), function(d) {
    ok <- d$successful
    data.frame(n = d$Omega_size[1], Variables = d$variables[1],
      "Mean rows" = round(mean(d$dominance_constraints), 1),
      Success = paste0(sum(ok), "/", nrow(d)),
      "LP time (s)" = bf_pair(d$seconds_three_LPs[ok]),
      "Peak RSS (MiB)" = bf_pair(d$peak_RSS_MiB, 1),
      "Focals L1/HD/CW" = if (any(ok)) sprintf("%.1f/%.1f/%.1f",
        mean(d$focal_L1[ok]), mean(d$focal_HD[ok]), mean(d$focal_CW[ok])) else "--",
      check.names = FALSE)
  }))
  scaling_table_r
}

summarise_boundary <- function(boundary_r) {
  data.frame(
    n = boundary_r$Omega_size, Variables = boundary_r$variables,
    Constraints = boundary_r$dominance_constraints,
    Nonzeros = boundary_r$dominance_nonzeros,
    "Matrix estimate (MiB)" = round(boundary_r$estimated_matrix_MiB, 1),
    "LP time (s)" = ifelse(is.na(boundary_r$seconds_three_LPs), "--",
                           sprintf("%.1f", boundary_r$seconds_three_LPs)),
    "Observed wall (s)" = round(boundary_r$wall_seconds, 1),
    "Peak RSS (MiB)" = round(boundary_r$peak_RSS_MiB, 1),
    "Peak commit (MiB)" = round(boundary_r$peak_commit_MiB, 1),
    Outcome = boundary_r$outcome, check.names = FALSE)
}
```

## 6. Display the results

The measured tables below are included directly in this document. They
use R 4.3.2, lpSolve 5.6.23, rcdd 1.6 and highs 1.9.0-1 on Windows.
Times and memory depend on the computer and its load. Fixed seeds
reproduce the inputs; a different solver version can choose another
point of the same optimal face.

Set run_benchmarks to TRUE in the execution block below to replace these
tables with fresh measurements. The additional packages for this step
are highs, Matrix, callr and ps. The small worked examples need only the
packages loaded earlier.

<details>
<summary>
Recorded measurements used for the displayed tables
</summary>

``` r
recorded <- list(
  small = read.table(header = TRUE, sep = "|", check.names = FALSE, text = "
Example|Task|Success|Time (s)|Peak RSS (MiB)|Peak commit (MiB)
Original ternary|LP_UAP|5/5|0.147 [0.240]|75.3 [75.4]|97.4 [97.5]
Original ternary|vertices|5/5|0.012 [0.033]|75.3 [75.3]|97.5 [97.5]
Modified ternary|LP_UAP|5/5|0.141 [0.170]|75.3 [75.3]|97.4 [97.5]
Modified ternary|vertices|5/5|0.009 [0.012]|75.3 [75.4]|97.4 [97.5]
Quaternary|LP_UAP|5/5|0.105 [0.121]|75.3 [75.4]|97.5 [97.7]
Quaternary|vertices|5/5|4.951 [6.090]|107.5 [107.6]|131.0 [131.1]"),
  scaling = read.table(header = TRUE, sep = "|", check.names = FALSE, text = "
n|Variables|Mean rows|Success|LP time (s)|Peak RSS (MiB)|Focals L1/HD/CW
5|32|27.4|20/20|0.021 [0.042]|175.0 [176.7]|21.8/21.5/20.3
6|64|57.5|20/20|0.014 [0.016]|174.9 [176.7]|33.5/32.0/30.2
7|128|105.7|20/20|0.023 [0.037]|174.8 [175.0]|48.4/44.9/40.2
8|256|226.8|20/20|0.065 [0.080]|176.4 [178.0]|71.7/61.8/53.9
9|512|356.2|20/20|0.115 [0.280]|179.6 [185.9]|89.3/77.0/64.2
10|1024|662.5|20/20|0.586 [1.277]|193.0 [204.2]|110.0/90.6/76.5"),
  boundary = read.table(header = TRUE, sep = "|", check.names = FALSE, text = "
n|Variables|Constraints|Nonzeros|Matrix estimate (MiB)|LP time (s)|Observed wall (s)|Peak RSS (MiB)|Peak commit (MiB)|Outcome
11|2048|2046|175098|2.0|6.0|7.5|249.8|377.7|solved
12|4096|2214|283516|3.3|11.3|12.9|287.9|463.9|solved
13|8192|2790|579286|6.7|25.1|27.5|316.3|522.4|solved
14|16384|16382|4766584|54.6|--|4.4|674.3|1576.1|process_memory_limit
15|32768|9883|4351334|49.9|--|5.9|632.2|1248.9|process_memory_limit
20|1048576|275642|709297932|8121.3|--|2.7|262.2|286.9|matrix_limit"))
```

</details>

Before a full rerun, we check that the sparse and dense LPs give the
same objective values on small inputs. This also checks the subset-sum
transform and the special treatment of the empty set. Different mass
vectors can be equally optimal.

<details>
<summary>
Small-case consistency checks
</summary>

``` r
bf_validate <- function() {
  for (n in 3:5) {
    input <- bf_generate(n, 700L + n)
    masks <- 0:(2^n - 1L)
    dense <- outer(masks, masks, function(a, b) bitwAnd(a, b) == b) * 1
    stopifnot(max(abs(bf_zeta(input$counts, n) -
                        as.vector(dense %*% input$counts))) == 0)
    rows <- bf_keep_rows(input$successes, n)
    sparse <- bf_sparse_matrix(rows, input$card, 2^n)
    stopifnot(identical(unname(as.matrix(sparse)),
                        unname(dense[rows + 1L, , drop = FALSE])))
    A <- rbind(sparse, Matrix::Matrix(1, 1, 2^n, sparse = TRUE))
    for (obj in bf_objectives(input$card, n)) {
      fit <- highs::highs_solve(
        L = obj, lower = 0, upper = c(0, rep(Inf, 2^n - 1L)),
        A = A, lhs = c(input$g[rows + 1L], 1),
        rhs = c(rep(Inf, length(rows)), 1),
        control = highs::highs_control(threads = 1L,
          primal_feasibility_tolerance = 1e-9,
          dual_feasibility_tolerance = 1e-9))
      reference <- lpSolve::lp("min", obj,
        rbind(dense, rep(1, 2^n), c(1, rep(0, 2^n - 1L))),
        c(rep(">=", 2^n), "=", "="), c(input$g, 1, 0))
      stopifnot(fit$status_message == "Optimal", reference$status == 0L,
                bf_violation(fit$primal_solution, input$g, n) < 1e-7,
                abs(fit$objective_value - reference$objval) < 1e-6)
    }
  }
  stopifnot(identical(bf_keep_rows(integer(8L), 3L), c(1L, 2L, 4L)))
  TRUE
}
```

</details>

``` r
# Change FALSE to TRUE to run all benchmarks.
# The option also allows a full run without editing this line.
run_benchmarks <- getOption("bf.run_benchmarks", FALSE)

if (run_benchmarks) {
  # Install once if needed: install.packages(c("highs", "Matrix", "callr", "ps"))
  library(highs)
  library(Matrix)
  library(callr)
  library(ps)
  stopifnot(bf_validate())
  small_r <- benchmark_small(repetitions = 5)
  scaling_r <- benchmark_frames(5:10, repetitions = 20)
  boundary_r <- benchmark_frames(c(11:15, 20), repetitions = 1,
                                  exploratory = TRUE)
  tables <- list(small = summarise_small(small_r),
                 scaling = summarise_scaling(scaling_r),
                 boundary = summarise_boundary(boundary_r))
} else {
  tables <- recorded
}
```

### Small examples

``` r
knitr::kable(tables$small, row.names = FALSE,
  caption = "Five fresh processes per task; median [maximum] time and peak memory.")
```

| Example | Task | Success | Time (s) | Peak RSS (MiB) | Peak commit (MiB) |
|:---|:---|:---|:---|:---|:---|
| Original ternary | LP_UAP | 5/5 | 0.147 \[0.240\] | 75.3 \[75.4\] | 97.4 \[97.5\] |
| Original ternary | vertices | 5/5 | 0.012 \[0.033\] | 75.3 \[75.3\] | 97.5 \[97.5\] |
| Modified ternary | LP_UAP | 5/5 | 0.141 \[0.170\] | 75.3 \[75.3\] | 97.4 \[97.5\] |
| Modified ternary | vertices | 5/5 | 0.009 \[0.012\] | 75.3 \[75.4\] | 97.4 \[97.5\] |
| Quaternary | LP_UAP | 5/5 | 0.105 \[0.121\] | 75.3 \[75.4\] | 97.5 \[97.7\] |
| Quaternary | vertices | 5/5 | 4.951 \[6.090\] | 107.5 \[107.6\] | 131.0 \[131.1\] |

Five fresh processes per task; median \[maximum\] time and peak memory.

### Twenty random inputs per frame size

``` r
knitr::kable(tables$scaling, row.names = FALSE,
  caption = "Twenty instances per frame size; median [maximum] time and peak RSS.")
```

| n | Variables | Mean rows | Success | LP time (s) | Peak RSS (MiB) | Focals L1/HD/CW |
|---:|---:|---:|:---|:---|:---|:---|
| 5 | 32 | 27.4 | 20/20 | 0.021 \[0.042\] | 175.0 \[176.7\] | 21.8/21.5/20.3 |
| 6 | 64 | 57.5 | 20/20 | 0.014 \[0.016\] | 174.9 \[176.7\] | 33.5/32.0/30.2 |
| 7 | 128 | 105.7 | 20/20 | 0.023 \[0.037\] | 174.8 \[175.0\] | 48.4/44.9/40.2 |
| 8 | 256 | 226.8 | 20/20 | 0.065 \[0.080\] | 176.4 \[178.0\] | 71.7/61.8/53.9 |
| 9 | 512 | 356.2 | 20/20 | 0.115 \[0.280\] | 179.6 \[185.9\] | 89.3/77.0/64.2 |
| 10 | 1024 | 662.5 | 20/20 | 0.586 \[1.277\] | 193.0 \[204.2\] | 110.0/90.6/76.5 |

Twenty instances per frame size; median \[maximum\] time and peak RSS.

### Exploratory larger frames

``` r
knitr::kable(tables$boundary, row.names = FALSE,
  caption = "One instance per size. Memory guards and pre-allocation skips are reported separately.")
```

| n | Variables | Constraints | Nonzeros | Matrix estimate (MiB) | LP time (s) | Observed wall (s) | Peak RSS (MiB) | Peak commit (MiB) | Outcome |
|---:|---:|---:|---:|---:|:---|---:|---:|---:|:---|
| 11 | 2048 | 2046 | 175098 | 2.0 | 6.0 | 7.5 | 249.8 | 377.7 | solved |
| 12 | 4096 | 2214 | 283516 | 3.3 | 11.3 | 12.9 | 287.9 | 463.9 | solved |
| 13 | 8192 | 2790 | 579286 | 6.7 | 25.1 | 27.5 | 316.3 | 522.4 | solved |
| 14 | 16384 | 16382 | 4766584 | 54.6 | – | 4.4 | 674.3 | 1576.1 | process_memory_limit |
| 15 | 32768 | 9883 | 4351334 | 49.9 | – | 5.9 | 632.2 | 1248.9 | process_memory_limit |
| 20 | 1048576 | 275642 | 709297932 | 8121.3 | – | 2.7 | 262.2 | 286.9 | matrix_limit |

One instance per size. Memory guards and pre-allocation skips are
reported separately.

In the recorded run, all 120 scaling instances were solved and
validated. Three different mass vectors were obtained in 117 cases; in
three cases some objectives selected the same vector. The larger probes
completed through n = 13. The n = 14 and n = 15 probes reached the
process-memory guard. For n = 20, the estimated storage of one dominance
matrix alone was 7.93 GiB, so the matrix was not allocated.

These observations concern the recorded inputs and resource limits, not
an intrinsic limit on frame size. In particular, a guarded run does not
imply that a dominating belief function does not exist. On rerunning the
experiment, read the outcome columns rather than assuming the same
stopping points.

``` r
cat("```text\n",
    paste(trimws(capture.output(sessionInfo()), which = "right"), collapse = "\n"),
    "\n```\n", sep = "")
```

``` text
R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26200)

Matrix products: default


locale:
[1] C
system code page: 65001

time zone: Europe/Prague
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] rcdd_1.6       knitr_1.50     lpSolve_5.6.23

loaded via a namespace (and not attached):
 [1] compiler_4.3.2  fastmap_1.2.0   cli_3.6.5       tools_4.3.2
 [5] htmltools_0.5.9 yaml_2.3.11     rmarkdown_2.30  xfun_0.54
 [9] digest_0.6.39   rlang_1.1.6     evaluate_1.0.5
```
