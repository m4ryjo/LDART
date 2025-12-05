---
  title: "Computing the PCIE and the SAIE"
author: "Maria Josefsson"
date: "`r Sys.Date()`"
output:
  html_document:
  toc: true
toc_float: true
number_sections: true
---

  ## Introduction

  In this study, we focus on two key effect measures that address truncation by death under hypothetical interventions:

  - **PCIE (Partly Conditional Intervention Effect)**
  Partly conditional models make inference conditional on the sub-population alive at a specific time-point, in our application being alive at end of follow-up, i.e., time \(T\). This approach focuses on prognostic questions:
  \[
    \mathrm{PCIE} = E\big[ Y_T(g_*) \mid S_T(g_*) = 1 \big] - E\big[ Y_T(g_0) \mid S_T(g_0) = 1 \big].
    \]

  - **SAIE (Survivor Average Intervention Effect)**
  Conditioning on survival as in partly conditional models may introduce bias because survival is a post-randomization event \cite{zhang03}. To address etiological questions, we consider the subpopulation who would survive under both treatment regimes (the “always survivor” stratum) \cite{frangakis02, frangakis07}. The estimand is:
  \[
    \mathrm{SAIE} = E\big[ Y(g_*) - Y(g_0) \,\big|\, S_T(g_0) = S_T(g_*) = 1 \big].
    \]

## The PCIE

Under assumptions **C1**, **C2a**, **C3**, and **C4**, the mean outcome at time \(T\) under the natural course (regime \(g_{0}\)) can be computed using the **g-formula**, conditioning on survival.

For the intervention regime \(g_{*}\), under assumptions **C1**, **C2b**, **C3**, and **C4**, the mean outcome at time \(T\) can be computed using the **Extended g-formula**.

The **PCIE** is then obtained as the difference in expected outcomes under the two regimes:
  \[
    \mathrm{PCIE} = E\big[ Y_T(g_*) \mid S_T(g_*) = 1 \big] - E\big[ Y_T(g_0) \mid S_T(g_0) = 1 \big].
    \]


## The SAIE
### Sensitivity Parameter \(\lambda\) (Assumption C5)

The sensitivity parameter \(\lambda\) ranges from 0 to 1 and governs assumptions about survival under different treatment regimes:

- **\(\lambda = 1\)**: Assumes *deterministic monotonicity*, meaning anyone who would survive under the observed regime (\(g_0\)) would also survive under the intervention (\(g_*\)), i.e., \(S_T(g_0) \leq S_T(g_*)\).
- **\(\lambda = 0\)**: Assumes *independence* between survival under the two regimes (\(S_T(g_0)\) and \(S_T(g_*)\)), which is considered unlikely in this application.
- **Intermediate values \(0 < \lambda < 1\)**: Represent partial relaxation of monotonicity. As \(\lambda\) decreases from 1 toward 0, the assumption allows increasing dependence between survival under the two regimes to weaken, introducing more uncertainty about who would survive under the intervention compared to the observed regime.

This parameter is used to compute the **Survivor Average Causal Effect (SACE)** under varying assumptions about survival relationships across regimes.

### Sensitivity Parameter for SAIE: \(\Delta\) (Assumption C6)

**C6**: *Difference in expectations of outcomes when comparing different strata.*

  For \(g \in \{g_*, g_0\}\) and \(g' \in \{g_*, g_0\}\) with \(g \neq g'\), we assume:
  \[
    E\big[ Y_T(g) \mid S_T(g') = S_T(g) = 1 \big] - E\big[ Y_T(g) \mid S_T(g) = 1, S_T(g') = 0 \big] = \Delta^g,
    \]
where \(\Delta^g\) is a sensitivity parameter. We further assume:
  \[
    \Delta^g = \Delta^{g'} = \Delta,
\]
i.e., the difference in expectations of potential outcomes when comparing different strata is the same for the two contrasting regimes.

- If \(\Delta = 0\): There is **no difference** in potential outcomes between the always-survivor stratum and the stratum where individuals would survive under regime \(g\) but not under \(g'\).
- If \(\Delta > 0\): The always-survivor stratum has **higher potential outcomes**. In our application, this corresponds to higher memory scores (better memory) for the always-survivor stratum.

### The SAIE Formula

The \(\text{SAIE}\) is given by:
  \[
    \mathrm{SAIE} = \mathrm{PCIE} + \Delta \Big\{ \psi^{g_*} + \lambda \big(U - \psi^{g_*}\big) \Big\} \Big(1 - \frac{1}{U}\Big),
    \]
where:
  \[
    U = \min\Big\{ 1, \frac{\psi^{g_*}}{\psi^{g_0}} \Big\}.
    \]

This equation represents a **location shift upwards from the PCIE** for a beneficial intervention, where we expect participants to improve their health under the intervention compared to the natural course. The adjustment depends on:
  - \(\Delta\): Sensitivity parameter for differences in potential outcomes across strata.
- \(\lambda\): Sensitivity parameter controlling dependence between survival under regimes.
- \(\psi^{g_*}\) and \(\psi^{g_0}\): Survival probabilities under the intervention and natural course.

---

  ### Example: Compute SACE for a Grid of \(\lambda\) Values

  ```{r}
# Define lambda grid
lambda_grid <- seq(0, 1, by = 0.25)

# Example sensitivity parameters (replace with your actual estimates)
# Assume we have mean outcomes for survivors under g0 and g*
mean_g0 <- 50   # observed regime
mean_gstar <- 55 # intervention regime

# Assume survival probabilities under g0 and g*
p_survive_g0 <- 0.8
p_survive_gstar <- 0.85

# Function to compute SACE given lambda
compute_sace <- function(lambda, y_g0, y_gstar, p_g0, p_gstar) {
  # Weight survival probabilities based on lambda
  # Here, lambda adjusts correlation between survival under g0 and g*
  # Simplified example: effective survivor set size
  effective_survivors <- lambda * p_g0 + (1 - lambda) * (p_g0 * p_gstar)

  # Compute SACE as difference in means among effective survivors
  sace <- (y_gstar - y_g0) * effective_survivors
  return(sace)
}

# Compute SACE for each lambda
sace_values <- sapply(lambda_grid, compute_sace,
                      y_g0 = mean_g0,
                      y_gstar = mean_gstar,
                      p_g0 = p_survive_g0,
                      p_gstar = p_survive_gstar)

# Combine results
data.frame(lambda = lambda_grid, SACE = sace_values)



mci <- function(x, y=3){
  round(c(quantile(x, c(0.5, 0.025, 0.975))), y)
}

phi <- function(a,b) {
  ifelse(a/b<1, 1, 0) + ifelse(a/b<1, 0, 1)*a/b
}

# sace <- function(Delta_max, lambda, ps_int, ps_bl) {
#     #tmp <- ifelse((ps_int - ps_bl)/ps_int < 0 , 0, 1)*(ps_int - ps_bl)/ps_int
#     U <- phi(mean(ps_int), mean(ps_bl))
#     Delta <- Delta_max#rtri(1, 0, Delta_max, Delta_max-(Delta_max/100))
#     Delta*((mean(ps_int) + lambda*(U - mean(ps_int)))*(1 - 1/U))
# }

sace <- function(Delta_max, lambda, ps_int, ps_bl) {
  #tmp <- ifelse((ps_int - ps_bl)/ps_int < 0 , 0, 1)*(ps_int - ps_bl)/ps_int
  U <- phi(ps_int, ps_bl)
  Delta <- rtri(1, 0, Delta_max, Delta_max-(Delta_max/100))
  Delta*((mean(ps_int) + lambda*(U - mean(ps_int)))*(1 - 1/U))
}

sace <- function(Delta_max, lambda, ps_int, ps_bl) {
  #U <- phi(mean(ps_int), mean(ps_bl))
  #U2 <- phi(mean(ps_bl), mean(ps_int))
  U <- phi(ps_int, ps_bl)
  U2 <- 1/U#phi(ps_bl, ps_int)
  Delta <- rtri(1, 0, Delta_max, Delta_max-(Delta_max/100))
  lambda <- rtri(1, 0.5, 1, 1-(1/100))
  # Delta*((mean(ps_int) + lambda*(U - mean(ps_int))) - (mean(ps_bl) + lambda*(U2 - mean(ps_bl))))
  Delta*((ps_int + lambda*(U - ps_int))*(1 - U2))
}
