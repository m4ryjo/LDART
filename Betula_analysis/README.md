## Introduction

In this study, we focus on two key effect measures that address
truncation by death under hypothetical interventions:

-   **PCIE (Partly Conditional Intervention Effect)** Partly conditional
    models make inference conditional on the sub-population alive at a
    specific time-point, in our application being alive at end of
    follow-up, i.e., time *T*. This approach focuses on prognostic
    questions:
    PCIE = *E*\[*Y*<sub>*T*</sub>(*g*<sub>\*</sub>) ∣ *S*<sub>*T*</sub>(*g*<sub>\*</sub>) = 1\] − *E*\[*Y*<sub>*T*</sub>(*g*<sub>0</sub>) ∣ *S*<sub>*T*</sub>(*g*<sub>0</sub>) = 1\].

-   **SAIE (Survivor Average Intervention Effect)** Conditioning on
    survival as in partly conditional models may introduce bias because
    survival is a post-randomization event . To address etiological
    questions, we consider the subpopulation who would survive under
    both treatment regimes (the “always survivor” stratum) . The
    estimand is:
    SAIE = *E*\[*Y*(*g*<sub>\*</sub>) − *Y*(*g*<sub>0</sub>) | *S*<sub>*T*</sub>(*g*<sub>0</sub>) = *S*<sub>*T*</sub>(*g*<sub>\*</sub>) = 1\].

## The PCIE

Under assumptions **C1**, **C2a**, **C3**, and **C4**, the mean outcome
at time *T* under the natural course (regime *g*<sub>0</sub>) can be
computed using the **g-formula**, conditioning on survival.

For the intervention regime *g*<sub>\*</sub>, under assumptions **C1**,
**C2b**, **C3**, and **C4**, the mean outcome at time *T* can be
computed using the **Extended g-formula**.

Both the g-formula and the Extended g-formula can be estimated using the
functions provided in the GcompBART package, which implements Bayesian
machine learning methods for flexible estimation under hypothetical
interventions.

The **PCIE** is then obtained as the difference in expected outcomes
under the two regimes:
PCIE = *E*\[*Y*<sub>*T*</sub>(*g*<sub>\*</sub>) ∣ *S*<sub>*T*</sub>(*g*<sub>\*</sub>) = 1\] − *E*\[*Y*<sub>*T*</sub>(*g*<sub>0</sub>) ∣ *S*<sub>*T*</sub>(*g*<sub>0</sub>) = 1\].

## The SAIE

### Sensitivity Parameter *λ* (Assumption C5)

The sensitivity parameter *λ* ranges from 0 to 1 and governs assumptions
about survival under different treatment regimes:

-   ***λ* = 1**: Assumes *deterministic monotonicity*, meaning anyone
    who would survive under the observed regime (*g*<sub>0</sub>) would
    also survive under the intervention (*g*<sub>\*</sub>), i.e.,
    *S*<sub>*T*</sub>(*g*<sub>0</sub>) ≤ *S*<sub>*T*</sub>(*g*<sub>\*</sub>).
-   ***λ* = 0**: Assumes *independence* between survival under the two
    regimes (*S*<sub>*T*</sub>(*g*<sub>0</sub>) and
    *S*<sub>*T*</sub>(*g*<sub>\*</sub>)), which is considered unlikely
    in this application.
-   **Intermediate values 0 &lt; *λ* &lt; 1**: Represent partial
    relaxation of monotonicity. As *λ* decreases from 1 toward 0, the
    assumption allows increasing dependence between survival under the
    two regimes to weaken, introducing more uncertainty about who would
    survive under the intervention compared to the observed regime.

This parameter is used to compute the **Survivor Average Causal Effect
(SACE)** under varying assumptions about survival relationships across
regimes.

### Sensitivity Parameter for SAIE: *Δ* (Assumption C6)

**C6**: *Difference in expectations of outcomes when comparing different
strata.*

For *g* ∈ {*g*<sub>\*</sub>, *g*<sub>0</sub>} and
*g*′ ∈ {*g*<sub>\*</sub>, *g*<sub>0</sub>} with *g* ≠ *g*′, we assume:
*E*\[*Y*<sub>*T*</sub>(*g*) ∣ *S*<sub>*T*</sub>(*g*′) = *S*<sub>*T*</sub>(*g*) = 1\] − *E*\[*Y*<sub>*T*</sub>(*g*) ∣ *S*<sub>*T*</sub>(*g*) = 1, *S*<sub>*T*</sub>(*g*′) = 0\] = *Δ*<sup>*g*</sup>,
where *Δ*<sup>*g*</sup> is a sensitivity parameter. We further assume:
*Δ*<sup>*g*</sup> = *Δ*<sup>*g*′</sup> = *Δ*,
i.e., the difference in expectations of potential outcomes when
comparing different strata is the same for the two contrasting regimes.

-   If *Δ* = 0: There is **no difference** in potential outcomes between
    the always-survivor stratum and the stratum where individuals would
    survive under regime *g* but not under *g*′.
-   If *Δ* &gt; 0: The always-survivor stratum has **higher potential
    outcomes**. In our application, this corresponds to higher memory
    scores (better memory) for the always-survivor stratum.

### The SAIE Formula

The SAIE is given by:
$$
    \mathrm{SAIE} = \mathrm{PCIE} + \Delta \Big\\ \psi^{g\_\*} + \lambda \big(U - \psi^{g\_\*}\big) \Big\\ \Big(1 - \frac{1}{U}\Big),
    $$
where:
$$
    U = \min\Big\\ 1, \frac{\psi^{g\_\*}}{\psi^{g\_0}} \Big\\.
    $$

This equation represents a **location shift upwards from the PCIE** for
a beneficial intervention, where we expect participants to improve their
health under the intervention compared to the natural course. The
adjustment depends on: - *Δ*: Sensitivity parameter for differences in
potential outcomes across strata. - *λ*: Sensitivity parameter
controlling dependence between survival under regimes. -
*ψ*<sup>*g*<sub>\*</sub></sup> and *ψ*<sup>*g*<sub>0</sub></sup>:
Survival probabilities under the intervention and natural course.

*Example: Workflow for estimating PCIE using GcompBART.*

    # Example: Estimating PCIE using GcompBART
    # 1. Fit Bayesian models for natural course
    BM <- BMfits(your_data,
                 var.type = your_vartype_nc,
                 opts = your_opts,
                 tgroup = your_tgroup,
                 base_hypers = your_basehypers)

    # 2. Compute g-formula under natural course
    out_nc <- gcompbart(your_data,
                 var.type = your_vartype_nc,
                 opts = your_opts,
                 tgroup = your_tgroup,
                 J = 10000,
                 BModels = BM)

    # 3. Compute extended g-formula under intervention regime
    out_int <- gcompbart(your_data,
                 var.type = your_vartype_int,
                 random.regime = rep("triangular", 3),
                 param = list(delta, delta, delta),
                 nat_value = TRUE,
                 above = TRUE,
                 threshold = sb_threshold,
                 incremental = TRUE,
                 opts = your_opts,
                 tgroup = your_tgroup,
                 J = 10000,
                 BModels = BM)

### Computing PCIE

Compute PCIE as difference between intervention and natural course
estimates using the extracted predicted means (y\_hat).

    pcie <- out_int$y_hat - out_nc$y_hat

### Computing SAIE

The Survivor Average Intervention Effect (SAIE) can be computed by
combining the PCIE with an adjustment term that accounts for the “always
survivor” stratum. Below is an illustrative example using a helper
function sace():

    ps_nc <- out_nc$s_hat
    ps_int <- out_int$s_hat

    saie <- pcie + sace(Delta, lambda, ps_int, ps_nc)
