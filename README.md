# Bayesian Pathway Analysis over Brain Network Mediators for Survival Data

A Bayesian framework for pathway (mediation) analysis in which the mediator is a
subject-specific **brain structural connectivity matrix** and the outcome is a
(possibly right-censored) **survival time**. The model jointly estimates

1. the effect of an exposure on the brain network (exposure → mediator), and
2. the effect of the brain network on time-to-event (mediator → outcome),

so that indirect (pathway) effects can be decomposed at the level of latent
sub-networks and individual nodes.

The reference application in this repository is the **ADNI** study, using APOE4
status as the exposure, structural connectivity derived from diffusion MRI as
the mediator, and time to Alzheimer's disease onset as the survival outcome.

---

## Model overview

Let `Z_i` be a scalar exposure, `X_i` a vector of covariates, `A_i` a
`p × p` symmetric connectivity matrix (mediator), and `T_i` a survival time with
censoring indicator `δ_i`. The model links these through two regressions that
share a low-rank sub-network representation:

- **Exposure → mediator (matrix response regression):**
  `A_i = Σ_h (a2h' X_i) · a1h a1h' + Σ_{k=1}^{K1} η_k · Z_i · α_k α_k' + E_i`
  The first term absorbs covariate effects on connectivity; the second term
  encodes the exposure's effect through `K1` latent sub-networks defined by
  loading vectors `α_k`, with sub-network effects `η_k`.

- **Mediator → outcome (AFT survival regression):**
  `log T_i = X_i' β_x + Z_i β_z + Σ_{k=1}^{K2} ω_k · β_k' A_i β_k + ε_i`
  The quadratic forms `β_k' A_i β_k` summarize `K2` latent sub-networks of the
  connectome that drive survival, with sub-network effects `ω_k`. Censored
  observations are imputed from a truncated normal in the AFT likelihood.

- **Sparsity and selection:** spike-and-slab priors on the node-level loadings
  `α_k` and `β_k` (latent indicators `t_{k,j}`, `v_{k,j}`) perform
  node selection within each sub-network. Bayesian lasso / inverse-Gaussian
  hyperpriors shrink the sub-network effects `η` and `ω`.

- **Indirect pathway effects** are obtained as posterior functionals of
  `(α, η)` (exposure → mediator) combined with `(β, ω)` (mediator → outcome).

Posterior inference is carried out by a **Gibbs sampler** whose full
conditionals are implemented in `mediation_ADNI.R`.

---

## Repository contents

| File | Description |
|---|---|
| `mediation_ADNI.R` | Full R implementation: Gibbs updates for all model blocks plus the top-level `mcmc()` driver and the ADNI analysis entry point. |
| `README.md` | This file. |

The script is self-contained: each MCMC update is a standalone function, and
`mcmc()` wires them together.

### Main components of `mediation_ADNI.R`

- `update_betax`, `update_betaz` — covariate and exposure coefficients in the
  AFT outcome model.
- `update_a1h`, `update_a2h` — rank-`h` covariate adjustment of the
  connectivity matrix.
- `update_alpha_k`, `update_eta` — exposure → mediator sub-network loadings and
  effects.
- `update_beta_k`, `update_omega` — mediator → outcome sub-network loadings
  and effects.
- `update_v`, `update_t` — spike-and-slab selection indicators for node-level
  sparsity.
- `update_tau_eta`, `update_tau_omega` — Bayesian lasso scale updates.
- `update_A`, `update_y`, `update_sigma0`, `update_sigmae` — latent
  connectivity, imputed log-survival times, and error variances.
- `mcmc(...)` — Gibbs driver returning posterior draws of all parameters.

---

## Requirements

- **R** (≥ 4.0 recommended)
- R packages:
  - `mvtnorm`
  - `invgamma`
  - `SuppDists`
  - `tensorr`
  - `Matrix`
  - `truncnorm`

Install with:

```r
install.packages(c("mvtnorm", "invgamma", "SuppDists",
                   "tensorr", "Matrix", "truncnorm"))
```

---

## Input data format

`mediation_ADNI.R` expects two CSV files in the working directory (paths can be
edited at the top of the script):

- **`matrix.csv`** — a `(p · n) × (p + 1)` stack of `p × p` structural
  connectivity matrices, one block per subject. The first column is a row
  identifier that is dropped on load. With the ADNI default `p = 83` (e.g.
  Desikan-Killiany cortical/subcortical parcellation), each subject contributes
  an 83-row block.
- **`demo.csv`** — subject-level table containing
  - `outcome` — observed time to AD onset (the script models `log(outcome)`),
  - `event` — censoring indicator (`1` = event observed, `0` = censored),
  - `APOE4` — exposure (APOE-ε4 status), and
  - three covariate columns at positions 8–10 (e.g. age, sex, education).

The row order of `demo.csv` must correspond to the block order of
`matrix.csv`.

> **Note.** The script uses an intercept-augmented design matrix
> `x <- cbind(1, x1)`. Ensure the covariate matrix is assigned to `x1` before
> running, or adapt this line to your data.

---

## Usage

After preparing `matrix.csv` and `demo.csv`:

```r
source("mediation_ADNI.R")

result <- mcmc(
  xtrain  = x,                        # covariates (with intercept)
  ztrain  = APOE4,                    # exposure
  Atrain  = structural_connectivity,  # stacked p*n x p connectivity
  ytrain  = log_outcome,              # log survival time
  censor  = censor_status,            # 1 = event, 0 = censored
  numite  = 10000,                    # MCMC iterations
  burn_in = 2000,
  k1      = 3,                        # # sub-networks: exposure -> mediator
  k2      = 3,                        # # sub-networks: mediator -> outcome
  h       = 2,                        # rank of covariate adjustment
  p       = 83                        # # ROIs (mediator dimension)
)
```

`mcmc()` returns a list of posterior draws including:

1. `pre_train` — predicted log survival times per iteration,
2. `A_new` — `α_k` (exposure → mediator loadings) per iteration,
3. `B_new` — `β_k` (mediator → outcome loadings) per iteration,
4. `v_new`, `t_new` — selection indicators,
5. `eta_train`, `omega_train` — sub-network effects,
6. `a1h_new`, `a2h_new` — covariate-adjustment factors,
7. `betax_new`, `betaz_new` — AFT coefficients,
8. `sigma0_new`, `sigmae_new` — residual variances.

Discard the first `burn_in` iterations before computing posterior summaries.
Indirect (pathway) effects can then be formed from paired draws of
`(α_k, η_k)` and `(β_k, ω_k)`.

---

## Tuning and practical notes

- **`k1`, `k2`** control the number of latent sub-networks driving the
  exposure → mediator and mediator → outcome pathways, respectively. Start
  small (2–4) and assess via posterior predictive checks.
- **`h`** is the rank of the covariate adjustment on the connectome. `h = 1`
  or `2` is usually sufficient.
- **`v0`, `v1`, `v2`** are the spike and slab scales of the sparsity prior.
  The defaults (`v0 = sqrt(0.1)`, `v1 = v2 = 2`) give a clear
  "near-zero vs. nonzero" contrast; tighten `v0` for stronger sparsity.
- The sampler is written for clarity rather than speed; for large `p` or large
  `n` consider vectorizing the inner loops or porting hot functions to C++
  (e.g. via `Rcpp`).

---

## Data availability

The ADNI dataset used in the reference analysis is **not** redistributed in
this repository. Access can be requested from the
Alzheimer's Disease Neuroimaging Initiative at
<https://adni.loni.usc.edu/>. Users must substitute their own
`matrix.csv` / `demo.csv` following the format above.

---

## Citation

If you use this code, please cite the accompanying paper. The methodology and
code were originally distributed as supplemental files (`ujae132`); the
reference PDF is available alongside the journal article.

---

## License

No license file is currently included. Until a license is added, the default
copyright applies and reuse requires the authors' permission. If you intend to
reuse the code, please open an issue to request a license.
