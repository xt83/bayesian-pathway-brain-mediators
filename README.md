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


## Citation

If you use this code, please cite the accompanying paper. The methodology and
code were originally distributed as supplemental files (`ujae132`); the
reference PDF is available alongside the journal article.

---

## License

No license file is currently included. Until a license is added, the default
copyright applies and reuse requires the authors' permission. If you intend to
reuse the code, please open an issue to request a license.
