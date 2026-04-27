# 2cT DEB model

This folder contains the STAN code and related R scripts to fit the Dynamic Energy Budget (DEB) model with two different temperature correction factors (2cT): a unimodal (Sharpe‑Schoolfield) correction for assimilation and an exponential (Arrhenius) correction for all other energy fluxes.

The model is used to reproduce the decoupling between growth and standard metabolic rate (SMR) observed in juvenile European seabass (*Dicentrarchus labrax*) under a natural thermal gradient (15.7–30.7 °C).

---

## Contents

| File | Description |
|------|-------------|
| `DEB_O2.stan` | STAN model code implementing the 2cT DEB framework. Includes the ODE system for state variables (reserve *E*, structure *V*, maturity *E_H*), temperature corrections, and the calculation of oxygen consumption. |
| `input.RData` | Processed dataset (length, weight, oxygen consumption) for 136 individual fish, along with all necessary constants and priors. |
| `out.RDS` | (Generated) Fitted model object. |
| `out_dignosi.txt` | (Generated) STAN diagnostic output (convergence, divergences, etc.). |

---

## Requirements

- R (≥ 4.1.3)
- `cmdstanr` package (see [installation guide](https://mc-stan.org/cmdstanr/))
- `deSolve` (for post‑fitting simulations if needed)

---

## Workflow

1. **Prepare environment**  
   Open R and set working directory to this folder.

2. **Load data and libraries**  
   Load `input.RData` which contains:
   - Observations (`obs`, `obs0`)
   - DEB fixed parameters (`f`, `E_G`, `kap`, ...)
   - Temperature parameters (sinusoidal fits per group)
   - Priors for all estimated parameters

3. **Compile the STAN model**  
   The script compiles `DEB_O2.stan` using `cmdstan_model()` with thread support for parallelisation.

4. **Run the inference**  
   The model samples the posterior using NUTS (2 chains, 1000 warmup + 1000 sampling iterations, 10 threads per chain).  
   *Note:* Running time is ~31.5 hours on a 40‑core server.

5. **Save outputs**  
   - `out.RDS` – fitted model object
   - `out_dignosi.txt` – diagnostic messages

---

## Key model features

- **Two temperature corrections**  
  - `cT1` (assimilation): unimodal (Sharpe‑Schoolfield)  
  - `cT2` (other fluxes): Arrhenius equation

- **Hierarchical structure**  
  - Individual‑level assimilation compound parameter `∀` (product of max assimilation and functional response)  
  - Common parameters for conductance, temperature coefficients, and error SDs

- **State variables**  
  - Energy reserve (*E*), structural volume (*V*), maturity (*E_H*)

- **Observables**  
  - Total length (`L_w`), wet weight (`W_w`), oxygen consumption rate (`J_O2`)

---

## Customisation

To adapt the model to your own data:

1. Replace `input.RData` with your own dataset (same variable names and structure).
2. Adjust priors and fixed parameters in the data block if necessary.
3. Modify the STAN model if you need alternative functional forms for temperature corrections.

---

## Citation

If you use this code, please cite the main paper:

Moro‑Martinez, I., Campos‑Candela, A., Tomàs‑Ferrer, J., & Palmer, M. (2026). Divergent thermal responses of energy acquisition and mobilization can explain decoupling between growth and metabolic rate in ectotherms. *Ecological Modelling*, submitted.

For questions: [palmer@imedea.uib-csic.es](mailto:palmer@imedea.uib-csic.es)
