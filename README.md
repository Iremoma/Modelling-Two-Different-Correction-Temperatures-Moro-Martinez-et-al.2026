# Divergent thermal responses of energy acquisition and mobilization Moro-Martinez et al., 2026

This repository contains all the code related to our work entitled **Divergent thermal responses of energy acquisition and mobilization explain decoupling between growth and metabolic rate in ectotherms**, by Irene Moro-Martinez, Andrea Campos-Candela, Joaquim Tomas-Ferrer and Miquel Palmer.
The work presented here is a contribution towards gaining greater ubderstanding of the thermal decoupling among physiological processes in Dynamic Energy Budget (**DEB**) models. [More about DEB theory [here](https://debportal.debtheory.org/docs/)]. For instance, we present a modified DEB model in which temperature dependence of assimilation energy flux is described using a unimodal Sharpe–Schoolfield curve, whereas temperature dependence of all other energy fluxes follow an Arrhenius (exponential) scaling.  We present our analytical strategy based on Bayesian inference developed in **STAN**. [More about STAN [here](https://mc-stan.org/)].

This repository is licensed under a [MIT License Copyright (c) 2026 Irene Moro-Martinez](./LICENSE). 

# Requirements
Install: 
- [R](https://www.r-project.org/about.html) (4.1.3 or higher)
- CmdStanR (Command Stan R) in R: CmdStanR is a lightweight interface to Stan for R users. If you are using Stan in R for the first time, or you have no experience with this package, we encourage you to visit [Getting started with CmdStanR]( https://mc-stan.org/cmdstanr/articles/cmdstanr.html)  and follow the recommended steps to install the package.

Recommended properties of the computer, and those that characterized the computer used for running these models, are: 
- Intel(R) Xeon(R)  CPU E5-2620 v4 @ 2.10GHz (40 cores and 64 GB RAM)

# Contents
We provide the R scripts necessary to reproduce the data exploration and analysis of the model presented on the paper Moro-Martinez _et al._, XXXX. **Note that the scripts are self-explanatory and have a lot of notations to guide the user on the different steps.** By running the code you will be able to reproduce the estimation of temperature parameters and explore the results.

All the code is organized in the following folder:

### [2cT DEB model]([https://github.com/Iremoma/Modelling-Two-Different-Correction-Temperatures-Moro-Martinez-et-al.2026/tree/main/2cT DEB model](https://github.com/Iremoma/Modelling-Two-Different-Correction-Temperatures-Moro-Martinez-et-al.2026/tree/main/2cT%20DEB%20model))

In this section you will find the STAN routine developed for estimating the DEB model parameters at the individual level for 136 juvenile European seabass (*Dicentrarchus labrax*), for which repeated measurements of length, wet weight and oxygen consumption were available. We provide:

- The input data file `input.RData` which the pre‑processed dataset used for the analysis.
- The STAN model (`DEB_O2.stan`) that implements the modified DEB framework with two temperature correction factors
- R code for running the Bayesian inference and saving the results
  
# Citation
If you use this code, please consider citing our work:

Moro-Martinez, I.; Campos-Candela, A.; Tomas-Ferrer, J.; and Palmer, M.; (submitted) Divergent thermal responses of energy acquisition and mobilization explain decoupling between growth and metabolic rate in ectotherms _Accepted in XXXX. Pending publication_.

Any questions: palmer@imedea.uib-csic.es , irenejues@gmail.com
