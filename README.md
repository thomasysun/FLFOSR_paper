# Fast longitudinal function-on-scalar regression (FLFOSR)

Paper: "Ultra-efficient MCMC for Bayesian longitudinal functional data analysis"


### Code

The code folder contains necessary files to fit the FLFOSR model and reproduce the results in the paper

* `./code/flfosr_act.R` produces results in the application section.
* `./code/flfosr_sims.R` produces results in the simulations section.
* `./code/flfosr1.R` code for fitting FLFOSR model.
* `./code/gen_flfosr.R` code that generates simulated longitudinal functional data for simulations.
* `./code/helper_functions.R` miscellaneous helper functions.
* `./code/process_nhanes_data_ts.R` loads and processes the NHANES physical activity data for application and plots.
* `./other_model_functions` contains code used to run models from competing methods, obtained from their respective repositories. Please refer to their documentation for more information.


#### Note on competing methods

The paper contains results from four competing methods:

* Gibbs sampler from `refund` package
* Variational Bayes from `refund` package
* Frequentist fixed effects inference with exchangeable correlations from Li et al. (2022) https://pubmed.ncbi.nlm.nih.gov/35491388/
* "FUI": Fast Univariate Inference from Cui et al. (2022) https://pubmed.ncbi.nlm.nih.gov/35712524/

The computation speed of the competing methods may be much slower, namely the Gibbs sampler from `refund`. Thus, these results are not run by default in the corresponding files. However, they may be included by setting `runrefund <- TRUE` inside the script.
