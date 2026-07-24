# README

Script to run cost-effectiveness analysis comparing IGRA and C-Tb testing for latent TB infection (LTBI), used for the [article](https://docs.google.com/document/d/1Neh9EVPZR6WmqFKuLN1wQFS35cu_wnxDa80BpOQXrgI/). The model is from [Vesga 2025], adapted for C-Tb.

## THIS MODEL WAS REVISED AND TREATMENT COSTS AND OTHER PARAMETERS CHANGE. THE UPDATED VERSION IS AT THE REPO [cTB_cea_UK2026](https://github.com/juanvesga/cTB_cea_UK2026)

# How to run?

Download the whole package, install libraries. In the file `14-run-incremental-case.R` the sensitivity analysis parameters can be adjusted as follows:

-   `uptake_val`: minimum value for uptake (1=100%), which is also the minimum for the return rate (r_ret) and treatment acceptance rate (r_acc)
-   `increm_val`: increments in r_ret and r_acc rates (0.01=1%)
-   `price_vals`: range of unit prices to include in the analysis
-   `pars_true_pos_ctb`: true positivity rate for C-Tb. By default it is set to the same value as QFT (stored in `model_parameters$model_pars$true_pos_ctb`). Enter a value manually to change it. In the 'higher sensitivity/specificity' scenario, this parameter was set to 1, which means equivalent performance to T-SPOT.TB.
-   `pars_true_neg_ctb`: true negativity rate for C-Tb. Same as for previous parameter.
-   `neg_test_diff_return`: scaling factor for return rate of negative result (C-Tb) patients. Set to 1 (=100%) by default. Values lower than 1 mean that those with a negative result (no induration) less likely to return for a second appointment.

To run the analysis enter `source("run.R")`.

Results will be saved in an appropriately named subfolder in `results/`. Eg. if the minimum uptake rate was set to 0.45 (45%), the increment to 0.02 (2% steps), the minimum price to 15, the maximum price to 25, the price increments to 2.5 and the scaling factor for negative return rate patients to 0.5 (50% of positive patients' return rate), then the subfolder will be named `uptakeval0.45incr0.02_pricemin15max25incr2.5_negretrate0.5`.

The results used for plots in the article are in `all_results_single_incremental_case.csv`.
