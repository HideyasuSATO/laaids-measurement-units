# Household Unit Values and LA/AIDS

Replication code for **Household Unit Values, Measurement Units, and the Current-Share Stone Index in LA/AIDS Estimation**, by Hideyasu Sato.

`laaids_public_analysis.R` builds the Malawi analysis from the original survey files and runs the exact-AIDS Monte Carlo experiments. It produces numerical table inputs, figures, construction diagnostics, and execution metadata. The empirical and simulation analyses can run separately.

## Files

- `laaids_public_analysis.R`: the complete analysis.
- `food_group_mapping.csv`: the item classifications and inclusion rules used in the paper.
- `.gitignore`: excludes survey data and generated files.
- `LICENSE`: the MIT License for the author-provided code, documentation, and food classification rules.

The ignore rules allow only these public files, this README, and an optional `CITATION.cff` file, including when a custom output directory is used.

## Data

Obtain the Stata files for the [Malawi Fifth Integrated Household Survey 2019–2020](https://microdata.worldbank.org/index.php/catalog/3818), reference **MWI_2019_IHS-V_v06_M**, from the World Bank Microdata Library. Version 06 includes the associated Market Survey. Follow the catalog's **Get Microdata** link and access requirements. Data access and redistribution are governed by the provider's conditions; download the files directly rather than obtaining survey microdata from this repository.

Place these files together in `data/`, or set `data_dir` to their directory. Keep the filenames and capitalization unchanged.

```text
hh_mod_a_filt.dta
HH_MOD_META.dta
HH_MOD_B.dta
HH_MOD_C.dta
HH_MOD_G1.dta
ihs5_consumption_aggregate.dta
householdgeovariables_ihs5.dta
ihs_foodconversion_factor_2020.dta
mrk_mod_a.dta
mrk_mod_d.dta
mrk_mod_otherspecify.dta
```

All eleven files are required for the empirical run. No survey data are needed for `monte_carlo` mode. Intermediate household and item files are generated locally and should remain outside the public repository.

## Software

Use R with these packages installed:

```r
install.packages(c(
  "haven", "dplyr", "tidyr", "stringr", "lubridate", "readr",
  "purrr", "tibble", "rlang", "broom", "ggplot2", "forcats", "scales"
))
```

The analysis checks its dependencies and does not install packages automatically. Simulation mode uses `readr`, `tibble`, and `ggplot2`; estimation uses base R's `optim`. Stata is not required.

The saved paper runs used R 4.6.1 on Windows, with dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0, purrr 1.2.2, stringr 1.6.0, tibble 3.3.1, rlang 1.3.0, broom 1.0.13, ggplot2 4.0.3, forcats 1.0.1, and scales 1.4.0. Each new run records its own `sessionInfo()`; numerical or graphical differences can depend on software versions.

## Run

Edit the `settings` list at the start of `laaids_public_analysis.R`:

```r
settings <- list(
  mode = "all",
  data_dir = "data",
  output_dir = "outputs",
  mapping_file = "food_group_mapping.csv",
  bootstrap_reps = 999L,
  mc_reps = 500L,
  mc_households = 600L,
  save_png = TRUE
)
```

Paths are relative to the script directory unless absolute paths are supplied. Choose `"empirical"`, `"monte_carlo"`, or `"all"`. The default settings run both analyses at the paper's replication counts. Use a new, empty output directory for every run; the script refuses to reuse a nonempty directory.

From the repository directory, run:

```sh
Rscript laaids_public_analysis.R
```

Or, in R:

```r
source("laaids_public_analysis.R", encoding = "UTF-8")
```

The supplied settings also read `LAAIDS_MODE`, `LAAIDS_DATA_DIR`, `LAAIDS_OUTPUT_DIR`, `LAAIDS_BOOT_REPS`, `LAAIDS_MC_REPS`, and `LAAIDS_MC_HOUSEHOLDS`. These environment variables let you choose a run without editing the script. For example, use a separate destination for a short empirical check:

```r
Sys.setenv(
  LAAIDS_MODE = "empirical",
  LAAIDS_OUTPUT_DIR = "outputs_smoke_empirical",
  LAAIDS_BOOT_REPS = "3"
)
source("laaids_public_analysis.R", encoding = "UTF-8")
```

For a short simulation check:

```r
Sys.setenv(
  LAAIDS_MODE = "monte_carlo",
  LAAIDS_OUTPUT_DIR = "outputs_smoke_mc",
  LAAIDS_MC_REPS = "2"
)
source("laaids_public_analysis.R", encoding = "UTF-8")
```

Small runs check execution; they do not reproduce the paper's intervals or simulation summaries. Clear these overrides or start a fresh R session before a final run. A complete run uses 999 paired bootstrap replications and 10 simulation scenarios with 500 replications of 600 households each. The simulation estimates 13 specifications per replication.

## Outputs

All output paths below are relative to `output_dir`.

| Directory | Contents |
|---|---|
| `processed/input/` | Survey controls, item data, and market-price intermediates |
| `empirical/` | Baseline results, bootstrap results, data-construction manifests, and local analysis files |
| `robustness/` | Eleven scenario results and combined sensitivity summaries |
| `figures/` | Main Figures 1–5 when their corresponding analyses run, and Supplementary Figures S.1–S.6 |
| `diagnostics/empirical/` | Identity, sample, restriction, and bootstrap checks |
| `monte_carlo/` | Simulation results, manifests, diagnostics, and Figure 1 |

Figures are saved as PDF and, when enabled, PNG. Figure 1 is also retained in `monte_carlo/`. Tables are exported as numerical CSV inputs; the manuscript's typeset tables are not compiled by this script. `run_settings.csv`, `sessionInfo.txt`, and `output_manifest.csv` record the run configuration, software environment, and generated-file inventory.

| Paper item | Principal output |
|---|---|
| Main Figure 1; Supplementary Tables S.2–S.5 | `monte_carlo/fig1.pdf`; `mc_replication_level_key_comparisons_with_mcse.csv`, `mc_theoretical_channels_summary_by_scenario.csv`, `mc_beta_summary_with_mcse.csv`, and `mc_invariant_index_unit_scale_sensitivity_green_alston.csv` in `monte_carlo/` |
| Main Table 1 | `empirical/table_sample_flow.csv` and `table_market_price_matching_main6_cov90.csv` |
| Main Table 2; Figure 4 | `empirical/table_index_sensitivity.csv`, `empirical_bootstrap_index_sensitivity_with_ci.csv`; `figures/fig4.pdf` |
| Main Table 3; Figure 5; Supplementary Table S.8 | `empirical/table_expenditure_response_levels_current_stone.csv`, `empirical_bootstrap_expenditure_levels_with_ci.csv`; `figures/fig5.pdf` |
| Main Tables 4–5; Supplementary Tables S.9–S.12 | `empirical/table_price_elasticity_levels_current_stone.csv`, `table_price_elasticity_cells_current_stone.csv`, `table_price_elasticity_sensitivity_current_stone.csv`, and `empirical_bootstrap_price_*_with_ci.csv` |
| Main Figures 2(a), 2(b), and 3 | `figures/fig2a.pdf`, `fig2b.pdf`, `fig3.pdf`; `empirical/table_unit_value_market_benchmark_main6_cov90.csv` |
| Supplementary Tables S.6–S.7 | `empirical/conversion_factor_sources_detailed.csv` and `table_group_coverage_main6.csv` |
| Supplementary Figures S.1–S.6 | `figures/figS1.pdf` through `figS6.pdf` |
| Supplementary Tables S.13–S.15 | `robustness/robustness_index_sensitivity.csv`, `robustness_expenditure_sensitivity.csv`, and `robustness_price_elasticity_sensitivity.csv` |
| Supplementary Tables S.1 and S.16–S.19 | Scenario, true-parameter, DGP-equation, model-specification, and optimizer manifests in `monte_carlo/` |

## Interpretation and reproduction checks

The empirical baseline uses purchased foods in six groups: cereals; roots and tubers; pulses and nuts; vegetables; animal protein; and oils, sugar, and condiments. Household conversion coverage must be at least 90%. Official conversions and documented exact/manual rules are used; the universal rough liquid fallback is included only in a robustness comparison. The code retains the paper's item and group winsorization, price-fill hierarchies, and fixed pooled item weights for market-price aggregation.

The main contrast holds market-price regressors fixed and changes raw versus standardized unit values inside the price index. Five share equations use household-weighted least squares; the vegetables equation is recovered by adding-up. Homogeneity is imposed through relative prices; symmetry and curvature are unrestricted. Zero shares are retained. Expenditure elasticities are conditional on purchased expenditure within the included groups. Price outputs distinguish the one-for-one LA-prime, partial-regressor, and Green–Alston conventions.

Bootstrap draws resample EA clusters within district strata, using the same draw for the raw and standardized specifications. Index reference means and weights, models, and elasticities are recomputed. Conversion factors, item/group winsorization thresholds, hierarchical fills, and constructed household variables stay fixed. The intervals therefore condition on the constructed household–group file. The empirical seed is 20260902; simulations retain their scenario-specific seeds and sequential random-number streams.

Useful checks for the paper settings are:

| Quantity | Paper value |
|---|---:|
| Household source sample | 11,434 |
| Final six-group estimation sample | 10,040 |
| Bootstrap design | 717 EA clusters; 32 district strata; 999 successful replications |
| Matched item observations in the market benchmark | 101,051 |
| Expenditure-weighted log correlation, raw / standardized | 0.192949 / 0.819103 |
| Centered mean absolute index gap, current Stone / corrected Stone / Törnqvist / base share | 0.627707 / 0.493293 / 0.402148 / 0.339335 |
| Mean absolute conventional expenditure-elasticity difference | 0.307047 |

Inspect the bootstrap design and success diagnostics, robustness scenario statuses, and Monte Carlo convergence outputs. The paper's simulation run had 65,000 converged model fits and no failed replications. These checks concern the maintained specification and its sensitivity; they do not establish that either empirical elasticity estimate is closest to an unknown true value.

## License

The author-provided analysis code, documentation, and food classification rules are available under the [MIT License](LICENSE), copyright (c) 2026 Hideyasu Sato.

Survey data are not distributed with this repository and are not covered by this license. Obtain them from the provider and follow the applicable access and use conditions. Third-party R packages retain their own licenses.
