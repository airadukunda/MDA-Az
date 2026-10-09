# Azithromycin MDA: age-structured AMR transmission models

R scripts for modelling macrolide-resistant carriage in *Escherichia coli*, *Streptococcus pneumoniae*, and *Staphylococcus aureus* following azithromycin mass drug administration (MDA). The model structure follows the mixed-carriage framework of [Davies et al](https://www.nature.com/articles/s41559-018-0786-x), with age-specific demography and contacts informed by Tanzania (and soon, Malawi and Niger). Main scenarios compare no MDA with annual, twice-yearly, and quarterly MDA, with treatment ending after year 10 in the 20-year simulations.

## Repository files

| File(s) | Purpose |
|---|---|
| `davies_2_8_tidy_{ecoli,spneumo,saureus}.R` | Define pathogen-specific models and run baseline equilibrium, MDA scenarios, summaries and figures. S. pneumoniae also generates the Malawi MORDOR comparison. |
| `{ecoli,spneumo,saureus}_oat_sensitivity_recalibrated_common_clearance.R` | Define recalibrated one-at-a-time (OAT) sensitivity cases. Re-fit transmission and baseline antibiotic pressure for each case; modify all carriage-clearance rates together for the clearance sensitivity. |
| `{ecoli,spneumo,saureus}_oat_worker.R` | Run one OAT sensitivity case per `SLURM_ARRAY_TASK_ID`, saving an `.rds` record. |
| `{ecoli,spneumo,saureus}_oat_combine.R` | Assemble the expected successful worker outputs, validate baseline-fit diagnostics and write OAT summaries / tornado CSVs. |
| `make_combined_report_figures.R` | Make the three-pathogen 20-year resistance time series and the six-panel annual/biannual **absolute-prevalence** sensitivity figure from saved CSVs. |
| `plot_oat_tornado_from_existing_csvs.R` | Make annual/biannual tornado plots for **percentage-point increase versus matched no-MDA**, from precomputed tornado CSVs without solving ODEs. |
| `fitness_cost_half_life_analysis.R` | Optional, computationally intensive fitness-cost analysis for *E. coli* and *S. pneumoniae* only, with post-MDA half-life estimates. |

All scripts are at repository root because they currently source one another by relative filename. Run from that directory.

## Required inputs

Place these four inputs in the repository root, or modify each main model's `config$data_dir` and its `input_files` list:

```text
Population_Afro_2023_1yearage.csv
3.U.1.Birth_1year_Afro.csv
3.U.1.AFRO_mortality_by_age_group_1yearage.csv
3.U.1.contact_Tanzania_1y.csv
```

Install R packages `deSolve`, `dplyr`, `ggplot2`, `purrr`, `readr`, `scales`, `tidyr`, `tibble`, and `stringr`. 

The long equilibrium calculations, especially for *S. aureus*, can be computationally expensive.

## Workflow

```text
   Demographic and contact inputs
                |
                v
       Pathogen model functions
          /              \
         v                v
 Main scenario runs    Recalibrated OAT cases
         |                |
         v                v
  Scenario CSVs      SLURM case RDS files
         |                |
         |                v
         |          Combine OAT results
         |                |
         v                v
  Combined trajectories and sensitivity figures

Optional: fitness-cost analysis -> persistence figures
```

### 1. Run the main models

```bash
Rscript davies_2_8_tidy_ecoli.R
Rscript davies_2_8_tidy_spneumo.R
Rscript davies_2_8_tidy_saureus.R
```

They write to `outputs_tanzania_mda/`, `outputs_tanzania_spneumo_mda/`, and `outputs_tanzania_mda_saureus/`, respectively, including `scenario_time_series.csv` and age-stratified outputs.

### 2. Run the OAT sensitivity analyses

The three `*_oat_sensitivity_recalibrated_common_clearance.R` files define the sensitivity cases and can be run directly, **or** used with the worker scripts to parallelise one case per array task. Do not run both approaches into the same output directory. The worker approach is preferred for the computationally intensive runs.

To run the first case of each worker locally (expensive):

```bash
SLURM_ARRAY_TASK_ID=1 Rscript ecoli_oat_worker.R
SLURM_ARRAY_TASK_ID=1 Rscript spneumo_oat_worker.R
SLURM_ARRAY_TASK_ID=1 Rscript saureus_oat_worker.R
```

On SLURM, submit each worker with `--array=1-N`, where `N = nrow(oat_specs)` **for that pathogen**. Each completed task writes `.../sensitivity_oat_<pathogen>_recalibrated/slurm_cases/case_###_<case_id>.rds`; failed tasks produce `_ERROR.rds`.

For each sensitivity case, the scripts calibrate transmission and background antibiotic pressure to recover the reference no-MDA carriage and resistance prevalence. The latest cluster copies reject calibration fits more than 1 percentage point from either reference target. Use the resulting `baseline_*_error_pp` columns as a quality check.

**Scope:** only annual, biannual and no-MDA sensitivity trajectories are included in the sensitivity analysis

### 3. Combine array results

Only after *all* expected worker cases have completed:

```bash
Rscript ecoli_oat_combine.R
Rscript spneumo_oat_combine.R
Rscript saureus_oat_combine.R
```

The combine scripts read the expected success filenames, reject missing results and unacceptable calibration errors, and ignore unrelated or failed records. They write `{pathogen}_oat_sensitivity_summary.csv`, time-series CSVs, and `{pathogen}_oat_tornado_{annual,biannual}_10y_all_ages.csv` in the corresponding recalibrated sensitivity directory.

### 4. Make combined figures without rerunning the models

```bash
Rscript make_combined_report_figures.R
```

This reads the three main scenario CSVs and three **recalibrated** sensitivity-summary CSVs. The 20-year legend order is annual, biannual, quarterly, no MDA. The six-panel OAT figure uses **absolute prevalence among carriers** at ten years.

To instead plot the **increase relative to matched no MDA** using already-processed annual/biannual tornado CSVs:

```bash
OAT_CSV_DIR=. Rscript plot_oat_tornado_from_existing_csvs.R
```

`OAT_CSV_DIR` should contain only one version of each input CSV. The plot script stops when it finds ambiguous duplicates.

### 5. Optional post-MDA persistence analysis

```bash
Rscript fitness_cost_half_life_analysis.R all       # expensive; runs two pathogens
Rscript fitness_cost_half_life_analysis.R combine   # fast; uses completed CSVs
```

For two separate cluster processes, run `ecoli` and `spneumo` as the argument, then `combine`. Default: biannual MDA for ten years, then thirty years of follow-up. The figure intentionally excludes *S. aureus* because resistance often does not halve within long follow-up. Results marked `censored = TRUE` are **lower bounds on half-life**, not observed crossing times.
