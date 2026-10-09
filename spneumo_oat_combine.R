Sys.setenv(RUN_ALL_CASES = "FALSE")

source("spneumo_oat_sensitivity_recalibrated_common_clearance.R")

case_out_dir <- file.path(config$output_dir, "slurm_cases")

case_files <- list.files(
  case_out_dir,
  pattern = "^case_.*\\.rds$",
  full.names = TRUE
)

if (length(case_files) == 0) {
  stop("No case files found in: ", case_out_dir)
}

oat_results <- lapply(case_files, readRDS)

sensitivity_summary <- purrr::map_dfr(oat_results, "case_summary")
sensitivity_time_series_all <- purrr::map_dfr(oat_results, "time_series_all")
sensitivity_time_series_by_age <- purrr::map_dfr(oat_results, "time_series_by_age")

readr::write_csv(
  sensitivity_summary,
  output_path("spneumo_oat_sensitivity_summary.csv")
)

readr::write_csv(
  sensitivity_time_series_all,
  output_path("spneumo_oat_time_series_all_ages.csv")
)

readr::write_csv(
  sensitivity_time_series_by_age,
  output_path("spneumo_oat_time_series_by_age.csv")
)

# Recreate your tornado plots/tables here.
# If your original script has plotting functions already defined above,
# you can call them exactly as before.

tornado_biannual_10y_all <- make_tornado_data(
  sensitivity_summary,
  scenario_name = "Biannual MDA",
  horizon = 10,
  outcome_group_name = "all_ages"
)

tornado_annual_10y_all <- make_tornado_data(
  sensitivity_summary,
  scenario_name = "Annual MDA",
  horizon = 10,
  outcome_group_name = "all_ages"
)

readr::write_csv(
  tornado_biannual_10y_all,
  output_path("spneumo_oat_tornado_biannual_10y_all_ages.csv")
)

readr::write_csv(
  tornado_annual_10y_all,
  output_path("spneumo_oat_tornado_annual_10y_all_ages.csv")
)

saveRDS(
  list(
    oat_specs = oat_specs,
    sensitivity_summary = sensitivity_summary,
    sensitivity_time_series_all = sensitivity_time_series_all,
    sensitivity_time_series_by_age = sensitivity_time_series_by_age,
    tornado_biannual_10y_all = tornado_biannual_10y_all,
    tornado_annual_10y_all = tornado_annual_10y_all
  ),
  file = output_path("spneumo_oat_sensitivity_outputs.rds")
)

message("Combined ", length(case_files), " case files.")
message("Outputs written to: ", normalizePath(config$output_dir, mustWork = FALSE))
