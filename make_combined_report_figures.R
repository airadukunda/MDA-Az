################################################################################
# Combined report figures for azithromycin MDA AMR models
#
# This script makes:
#   1. A single 3-panel 20-year resistance time-series figure for E. coli,
#      S. pneumoniae, and S. aureus.
#   2. A single 6-panel sensitivity figure showing absolute 10-year resistance
#      prevalence under annual and biannual MDA for all three pathogens.
#
# It uses CSV outputs that are already written by the main model scripts and
# OAT sensitivity scripts. Run the pathogen main scripts and sensitivity scripts
# before running this script.
################################################################################

# -----------------------------------------------------------------------------
# 0. Packages and paths
# -----------------------------------------------------------------------------

required_packages <- c("dplyr", "readr", "ggplot2", "purrr", "tidyr", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install missing packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

# Set this to the folder containing the output folders from the three pathogen scripts.
# Usually this can stay as "." if you run from the project root.
project_dir <- "."

combined_output_dir <- file.path(project_dir, "outputs_combined_report_figures")
dir.create(combined_output_dir, recursive = TRUE, showWarnings = FALSE)

output_path <- function(filename) {
  file.path(combined_output_dir, filename)
}

pathogen_config <- tibble::tibble(
  pathogen = c("ecoli", "spneumo", "saureus"),
  pathogen_label = c("E. coli", "S. pneumoniae", "S. aureus"),
  main_output_dir = c(
    "outputs_tanzania_mda",
    "outputs_tanzania_spneumo_mda",
    "outputs_tanzania_mda_saureus"
  ),
  sensitivity_output_dir = c(
    "outputs_tanzania_mda/sensitivity_oat_ecoli",
    "outputs_tanzania_spneumo_mda/sensitivity_oat_spneumo",
    "outputs_tanzania_mda_saureus/sensitivity_oat_saureus"
  ),
  sensitivity_summary_file = c(
    "ecoli_oat_sensitivity_summary.csv",
    "spneumo_oat_sensitivity_summary.csv",
    "saureus_oat_sensitivity_summary.csv"
  )
)

# Keep scenario colours consistent with the individual pathogen plots:
# Annual = red, Biannual = green, No MDA = blue.
scenario_levels <- c("Annual MDA", "Biannual MDA", "Quarterly MDA", "No MDA")
scenario_colours <- c(
  "Annual MDA" = "#F8766D",
  "Biannual MDA" = "#00BA38",
  "No MDA" = "#619CFF",
  "Quarterly MDA" = "#C77CFF"
)

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------

check_file_exists <- function(path, hint = NULL) {
  if (!file.exists(path)) {
    msg <- paste0("File not found: ", path)
    if (!is.null(hint)) {
      msg <- paste0(msg, "\n", hint)
    }
    stop(msg, call. = FALSE)
  }
  invisible(path)
}

standardise_scenario <- function(x) {
  dplyr::case_when(
    x == "Annual MDA" ~ "Annual MDA",
    x == "Biannual MDA" ~ "Biannual MDA",
    x == "Quarterly MDA" ~ "Quarterly MDA",
    x == "No MDA" ~ "No MDA",
    TRUE ~ x
  )
}

read_pathogen_time_series <- function(pathogen, pathogen_label, main_output_dir) {
  csv_path <- file.path(project_dir, main_output_dir, "scenario_time_series.csv")
  check_file_exists(
    csv_path,
    hint = paste0(
      "Run the main model script for ", pathogen_label,
      " first, then check that scenario_time_series.csv was written."
    )
  )

  dat <- readr::read_csv(csv_path, show_col_types = FALSE)

  required_cols <- c("scenario", "horizon_years", "time_years", "resistance_prevalence")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", csv_path, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  dat |>
    dplyr::mutate(
      pathogen = pathogen,
      pathogen_label = pathogen_label,
      scenario = standardise_scenario(scenario),
      scenario = factor(scenario, levels = scenario_levels),
      resistant_percent = resistance_prevalence
    ) |>
    dplyr::filter(
      horizon_years == 20,
      scenario %in% scenario_levels
    )
}

read_sensitivity_summary <- function(pathogen,
                                     pathogen_label,
                                     sensitivity_output_dir,
                                     sensitivity_summary_file) {
  csv_path <- file.path(project_dir, sensitivity_output_dir, sensitivity_summary_file)
  check_file_exists(
    csv_path,
    hint = paste0(
      "Run the OAT sensitivity script for ", pathogen_label,
      " first, then check that ", sensitivity_summary_file, " was written."
    )
  )

  dat <- readr::read_csv(csv_path, show_col_types = FALSE)

  # The S. pneumoniae sensitivity script may call the absolute outcome
  # resistant_among_carried, whereas E. coli/S. aureus usually call it
  # resistance_prevalence. Standardise to resistance_prevalence.
  if (!"resistance_prevalence" %in% names(dat) && "resistant_among_carried" %in% names(dat)) {
    dat <- dat |>
      dplyr::mutate(resistance_prevalence = resistant_among_carried)
  }

  required_cols <- c(
    "parameter", "scenario", "horizon_years", "outcome_group",
    "resistance_prevalence"
  )
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", csv_path, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  dat |>
    dplyr::mutate(
      pathogen = pathogen,
      pathogen_label = pathogen_label,
      scenario = standardise_scenario(scenario),
      scenario = factor(scenario, levels = scenario_levels)
    )
}

make_absolute_tornado_data <- function(sensitivity_summary,
                                       scenario_name,
                                       horizon = 10,
                                       outcome_group_name = "all_ages",
                                       outcome = "resistance_prevalence") {
  primary <- sensitivity_summary |>
    dplyr::filter(
      horizon_years == horizon,
      scenario == scenario_name,
      outcome_group == outcome_group_name,
      !is.na(.data[[outcome]])
    )

  baseline_value <- primary |>
    dplyr::filter(parameter == "baseline") |>
    dplyr::slice(1) |>
    dplyr::pull(.data[[outcome]])

  if (length(baseline_value) == 0) {
    warning("No baseline value found for ", scenario_name, "; returning empty tornado data.")
    return(tibble::tibble())
  }

  primary |>
    dplyr::filter(parameter != "baseline") |>
    dplyr::group_by(pathogen, pathogen_label, scenario, parameter) |>
    dplyr::summarise(
      min_effect = min(.data[[outcome]], na.rm = TRUE),
      max_effect = max(.data[[outcome]], na.rm = TRUE),
      range = max_effect - min_effect,
      baseline_effect = baseline_value,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      scenario_label = as.character(scenario),
      parameter_label = dplyr::recode(
        parameter,
        beta.S_start_multiplier = "Transmission",
        beta.S_multiplier = "Transmission",
        macrolide_use_multiplier = "Background macrolide use",
        background_antibiotic_multiplier = "Background antibiotic use",
        mda_effect_multiplier = "MDA treatment effect",
        mda_cov = "MDA coverage",
        mda_duration = "MDA duration",
        c = "Resistance fitness cost",
        k = "Co-colonisation efficiency",
        u.S_monthly_multiplier = "Sensitive carriage clearance",
        u.R_monthly_multiplier = "Resistant carriage clearance",
        u.C_monthly_multiplier = "Mixed-carriage clearance",
        u.S_multiplier = "Sensitive carriage clearance",
        u.R_multiplier = "Resistant carriage clearance",
        u.C_multiplier = "Mixed-carriage clearance",
        clearance_multiplier = "Carriage clearance",
        background_antibiotic_use_multiplier = "Background antibiotic use",
        background_macrolide_use_multiplier = "Background macrolide use",
        beta_multiplier = "Transmission",
        .default = parameter
      )
    )
}

# -----------------------------------------------------------------------------
# 2. Combined 20-year time-series figure for all three pathogens
# -----------------------------------------------------------------------------

combined_time_series_20y <- purrr::pmap_dfr(
  pathogen_config,
  function(pathogen, pathogen_label, main_output_dir,
           sensitivity_output_dir, sensitivity_summary_file) {
    read_pathogen_time_series(
      pathogen = pathogen,
      pathogen_label = pathogen_label,
      main_output_dir = main_output_dir
    )
  }
)

# Order panels in the desired report order.
combined_time_series_20y <- combined_time_series_20y |>
  dplyr::mutate(
    pathogen_label = factor(
      pathogen_label,
      levels = c("E. coli", "S. pneumoniae", "S. aureus")
    )
  )

plot_combined_20y_time_series <- function(dat,
                                          ncol = 1,
                                          free_y = TRUE) {
  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = time_years,
      y = resistant_percent,
      colour = scenario
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(
      values = c(
        "Annual MDA" = "#F8766D",
        "Biannual MDA" = "#00BA38",
        "Quarterly MDA" = "#C77CFF",
        "No MDA" = "#619CFF"
      ),
      breaks = c("Annual MDA", "Biannual MDA", "Quarterly MDA", "No MDA")
    ) +
    ggplot2::facet_wrap(
      ~pathogen_label,
      ncol = ncol,
      scales = if (free_y) "free_y" else "fixed"
    ) +
    ggplot2::scale_colour_manual(
      values = scenario_colours,
      breaks = scenario_levels,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Macrolide-resistance prevalence over 20 years",
      subtitle = "",
      x = "Time since scenario start (years)",
      y = "Macrolide-resistant among carriers/colonised (%)",
      colour = "Scenario"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

combined_time_series_plot <- plot_combined_20y_time_series(
  combined_time_series_20y,
  ncol = 1,
  free_y = TRUE
)

combined_time_series_plot_3col <- plot_combined_20y_time_series(
  combined_time_series_20y,
  ncol = 3,
  free_y = TRUE
)

readr::write_csv(
  combined_time_series_20y,
  output_path("combined_20y_time_series_data.csv")
)

ggplot2::ggsave(
  output_path("combined_20y_resistance_time_series_three_pathogens.png"),
  combined_time_series_plot,
  width = 8,
  height = 9,
  dpi = 300
)

ggplot2::ggsave(
  output_path("combined_20y_resistance_time_series_three_pathogens_3col.png"),
  combined_time_series_plot_3col,
  width = 12,
  height = 4.5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# 3. Combined 6-panel sensitivity tornado figure
# -----------------------------------------------------------------------------

combined_sensitivity_summary <- purrr::pmap_dfr(
  pathogen_config,
  function(pathogen, pathogen_label, main_output_dir,
           sensitivity_output_dir, sensitivity_summary_file) {
    read_sensitivity_summary(
      pathogen = pathogen,
      pathogen_label = pathogen_label,
      sensitivity_output_dir = sensitivity_output_dir,
      sensitivity_summary_file = sensitivity_summary_file
    )
  }
)

absolute_tornado_combined <- purrr::map_dfr(
  c("Annual MDA", "Biannual MDA"),
  function(scenario_name) {
    combined_sensitivity_summary |>
      dplyr::group_split(pathogen, pathogen_label) |>
      purrr::map_dfr(function(dat) {
        make_absolute_tornado_data(
          sensitivity_summary = dat,
          scenario_name = scenario_name,
          horizon = 10,
          outcome_group_name = "all_ages",
          outcome = "resistance_prevalence"
        )
      })
  }
) |>
  dplyr::mutate(
    pathogen_label = factor(
      pathogen_label,
      levels = c("E. coli", "S. pneumoniae", "S. aureus")
    ),
    scenario_label = factor(
      scenario_label,
      levels = c("Annual MDA", "Biannual MDA")
    )
  )

parameters_to_exclude <- c(
  "Resistant carriage clearance",
  "Carriage clearance",
  "theta"
)

absolute_tornado_combined <- absolute_tornado_combined |>
  dplyr::filter(!parameter_label %in% parameters_to_exclude)

# Keep a consistent y-axis order within each panel by ordering parameters by the
# average range across all six panels. This makes the figure easier to scan.
parameter_order <- absolute_tornado_combined |>
  dplyr::group_by(parameter_label) |>
  dplyr::summarise(mean_range = mean(abs(range), na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(mean_range) |>
  dplyr::pull(parameter_label)

absolute_tornado_combined <- absolute_tornado_combined |>
  dplyr::mutate(
    parameter_label = factor(parameter_label, levels = parameter_order)
  )

plot_combined_absolute_tornado <- function(tornado_data) {
  ggplot2::ggplot(
    tornado_data,
    ggplot2::aes(y = parameter_label)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = min_effect,
        xend = max_effect,
        yend = parameter_label
      ),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = min_effect),
      size = 1.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = max_effect),
      size = 1.8
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = baseline_effect),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::facet_grid(
      pathogen_label ~ scenario_label,
      scales = "free",
      space = "free_y"
    ) +
    ggplot2::labs(
      title = "One-at-a-time sensitivity analysis on resistance prevalence after 10 years",
      subtitle = "",
      x = "Macrolide-resistant among carriers/colonised (%)",
      y = "Parameter varied"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(1.0, "lines")
    )
}

combined_sensitivity_tornado_plot <- plot_combined_absolute_tornado(
  absolute_tornado_combined
)

readr::write_csv(
  combined_sensitivity_summary,
  output_path("combined_oat_sensitivity_summary.csv")
)

readr::write_csv(
  absolute_tornado_combined,
  output_path("combined_oat_absolute_tornado_annual_biannual_10y.csv")
)

ggplot2::ggsave(
  output_path("combined_oat_absolute_tornado_annual_biannual_10y.png"),
  combined_sensitivity_tornado_plot,
  width = 11,
  height = 8.5,
  dpi = 300
)

message("Combined report figures written to: ", normalizePath(combined_output_dir, mustWork = FALSE))
