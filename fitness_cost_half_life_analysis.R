# =============================================================================
# Fitness cost vs half-life of excess macrolide resistance after MDA cessation
# =============================================================================
#
# Purpose:
#   For each pathogen, vary the resistance fitness cost c, recalibrate baseline
#   transmission and background antibiotic pressure, run a 10-year MDA programme
#   followed by no MDA, and estimate the time after stopping for excess
#   resistance to fall by half.
#
# Definition of half-life used here:
#   Let R0    = no-MDA baseline resistance among carriers/colonised.
#       Rstop = resistance at the time MDA stops.
#   The excess at stopping is Rstop - R0.
#   The half-life is the first time after stopping when resistance falls to:
#       R0 + 0.5 * (Rstop - R0)
#
# Usage:
#   Rscript fitness_cost_half_life_analysis.R all
#   Rscript fitness_cost_half_life_analysis.R ecoli
#   Rscript fitness_cost_half_life_analysis.R spneumo
#   Rscript fitness_cost_half_life_analysis.R saureus
#   Rscript fitness_cost_half_life_analysis.R combine
#
# Requirements:
#   Run from the directory containing these scripts:
#     ecoli_oat_sensitivity_recalibrated_common_clearance.R
#     spneumo_oat_sensitivity_recalibrated_common_clearance.R
#     saureus_oat_sensitivity_recalibrated_common_clearance.R
#   and the corresponding main model scripts sourced by them.
#
#   The recalibrated sensitivity scripts must respect RUN_ALL_CASES=FALSE.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

get_this_script <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 0) {
    return("fitness_cost_half_life_analysis.R")
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
}

mode <- commandArgs(trailingOnly = TRUE)[1] %||% "all"

output_dir <- Sys.getenv(
  "HALF_LIFE_OUTPUT_DIR",
  unset = "outputs_fitness_cost_half_life"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Scenario defaults. Edit these if you want annual MDA or longer follow-up.
mda_frequency_per_year <- as.numeric(Sys.getenv("HALF_LIFE_MDA_FREQ", unset = "2"))
mda_years <- as.numeric(Sys.getenv("HALF_LIFE_MDA_YEARS", unset = "10"))
horizon_years <- as.numeric(Sys.getenv("HALF_LIFE_HORIZON_YEARS", unset = "40"))

if (horizon_years <= mda_years) {
  stop("horizon_years must be greater than mda_years.")
}

pathogen_specs <- tibble::tribble(
  ~pathogen,  ~pathogen_label,     ~sensitivity_script,                                                   ~cost_grid,
  "ecoli",    "E. coli",          "ecoli_oat_sensitivity_recalibrated_common_clearance.R",              list(seq(0.00, 0.30, by = 0.025)),
  "spneumo",  "S. pneumoniae",    "spneumo_oat_sensitivity_recalibrated_common_clearance.R",            list(seq(0.00, 0.20, by = 0.020)),
  "saureus",  "S. aureus",        "saureus_oat_sensitivity_recalibrated_common_clearance.R",            list(seq(0.00, 0.25, by = 0.025))
)

extract_resistance_column <- function(summary_df) {
  if ("resistance_prevalence" %in% names(summary_df)) {
    return("resistance_prevalence")
  }
  if ("resistant_among_carried" %in% names(summary_df)) {
    return("resistant_among_carried")
  }
  stop("Could not identify resistance column in model summary.")
}

interp_at_time <- function(df, time_col, value_col, t) {
  df <- df |>
    dplyr::arrange(.data[[time_col]]) |>
    dplyr::filter(is.finite(.data[[time_col]]), is.finite(.data[[value_col]]))

  if (nrow(df) == 0) return(NA_real_)

  if (t <= min(df[[time_col]])) return(df[[value_col]][which.min(df[[time_col]])])
  if (t >= max(df[[time_col]])) return(df[[value_col]][which.max(df[[time_col]])])

  stats::approx(
    x = df[[time_col]],
    y = df[[value_col]],
    xout = t,
    rule = 2
  )$y
}

find_resistance_half_life <- function(time_series,
                                      stop_year,
                                      baseline_resistance,
                                      time_col = "time_years",
                                      resistance_col = "resistance_prevalence") {
  ts <- time_series |>
    dplyr::arrange(.data[[time_col]]) |>
    dplyr::filter(is.finite(.data[[time_col]]), is.finite(.data[[resistance_col]]))

  r_stop <- interp_at_time(ts, time_col, resistance_col, stop_year)

  if (!is.finite(r_stop) || !is.finite(baseline_resistance)) {
    return(tibble(
      resistance_at_stop = r_stop,
      baseline_resistance = baseline_resistance,
      excess_at_stop = NA_real_,
      half_resistance_threshold = NA_real_,
      half_life_years = NA_real_,
      half_life_absolute_time_years = NA_real_,
      censored = TRUE
    ))
  }

  excess <- r_stop - baseline_resistance

  # If there is no excess resistance at stopping, define half-life as 0.
  if (excess <= 0) {
    return(tibble(
      resistance_at_stop = r_stop,
      baseline_resistance = baseline_resistance,
      excess_at_stop = excess,
      half_resistance_threshold = baseline_resistance,
      half_life_years = 0,
      half_life_absolute_time_years = stop_year,
      censored = FALSE
    ))
  }

  threshold <- baseline_resistance + 0.5 * excess

  post <- ts |>
    dplyr::filter(.data[[time_col]] >= stop_year)

  if (nrow(post) == 0) {
    return(tibble(
      resistance_at_stop = r_stop,
      baseline_resistance = baseline_resistance,
      excess_at_stop = excess,
      half_resistance_threshold = threshold,
      half_life_years = NA_real_,
      half_life_absolute_time_years = NA_real_,
      censored = TRUE
    ))
  }

  # Add exact stopping-time point if it is not already in the output.
  if (!any(abs(post[[time_col]] - stop_year) < 1e-8)) {
    post <- dplyr::bind_rows(
      tibble::tibble(
        !!time_col := stop_year,
        !!resistance_col := r_stop
      ),
      post |>
        dplyr::select(dplyr::all_of(c(time_col, resistance_col)))
    ) |>
      dplyr::arrange(.data[[time_col]])
  } else {
    post <- post |>
      dplyr::select(dplyr::all_of(c(time_col, resistance_col)))
  }

  crossed <- which(post[[resistance_col]] <= threshold)
  crossed <- crossed[post[[time_col]][crossed] >= stop_year]

  if (length(crossed) == 0) {
    max_followup <- max(post[[time_col]], na.rm = TRUE) - stop_year
    return(tibble(
      resistance_at_stop = r_stop,
      baseline_resistance = baseline_resistance,
      excess_at_stop = excess,
      half_resistance_threshold = threshold,
      half_life_years = max_followup,
      half_life_absolute_time_years = max(post[[time_col]], na.rm = TRUE),
      censored = TRUE
    ))
  }

  i <- crossed[[1]]

  if (i == 1) {
    t_cross <- post[[time_col]][i]
  } else {
    t0 <- post[[time_col]][i - 1]
    t1 <- post[[time_col]][i]
    y0 <- post[[resistance_col]][i - 1]
    y1 <- post[[resistance_col]][i]

    if (isTRUE(all.equal(y0, y1))) {
      t_cross <- t1
    } else {
      frac <- (threshold - y0) / (y1 - y0)
      frac <- max(0, min(1, frac))
      t_cross <- t0 + frac * (t1 - t0)
    }
  }

  tibble(
    resistance_at_stop = r_stop,
    baseline_resistance = baseline_resistance,
    excess_at_stop = excess,
    half_resistance_threshold = threshold,
    half_life_years = t_cross - stop_year,
    half_life_absolute_time_years = t_cross,
    censored = FALSE
  )
}

run_pathogen <- function(pathogen_name) {
  spec <- pathogen_specs |>
    dplyr::filter(pathogen == pathogen_name) |>
    dplyr::slice(1)

  if (nrow(spec) != 1) {
    stop("Unknown pathogen: ", pathogen_name)
  }

  message("\n=== Running half-life analysis for ", spec$pathogen_label, " ===")

  if (!file.exists(spec$sensitivity_script)) {
    stop("Cannot find sensitivity script: ", spec$sensitivity_script)
  }

  # The sourced sensitivity script should define:
  #   baseline_parameters_reference
  #   calibrate_baseline_for_case()
  #   run_scenario()
  #   summarise_model_output()
  #   indices, config
  Sys.setenv(RUN_ALL_CASES = "FALSE")
  source(spec$sensitivity_script)

  cost_values <- sort(unique(c(unlist(spec$cost_grid), baseline_parameters_reference$c)))

  half_life_results <- list()
  time_series_results <- list()

  for (i in seq_along(cost_values)) {
    c_value <- cost_values[[i]]

    message(
      "  ", spec$pathogen_label,
      ": fitness cost c = ", signif(c_value, 4),
      " (", i, "/", length(cost_values), ")"
    )

    bp_case <- baseline_parameters_reference
    bp_case$c <- c_value

    calibrated_case <- calibrate_baseline_for_case(bp_case)

    scenario_result <- run_scenario(
      name = "Biannual MDA for 10 years then stop",
      horizon_years = horizon_years,
      base_state = calibrated_case$equilibrium_state,
      base_parameters = calibrated_case$parameters,
      frequency_per_year = mda_frequency_per_year,
      mda_years = mda_years
    )

    summary <- summarise_model_output(
      out = scenario_result$output,
      indices = indices,
      days_per_year = config$days_per_year
    )

    resistance_col <- extract_resistance_column(summary)

    ts <- summary |>
      dplyr::mutate(
        pathogen = spec$pathogen,
        pathogen_label = spec$pathogen_label,
        fitness_cost = c_value,
        mda_frequency_per_year = mda_frequency_per_year,
        mda_years = mda_years,
        horizon_years = horizon_years,
        resistance = .data[[resistance_col]]
      )

    baseline_resistance <- calibrated_case$calibration_diagnostics$achieved_resistance

    half_life <- find_resistance_half_life(
      time_series = ts,
      stop_year = mda_years,
      baseline_resistance = baseline_resistance,
      time_col = "time_years",
      resistance_col = "resistance"
    ) |>
      dplyr::mutate(
        pathogen = spec$pathogen,
        pathogen_label = spec$pathogen_label,
        fitness_cost = c_value,
        mda_frequency_per_year = mda_frequency_per_year,
        mda_years = mda_years,
        horizon_years = horizon_years,
        calibrated_beta = calibrated_case$calibration_diagnostics$calibrated_beta,
        calibrated_background_pressure = calibrated_case$calibration_diagnostics$calibrated_background_pressure,
        achieved_baseline_carriage = calibrated_case$calibration_diagnostics$achieved_carriage,
        achieved_baseline_resistance = calibrated_case$calibration_diagnostics$achieved_resistance,
        baseline_carriage_error_pp = calibrated_case$calibration_diagnostics$carriage_error_pp,
        baseline_resistance_error_pp = calibrated_case$calibration_diagnostics$resistance_error_pp,
        .before = 1
      )

    half_life_results[[i]] <- half_life
    time_series_results[[i]] <- ts
  }

  pathogen_half_life <- dplyr::bind_rows(half_life_results)
  pathogen_time_series <- dplyr::bind_rows(time_series_results)

  readr::write_csv(
    pathogen_half_life,
    file.path(output_dir, paste0(spec$pathogen, "_fitness_cost_half_life.csv"))
  )

  readr::write_csv(
    pathogen_time_series,
    file.path(output_dir, paste0(spec$pathogen, "_fitness_cost_half_life_time_series.csv"))
  )

  message("Saved results for ", spec$pathogen_label)
}

combine_outputs <- function() {
  half_life_files <- file.path(
    output_dir,
    paste0(pathogen_specs$pathogen, "_fitness_cost_half_life.csv")
  )

  missing_files <- half_life_files[!file.exists(half_life_files)]
  if (length(missing_files) > 0) {
    stop("Missing half-life files:\n", paste(missing_files, collapse = "\n"))
  }

  half_life <- purrr::map_dfr(half_life_files, readr::read_csv, show_col_types = FALSE) |>
    dplyr::mutate(
      pathogen_label = factor(
        pathogen_label,
        levels = c("E. coli", "S. pneumoniae", "S. aureus")
      ),
      half_life_plot = half_life_years,
      half_life_label = dplyr::if_else(
        censored,
        paste0(">", signif(half_life_years, 3)),
        as.character(signif(half_life_years, 3))
      )
    )

  readr::write_csv(
    half_life,
    file.path(output_dir, "combined_fitness_cost_half_life.csv")
  )

  p <- ggplot2::ggplot(
    half_life,
    ggplot2::aes(
      x = fitness_cost,
      y = half_life_plot,
      colour = pathogen_label,
      group = pathogen_label
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(shape = censored), size = 2.4, na.rm = TRUE) +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 17),
      labels = c(`FALSE` = "Reached half level", `TRUE` = "Not reached by end of follow-up"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Persistence of excess macrolide resistance after MDA stops",
      x = "Fitness cost of resistance, c",
      y = "Half-life of excess resistance after stopping MDA (years)",
      colour = "Pathogen"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(
    filename = file.path(output_dir, "fitness_cost_half_life_three_pathogens.png"),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )

  p_facet <- ggplot2::ggplot(
    half_life,
    ggplot2::aes(
      x = fitness_cost,
      y = half_life_plot
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(shape = censored), size = 2.2, na.rm = TRUE) +
    ggplot2::facet_wrap(~ pathogen_label, scales = "free_x") +
    ggplot2::scale_shape_manual(
      values = c(`FALSE` = 16, `TRUE` = 17),
      labels = c(`FALSE` = "Reached half level", `TRUE` = "Not reached by end of follow-up"),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Persistence of excess macrolide resistance after MDA stops",
      x = "Fitness cost of resistance, c",
      y = "Half-life of excess resistance after stopping MDA (years)"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(
    filename = file.path(output_dir, "fitness_cost_half_life_three_pathogens_faceted.png"),
    plot = p_facet,
    width = 9,
    height = 4.5,
    dpi = 300
  )

  message("Combined outputs written to: ", normalizePath(output_dir, mustWork = FALSE))
}

if (mode == "all") {
  # Run each pathogen in a separate R process to avoid object/function name clashes
  # between the model scripts.
  script <- get_this_script()
  rscript <- file.path(R.home("bin"), "Rscript")

  for (p in pathogen_specs$pathogen) {
    message("Launching worker for ", p)
    status <- system2(rscript, args = c(script, p))
    if (!identical(status, 0L)) {
      stop("Worker failed for pathogen: ", p)
    }
  }

  status <- system2(rscript, args = c(script, "combine"))
  if (!identical(status, 0L)) {
    stop("Combine step failed.")
  }
} else if (mode %in% pathogen_specs$pathogen) {
  run_pathogen(mode)
} else if (mode == "combine") {
  combine_outputs()
} else {
  stop("Unknown mode: ", mode, ". Use all, combine, ecoli, spneumo, or saureus.")
}
