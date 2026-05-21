# =============================================================================
# One-at-a-time sensitivity analysis for S. pneumoniae MDA model
# =============================================================================
# -----------------------------------------------------------------------------
# 0. Load calibrated model
# -----------------------------------------------------------------------------

source("davies_2_8_tidy_spneumo.R")

# Put sensitivity outputs in a subfolder.
base_output_dir <- config$output_dir
config$output_dir <- file.path(base_output_dir, "sensitivity_oat_spneumo")
make_dir(config$output_dir)

message("Sensitivity outputs will be written to: ", normalizePath(config$output_dir))

# Save baseline objects so that each sensitivity case starts from the same point.
baseline_parameters_reference <- baseline_parameters
config_reference <- config

# -----------------------------------------------------------------------------
# 1. Sensitivity specifications
# -----------------------------------------------------------------------------

make_oat_specs_spneumo <- function(bp) {
  tibble::tribble(
    ~case_id, ~parameter, ~level, ~value, ~description,
    "baseline", "baseline", "baseline", NA_real_,
    "Calibrated baseline parameter set",
    "c_low", "c", "low", 0.02,
    "Low resistance fitness cost",
    "c_high", "c", "high", 0.15,
    "High resistance fitness cost",
    "k_low", "k", "low", 0.20,
    "Low co-colonisation efficiency",
    "k_high", "k", "high", 1.00,
    "High co-colonisation efficiency",
    "mda_cov_low", "mda_cov", "low", 0.70,
    "Lower MDA coverage",
    "mda_cov_high", "mda_cov", "high", 0.95,
    "Higher MDA coverage",
    "mda_effect_low", "mda_effect_multiplier", "low", 0.50,
    "Half baseline MDA treatment effect",
    "mda_effect_high", "mda_effect_multiplier", "high", 2.00,
    "Double baseline MDA treatment effect",
    "mda_duration_short", "mda_duration", "low", 14,
    "Shorter MDA campaign window",
    "mda_duration_long", "mda_duration", "high", 60,
    "Longer MDA campaign window",
    "macrolide_use_low", "macrolide_use_multiplier", "low", 0.50,
    "Half baseline background macrolide pressure",
    "macrolide_use_high", "macrolide_use_multiplier", "high", 1.50,
    "One and a half times baseline background macrolide pressure",
    "beta_low", "beta.S_start_multiplier", "low", 0.75,
    "Lower transmission coefficient",
    "beta_high", "beta.S_start_multiplier", "high", 1.25,
    "Higher transmission coefficient"
  )
}

oat_specs <- make_oat_specs_spneumo(baseline_parameters_reference)

# -----------------------------------------------------------------------------
# 2. Apply one sensitivity change to baseline_parameters
# -----------------------------------------------------------------------------

apply_spneumo_sensitivity_change <- function(bp, parameter, value) {
  if (parameter == "baseline") {
    return(bp)
  }

  if (parameter == "c") {
    bp$c <- value
    return(bp)
  }

  if (parameter == "k") {
    bp$k <- value
    return(bp)
  }

  if (parameter == "mda_cov") {
    bp$mda_cov <- value
    return(bp)
  }

  if (parameter == "mda_duration") {
    bp$mda_duration <- value
    return(bp)
  }

  if (parameter == "macrolide_use_multiplier") {
    bp$macrolide_use_ddd_per_1000_per_day <-
      baseline_parameters_reference$macrolide_use_ddd_per_1000_per_day * value
    return(bp)
  }

  if (parameter == "beta.S_start_multiplier") {
    bp$beta.S_start <- baseline_parameters_reference$beta.S_start * value
    return(bp)
  }

  if (parameter == "mda_effect_multiplier") {
    # Supports either the old formulation:
    #   a, a.C
    # or the newer realistic-coverage formulation:
    #   mda_p_clear_S, mda_p_select_C
    if ("mda_p_clear_S" %in% names(bp) && "mda_p_select_C" %in% names(bp)) {
      bp$mda_p_clear_S <- baseline_parameters_reference$mda_p_clear_S * value
      bp$mda_p_select_C <- baseline_parameters_reference$mda_p_select_C * value
    } else {
      bp$a <- baseline_parameters_reference$a * value
      bp$a.C <- baseline_parameters_reference$a.C * value
    }

    return(bp)
  }

  stop("Unknown sensitivity parameter: ", parameter)
}

# -----------------------------------------------------------------------------
# 3. Scenario set for sensitivity analysis
# -----------------------------------------------------------------------------

make_sensitivity_scenarios <- function() {
  tibble::tribble(
    ~scenario, ~horizon_years, ~frequency_per_year, ~mda_years,
    "No MDA", 10, NA_real_, NA_real_,
    "Annual MDA", 10, 1, NA_real_,
    "Biannual MDA", 10, 2, NA_real_,
    "Annual MDA stopped after 5 years", 10, 1, 5,
    "Biannual MDA stopped after 5 years", 10, 2, 5,
    "No MDA", 5.5, NA_real_, NA_real_,
    "MORDOR-style biannual MDA for 2 years", 5.5, 2, 2
  )
}

sensitivity_scenarios <- make_sensitivity_scenarios()

# -----------------------------------------------------------------------------
# 4. Helpers for summarising outcomes
# -----------------------------------------------------------------------------

summarise_age_group_endpoint <- function(time_series_by_age,
                                         age_bands,
                                         group_name) {
  time_series_by_age |>
    dplyr::filter(age_band %in% age_bands) |>
    dplyr::group_by(
      scenario,
      horizon_years,
      frequency_per_year,
      mda_years,
      time_days,
      time_years
    ) |>
    dplyr::summarise(
      total_population = sum(total_population, na.rm = TRUE),
      colonised = sum(colonised, na.rm = TRUE),
      resistant = sum(resistant, na.rm = TRUE),
      carriage_prevalence = dplyr::if_else(
        total_population > 0,
        100 * colonised / total_population,
        NA_real_
      ),
      resistant_among_carried = dplyr::if_else(
        colonised > 0,
        100 * resistant / colonised,
        NA_real_
      ),
      .groups = "drop"
    ) |>
    dplyr::group_by(horizon_years, scenario) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup() |>
    dplyr::mutate(outcome_group = group_name)
}

add_no_mda_comparison_age_group <- function(group_endpoints) {
  no_mda <- group_endpoints |>
    dplyr::filter(scenario == "No MDA") |>
    dplyr::select(
      horizon_years,
      outcome_group,
      no_mda_resistant_among_carried = resistant_among_carried,
      no_mda_carriage_prevalence = carriage_prevalence
    )

  group_endpoints |>
    dplyr::left_join(
      no_mda,
      by = c("horizon_years", "outcome_group")
    ) |>
    dplyr::mutate(
      resistance_difference_percentage_points =
        resistant_among_carried - no_mda_resistant_among_carried,
      carriage_difference_percentage_points =
        carriage_prevalence - no_mda_carriage_prevalence
    )
}

extract_case_outputs <- function(case_metadata,
                                 time_series_all,
                                 time_series_by_age) {
  endpoints_all <- endpoint_summary(time_series_all)
  comparison_all <- compare_with_no_mda(endpoints_all) |>
    dplyr::mutate(outcome_group = "all_ages")

  under5_endpoints <- summarise_age_group_endpoint(
    time_series_by_age,
    age_bands = c("0", "1-4"),
    group_name = "under5"
  ) |>
    add_no_mda_comparison_age_group()

  adult_endpoints <- summarise_age_group_endpoint(
    time_series_by_age,
    age_bands = c("18-64", "65+"),
    group_name = "adults_18plus"
  ) |>
    add_no_mda_comparison_age_group()

  comparison_age_groups <- dplyr::bind_rows(
    under5_endpoints,
    adult_endpoints
  )

  all_age_summary <- comparison_all |>
    dplyr::select(
      horizon_years,
      scenario,
      outcome_group,
      resistance_prevalence,
      colonisation_prevalence,
      resistance_difference_percentage_points
    ) |>
    dplyr::rename(
      resistant_among_carried = resistance_prevalence,
      carriage_prevalence = colonisation_prevalence
    )

  age_group_summary <- comparison_age_groups |>
    dplyr::select(
      horizon_years,
      scenario,
      outcome_group,
      resistant_among_carried,
      carriage_prevalence,
      resistance_difference_percentage_points
    )

  dplyr::bind_rows(
    all_age_summary,
    age_group_summary
  ) |>
    dplyr::mutate(
      case_id = case_metadata$case_id,
      parameter = case_metadata$parameter,
      level = case_metadata$level,
      value = case_metadata$value,
      description = case_metadata$description,
      .before = 1
    )
}

# -----------------------------------------------------------------------------
# 5. Run one sensitivity case
# -----------------------------------------------------------------------------

run_spneumo_oat_case <- function(case_row) {
  case_row <- as.list(case_row)

  message(
    "Running sensitivity case: ",
    case_row$case_id,
    " (", case_row$description, ")"
  )

  # Modify baseline_parameters for this case.
  bp_case <- apply_spneumo_sensitivity_change(
    bp = baseline_parameters_reference,
    parameter = case_row$parameter,
    value = case_row$value
  )

  # make_parameters() reads the global baseline_parameters object, so update it
  # temporarily for this case.
  baseline_parameters <<- bp_case

  parameters_case <- make_parameters(inputs, indices, ageing)

  # Recalculate the no-MDA equilibrium for this parameter set.
  equilibrium_parameters_case <- make_no_mda_parameters(parameters_case)

  equilibrium_output_case <- solve_model(
    times = make_times(config$equilibrium_years),
    state = initial_state,
    parameters = equilibrium_parameters_case
  )

  equilibrium_state_case <- as.numeric(
    equilibrium_output_case[nrow(equilibrium_output_case), -1]
  )
  names(equilibrium_state_case) <- make_state_names(inputs$age_groups)

  # Run compact scenario set.
  case_results <- purrr::pmap(
    sensitivity_scenarios,
    function(scenario, horizon_years, frequency_per_year, mda_years) {
      run_scenario(
        name = scenario,
        horizon_years = horizon_years,
        base_state = equilibrium_state_case,
        base_parameters = parameters_case,
        frequency_per_year = frequency_per_year,
        mda_years = if (is.na(mda_years)) NULL else mda_years
      )
    }
  )

  time_series_all <- purrr::map_dfr(case_results, function(result) {
    summarise_model_output(
      out = result$output,
      indices = indices,
      days_per_year = config$days_per_year
    ) |>
      dplyr::mutate(
        scenario = result$scenario,
        horizon_years = result$horizon_years,
        frequency_per_year = result$frequency_per_year,
        mda_years = result$mda_years %||% result$horizon_years,
        .before = 1
      )
  })

  time_series_by_age_case <- purrr::map_dfr(case_results, function(result) {
    summarise_by_age_band(
      out = result$output,
      indices = indices,
      age_groups = inputs$age_groups,
      days_per_year = config$days_per_year
    ) |>
      dplyr::mutate(
        scenario = result$scenario,
        horizon_years = result$horizon_years,
        frequency_per_year = result$frequency_per_year,
        mda_years = result$mda_years %||% result$horizon_years,
        .before = 1
      )
  })

  case_summary <- extract_case_outputs(
    case_metadata = case_row,
    time_series_all = time_series_all,
    time_series_by_age = time_series_by_age_case
  )

  list(
    case_summary = case_summary,
    time_series_all = time_series_all |>
      dplyr::mutate(
        case_id = case_row$case_id,
        parameter = case_row$parameter,
        level = case_row$level,
        value = case_row$value,
        .before = 1
      ),
    time_series_by_age = time_series_by_age_case |>
      dplyr::mutate(
        case_id = case_row$case_id,
        parameter = case_row$parameter,
        level = case_row$level,
        value = case_row$value,
        .before = 1
      )
  )
}

# -----------------------------------------------------------------------------
# 6. Run all OAT cases
# -----------------------------------------------------------------------------

oat_results <- purrr::map(
  seq_len(nrow(oat_specs)),
  function(i) {
    tryCatch(
      run_spneumo_oat_case(oat_specs[i, ]),
      error = function(e) {
        warning(
          "Sensitivity case failed: ",
          oat_specs$case_id[i],
          " -- ",
          conditionMessage(e)
        )

        list(
          case_summary = tibble::tibble(
            case_id = oat_specs$case_id[i],
            parameter = oat_specs$parameter[i],
            level = oat_specs$level[i],
            value = oat_specs$value[i],
            description = oat_specs$description[i],
            horizon_years = NA_real_,
            scenario = NA_character_,
            outcome_group = NA_character_,
            resistant_among_carried = NA_real_,
            carriage_prevalence = NA_real_,
            resistance_difference_percentage_points = NA_real_
          ),
          time_series_all = tibble::tibble(),
          time_series_by_age = tibble::tibble()
        )
      }
    )
  }
)

# Restore baseline parameters.
baseline_parameters <<- baseline_parameters_reference

sensitivity_summary <- purrr::map_dfr(oat_results, "case_summary")
sensitivity_time_series_all <- purrr::map_dfr(oat_results, "time_series_all")
sensitivity_time_series_by_age <- purrr::map_dfr(oat_results, "time_series_by_age")

# -----------------------------------------------------------------------------
# 7. Tornado plot helpers
# -----------------------------------------------------------------------------

make_tornado_data <- function(sensitivity_summary,
                              scenario_name = "Biannual MDA",
                              horizon = 10,
                              outcome_group_name = "all_ages") {
  primary <- sensitivity_summary |>
    dplyr::filter(
      horizon_years == horizon,
      scenario == scenario_name,
      outcome_group == outcome_group_name,
      !is.na(resistance_difference_percentage_points)
    )

  baseline_value <- primary |>
    dplyr::filter(parameter == "baseline") |>
    dplyr::slice(1) |>
    dplyr::pull(resistance_difference_percentage_points)

  primary |>
    dplyr::filter(parameter != "baseline") |>
    dplyr::group_by(parameter) |>
    dplyr::summarise(
      min_effect = min(resistance_difference_percentage_points, na.rm = TRUE),
      max_effect = max(resistance_difference_percentage_points, na.rm = TRUE),
      range = max_effect - min_effect,
      baseline_effect = baseline_value,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(abs(range)))
}

plot_oat_tornado <- function(tornado_data,
                             title = "One-at-a-time sensitivity analysis",
                             x_label = "10-year increase vs no MDA (percentage points)") {
  ggplot2::ggplot(
    tornado_data,
    ggplot2::aes(
      y = reorder(parameter, abs(range))
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = min_effect,
        xend = max_effect,
        yend = reorder(parameter, abs(range))
      ),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = min_effect),
      size = 2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = max_effect),
      size = 2
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = baseline_effect),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = "Parameter varied"
    ) +
    ggplot2::theme_classic(base_size = 12)
}

plot_oat_endpoint_points <- function(sensitivity_summary,
                                     scenario_name = "Biannual MDA",
                                     horizon = 10,
                                     outcome_group_name = "all_ages") {
  sensitivity_summary |>
    dplyr::filter(
      horizon_years == horizon,
      scenario == scenario_name,
      outcome_group == outcome_group_name,
      !is.na(resistance_difference_percentage_points)
    ) |>
    dplyr::mutate(
      case_label = dplyr::if_else(
        parameter == "baseline",
        "baseline",
        paste0(parameter, " / ", level)
      )
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = resistance_difference_percentage_points,
        y = reorder(case_label, resistance_difference_percentage_points)
      )
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::labs(
      title = paste0(
        "S. pneumoniae OAT sensitivity: ",
        scenario_name,
        ", ",
        horizon,
        " years"
      ),
      x = "Increase in resistant among carriers vs no MDA (percentage points)",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 12)
}


# Absolute-outcome tornado data.
# This differs from make_tornado_data(), which uses the MDA-minus-no-MDA
# difference as the outcome. Here the outcome is the scenario-specific
# endpoint resistance prevalence itself.
make_tornado_data_absolute <- function(sensitivity_summary,
                                       scenario_name = "Biannual MDA",
                                       horizon = 10,
                                       outcome_group_name = "all_ages",
                                       outcome_col = "resistant_among_carried") {
  primary <- sensitivity_summary |>
    dplyr::filter(
      horizon_years == horizon,
      scenario == scenario_name,
      outcome_group == outcome_group_name,
      !is.na(.data[[outcome_col]])
    )

  baseline_value <- primary |>
    dplyr::filter(parameter == "baseline") |>
    dplyr::slice(1) |>
    dplyr::pull(.data[[outcome_col]])

  primary |>
    dplyr::filter(parameter != "baseline") |>
    dplyr::group_by(parameter) |>
    dplyr::summarise(
      min_effect = min(.data[[outcome_col]], na.rm = TRUE),
      max_effect = max(.data[[outcome_col]], na.rm = TRUE),
      range = max_effect - min_effect,
      baseline_effect = baseline_value,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      scenario = scenario_name,
      horizon_years = horizon,
      outcome_group = outcome_group_name,
      outcome = outcome_col,
      .before = 1
    ) |>
    dplyr::arrange(dplyr::desc(abs(range)))
}

make_tornado_data_absolute_multi <- function(sensitivity_summary,
                                             scenario_names = c(
                                               "No MDA",
                                               "Annual MDA",
                                               "Biannual MDA"
                                             ),
                                             horizon = 10,
                                             outcome_group_name = "all_ages",
                                             outcome_col = "resistant_among_carried") {
  purrr::map_dfr(
    scenario_names,
    function(scenario_name) {
      make_tornado_data_absolute(
        sensitivity_summary = sensitivity_summary,
        scenario_name = scenario_name,
        horizon = horizon,
        outcome_group_name = outcome_group_name,
        outcome_col = outcome_col
      )
    }
  )
}

plot_oat_tornado_absolute <- function(tornado_data,
                                      title = "One-at-a-time sensitivity analysis",
                                      x_label = "Endpoint resistance prevalence (%)") {
  baseline_lines <- tornado_data |>
    dplyr::distinct(baseline_effect)

  ggplot2::ggplot(
    tornado_data,
    ggplot2::aes(
      y = reorder(parameter, abs(range))
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = min_effect,
        xend = max_effect,
        yend = reorder(parameter, abs(range))
      ),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = min_effect),
      size = 2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = max_effect),
      size = 2
    ) +
    ggplot2::geom_vline(
      data = baseline_lines,
      ggplot2::aes(xintercept = baseline_effect),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = "Parameter varied"
    ) +
    ggplot2::theme_classic(base_size = 12)
}

plot_oat_tornado_absolute_faceted <- function(tornado_data,
                                              title = "S. pneumoniae sensitivity: absolute endpoint resistance",
                                              x_label = "10-year resistant among carriers (%)") {
  baseline_lines <- tornado_data |>
    dplyr::distinct(scenario, baseline_effect)

  ggplot2::ggplot(
    tornado_data,
    ggplot2::aes(
      y = reorder(parameter, abs(range))
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = min_effect,
        xend = max_effect,
        yend = reorder(parameter, abs(range))
      ),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = min_effect),
      size = 2
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = max_effect),
      size = 2
    ) +
    ggplot2::geom_vline(
      data = baseline_lines,
      ggplot2::aes(xintercept = baseline_effect),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::facet_wrap(~scenario, ncol = 1) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = "Parameter varied"
    ) +
    ggplot2::theme_classic(base_size = 12)
}

# -----------------------------------------------------------------------------
# 8. Save outputs
# -----------------------------------------------------------------------------

readr::write_csv(
  oat_specs,
  output_path("spneumo_oat_sensitivity_specs.csv")
)

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


# Absolute endpoint resistance prevalence tornado data for each scenario.
absolute_tornado_no_mda_10y_all <- make_tornado_data_absolute(
  sensitivity_summary,
  scenario_name = "No MDA",
  horizon = 10,
  outcome_group_name = "all_ages",
  outcome_col = "resistant_among_carried"
)

absolute_tornado_annual_10y_all <- make_tornado_data_absolute(
  sensitivity_summary,
  scenario_name = "Annual MDA",
  horizon = 10,
  outcome_group_name = "all_ages",
  outcome_col = "resistant_among_carried"
)

absolute_tornado_biannual_10y_all <- make_tornado_data_absolute(
  sensitivity_summary,
  scenario_name = "Biannual MDA",
  horizon = 10,
  outcome_group_name = "all_ages",
  outcome_col = "resistant_among_carried"
)

absolute_tornado_10y_all_by_scenario <- make_tornado_data_absolute_multi(
  sensitivity_summary,
  scenario_names = c("No MDA", "Annual MDA", "Biannual MDA"),
  horizon = 10,
  outcome_group_name = "all_ages",
  outcome_col = "resistant_among_carried"
)

readr::write_csv(
  absolute_tornado_no_mda_10y_all,
  output_path("spneumo_oat_tornado_absolute_no_mda_10y_all_ages.csv")
)

readr::write_csv(
  absolute_tornado_annual_10y_all,
  output_path("spneumo_oat_tornado_absolute_annual_10y_all_ages.csv")
)

readr::write_csv(
  absolute_tornado_biannual_10y_all,
  output_path("spneumo_oat_tornado_absolute_biannual_10y_all_ages.csv")
)

readr::write_csv(
  absolute_tornado_10y_all_by_scenario,
  output_path("spneumo_oat_tornado_absolute_10y_all_ages_by_scenario.csv")
)

tornado_biannual_plot <- plot_oat_tornado(
  tornado_biannual_10y_all,
  title = "S. pneumoniae sensitivity: biannual MDA at 10 years"
)

tornado_annual_plot <- plot_oat_tornado(
  tornado_annual_10y_all,
  title = "S. pneumoniae sensitivity: annual MDA at 10 years"
)

endpoint_points_plot <- plot_oat_endpoint_points(
  sensitivity_summary,
  scenario_name = "Biannual MDA",
  horizon = 10,
  outcome_group_name = "all_ages"
)


absolute_tornado_no_mda_plot <- plot_oat_tornado_absolute(
  absolute_tornado_no_mda_10y_all,
  title = "S. pneumoniae sensitivity: no MDA at 10 years",
  x_label = "10-year resistant among carriers (%)"
)

absolute_tornado_annual_plot <- plot_oat_tornado_absolute(
  absolute_tornado_annual_10y_all,
  title = "S. pneumoniae sensitivity: annual MDA at 10 years",
  x_label = "10-year resistant among carriers (%)"
)

absolute_tornado_biannual_plot <- plot_oat_tornado_absolute(
  absolute_tornado_biannual_10y_all,
  title = "S. pneumoniae sensitivity: biannual MDA at 10 years",
  x_label = "10-year resistant among carriers (%)"
)

absolute_tornado_faceted_plot <- plot_oat_tornado_absolute_faceted(
  absolute_tornado_10y_all_by_scenario,
  title = "S. pneumoniae sensitivity: endpoint resistance by scenario",
  x_label = "10-year resistant among carriers (%)"
)

save_plot(
  tornado_biannual_plot,
  "spneumo_oat_tornado_biannual_10y_all_ages.png",
  width = 8,
  height = 5
)

save_plot(
  tornado_annual_plot,
  "spneumo_oat_tornado_annual_10y_all_ages.png",
  width = 8,
  height = 5
)

save_plot(
  endpoint_points_plot,
  "spneumo_oat_endpoint_points_biannual_10y_all_ages.png",
  width = 8,
  height = 6
)


save_plot(
  absolute_tornado_no_mda_plot,
  "spneumo_oat_tornado_absolute_no_mda_10y_all_ages.png",
  width = 8,
  height = 5
)

save_plot(
  absolute_tornado_annual_plot,
  "spneumo_oat_tornado_absolute_annual_10y_all_ages.png",
  width = 8,
  height = 5
)

save_plot(
  absolute_tornado_biannual_plot,
  "spneumo_oat_tornado_absolute_biannual_10y_all_ages.png",
  width = 8,
  height = 5
)

save_plot(
  absolute_tornado_faceted_plot,
  "spneumo_oat_tornado_absolute_10y_all_ages_by_scenario.png",
  width = 8,
  height = 10
)

saveRDS(
  list(
    oat_specs = oat_specs,
    sensitivity_summary = sensitivity_summary,
    sensitivity_time_series_all = sensitivity_time_series_all,
    sensitivity_time_series_by_age = sensitivity_time_series_by_age,
    tornado_biannual_10y_all = tornado_biannual_10y_all,
    tornado_annual_10y_all = tornado_annual_10y_all,
    absolute_tornado_no_mda_10y_all = absolute_tornado_no_mda_10y_all,
    absolute_tornado_annual_10y_all = absolute_tornado_annual_10y_all,
    absolute_tornado_biannual_10y_all = absolute_tornado_biannual_10y_all,
    absolute_tornado_10y_all_by_scenario = absolute_tornado_10y_all_by_scenario
  ),
  file = output_path("spneumo_oat_sensitivity_outputs.rds")
)

message("OAT sensitivity analysis complete.")
message("Outputs written to: ", normalizePath(config$output_dir, mustWork = FALSE))

# -----------------------------------------------------------------------------
# 9. Parameter-range summary table for report
# -----------------------------------------------------------------------------

make_parameter_range_table <- function(oat_specs,
                                       baseline_parameters_reference) {
  baseline_row <- tibble::tibble(
    parameter = c(
      "c",
      "k",
      "mda_cov",
      "mda_duration",
      "macrolide_use_ddd_per_1000_per_day",
      "beta.S_start",
      if ("mda_p_clear_S" %in% names(baseline_parameters_reference)) {
        "mda_p_clear_S"
      } else {
        "a"
      },
      if ("mda_p_select_C" %in% names(baseline_parameters_reference)) {
        "mda_p_select_C"
      } else {
        "a.C"
      }
    ),
    baseline_value = c(
      baseline_parameters_reference$c,
      baseline_parameters_reference$k,
      baseline_parameters_reference$mda_cov,
      baseline_parameters_reference$mda_duration,
      baseline_parameters_reference$macrolide_use_ddd_per_1000_per_day,
      baseline_parameters_reference$beta.S_start,
      if ("mda_p_clear_S" %in% names(baseline_parameters_reference)) {
        baseline_parameters_reference$mda_p_clear_S
      } else {
        baseline_parameters_reference$a
      },
      if ("mda_p_select_C" %in% names(baseline_parameters_reference)) {
        baseline_parameters_reference$mda_p_select_C
      } else {
        baseline_parameters_reference$a.C
      }
    )
  )

  low_high <- oat_specs |>
    dplyr::filter(parameter != "baseline") |>
    dplyr::group_by(parameter) |>
    dplyr::summarise(
      low_value = suppressWarnings(min(value, na.rm = TRUE)),
      high_value = suppressWarnings(max(value, na.rm = TRUE)),
      .groups = "drop"
    )

  baseline_row |>
    dplyr::left_join(low_high, by = "parameter") |>
    dplyr::mutate(
      parameter_description = dplyr::case_when(
        parameter == "c" ~ "Fitness cost of resistance",
        parameter == "k" ~ "Co-colonisation efficiency",
        parameter == "mda_cov" ~ "MDA coverage",
        parameter == "mda_duration" ~ "MDA campaign duration (days)",
        parameter == "macrolide_use_ddd_per_1000_per_day" ~
          "Background macrolide use",
        parameter == "beta.S_start" ~ "Transmission coefficient",
        parameter %in% c("a", "mda_p_clear_S") ~
          "MDA treatment effect on susceptible carriage",
        parameter %in% c("a.C", "mda_p_select_C") ~
          "MDA treatment effect in mixed carriage",
        TRUE ~ parameter
      )
    ) |>
    dplyr::select(
      parameter,
      parameter_description,
      baseline_value,
      low_value,
      high_value
    ) |>
    dplyr::arrange(parameter)
}

parameter_range_table <- make_parameter_range_table(
  oat_specs = oat_specs,
  baseline_parameters_reference = baseline_parameters_reference
)

# Save CSV version.
readr::write_csv(
  parameter_range_table,
  output_path("spneumo_oat_parameter_ranges.csv")
)

# Print nicely in console.
print(parameter_range_table)

# Create markdown-style table for report copy-paste.
parameter_range_table_markdown <- parameter_range_table |>
  dplyr::mutate(
    baseline_value = signif(baseline_value, 3),
    low_value = signif(low_value, 3),
    high_value = signif(high_value, 3)
  ) |>
  dplyr::rename(
    Parameter = parameter,
    Description = parameter_description,
    Baseline = baseline_value,
    Low = low_value,
    High = high_value
  )

message("\n\nSensitivity-analysis parameter ranges:\n")

print(parameter_range_table_markdown)

# Optional: save a simple markdown-ready text file.
markdown_lines <- c(
  "| Parameter | Description | Baseline | Low | High |",
  "|---|---|---:|---:|---:|"
)

markdown_rows <- apply(
  parameter_range_table_markdown,
  1,
  function(row) {
    paste0(
      "| ",
      row["Parameter"], " | ",
      row["Description"], " | ",
      row["Baseline"], " | ",
      row["Low"], " | ",
      row["High"], " |"
    )
  }
)

writeLines(
  c(markdown_lines, markdown_rows),
  con = output_path("spneumo_oat_parameter_ranges.md")
)
