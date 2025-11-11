#!/usr/bin/env Rscript

# Advanced Bacterial Dynamics Model — Persistence + Proportion Visuals (fixed)
# ---------------------------------------------------------------------------

# Dependencies -------------------------------------------------------------####
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
suppressWarnings(suppressMessages({
  tryCatch({
    pacman::p_load(deSolve, viridis, tidyverse, data.table, patchwork)
  }, error = function(e) {
    pkgs <- c("deSolve","viridis","tidyverse","data.table","patchwork",
              "dplyr","ggplot2","readr","tibble","stringr")
    for (pkg in pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        tryCatch({
          install.packages(pkg, repos = "http://cran.us.r-project.org")
          library(pkg, character.only = TRUE)
        }, silent = TRUE)
      } else library(pkg, character.only = TRUE)
    }
  })
}))

# Config -------------------------------------------------------------------####
CONFIG <- list(
  COUNTRY = "Afghanistan",
  YEARS = 2000:2023,
  OUTPUT_DIR = file.path(getwd(), "outputs"),
  DEBUG_MODE = TRUE,
  SIMULATION_DURATION = 365.25 * 20,  # 20y
  RANDOM_SEED = 42
)
dir.create(CONFIG$OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Logging ------------------------------------------------------------------####
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s: %s\n", level, timestamp, message))
}

# Data Loading -------------------------------------------------------------####
load_demographic_data <- function(country = CONFIG$COUNTRY,
                                  target_year = max(CONFIG$YEARS)) {
  log_message(sprintf("Loading data for %s in year %d", country, target_year))
  
  population_path <- file.path(getwd(), "Population_emro_2023_1yearage.csv")
  birth_path      <- file.path(getwd(), "3.U.1.Birth_1year_emro.csv")
  contact_path    <- file.path(getwd(), "3.U.1.contact_Pakistan_1y.csv")
  
  tryCatch({
    population_data <- readr::read_csv(population_path, show_col_types = FALSE) %>%
      dplyr::filter(Country == !!country, Year == target_year) %>%
      dplyr::mutate(
        Population_age    = Population_age * 1000,
        Annual_population = Annual_population * 1000
      ) %>%
      dplyr::arrange(Age_Category)
    
    n_age <- nrow(population_data)
    if (n_age == 0) stop("No population rows after filtering—check Country/Year.")
    
    birth_data <- readr::read_csv(birth_path, show_col_types = FALSE) %>%
      dplyr::filter(Country == !!country, Year == target_year) %>%
      dplyr::arrange(Age)
    
    births_vector <- suppressWarnings(as.numeric(birth_data$Birth))
    if (length(births_vector) != n_age) {
      stop(sprintf("Birth vector length (%d) != n_age (%d).", length(births_vector), n_age))
    }
    births_vector[is.na(births_vector)] <- 0
    
    contact_matrix_raw <- utils::read.csv(contact_path, header = FALSE, check.names = FALSE)
    contact_matrix_char <- as.matrix(contact_matrix_raw[-1, -1, drop = FALSE])
    contact_matrix_num <- matrix(
      as.numeric(gsub(",", "", contact_matrix_char)),
      nrow = nrow(contact_matrix_char),
      ncol = ncol(contact_matrix_char)
    )
    if (anyNA(contact_matrix_num)) {
      bad <- which(is.na(contact_matrix_num), arr.ind = TRUE)
      stop(sprintf("Contact matrix has non-numeric cells, e.g. [%s]",
                   paste(utils::head(apply(bad, 1, paste, collapse=","), 6), collapse="; ")))
    }
    contact_matrix <- contact_matrix_num
    storage.mode(contact_matrix) <- "double"
    
    if (!all(dim(contact_matrix) == c(n_age, n_age))) {
      if (nrow(contact_matrix) >= n_age && ncol(contact_matrix) >= n_age) {
        contact_matrix <- contact_matrix[1:n_age, 1:n_age, drop = FALSE]
      } else {
        pad_val <- mean(contact_matrix)
        pad_mat <- matrix(pad_val, nrow = n_age, ncol = n_age)
        pad_mat[seq_len(nrow(contact_matrix)), seq_len(ncol(contact_matrix))] <- contact_matrix
        contact_matrix <- pad_mat
      }
    }
    
    log_message("Demographic data loaded successfully")
    list(population = population_data,
         births_vector = births_vector,
         contact_matrix = contact_matrix,
         n_age = n_age)
  }, error = function(e) {
    log_message(sprintf("Error loading demographic data: %s", e$message), "ERROR")
    NULL
  })
}
#------------------------------------------------------------------------####
safe_write_csv <- function(df, path, retries = 6, sleep_sec = 0.7, fallback_timestamp = TRUE) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  for (i in seq_len(retries)) {
    ok <- tryCatch({
      tmp <- tempfile(pattern = "tmp_csv_", fileext = ".csv")
      readr::write_csv(df, tmp)
      if (file.exists(path)) suppressWarnings(file.remove(path))
      file.rename(tmp, path)
    }, error = function(e) FALSE)

    if (ok) {
      log_message(paste("Wrote:", normalizePath(path, mustWork = FALSE)))
      return(invisible(TRUE))
    }
    Sys.sleep(sleep_sec)
  }

  if (fallback_timestamp) {
    alt <- file.path(
      dirname(path),
      paste0(tools::file_path_sans_ext(basename(path)), "_",
             format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv")
    )
    readr::write_csv(df, alt)
    log_message(paste("Target locked; wrote to:", normalizePath(alt, mustWork = FALSE)), "WARNING")
    return(invisible(FALSE))
  } else {
    stop(sprintf("Failed to write after %d retries: %s", retries, path))
  }
}


# Parameters ---------------------------------------------------------------####
prepare_model_parameters <- function(demographic_data) {
  set.seed(CONFIG$RANDOM_SEED)
  list(
    beta_sensitive       = 10 * 12 / 365.25,
    beta_resistant       =  8 * 12 / 365.25,
    clearance_sensitive  =  0.2 * 12 / 365.25,
    clearance_resistant  =  0.15 * 12 / 365.25,
    co_colonization_rate = 0.5,
    resistance_cost      = 0.1,
    base_mortality       = 0.001
  )
}

# ODEs --------------------------------------------------------------------####
create_bacterial_odes <- function(n_age, contact_matrix, births_vector) {
  function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      X  <- state[1:n_age]
      S  <- state[(n_age + 1):(2*n_age)]
      R  <- state[(2*n_age + 1):(3*n_age)]
      Sr <- state[(3*n_age + 1):(4*n_age)]
      Rs <- state[(4*n_age + 1):(5*n_age)]
      
      total_population <- X + S + R + Sr + Rs
      denom <- pmax(total_population, 1e-10)
      
      vS <- (S + Sr) / denom
      vR <- (R + Rs) / denom
      
      contact_weighted_S <- contact_matrix %*% vS
      contact_weighted_R <- contact_matrix %*% vR
      
      lambda_S <- as.numeric(beta_sensitive * contact_weighted_S)
      lambda_R <- as.numeric(beta_resistant * contact_weighted_R)
      
      dX  <- births_vector - lambda_S * X - lambda_R * X
      dS  <- lambda_S * X - clearance_sensitive * S - co_colonization_rate * lambda_R * S
      dR  <- lambda_R * X - clearance_resistant * R - co_colonization_rate * lambda_S * R
      dSr <- co_colonization_rate * lambda_R * S - clearance_sensitive * Sr
      dRs <- co_colonization_rate * lambda_S * R - clearance_resistant * Rs
      
      list(c(dX, dS, dR, dSr, dRs))
    })
  }
}

# Run model ---------------------------------------------------------------####
run_bacterial_dynamics_model <- function() {
  demographic_data <- load_demographic_data()
  if (is.null(demographic_data)) stop("Failed to load demographic data")
  
  parameters <- prepare_model_parameters(demographic_data)
  n_age <- demographic_data$n_age
  
  initial_state <- c(
    X  = demographic_data$population$Population_age * 0.95,
    S  = demographic_data$population$Population_age * 0.025,
    R  = demographic_data$population$Population_age * 0.025,
    Sr = rep(0, n_age),
    Rs = rep(0, n_age)
  )
  
  stopifnot(
    is.numeric(demographic_data$contact_matrix),
    all(dim(demographic_data$contact_matrix) == c(n_age, n_age)),
    is.numeric(demographic_data$births_vector),
    length(demographic_data$births_vector) == n_age,
    is.numeric(initial_state),
    length(initial_state) == 5 * n_age
  )
  
  time_vector <- seq(0, CONFIG$SIMULATION_DURATION, by = 1)
  
  tryCatch({
    ode_result <- deSolve::ode(
      y     = initial_state,
      times = time_vector,
      func  = create_bacterial_odes(n_age,
                                    demographic_data$contact_matrix, demographic_data$births_vector),
      parms = parameters
    )
    list(results = as.data.frame(ode_result),
         parameters = parameters,
         demographic_data = demographic_data)
  }, error = function(e) {
    log_message(sprintf("Error in ODE solving: %s", e$message), "ERROR"); NULL
  })
}

# ---- Tidy & summaries ---------------------------------------------------####
tidy_results <- function(model_results) {
  model_results$results %>%
    tidyr::pivot_longer(cols = -time, names_to = "compartment", values_to = "value") %>%
    dplyr::mutate(
      group = stringr::str_extract(compartment, "^[A-Za-z]+"),
      age   = suppressWarnings(as.numeric(stringr::str_extract(compartment, "\\d+")))
    ) %>% dplyr::filter(!is.na(age))
}

add_total_and_props <- function(results_long) {
  totals <- results_long %>% dplyr::group_by(time) %>%
    dplyr::summarise(total = sum(value), .groups = "drop")
  
  by_group <- results_long %>%
    dplyr::group_by(time, group) %>%
    dplyr::summarise(total_group = sum(value), .groups = "drop") %>%
    dplyr::left_join(totals, by = "time") %>%
    dplyr::mutate(prop_group = dplyr::if_else(total > 0, total_group/total, 0))
  
  prev <- results_long %>%
    dplyr::mutate(strain = dplyr::case_when(
      group %in% c("S","Sr") ~ "Sensitive",
      group %in% c("R","Rs") ~ "Resistant",
      TRUE                   ~ "Uncolonized"
    )) %>%
    dplyr::group_by(time, strain) %>%
    dplyr::summarise(total_strain = sum(value), .groups = "drop") %>%
    dplyr::left_join(totals, by = "time") %>%
    dplyr::mutate(prop_strain = dplyr::if_else(total > 0, total_strain/total, 0))
  
  list(totals = totals, by_group = by_group, prev = prev)
}

# --------- Plot helpers --------------------------------------------------####
theme_clean <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold"))
}

plot_total_counts <- function(results_long) {
  results_long %>% dplyr::group_by(time, group) %>%
    dplyr::summarise(total_value = sum(value), .groups = "drop") %>%
    ggplot2::ggplot(ggplot2::aes(time, total_value, color = group)) +
    ggplot2::geom_line(linewidth = 1.05, lineend = "round") +
    ggplot2::labs(title = "Total Bacterial Population Dynamics",
                  x = "Time (Days)", y = "Total Population") +
    viridis::scale_color_viridis(discrete = TRUE, end = 0.9) +
    theme_clean()
}

plot_total_counts_zoom <- function(results_long, max_days = 1000) {
  plot_total_counts(results_long) +
    ggplot2::coord_cartesian(xlim = c(0, max_days)) +
    ggplot2::labs(title = sprintf("Total Dynamics (First %d Days)", max_days))
}

plot_strain_prevalence <- function(prev_df) {
  ggplot2::ggplot(prev_df, ggplot2::aes(time, prop_strain, fill = strain)) +
    ggplot2::geom_area(alpha = 0.95) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    viridis::scale_fill_viridis(discrete = TRUE, end = 0.9) +
    ggplot2::labs(title = "Strain Prevalence Over Time (Proportions)",
                  x = "Time (Days)", y = "Population Share") +
    theme_clean()
}

plot_group_proportions <- function(by_group_df) {
  ggplot2::ggplot(by_group_df, ggplot2::aes(time, prop_group, color = group)) +
    ggplot2::geom_line(linewidth = 1.05) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    viridis::scale_color_viridis(discrete = TRUE, end = 0.9) +
    ggplot2::labs(title = "Compartment Proportions Over Time",
                  x = "Time (Days)", y = "Share of Total Population") +
    theme_clean()
}

# age prevalence at times --------------------------------------------------####
plot_age_profile_at_times <- function(results_long,
                                      times = c(365, 5*365, 10*365, 20*365)) {
  # Build strain-wise totals per (time, age)
  df <- results_long %>%
    dplyr::filter(group %in% c("X","S","R","Sr","Rs")) %>%
    dplyr::mutate(strain = dplyr::case_when(
      group %in% c("S","Sr") ~ "Sensitive",
      group %in% c("R","Rs") ~ "Resistant",
      TRUE                   ~ "Uncolonized"
    )) %>%
    dplyr::group_by(time, age, strain) %>%
    dplyr::summarise(value = sum(value), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = strain, values_from = value, values_fill = 0) %>%
    dplyr::mutate(
      N       = `Sensitive` + `Resistant` + `Uncolonized`,
      prop_S  = dplyr::if_else(N > 0, `Sensitive` / N, 0),
      prop_R  = dplyr::if_else(N > 0, `Resistant` / N, 0)
    )
  
  # nearest simulated times to requested "times"
  tvals <- sort(unique(df$time))
  nearest <- sapply(times, function(tt) tvals[which.min(abs(tvals - tt))])
  lab_map <- tibble::tibble(
    time = as.double(nearest),
    time_label = paste0("t≈", round(nearest / 365.25, 1), " y")
  )
  
  df_plot <- df %>%
    dplyr::filter(time %in% nearest) %>%
    dplyr::select(time, age, prop_S, prop_R) %>%
    tidyr::pivot_longer(cols = c(prop_S, prop_R),
                        names_to = "strain", values_to = "prop") %>%
    dplyr::mutate(strain = dplyr::recode(strain, prop_S = "Sensitive", prop_R = "Resistant")) %>%
    dplyr::left_join(lab_map, by = "time")
  
  ggplot2::ggplot(df_plot, ggplot2::aes(age, prop, color = strain, fill = strain)) +
    ggplot2::geom_area(position = "identity", alpha = 0.55) +
    ggplot2::facet_wrap(~ time_label, ncol = 2, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    viridis::scale_color_viridis(discrete = TRUE, end = 0.9) +
    viridis::scale_fill_viridis(discrete = TRUE, end = 0.9) +
    ggplot2::labs(title = "Age-Specific Prevalence at Selected Times",
                  x = "Age (Index)", y = "Prevalence (Proportion)") +
    theme_clean()
}

# Save wrapper ------------------------------------------------------------####
#save_plot <- function(plot, filename, width, height, dpi = 150) {
#  path <- file.path(CONFIG$OUTPUT_DIR, filename)
 # ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
 # log_message(paste("Saved:", normalizePath(path, mustWork = FALSE)))
#}

# Save wrapper ------------------------------------------------------------####
save_plot <- function(plot, filename, width, height, dpi = 150) {
  path <- file.path(CONFIG$OUTPUT_DIR, filename)
  
  # ✅ Ensure white background in the saved PNG
  plot <- plot + ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA)
  )
  
  ggplot2::ggsave(
    path,
    plot   = plot,
    width  = width,
    height = height,
    dpi    = dpi,
    bg     = "white"  # ensures PNG background itself is white (not transparent)
  )
  
  log_message(paste("Saved:", normalizePath(path, mustWork = FALSE)))
}


# Plot pipeline -----------------------------------------------------------####

plot_bacterial_dynamics <- function(model_results) {
  if (is.null(model_results)) { log_message("No model results to plot", "WARNING"); return(NULL) }
  
  results_long <- tidy_results(model_results)
  summaries <- add_total_and_props(results_long)
  
  p_counts_full <- plot_total_counts(results_long)
  p_counts_zoom <- plot_total_counts_zoom(results_long, max_days = 1000)
  p_prev        <- plot_strain_prevalence(summaries$prev)
  p_props       <- plot_group_proportions(summaries$by_group)
  p_age_times   <- plot_age_profile_at_times(results_long)
  
  save_plot(p_counts_full, "01_total_counts_full.png", 10, 6)
  save_plot(p_counts_zoom, "02_total_counts_first1000d.png", 10, 6)
  save_plot(p_prev,        "03_strain_prevalence_proportions.png", 10, 6)
  save_plot(p_props,       "04_compartment_proportions.png", 10, 6)
  save_plot(p_age_times,   "05_age_prevalence_selected_times.png", 12, 8)
  
  (p_counts_zoom / p_prev) | (p_props / p_age_times)
}

# Main ---------------------------------------------------------------------####
main <- function() {
  tryCatch({
    model_results <- run_bacterial_dynamics_model()
    if (!is.null(model_results)) {
      plot_bacterial_dynamics(model_results)
      
    #  readr::write_csv(model_results$results,
     #                  file.path(CONFIG$OUTPUT_DIR, "bacterial_dynamics_results.csv"))
      
     # params_tbl <- tibble::enframe(model_results$parameters, name = "parameter", value = "value")
     # readr::write_csv(params_tbl,
      #                 file.path(CONFIG$OUTPUT_DIR, "model_parameters.csv"))
      
      #results_long <- tidy_results(model_results)
      #summaries <- add_total_and_props(results_long)
      #readr::write_csv(summaries$prev,
       #                file.path(CONFIG$OUTPUT_DIR, "summary_prevalence_by_strain.csv"))
      #readr::write_csv(summaries$by_group,
       #                file.path(CONFIG$OUTPUT_DIR, "summary_proportions_by_compartment.csv"))
      #readr::write_csv(summaries$totals,
       #                file.path(CONFIG$OUTPUT_DIR, "summary_totals.csv"))
      
      safe_write_csv(model_results$results,
                     file.path(CONFIG$OUTPUT_DIR, "bacterial_dynamics_results.csv"))
      
      params_tbl <- tibble::enframe(model_results$parameters, name = "parameter", value = "value")
      safe_write_csv(params_tbl,
                     file.path(CONFIG$OUTPUT_DIR, "model_parameters.csv"))
      
      results_long <- tidy_results(model_results)
      summaries <- add_total_and_props(results_long)
      safe_write_csv(summaries$prev,
                     file.path(CONFIG$OUTPUT_DIR, "summary_prevalence_by_strain.csv"))
      safe_write_csv(summaries$by_group,
                     file.path(CONFIG$OUTPUT_DIR, "summary_proportions_by_compartment.csv"))
      safe_write_csv(summaries$totals,
                     file.path(CONFIG$OUTPUT_DIR, "summary_totals.csv"))
      
    }
    model_results
  }, error = function(e) {
    log_message(sprintf("Critical error in model execution: %s", e$message), "CRITICAL")
    traceback(); NULL
  })
}

# Execute ------------------------------------------------------------------####
results <- main()
if (!is.null(results)) {
  log_message("Model execution completed successfully", "SUCCESS")
  log_message(sprintf("Outputs: %s", normalizePath(CONFIG$OUTPUT_DIR, mustWork = FALSE)))
} else {
  log_message("Model execution failed", "FAILURE")
}

