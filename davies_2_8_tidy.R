################################################################################
# Azithromycin MDA and antimicrobial-resistant E. coli in Tanzania
# Tidied analysis script
#
# This version keeps the core model, scenario runs, summary tables, and key plots.
# It removes the early exploratory data-processing section and the later
# sensitivity-analysis sections from the original script.
#
# How to use:
#   1. Put this script in the same folder as the input CSV files, or edit
#      config$data_dir below.
#   2. Check the assumptions in config and baseline_parameters.
#   3. Run the script. Tables and figures are written to config$output_dir.
################################################################################

# -----------------------------------------------------------------------------
# 0. User settings
# -----------------------------------------------------------------------------

config <- list(
  country = "United Republic of Tanzania",
  country_label = "Tanzania",
  year = 2023,

  # Folder containing the demographic/contact CSV files.
  # Use "." if the CSV files are in the same folder as this script.
  data_dir = ".",

  # Where to write output CSV files and figures.
  output_dir = "outputs_tanzania_mda",
  save_plots = TRUE,

  # Input units. The original files used values in thousands.
  population_units = "thousands", # "thousands" or "persons"
  births_units = "thousands", # "thousands" or "persons"
  deaths_units = "thousands", # "thousands" or "persons"

  # The original contact matrix was divided by 5 before use.
  contact_scaling_divisor = 5,

  # Time settings.
  days_per_year = 365.25,
  time_step_days = 1,
  equilibrium_years = 110,
  scenario_horizons_years = c(1, 5, 10),
  long_horizon_years = 20,

  # Population model:
  #   "static"  balances births to total deaths, keeping population size stable.
  #   "dynamic" uses the age-specific birth rates from the birth file.
  population_model = "static",

  # Initial colonisation proportions, by age.
  # These must sum to 1.
  initial_uncolonised = 0.95,
  initial_sensitive = 0.025,
  initial_resistant = 0.025,
  initial_sensitive_treated = 0,
  initial_resistant_treated = 0
)

input_files <- list(
  population = "Population_Afro_2023_1yearage.csv",
  births = "3.U.1.Birth_1year_Afro.csv",
  mortality = "3.U.1.AFRO_mortality_by_age_group_1yearage.csv",
  contact_matrix = "3.U.1.contact_Tanzania_1y.csv"
)

# Baseline model parameters. Rates u.S, u.R, and u.C are supplied below as
# monthly rates and converted to daily rates in make_parameters().
baseline_parameters <- list(
  # Pathogen parameters.
  beta.S = 0.03,
  u.S_monthly = 1,
  u.R_monthly = 1,
  u.C_monthly = 1,
  k = 0.5,
  c = 0.20,

  # MDA azithromycin parameters.
  a = 0.16,
  a.C = 0.16,
  mda_duration = 30,
  mda_cov = 1,
  theta = 0.13,
  targeted_age_indices = 1:5, # age groups 0, 1, 2, 3, 4

  # Baseline antibiotic use. a.use_p is DDD/1000 inhabitants/day.
  a.use_p = 23.1,
  d = 7,

  # AMR-attributable mortality rate.
  amrd_rate = (27.3 / 100000 / 365) / (0.15 * 0.9),

  # Other parameters kept for transparency, but not used directly by the ODE.
  kappa = 0
)

# -----------------------------------------------------------------------------
# 1. Packages and small utilities
# -----------------------------------------------------------------------------

required_packages <- c(
  "deSolve", "dplyr", "ggplot2", "purrr", "readr", "scales", "tidyr", "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

make_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

input_path <- function(filename) {
  file.path(config$data_dir, filename)
}

output_path <- function(filename) {
  file.path(config$output_dir, filename)
}

read_input_csv <- function(filename) {
  path <- input_path(filename)
  if (!file.exists(path)) {
    stop("Missing input file: ", path, call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

save_plot <- function(plot, filename, width = 8, height = 5, dpi = 300) {
  if (!isTRUE(config$save_plots)) {
    return(invisible(NULL))
  }
  ggplot2::ggsave(
    filename = output_path(filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  invisible(NULL)
}

find_first_col <- function(data, candidates) {
  found <- intersect(candidates, names(data))
  if (length(found) == 0) {
    return(NA_character_)
  }
  found[[1]]
}

choose_col <- function(data, candidates, fallback_index, label) {
  found <- find_first_col(data, candidates)
  if (!is.na(found)) {
    return(found)
  }

  if (!is.na(fallback_index) && fallback_index <= ncol(data)) {
    warning(
      "Could not find a named ", label, " column. Using column ",
      fallback_index, " ('", names(data)[fallback_index], "').",
      call. = FALSE
    )
    return(names(data)[fallback_index])
  }

  stop("Could not identify ", label, " column.", call. = FALSE)
}

age_to_number <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9]", "", as.character(x))))
}

arrange_by_age <- function(data) {
  age_col <- find_first_col(data, c("Age", "Age_Category", "Age_Group", "AgeGroup"))
  if (is.na(age_col)) {
    return(data)
  }
  data[order(age_to_number(data[[age_col]])), , drop = FALSE]
}

filter_country_year <- function(data, country, year) {
  if (!all(c("Country", "Year") %in% names(data))) {
    stop("Input data must contain Country and Year columns.", call. = FALSE)
  }

  year_values <- suppressWarnings(as.integer(as.character(data$Year)))
  out <- data[data$Country == country & year_values == year, , drop = FALSE]

  if (nrow(out) == 0) {
    stop(
      "No rows found for country '", country, "' and year ", year, ".",
      call. = FALSE
    )
  }

  arrange_by_age(out)
}

convert_units_to_persons <- function(x, units) {
  units <- match.arg(units, c("thousands", "persons"))
  if (units == "thousands") {
    return(1000 * x)
  }
  x
}

# -----------------------------------------------------------------------------
# 2. Load Tanzania demographic and contact inputs
# -----------------------------------------------------------------------------

read_contact_matrix <- function(filename, n_age, age_groups) {
  path <- input_path(filename)
  if (!file.exists(path)) {
    stop("Missing contact matrix file: ", path, call. = FALSE)
  }

  raw <- utils::read.csv(path, check.names = FALSE)

  # If the CSV contains a row-label column, remove it.
  if (ncol(raw) == n_age + 1) {
    raw <- raw[, -1, drop = FALSE]
  }

  mat <- as.matrix(data.frame(lapply(raw, as.numeric), check.names = FALSE))

  if (anyNA(mat)) {
    stop("The contact matrix contains non-numeric values.", call. = FALSE)
  }

  if (!all(dim(mat) == c(n_age, n_age))) {
    stop(
      "Expected a ", n_age, " x ", n_age, " contact matrix, but found ",
      paste(dim(mat), collapse = " x "), ".",
      call. = FALSE
    )
  }

  mat <- mat / config$contact_scaling_divisor
  dimnames(mat) <- list(contactee = age_groups, contactor = age_groups)
  mat
}

load_model_inputs <- function() {
  population_raw <- read_input_csv(input_files$population)
  population <- filter_country_year(population_raw, config$country, config$year)

  pop_col <- choose_col(
    population,
    candidates = c("Population_age", "Population", "population"),
    fallback_index = 5,
    label = "population"
  )

  age_col <- find_first_col(population, c("Age", "Age_Category", "Age_Group", "AgeGroup"))
  age_groups <- if (!is.na(age_col)) {
    as.character(population[[age_col]])
  } else {
    as.character(seq_len(nrow(population)) - 1)
  }

  if (length(age_groups) > 0 && tail(age_groups, 1) == "100") {
    age_groups[length(age_groups)] <- "100+"
  }

  population$population_persons <- convert_units_to_persons(
    as.numeric(population[[pop_col]]),
    config$population_units
  )

  n_age <- nrow(population)
  if (n_age < 2) stop("Population input must contain at least two age groups.")

  births_raw <- read_input_csv(input_files$births)
  births <- filter_country_year(births_raw, config$country, config$year)
  if (nrow(births) < n_age) {
    stop("Birth input has fewer age groups than the population input.", call. = FALSE)
  }
  births <- births[seq_len(n_age), , drop = FALSE]
  births_col <- choose_col(
    births,
    candidates = c("Births", "births"),
    fallback_index = 5,
    label = "births"
  )
  births_count <- convert_units_to_persons(as.numeric(births[[births_col]]), config$births_units)
  birth_rate <- births_count / (population$population_persons * config$days_per_year)

  mortality_raw <- read_input_csv(input_files$mortality)
  mortality <- filter_country_year(mortality_raw, config$country, config$year)
  if (nrow(mortality) < n_age) {
    stop("Mortality input has fewer age groups than the population input.", call. = FALSE)
  }
  mortality <- mortality[seq_len(n_age), , drop = FALSE]
  deaths_col <- choose_col(
    mortality,
    candidates = c("Deaths", "deaths", "Death", "death", "Percentage"),
    fallback_index = 5,
    label = "mortality/deaths"
  )
  death_count <- convert_units_to_persons(as.numeric(mortality[[deaths_col]]), config$deaths_units)
  mortality_rate <- death_count / (population$population_persons * config$days_per_year)

  contact_matrix <- read_contact_matrix(input_files$contact_matrix, n_age, age_groups)

  list(
    population = population,
    births = births,
    mortality = mortality,
    contact_matrix = contact_matrix,
    population_vector = population$population_persons,
    birth_rate = birth_rate,
    mortality_rate = mortality_rate,
    age_groups = age_groups,
    n_age = n_age
  )
}

# -----------------------------------------------------------------------------
# 3. Model setup
# -----------------------------------------------------------------------------

make_ageing_matrix <- function(n_age, days_per_year) {
  ageing <- t(diff(diag(n_age), lag = 1)) / days_per_year
  cbind(ageing, rep(0, n_age))
}

make_indices <- function(n_age) {
  list(
    X = seq_len(n_age),
    S = n_age + seq_len(n_age),
    R = 2 * n_age + seq_len(n_age),
    Sr = 3 * n_age + seq_len(n_age),
    Rs = 4 * n_age + seq_len(n_age),
    D = 5 * n_age + seq_len(n_age),
    CumIncR = 6 * n_age + seq_len(n_age),
    AMRD = 7 * n_age + seq_len(n_age)
  )
}

make_state_names <- function(age_groups) {
  compartments <- c("X", "S", "R", "Sr", "Rs", "D", "CumIncR", "AMRD")
  unlist(lapply(compartments, function(comp) paste0(comp, "_", age_groups)))
}

make_initial_state <- function(population_vector, age_groups) {
  proportions <- c(
    config$initial_uncolonised,
    config$initial_sensitive,
    config$initial_resistant,
    config$initial_sensitive_treated,
    config$initial_resistant_treated
  )

  if (abs(sum(proportions) - 1) > 1e-8) {
    stop("Initial compartment proportions must sum to 1.", call. = FALSE)
  }

  state <- c(
    config$initial_uncolonised * population_vector,
    config$initial_sensitive * population_vector,
    config$initial_resistant * population_vector,
    config$initial_sensitive_treated * population_vector,
    config$initial_resistant_treated * population_vector,
    rep(0, length(population_vector)),
    rep(0, length(population_vector)),
    rep(0, length(population_vector))
  )

  names(state) <- make_state_names(age_groups)
  state
}

mda_active <- function(time, starts, duration) {
  length(starts) > 0 && any(time >= starts & time < starts + duration)
}

mda_schedule <- function(total_years, frequency_per_year, mda_years = NULL) {
  if (is.na(frequency_per_year) || frequency_per_year <= 0) {
    return(numeric(0))
  }

  active_years <- min(total_years, mda_years %||% total_years)
  if (active_years <= 0) {
    return(numeric(0))
  }

  interval_days <- config$days_per_year / frequency_per_year
  end_day <- active_years * config$days_per_year
  seq(0, end_day - 1e-8, by = interval_days)
}

make_parameters <- function(inputs, indices, ageing) {
  daily_u.S <- baseline_parameters$u.S_monthly * 12 / config$days_per_year
  daily_u.R <- baseline_parameters$u.R_monthly * 12 / config$days_per_year
  daily_u.C <- baseline_parameters$u.C_monthly * 12 / config$days_per_year

  tau <- (baseline_parameters$a.use_p / 1000) / baseline_parameters$d

  azt <- rep(0, inputs$n_age)
  targeted <- baseline_parameters$targeted_age_indices
  targeted <- targeted[targeted >= 1 & targeted <= inputs$n_age]
  azt[targeted] <- 1

  r_mda <- if (baseline_parameters$mda_cov >= 1) {
    Inf
  } else {
    -log(1 - baseline_parameters$mda_cov) / baseline_parameters$mda_duration
  }

  list(
    beta.S = baseline_parameters$beta.S,
    u.S = daily_u.S,
    u.R = daily_u.R,
    u.C = daily_u.C,
    k = baseline_parameters$k,
    c = baseline_parameters$c,
    a = baseline_parameters$a,
    a.C = baseline_parameters$a.C,
    mda_duration = baseline_parameters$mda_duration,
    mda_cov = baseline_parameters$mda_cov,
    theta = baseline_parameters$theta,
    tau = tau,
    azt = azt,
    r_mda = r_mda,
    amrd_rate = baseline_parameters$amrd_rate,
    kappa = baseline_parameters$kappa,
    use_mda = TRUE,
    mda_start_times = numeric(0),
    population_model = config$population_model,

    # Data inputs used by the ODE.
    idx = indices,
    n_age = inputs$n_age,
    m_contact = inputs$contact_matrix,
    ageing = ageing,
    birth_rate = inputs$birth_rate,
    mort = inputs$mortality_rate
  )
}

make_no_mda_parameters <- function(parameters) {
  parameters$a <- 0
  parameters$a.C <- 0
  parameters$theta <- 0
  parameters$use_mda <- FALSE
  parameters$mda_start_times <- numeric(0)
  parameters
}

# -----------------------------------------------------------------------------
# 4. ODE system and model runner
# -----------------------------------------------------------------------------

ecoli_odes <- function(t, state, parameters) {
  idx <- parameters$idx

  X <- state[idx$X]
  S <- state[idx$S]
  R <- state[idx$R]
  Sr <- state[idx$Sr]
  Rs <- state[idx$Rs]

  N <- X + S + R + Sr + Rs
  safe_N <- ifelse(N > 0, N, 1)

  S_total <- S + Sr
  R_total <- R + Rs

  lambda.S <- as.vector(parameters$beta.S * (parameters$m_contact %*% (S_total / safe_N)))
  beta.R <- parameters$beta.S * (1 - parameters$c)
  lambda.R <- as.vector(beta.R * (parameters$m_contact %*% (R_total / safe_N)))

  is_mda <- isTRUE(parameters$use_mda) &&
    mda_active(t, parameters$mda_start_times, parameters$mda_duration)

  mda_effect <- as.numeric(is_mda) * parameters$mda_cov * parameters$azt
  a_t <- parameters$tau + parameters$a * mda_effect
  a.C_t <- parameters$tau + parameters$a.C * mda_effect

  mort_eff <- parameters$mort
  if (is_mda) {
    under_five <- seq_len(min(5, length(mort_eff)))
    mort_eff[under_five] <- parameters$mort[under_five] * (1 - parameters$theta)
  }

  births <- rep(0, parameters$n_age)
  if (parameters$population_model == "dynamic") {
    births[1] <- sum(parameters$birth_rate * N)
  } else if (parameters$population_model == "static") {
    births[1] <- sum(mort_eff * N + parameters$amrd_rate * (R + Rs))
  } else {
    stop("population_model must be 'static' or 'dynamic'.")
  }

  age_X <- as.vector(parameters$ageing %*% X)
  age_S <- as.vector(parameters$ageing %*% S)
  age_R <- as.vector(parameters$ageing %*% R)
  age_Sr <- as.vector(parameters$ageing %*% Sr)
  age_Rs <- as.vector(parameters$ageing %*% Rs)

  dX <- births + (parameters$u.S + a_t) * S + parameters$u.R * R +
    parameters$u.C * (Sr + Rs) - (lambda.S + lambda.R) * X +
    age_X - mort_eff * X

  dS <- lambda.S * X - parameters$u.S * S - parameters$k * lambda.R * S -
    a_t * S + age_S - mort_eff * S

  dR <- lambda.R * X - parameters$u.R * R - parameters$k * lambda.S * R +
    a.C_t * (Sr + Rs) + age_R - (mort_eff + parameters$amrd_rate) * R

  dSr <- parameters$k * lambda.R * S - parameters$u.C * Sr - a.C_t * Sr +
    age_Sr - mort_eff * Sr

  dRs <- parameters$k * lambda.S * R - parameters$u.C * Rs - a.C_t * Rs +
    age_Rs - (mort_eff + parameters$amrd_rate) * Rs

  dD <- mort_eff * (X + S + R + Sr + Rs) + parameters$amrd_rate * (R + Rs)
  dCumIncR <- lambda.R * X + parameters$k * lambda.R * S
  dAMRD <- parameters$amrd_rate * (R + Rs)

  list(c(dX, dS, dR, dSr, dRs, dD, dCumIncR, dAMRD))
}

solve_model <- function(times, state, parameters) {
  as.data.frame(deSolve::ode(
    y = state,
    times = times,
    func = ecoli_odes,
    parms = parameters,
    method = "lsoda"
  ))
}

make_times <- function(years) {
  seq(0, years * config$days_per_year, by = config$time_step_days)
}

# -----------------------------------------------------------------------------
# 5. Summary metrics
# -----------------------------------------------------------------------------

summarise_model_output <- function(out, indices, days_per_year) {
  idx <- indices

  total_population <- rowSums(out[, c(idx$X, idx$S, idx$R, idx$Sr, idx$Rs) + 1])
  colonised <- rowSums(out[, c(idx$S, idx$R, idx$Sr, idx$Rs) + 1])
  sensitive <- rowSums(out[, c(idx$S, idx$Sr) + 1])
  resistant <- rowSums(out[, c(idx$R, idx$Rs) + 1])
  cumulative_deaths <- rowSums(out[, idx$D + 1])
  cumulative_resistant_incidence <- rowSums(out[, idx$CumIncR + 1])
  cumulative_amr_deaths <- rowSums(out[, idx$AMRD + 1])

  tibble::tibble(
    time_days = out[[1]],
    time_years = out[[1]] / days_per_year,
    total_population = total_population,
    colonised = colonised,
    sensitive = sensitive,
    resistant = resistant,
    resistance_prevalence = 100 * resistant / colonised,
    colonisation_prevalence = 100 * colonised / total_population,
    cumulative_deaths = cumulative_deaths,
    non_amr_deaths = cumulative_deaths - cumulative_amr_deaths,
    cumulative_amr_deaths = cumulative_amr_deaths,
    cumulative_resistant_incidence = cumulative_resistant_incidence,
    daily_deaths = c(0, diff(cumulative_deaths)),
    daily_amr_deaths = c(0, diff(cumulative_amr_deaths))
  )
}

endpoint_summary <- function(time_series) {
  time_series |>
    dplyr::group_by(horizon_years, scenario) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()
}

compare_with_no_mda <- function(endpoints) {
  baseline <- endpoints |>
    dplyr::filter(scenario == "No MDA") |>
    dplyr::select(
      horizon_years,
      no_mda_non_amr_deaths = non_amr_deaths,
      no_mda_amr_deaths = cumulative_amr_deaths,
      no_mda_total_deaths = cumulative_deaths,
      no_mda_resistance = resistance_prevalence,
      no_mda_resistant_incidence = cumulative_resistant_incidence
    )

  endpoints |>
    dplyr::left_join(baseline, by = "horizon_years") |>
    dplyr::mutate(
      non_amr_deaths_averted = no_mda_non_amr_deaths - non_amr_deaths,
      extra_amr_deaths = cumulative_amr_deaths - no_mda_amr_deaths,
      net_deaths_averted = non_amr_deaths_averted - extra_amr_deaths,
      total_deaths_averted = no_mda_total_deaths - cumulative_deaths,
      resistance_difference_percentage_points = resistance_prevalence - no_mda_resistance,
      extra_resistant_incidence = cumulative_resistant_incidence - no_mda_resistant_incidence
    )
}

# -----------------------------------------------------------------------------
# 6. Plotting helpers
# -----------------------------------------------------------------------------

plot_population_structure <- function(inputs) {
  population <- inputs$population
  age_groups <- inputs$age_groups
  population$age_group <- factor(age_groups, levels = age_groups)

  ggplot2::ggplot(population, ggplot2::aes(x = age_group, y = population_persons)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = " M")) +
    ggplot2::labs(
      title = paste0(config$country_label, " population structure (", config$year, ")"),
      x = "Age group",
      y = "Population"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
}

plot_contact_matrix <- function(inputs) {
  contact_df <- as.data.frame(as.table(inputs$contact_matrix), responseName = "contacts")
  names(contact_df) <- c("contactee", "contactor", "contacts")

  ggplot2::ggplot(contact_df, ggplot2::aes(x = contactor, y = contactee, fill = contacts)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    ggplot2::labs(
      title = paste0(config$country_label, " contact matrix"),
      x = "Contactor age",
      y = "Contactee age",
      fill = "Contacts"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 6)
    )
}

plot_resistance_time_series <- function(time_series) {
  time_series |>
    dplyr::filter(scenario %in% c("No MDA", "Annual MDA", "Biannual MDA")) |>
    ggplot2::ggplot(ggplot2::aes(
      x = time_years,
      y = resistance_prevalence,
      colour = scenario
    )) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~horizon_years, scales = "free_x", labeller = ggplot2::label_both) +
    ggplot2::labs(
      title = "Resistance prevalence over time",
      x = "Time (years)",
      y = "Resistant among colonised (%)",
      colour = "Scenario"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_net_deaths_averted <- function(comparison) {
  comparison |>
    dplyr::filter(scenario != "No MDA") |>
    ggplot2::ggplot(ggplot2::aes(
      x = factor(horizon_years),
      y = net_deaths_averted,
      fill = scenario
    )) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::labs(
      title = "Net deaths averted relative to no MDA",
      x = "Horizon (years)",
      y = "Net deaths averted",
      fill = "Scenario"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_resistance_endpoint <- function(endpoints) {
  endpoints |>
    dplyr::filter(scenario %in% c("No MDA", "Annual MDA", "Biannual MDA")) |>
    ggplot2::ggplot(ggplot2::aes(
      x = factor(horizon_years),
      y = resistance_prevalence,
      fill = scenario
    )) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(
      title = "Resistance prevalence at the end of each horizon",
      x = "Horizon (years)",
      y = "Resistant among colonised (%)",
      fill = "Scenario"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# 7. Run model scenarios
# -----------------------------------------------------------------------------

run_scenario <- function(name, horizon_years, base_state, base_parameters,
                         frequency_per_year = NA_real_, mda_years = NULL) {
  parameters <- base_parameters

  if (is.na(frequency_per_year)) {
    parameters <- make_no_mda_parameters(parameters)
  } else {
    parameters$use_mda <- TRUE
    parameters$mda_start_times <- mda_schedule(
      total_years = horizon_years,
      frequency_per_year = frequency_per_year,
      mda_years = mda_years
    )
  }

  times <- make_times(horizon_years)
  output <- solve_model(times, base_state, parameters)

  list(
    scenario = name,
    horizon_years = horizon_years,
    frequency_per_year = frequency_per_year,
    mda_years = mda_years,
    parameters = parameters,
    output = output
  )
}

run_main_scenarios <- function(state, parameters) {
  horizons <- c(config$scenario_horizons_years, config$long_horizon_years)

  specs <- tidyr::expand_grid(
    horizon_years = horizons,
    scenario = c("No MDA", "Annual MDA", "Biannual MDA")
  ) |>
    dplyr::mutate(
      frequency_per_year = dplyr::case_when(
        scenario == "Annual MDA" ~ 1,
        scenario == "Biannual MDA" ~ 2,
        TRUE ~ NA_real_
      ),
      # For the long horizon, keep MDA to the first 10 years.
      # For shorter horizons, MDA continues for the full horizon.
      mda_years = dplyr::if_else(
        horizon_years == config$long_horizon_years & scenario != "No MDA",
        10,
        NA_real_
      )
    )

  purrr::pmap(
    specs,
    function(horizon_years, scenario, frequency_per_year, mda_years) {
      run_scenario(
        name = scenario,
        horizon_years = horizon_years,
        base_state = state,
        base_parameters = parameters,
        frequency_per_year = frequency_per_year,
        mda_years = if (is.na(mda_years)) NULL else mda_years
      )
    }
  )
}

run_stop_scenarios <- function(state, parameters) {
  specs <- tidyr::expand_grid(
    frequency_per_year = c(1, 2),
    mda_years = c(5, 6, 7)
  ) |>
    dplyr::mutate(
      scenario = dplyr::case_when(
        frequency_per_year == 1 ~ paste0("Annual MDA stopped after ", mda_years, " years"),
        frequency_per_year == 2 ~ paste0("Biannual MDA stopped after ", mda_years, " years")
      ),
      horizon_years = 10
    )

  purrr::pmap(
    specs,
    function(frequency_per_year, mda_years, scenario, horizon_years) {
      run_scenario(
        name = scenario,
        horizon_years = horizon_years,
        base_state = state,
        base_parameters = parameters,
        frequency_per_year = frequency_per_year,
        mda_years = mda_years
      )
    }
  )
}

# -----------------------------------------------------------------------------
# 8. Execute analysis
# -----------------------------------------------------------------------------

make_dir(config$output_dir)

message("Loading model inputs...")
inputs <- load_model_inputs()

indices <- make_indices(inputs$n_age)
ageing <- make_ageing_matrix(inputs$n_age, config$days_per_year)
parameters <- make_parameters(inputs, indices, ageing)
initial_state <- make_initial_state(inputs$population_vector, inputs$age_groups)

message("Running no-MDA equilibrium burn-in...")
equilibrium_parameters <- make_no_mda_parameters(parameters)
equilibrium_times <- make_times(config$equilibrium_years)
equilibrium_output <- solve_model(equilibrium_times, initial_state, equilibrium_parameters)
equilibrium_state <- as.numeric(equilibrium_output[nrow(equilibrium_output), -1])
names(equilibrium_state) <- make_state_names(inputs$age_groups)

message("Running main scenarios...")
main_results <- run_main_scenarios(equilibrium_state, parameters)

message("Running stop-MDA scenarios...")
stop_results <- run_stop_scenarios(equilibrium_state, parameters)

all_results <- c(main_results, stop_results)

message("Summarising outputs...")
time_series <- purrr::map_dfr(all_results, function(result) {
  summarise_model_output(result$output, indices, config$days_per_year) |>
    dplyr::mutate(
      scenario = result$scenario,
      horizon_years = result$horizon_years,
      frequency_per_year = result$frequency_per_year,
      mda_years = result$mda_years %||% result$horizon_years,
      .before = 1
    )
})

endpoints <- endpoint_summary(time_series)
comparison <- compare_with_no_mda(endpoints)

equilibrium_summary <- summarise_model_output(
  equilibrium_output,
  indices,
  config$days_per_year
)

# -----------------------------------------------------------------------------
# 9. Write outputs
# -----------------------------------------------------------------------------

readr::write_csv(time_series, output_path("scenario_time_series.csv"))
readr::write_csv(endpoints, output_path("scenario_endpoints.csv"))
readr::write_csv(comparison, output_path("scenario_comparison_vs_no_mda.csv"))
readr::write_csv(equilibrium_summary, output_path("equilibrium_burn_in.csv"))

saveRDS(
  list(
    config = config,
    baseline_parameters = baseline_parameters,
    inputs = inputs,
    parameters = parameters,
    equilibrium_output = equilibrium_output,
    scenario_results = all_results,
    time_series = time_series,
    endpoints = endpoints,
    comparison = comparison
  ),
  file = output_path("model_results.rds")
)

# -----------------------------------------------------------------------------
# 10. Key figures
# -----------------------------------------------------------------------------

population_plot <- plot_population_structure(inputs)
contact_plot <- plot_contact_matrix(inputs)
resistance_time_plot <- plot_resistance_time_series(time_series)
net_deaths_plot <- plot_net_deaths_averted(comparison)
endpoint_resistance_plot <- plot_resistance_endpoint(endpoints)

save_plot(population_plot, "population_structure.png", width = 10, height = 5)
save_plot(contact_plot, "contact_matrix.png", width = 9, height = 8)
save_plot(resistance_time_plot, "resistance_prevalence_time_series.png", width = 11, height = 6)
save_plot(net_deaths_plot, "net_deaths_averted.png", width = 8, height = 5)
save_plot(endpoint_resistance_plot, "endpoint_resistance_prevalence.png", width = 8, height = 5)

# Show a compact endpoint table in the console.
print(
  comparison |>
    dplyr::select(
      horizon_years,
      scenario,
      resistance_prevalence,
      non_amr_deaths_averted,
      extra_amr_deaths,
      net_deaths_averted,
      total_deaths_averted
    ) |>
    dplyr::arrange(horizon_years, scenario)
)

message("Done. Outputs written to: ", normalizePath(config$output_dir, mustWork = FALSE))
