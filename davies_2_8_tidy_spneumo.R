################################################################################
# Azithromycin MDA and antimicrobial-resistant S. pneumoniae in Tanzania
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
  output_dir = "outputs_tanzania_spneumo_mda",
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
  population_model = "static"
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
  # Transmission: calibrate beta.S
  # Calibrated beta.S: 0.0110696. Initial guess was 0.03
  beta.S_start = 0.0110696,

  # Resistance biology.
  c = 0.08, # fitness cost of macrolide resistance; sensitivity 0.00-0.15
  k = 0.40, # relative susceptibility to mixed carriage; sensitivity 0.20-1.00

  # Azithromycin MDA effect on pneumococcal carriage.
  mda_duration = 30,
  mda_cov = 0.9,
  mda_p_clear_S = 0.14,
  mda_p_select_C = 0.14,

  # Mortality effect: off in carriage-only pneumococcal model.
  theta = 0,

  # Approximation for children aged 1-59 months using one-year age groups.
  targeted_age_indices = 1:5,

  # Background macrolide pressure.
  # Prefer calibration to baseline resistance over hard-coding.
  # baseline tau (macrolide use = 2) is 4e-4, Calibrated tau: 0.000886033
  macrolide_use_ddd_per_1000_per_day = 1.25 * 0.000886033 / 4e-4,
  treatment_duration_days = 5,

  # AMR-attributable pneumococcal mortality: off initially.
  amrd_rate = 0,
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
  age_col <- find_first_col(
    data, c("Age", "Age_Category", "Age_Group", "AgeGroup")
  )
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

check_mda_coverage <- function(mda_cov, mda_duration) {
  mda_cov_safe <- min(max(mda_cov, 0), 0.999999)
  r_mda <- -log1p(-mda_cov_safe) / mda_duration
  achieved_cov <- 1 - exp(-r_mda * mda_duration)

  tibble::tibble(
    requested_coverage = mda_cov,
    mda_duration_days = mda_duration,
    treatment_rate_per_day = r_mda,
    achieved_coverage = achieved_cov
  )
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

  age_col <- find_first_col(
    population, c("Age", "Age_Category", "Age_Group", "AgeGroup")
  )
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
    stop(
      "Birth input has fewer age groups than the population input.",
      call. = FALSE
    )
  }
  births <- births[seq_len(n_age), , drop = FALSE]
  births_col <- choose_col(
    births,
    candidates = c("Births", "births"),
    fallback_index = 5,
    label = "births"
  )
  births_count <- convert_units_to_persons(
    as.numeric(births[[births_col]]),
    config$births_units
  )
  birth_rate <- births_count /
    (population$population_persons * config$days_per_year)

  mortality_raw <- read_input_csv(input_files$mortality)
  mortality <- filter_country_year(mortality_raw, config$country, config$year)
  if (nrow(mortality) < n_age) {
    stop(
      "Mortality input has fewer age groups than the population input.",
      call. = FALSE
    )
  }
  mortality <- mortality[seq_len(n_age), , drop = FALSE]
  deaths_col <- choose_col(
    mortality,
    candidates = c("Deaths", "deaths", "Death", "death", "Percentage"),
    fallback_index = 5,
    label = "mortality/deaths"
  )
  death_count <- convert_units_to_persons(
    as.numeric(mortality[[deaths_col]]), config$deaths_units
  )
  mortality_rate <- death_count /
    (population$population_persons * config$days_per_year)

  contact_matrix <- read_contact_matrix(
    input_files$contact_matrix, n_age, age_groups
  )

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

make_initial_state_spneumo <- function(population_vector, age_groups, parameters) {
  age <- age_to_number(age_groups)

  carriage <- dplyr::case_when(
    age < 1 ~ 0.50,
    age < 5 ~ 0.45,
    age < 18 ~ 0.20,
    age < 65 ~ 0.08,
    TRUE ~ 0.10
  )

  resistant_fraction <- dplyr::case_when(
    age < 5 ~ 0.12,
    age < 18 ~ 0.10,
    TRUE ~ 0.08
  )

  initial_resistant <- carriage * resistant_fraction
  initial_sensitive <- carriage - initial_resistant
  initial_uncolonised <- 1 - carriage

  state <- c(
    initial_uncolonised * population_vector,
    initial_sensitive * population_vector,
    initial_resistant * population_vector,
    rep(0, length(population_vector)),
    rep(0, length(population_vector)),
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

mda_schedule <- function(total_years,
                         frequency_per_year,
                         mda_years = NULL) {
  if (is.na(frequency_per_year) || frequency_per_year <= 0) {
    return(numeric(0))
  }

  active_years <- min(
    total_years,
    mda_years %||% total_years
  )

  if (active_years <= 0) {
    return(numeric(0))
  }

  interval_days <- config$days_per_year / frequency_per_year

  # Examples:
  # annual for 2 years   -> 2 rounds
  # biannual for 2 years -> 4 rounds
  n_rounds <- floor(active_years * frequency_per_year)

  if (n_rounds <= 0) {
    return(numeric(0))
  }

  seq(
    from = 0,
    by = interval_days,
    length.out = n_rounds
  )
}

make_parameters <- function(inputs, indices, ageing) {
  age_pars <- make_spneumo_age_parameters(inputs$age_groups)

  daily_u.S <- age_pars$u.S_monthly_age * 12 / config$days_per_year
  daily_u.R <- age_pars$u.R_monthly_age * 12 / config$days_per_year
  daily_u.C <- age_pars$u.C_monthly_age * 12 / config$days_per_year

  scalar_tau <- (baseline_parameters$macrolide_use_ddd_per_1000_per_day / 1000) /
    baseline_parameters$treatment_duration_days

  tau <- scalar_tau * age_pars$tau_age_multiplier

  azt <- rep(0, inputs$n_age)
  targeted <- baseline_parameters$targeted_age_indices
  targeted <- targeted[targeted >= 1 & targeted <= inputs$n_age]
  azt[targeted] <- 1

  mda_cov_safe <- min(max(baseline_parameters$mda_cov, 0), 0.999999)

  r_mda <- -log1p(-mda_cov_safe) / baseline_parameters$mda_duration

  list(
    beta.S = baseline_parameters$beta.S_start,
    u.S = daily_u.S,
    u.R = daily_u.R,
    u.C = daily_u.C,
    k = baseline_parameters$k,
    c = baseline_parameters$c,
    mda_p_clear_S = baseline_parameters$mda_p_clear_S,
    mda_p_select_C = baseline_parameters$mda_p_select_C,
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
    mort = inputs$mortality_rate,
    infectiousness_age = age_pars$infectiousness_age,
    susceptibility_age = age_pars$susceptibility_age
  )
}

make_no_mda_parameters <- function(parameters) {
  parameters$theta <- 0
  parameters$use_mda <- FALSE
  parameters$mda_start_times <- numeric(0)
  parameters
}

make_spneumo_age_parameters <- function(age_groups) {
  age <- age_to_number(age_groups)

  susceptibility_age <- dplyr::case_when(
    age < 1 ~ 1.35,
    age < 5 ~ 1.25,
    age < 18 ~ 0.65,
    age < 65 ~ 0.25,
    TRUE ~ 0.5
  )

  infectiousness_age <- dplyr::case_when(
    age < 1 ~ 1.3,
    age < 5 ~ 1.3,
    age < 18 ~ 0.80,
    age < 65 ~ 0.35,
    TRUE ~ 0.40
  )

  u.S_monthly_age <- dplyr::case_when(
    age < 1 ~ 0.3,
    age < 5 ~ 0.55,
    age < 18 ~ 0.87,
    age < 65 ~ 0.80,
    TRUE ~ 0.70
  )

  u.R_multiplier_age <- dplyr::case_when(
    age < 1 ~ 1,
    age < 5 ~ 1,
    age < 18 ~ 1,
    age < 65 ~ 1,
    TRUE ~ 1
  )

  tau_age_multiplier <- dplyr::case_when(
    age < 1 ~ 2,
    age < 5 ~ 1.75,
    age < 18 ~ 1.2,
    age < 65 ~ 0.1,
    TRUE ~ 0.1
  )

  list(
    susceptibility_age = susceptibility_age,
    infectiousness_age = infectiousness_age,
    u.S_monthly_age = u.S_monthly_age,
    u.R_monthly_age = u.S_monthly_age * u.R_multiplier_age,
    u.C_monthly_age = u.S_monthly_age * 0.80,
    u.R_multiplier_age = u.R_multiplier_age,
    tau_age_multiplier = tau_age_multiplier
  )
}

make_scenario_times <- function(horizon_years, parameters) {
  base_times <- make_times(horizon_years)

  mda_times <- c(
    parameters$mda_start_times,
    parameters$mda_start_times + parameters$mda_duration
  )

  sort(unique(c(base_times, mda_times)))
}

# -----------------------------------------------------------------------------
# 4. ODE system and model runner
# -----------------------------------------------------------------------------

spneumo_odes <- function(t, state, parameters) {
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

  S_infectious <- parameters$infectiousness_age * (S_total / safe_N)
  R_infectious <- parameters$infectiousness_age * (R_total / safe_N)

  lambda.S <- parameters$susceptibility_age *
    as.vector(parameters$beta.S * (parameters$m_contact %*% S_infectious))

  beta.R <- parameters$beta.S * (1 - parameters$c)

  lambda.R <- parameters$susceptibility_age *
    as.vector(beta.R * (parameters$m_contact %*% R_infectious))

  is_mda <- isTRUE(parameters$use_mda) &&
    mda_active(t, parameters$mda_start_times, parameters$mda_duration)

  mda_treatment_rate <- if (is_mda) {
    parameters$r_mda * parameters$azt
  } else {
    rep(0, parameters$n_age)
  }

  a_t <- parameters$tau +
    parameters$mda_p_clear_S * mda_treatment_rate

  a.C_t <- parameters$tau +
    parameters$mda_p_select_C * mda_treatment_rate

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
    func = spneumo_odes,
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
    interval_deaths = c(0, diff(cumulative_deaths)),
    interval_amr_deaths = c(0, diff(cumulative_amr_deaths))
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

plot_macrolide_resistance_time_series <- function(time_series) {
  time_series |>
    dplyr::filter(scenario %in% c("No MDA", "Annual MDA", "Biannual MDA")) |>
    ggplot2::ggplot(ggplot2::aes(
      x = time_years,
      y = resistance_prevalence,
      colour = scenario
    )) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(
      ~horizon_years,
      scales = "free_x", labeller = ggplot2::label_both
    ) +
    ggplot2::labs(
      title = "Macrolide-resistance prevalence over time (S. pneumoniae)",
      x = "Time (years)",
      y = "Macrolide-resistant among carriers (%)",
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
      y = "Macrolide-resistant among carriers (%)",
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

  times <- make_scenario_times(horizon_years, parameters)
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
initial_state <- make_initial_state_spneumo(
  inputs$population_vector, inputs$age_groups, parameters
)

print(
  check_mda_coverage(
    mda_cov = baseline_parameters$mda_cov,
    mda_duration = baseline_parameters$mda_duration
  )
)

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
resistance_time_plot <- plot_macrolide_resistance_time_series(time_series)
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

# -----------------------------------------------------------------------------
# 11. Pneumococcal calibration outputs and age-band diagnostics
# -----------------------------------------------------------------------------

make_age_bands <- function(age_groups) {
  age <- age_to_number(age_groups)

  dplyr::case_when(
    age < 1 ~ "0",
    age < 5 ~ "1-4",
    age < 18 ~ "5-17",
    age < 65 ~ "18-64",
    TRUE ~ "65+"
  )
}

summarise_by_age_band <- function(out, indices, age_groups, days_per_year) {
  idx <- indices

  age_band_lookup <- tibble::tibble(
    age_index = seq_along(age_groups),
    age_group = age_groups,
    age_band = make_age_bands(age_groups)
  )

  purrr::map_dfr(unique(age_band_lookup$age_band), function(band) {
    age_idx <- age_band_lookup$age_index[age_band_lookup$age_band == band]

    X <- as.matrix(out[, idx$X[age_idx] + 1, drop = FALSE])
    S <- as.matrix(out[, idx$S[age_idx] + 1, drop = FALSE])
    R <- as.matrix(out[, idx$R[age_idx] + 1, drop = FALSE])
    Sr <- as.matrix(out[, idx$Sr[age_idx] + 1, drop = FALSE])
    Rs <- as.matrix(out[, idx$Rs[age_idx] + 1, drop = FALSE])

    D <- as.matrix(out[, idx$D[age_idx] + 1, drop = FALSE])
    CumIncR <- as.matrix(out[, idx$CumIncR[age_idx] + 1, drop = FALSE])
    AMRD <- as.matrix(out[, idx$AMRD[age_idx] + 1, drop = FALSE])

    total_population <- rowSums(X + S + R + Sr + Rs)
    colonised <- rowSums(S + R + Sr + Rs)
    sensitive <- rowSums(S + Sr)
    resistant <- rowSums(R + Rs)
    cumulative_deaths <- rowSums(D)
    cumulative_resistant_incidence <- rowSums(CumIncR)
    cumulative_amr_deaths <- rowSums(AMRD)

    tibble::tibble(
      time_days = out[[1]],
      time_years = out[[1]] / days_per_year,
      age_band = band,
      total_population = total_population,
      colonised = colonised,
      sensitive = sensitive,
      resistant = resistant,
      carriage_prevalence = 100 * colonised / total_population,
      resistant_among_carried = dplyr::if_else(
        colonised > 0,
        100 * resistant / colonised,
        NA_real_
      ),
      resistant_prevalence_population = 100 * resistant / total_population,
      cumulative_deaths = cumulative_deaths,
      non_amr_deaths = cumulative_deaths - cumulative_amr_deaths,
      cumulative_amr_deaths = cumulative_amr_deaths,
      cumulative_resistant_incidence = cumulative_resistant_incidence,
      interval_deaths = c(0, diff(cumulative_deaths)),
      interval_amr_deaths = c(0, diff(cumulative_amr_deaths))
    )
  }) |>
    dplyr::mutate(
      age_band = factor(age_band, levels = c("0", "1-4", "5-17", "18-64", "65+"))
    ) |>
    dplyr::arrange(age_band, time_days)
}

# Initial calibration targets.
carriage_targets <- tibble::tibble(
  age_band = factor(
    c("0", "1-4", "5-17", "18-64", "65+"),
    levels = c("0", "1-4", "5-17", "18-64", "65+")
  ),
  target_carriage = c(0.50, 0.45, 0.20, 0.08, 0.10),
  target_resistant_among_carried = c(0.21, 0.21, 0.18, 0.14, 0.14),
  weight_carriage = c(2, 2, 1, 1, 1),
  weight_resistance = c(2, 2, 1, 1, 1)
)

message("Creating age-band calibration summaries...")

equilibrium_age_summary <- summarise_by_age_band(
  out = equilibrium_output,
  indices = indices,
  age_groups = inputs$age_groups,
  days_per_year = config$days_per_year
)

time_series_by_age <- purrr::map_dfr(all_results, function(result) {
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

endpoints_by_age <- time_series_by_age |>
  dplyr::group_by(horizon_years, scenario, age_band) |>
  dplyr::slice_tail(n = 1) |>
  dplyr::ungroup()

calibration_check <- equilibrium_age_summary |>
  dplyr::group_by(age_band) |>
  dplyr::slice_tail(n = 1) |>
  dplyr::ungroup() |>
  dplyr::left_join(carriage_targets, by = "age_band") |>
  dplyr::mutate(
    model_carriage = carriage_prevalence / 100,
    model_resistant_among_carried = resistant_among_carried / 100,
    carriage_difference = model_carriage - target_carriage,
    resistance_difference = model_resistant_among_carried -
      target_resistant_among_carried,
    weighted_carriage_error = weight_carriage * carriage_difference^2,
    weighted_resistance_error = weight_resistance * resistance_difference^2
  )

calibration_score <- calibration_check |>
  dplyr::summarise(
    total_weighted_error = sum(
      weighted_carriage_error + weighted_resistance_error,
      na.rm = TRUE
    )
  ) |>
  dplyr::pull(total_weighted_error)

# Optional helpers for actual parameter calibration. These functions are not run
# automatically because each objective evaluation solves the ODE over the full
# equilibrium horizon.
run_no_mda_equilibrium <- function(parameters, initial_state, years = config$equilibrium_years) {
  p <- make_no_mda_parameters(parameters)

  solve_model(
    times = make_times(years),
    state = initial_state,
    parameters = p
  )
}

extract_final_age_summary <- function(out, indices, age_groups, days_per_year) {
  summarise_by_age_band(out, indices, age_groups, days_per_year) |>
    dplyr::group_by(age_band) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()
}

calibration_error <- function(model_summary, targets) {
  model_summary |>
    dplyr::left_join(targets, by = "age_band") |>
    dplyr::mutate(
      model_carriage = carriage_prevalence / 100,
      model_resistance = resistant_among_carried / 100,
      carriage_error = weight_carriage * (model_carriage - target_carriage)^2,
      resistance_error = weight_resistance *
        (model_resistance - target_resistant_among_carried)^2
    ) |>
    dplyr::summarise(error = sum(carriage_error + resistance_error, na.rm = TRUE)) |>
    dplyr::pull(error)
}

calibrate_beta <- function(parameters, initial_state, indices, age_groups,
                           targets, lower = 0.001, upper = 0.2) {
  objective <- function(beta) {
    p <- parameters
    p$beta.S <- beta

    out <- run_no_mda_equilibrium(p, initial_state)
    final_summary <- extract_final_age_summary(
      out = out,
      indices = indices,
      age_groups = age_groups,
      days_per_year = config$days_per_year
    )

    final_summary |>
      dplyr::left_join(targets, by = "age_band") |>
      dplyr::mutate(
        model_carriage = carriage_prevalence / 100,
        carriage_error = weight_carriage *
          (model_carriage - target_carriage)^2
      ) |>
      dplyr::summarise(error = sum(carriage_error, na.rm = TRUE)) |>
      dplyr::pull(error)
  }

  stats::optimize(objective, interval = c(lower, upper))$minimum
}

calibrate_tau <- function(parameters, initial_state, indices, age_groups,
                          targets, lower = 0, upper = 0.003) {
  objective <- function(tau) {
    p <- parameters
    p$tau <- tau

    out <- run_no_mda_equilibrium(p, initial_state)
    final_summary <- extract_final_age_summary(
      out = out,
      indices = indices,
      age_groups = age_groups,
      days_per_year = config$days_per_year
    )

    final_summary |>
      dplyr::left_join(targets, by = "age_band") |>
      dplyr::mutate(
        model_resistance = resistant_among_carried / 100,
        resistance_error = weight_resistance *
          (model_resistance - target_resistant_among_carried)^2
      ) |>
      dplyr::summarise(error = sum(resistance_error, na.rm = TRUE)) |>
      dplyr::pull(error)
  }

  stats::optimize(objective, interval = c(lower, upper))$minimum
}

# To run a staged calibration, uncomment this block. After updating beta.S and
# tau, rerun the equilibrium and scenario sections above so all outputs use the
# calibrated values.
# parameters$beta.S <- calibrate_beta(
#   parameters = parameters,
#   initial_state = initial_state,
#   indices = indices,
#   age_groups = inputs$age_groups,
#   targets = carriage_targets
# )
# message("Calibrated beta.S: ", signif(parameters$beta.S, 6))

# parameters$tau <- calibrate_tau(
#   parameters = parameters,
#   initial_state = initial_state,
#   indices = indices,
#   age_groups = inputs$age_groups,
#   targets = carriage_targets
# )
# message("Calibrated tau: ", signif(parameters$tau, 6))

make_mda_impact_by_age <- function(endpoints_by_age) {
  no_mda <- endpoints_by_age |>
    dplyr::filter(scenario == "No MDA") |>
    dplyr::select(
      horizon_years,
      age_band,
      no_mda_carriage = carriage_prevalence,
      no_mda_resistant_population = resistant_prevalence_population,
      no_mda_resistant_among_carried = resistant_among_carried
    )

  endpoints_by_age |>
    dplyr::filter(scenario != "No MDA") |>
    dplyr::left_join(no_mda, by = c("horizon_years", "age_band")) |>
    dplyr::mutate(
      carriage_difference_pp = carriage_prevalence - no_mda_carriage,
      resistant_population_difference_pp =
        resistant_prevalence_population - no_mda_resistant_population,
      resistant_among_carried_difference_pp =
        resistant_among_carried - no_mda_resistant_among_carried
    )
}

plot_observed_vs_modelled_carriage_by_age <- function(calibration_check) {
  plot_data <- calibration_check |>
    dplyr::select(age_band, model_carriage, target_carriage) |>
    tidyr::pivot_longer(
      cols = c(model_carriage, target_carriage),
      names_to = "series",
      values_to = "value"
    ) |>
    dplyr::mutate(
      series = dplyr::recode(
        series,
        model_carriage = "Modelled",
        target_carriage = "Observed / target"
      )
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = age_band, y = 100 * value, fill = series)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.7
    ) +
    ggplot2::labs(
      title = "Observed vs modelled pneumococcal carriage by age",
      x = "Age band",
      y = "Carriage prevalence (%)",
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_observed_vs_modelled_resistance_by_age <- function(calibration_check) {
  plot_data <- calibration_check |>
    dplyr::select(
      age_band,
      model_resistant_among_carried,
      target_resistant_among_carried
    ) |>
    tidyr::pivot_longer(
      cols = c(model_resistant_among_carried, target_resistant_among_carried),
      names_to = "series",
      values_to = "value"
    ) |>
    dplyr::mutate(
      series = dplyr::recode(
        series,
        model_resistant_among_carried = "Modelled",
        target_resistant_among_carried = "Observed / target"
      )
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = age_band, y = 100 * value, fill = series)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.7
    ) +
    ggplot2::labs(
      title = "Observed vs modelled macrolide resistance by age",
      x = "Age band",
      y = "Macrolide-resistant among carriers (%)",
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_no_mda_equilibrium_age_profile <- function(equilibrium_age_summary) {
  endpoint <- equilibrium_age_summary |>
    dplyr::group_by(age_band) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()

  plot_data <- endpoint |>
    dplyr::select(
      age_band,
      carriage_prevalence,
      resistant_prevalence_population,
      resistant_among_carried
    ) |>
    tidyr::pivot_longer(
      cols = c(
        carriage_prevalence,
        resistant_prevalence_population,
        resistant_among_carried
      ),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        carriage_prevalence = "Any pneumococcal carriage",
        resistant_prevalence_population = "Resistant carriage in population",
        resistant_among_carried = "Resistant among carriers"
      )
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = age_band, y = value, fill = metric)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.7
    ) +
    ggplot2::labs(
      title = "No-MDA equilibrium age profile",
      x = "Age band",
      y = "Percent",
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_mda_impact_by_age_group <- function(mda_impact_by_age,
                                         horizon = 10,
                                         metric = "resistant_among_carried_difference_pp") {
  metric_labels <- c(
    carriage_difference_pp =
      "Change in carriage prevalence vs no MDA (percentage points)",
    resistant_population_difference_pp =
      "Change in resistant carriage prevalence vs no MDA (percentage points)",
    resistant_among_carried_difference_pp =
      "Change in resistant among carriers vs no MDA (percentage points)"
  )

  y_label <- metric_labels[[metric]] %||% metric

  mda_impact_by_age |>
    dplyr::filter(
      horizon_years == horizon,
      scenario %in% c("Annual MDA", "Biannual MDA")
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = age_band,
        y = .data[[metric]],
        fill = scenario
      )
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.7
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::labs(
      title = paste0("MDA impact by age group at ", horizon, " years"),
      x = "Age band",
      y = y_label,
      fill = "Scenario"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_resistant_carriage_prevalence_over_time <- function(time_series_by_age,
                                                         horizon = 10) {
  time_series_by_age |>
    dplyr::filter(
      horizon_years == horizon,
      scenario %in% c("No MDA", "Annual MDA", "Biannual MDA")
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = time_years,
        y = resistant_prevalence_population,
        colour = scenario
      )
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~age_band, scales = "free_y") +
    ggplot2::labs(
      title = paste0("Resistant pneumococcal carriage prevalence over ", horizon, " years"),
      x = "Time since scenario start (years)",
      y = "Resistant carriage prevalence in population (%)",
      colour = "Scenario"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_resistant_fraction_among_carriers_over_time <- function(time_series_by_age,
                                                             horizon = 10) {
  time_series_by_age |>
    dplyr::filter(
      horizon_years == horizon,
      scenario %in% c("No MDA", "Annual MDA", "Biannual MDA")
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = time_years,
        y = resistant_among_carried,
        colour = scenario
      )
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~age_band, scales = "free_y") +
    ggplot2::labs(
      title = paste0("Macrolide-resistance prevalence by age group over ", horizon, " years (S. pneumoniae)"),
      x = "Time since scenario start (years)",
      y = "Macrolide-resistant among carriers (%)",
      colour = "Scenario"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

mda_impact_by_age <- make_mda_impact_by_age(endpoints_by_age)

observed_modelled_carriage_plot <-
  plot_observed_vs_modelled_carriage_by_age(calibration_check)

observed_modelled_resistance_plot <-
  plot_observed_vs_modelled_resistance_by_age(calibration_check)

no_mda_equilibrium_age_profile_plot <-
  plot_no_mda_equilibrium_age_profile(equilibrium_age_summary)

mda_impact_age_plot <-
  plot_mda_impact_by_age_group(
    mda_impact_by_age,
    horizon = 10,
    metric = "resistant_among_carried_difference_pp"
  )

resistant_carriage_time_plot <-
  plot_resistant_carriage_prevalence_over_time(
    time_series_by_age,
    horizon = 10
  )

resistant_fraction_time_plot <-
  plot_resistant_fraction_among_carriers_over_time(
    time_series_by_age,
    horizon = 10
  )

readr::write_csv(
  carriage_targets,
  output_path("calibration_targets_by_age_band.csv")
)

readr::write_csv(
  equilibrium_age_summary,
  output_path("equilibrium_burn_in_by_age_band.csv")
)

readr::write_csv(
  time_series_by_age,
  output_path("scenario_time_series_by_age_band.csv")
)

readr::write_csv(
  endpoints_by_age,
  output_path("scenario_endpoints_by_age_band.csv")
)

readr::write_csv(
  calibration_check,
  output_path("calibration_check_by_age_band.csv")
)

readr::write_csv(
  mda_impact_by_age,
  output_path("mda_impact_by_age_band.csv")
)

save_plot(
  observed_modelled_carriage_plot,
  "observed_vs_modelled_carriage_by_age.png",
  width = 8,
  height = 5
)

save_plot(
  observed_modelled_resistance_plot,
  "observed_vs_modelled_resistance_by_age.png",
  width = 8,
  height = 5
)

save_plot(
  no_mda_equilibrium_age_profile_plot,
  "no_mda_equilibrium_age_profile.png",
  width = 9,
  height = 5
)

save_plot(
  mda_impact_age_plot,
  "mda_impact_by_age_group.png",
  width = 9,
  height = 5
)

save_plot(
  resistant_carriage_time_plot,
  "resistant_carriage_prevalence_over_time_by_age.png",
  width = 10,
  height = 7
)

save_plot(
  resistant_fraction_time_plot,
  "resistant_fraction_among_carriers_over_time_by_age.png",
  width = 10,
  height = 7
)

saveRDS(
  list(
    carriage_targets = carriage_targets,
    equilibrium_age_summary = equilibrium_age_summary,
    time_series_by_age = time_series_by_age,
    endpoints_by_age = endpoints_by_age,
    calibration_check = calibration_check,
    calibration_score = calibration_score,
    mda_impact_by_age = mda_impact_by_age
  ),
  file = output_path("calibration_outputs.rds")
)

message("Calibration weighted error: ", signif(calibration_score, 4))
message("Done. Outputs written to: ", normalizePath(config$output_dir, mustWork = FALSE))

# -----------------------------------------------------------------------------
# MORDOR-style resistance calibration plot
# Twice-yearly MDA for 2 years, then 3.5 years of follow-up
# -----------------------------------------------------------------------------

message("Running MORDOR-style resistance calibration scenario...")

# Paper targets from the Malawi long-term follow-up paper.
# These are macrolide-resistant among pneumococcal carriage isolates.
mordor_malawi_resistance_targets <- tibble::tibble(
  target_group = c(
    "Azithromycin clusters",
    "Azithromycin clusters",
    "Azithromycin clusters",
    "Placebo clusters",
    "Placebo clusters",
    "Placebo clusters",
    "Non-MDA site",
    "Non-MDA site",
    "Non-MDA site"
  ),
  comparison_group = c(
    "azithromycin",
    "azithromycin",
    "azithromycin",
    "placebo",
    "placebo",
    "placebo",
    "non_mda",
    "non_mda",
    "non_mda"
  ),
  paper_timepoint = c(
    "Baseline",
    "6 months post-MDA",
    "3.5 years post-MDA",
    "Baseline",
    "6 months post-MDA",
    "3.5 years post-MDA",
    "Baseline",
    "6 months",
    "2.5 years"
  ),

  # The model time axis starts at the beginning of the MDA programme.
  # MDA lasts 2 years.
  # 6 months post-MDA is therefore 2.5 years from model start.
  # 3.5 years post-MDA is therefore 5.5 years from model start.
  time_years = c(
    0.0,
    2.5,
    5.5,
    0.0,
    2.5,
    5.5,
    0.0,
    2.5,
    5.5
  ),
  resistant_among_carried = c(
    0.217,
    0.319,
    0.321,
    0.210,
    0.250,
    0.309,
    0.169,
    0.165,
    0.165
  ),
  lower = c(
    0.142,
    0.242,
    0.259,
    0.135,
    0.177,
    0.254,
    0.128,
    0.133,
    0.125
  ),
  upper = c(
    0.317,
    0.408,
    0.390,
    0.311,
    0.341,
    0.371,
    0.218,
    0.203,
    0.214
  )
)

# Optional carriage targets from the same paper's long-term follow-up survey.
# These are not used in the resistance plot, but they are useful to save.
mordor_malawi_carriage_targets <- tibble::tibble(
  target_group = c(
    "Azithromycin clusters, MORDOR-exposed children",
    "Placebo clusters, MORDOR-exposed children",
    "Azithromycin clusters, MORDOR-unexposed children",
    "Placebo clusters, MORDOR-unexposed children"
  ),
  comparison_group = c(
    "azithromycin",
    "placebo",
    "azithromycin",
    "placebo"
  ),
  age_group = c(
    "4-9",
    "4-9",
    "1-3",
    "1-3"
  ),
  carriage_prevalence = c(
    145 / 222,
    161 / 222,
    191 / 226,
    181 / 222
  )
)

# Run the model version corresponding to the paper:
# twice-yearly MDA for 2 years, then no MDA until 5.5 years.
mordor_biannual_result <- run_scenario(
  name = "MORDOR-style biannual MDA for 2 years",
  horizon_years = 5.5,
  base_state = equilibrium_state,
  base_parameters = parameters,
  frequency_per_year = 2,
  mda_years = 2
)

# Also run a matched no-MDA comparator over the same 5.5-year horizon.
mordor_no_mda_result <- run_scenario(
  name = "No MDA, 5.5-year comparator",
  horizon_years = 5.5,
  base_state = equilibrium_state,
  base_parameters = parameters,
  frequency_per_year = NA_real_,
  mda_years = NULL
)

# Summarise both outputs by age band.
mordor_biannual_by_age <- summarise_by_age_band(
  out = mordor_biannual_result$output,
  indices = indices,
  age_groups = inputs$age_groups,
  days_per_year = config$days_per_year
) |>
  dplyr::mutate(scenario = "Model: biannual MDA for 2 years")

mordor_no_mda_by_age <- summarise_by_age_band(
  out = mordor_no_mda_result$output,
  indices = indices,
  age_groups = inputs$age_groups,
  days_per_year = config$days_per_year
) |>
  dplyr::mutate(scenario = "Model: no MDA")

mordor_by_age <- dplyr::bind_rows(
  mordor_biannual_by_age,
  mordor_no_mda_by_age
)

# Combine child age bands to approximate the paper's child isolate population.
# Your model age bands do not exactly match 1-3 and 4-9 years, so this uses
# the broad child bands 1-4 and 5-17. If you later create exact 1-3 and 4-9
# summaries, replace this function.
summarise_mordor_child_resistance <- function(age_summary) {
  age_summary |>
    dplyr::filter(age_band %in% c("1-4", "5-17")) |>
    dplyr::group_by(scenario, time_days, time_years) |>
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
    )
}

mordor_child_modelled <- summarise_mordor_child_resistance(mordor_by_age)

# Extract modelled values closest to the paper's target timepoints.
extract_model_at_target_times <- function(modelled_summary, target_times) {
  purrr::map_dfr(target_times, function(target_time) {
    modelled_summary |>
      dplyr::group_by(scenario) |>
      dplyr::slice_min(
        order_by = abs(time_years - target_time),
        n = 1,
        with_ties = FALSE
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(target_time_years = target_time)
  })
}

mordor_model_at_targets <- extract_model_at_target_times(
  modelled_summary = mordor_child_modelled,
  target_times = c(0, 2.5, 5.5)
)

# Use azithromycin-cluster paper targets as the main calibration target.
mordor_azithro_targets <- mordor_malawi_resistance_targets |>
  dplyr::filter(comparison_group == "azithromycin") |>
  dplyr::transmute(
    target_time_years = time_years,
    paper_timepoint = paper_timepoint,
    target_percent = 100 * resistant_among_carried,
    lower_percent = 100 * lower,
    upper_percent = 100 * upper
  )

# Build a comparison table.
mordor_resistance_calibration_table <- mordor_model_at_targets |>
  dplyr::filter(scenario == "Model: biannual MDA for 2 years") |>
  dplyr::select(
    scenario,
    target_time_years,
    model_time_years = time_years,
    model_resistant_among_carried = resistant_among_carried
  ) |>
  dplyr::left_join(
    mordor_azithro_targets,
    by = "target_time_years"
  ) |>
  dplyr::mutate(
    model_minus_target_pp = model_resistant_among_carried - target_percent,
    squared_error = model_minus_target_pp^2
  )

mordor_resistance_score <- mordor_resistance_calibration_table |>
  dplyr::summarise(
    score = sum(squared_error, na.rm = TRUE)
  ) |>
  dplyr::pull(score)

message(
  "MORDOR-style resistance calibration score: ",
  round(mordor_resistance_score, 3)
)

# Plot model trajectory against the paper's azithromycin-cluster targets.
plot_mordor_resistance_fit <- function(modelled_child_summary,
                                       azithro_targets) {
  modelled_plot <- modelled_child_summary |>
    dplyr::filter(
      scenario %in% c(
        "Model: biannual MDA for 2 years",
        "Model: no MDA"
      )
    ) |>
    dplyr::mutate(
      series = scenario,
      value = resistant_among_carried
    ) |>
    dplyr::select(time_years, value, series)

  target_plot <- azithro_targets |>
    dplyr::mutate(
      series = "Paper: azithromycin clusters",
      value = target_percent
    ) |>
    dplyr::select(
      time_years = target_time_years,
      value,
      lower_percent,
      upper_percent,
      paper_timepoint,
      series
    )

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = modelled_plot,
      ggplot2::aes(
        x = time_years,
        y = value,
        colour = series
      ),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      data = target_plot,
      ggplot2::aes(
        x = time_years,
        y = value,
        colour = series
      ),
      size = 2.8
    ) +
    ggplot2::geom_errorbar(
      data = target_plot,
      ggplot2::aes(
        x = time_years,
        ymin = lower_percent,
        ymax = upper_percent,
        colour = series
      ),
      width = 0.08,
      linewidth = 0.6
    ) +
    ggplot2::geom_vline(
      xintercept = 2,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::annotate(
      "text",
      x = 2,
      y = Inf,
      label = "MDA stops",
      vjust = 1.4,
      hjust = -0.1,
      size = 3.2
    ) +
    ggplot2::labs(
      title = "",
      subtitle = "",
      x = "Years since start of MDA programme",
      y = "Macrolide-resistant among carriers (%)",
      colour = NULL
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 5.5),
      ylim = c(0, 100)
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "Model: no MDA" = "#619CFF",
        "Paper: azithromycin clusters" = "black",
        "Model: biannual MDA for 2 years" = "#F8766D"
      )
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom"
    )
}

mordor_resistance_plot <- plot_mordor_resistance_fit(
  modelled_child_summary = mordor_child_modelled,
  azithro_targets = mordor_azithro_targets
)

save_plot(
  mordor_resistance_plot,
  "mordor_malawi_resistance_calibration.png",
  width = 8,
  height = 5
)

# Write useful diagnostic tables.
readr::write_csv(
  mordor_malawi_resistance_targets,
  output_path("mordor_malawi_resistance_targets.csv")
)

readr::write_csv(
  mordor_malawi_carriage_targets,
  output_path("mordor_malawi_carriage_targets.csv")
)

readr::write_csv(
  mordor_child_modelled,
  output_path("mordor_style_child_modelled_time_series.csv")
)

readr::write_csv(
  mordor_model_at_targets,
  output_path("mordor_style_model_at_target_times.csv")
)

readr::write_csv(
  mordor_resistance_calibration_table,
  output_path("mordor_resistance_calibration_table.csv")
)

saveRDS(
  list(
    mordor_malawi_resistance_targets = mordor_malawi_resistance_targets,
    mordor_malawi_carriage_targets = mordor_malawi_carriage_targets,
    mordor_biannual_result = mordor_biannual_result,
    mordor_no_mda_result = mordor_no_mda_result,
    mordor_by_age = mordor_by_age,
    mordor_child_modelled = mordor_child_modelled,
    mordor_model_at_targets = mordor_model_at_targets,
    mordor_resistance_calibration_table = mordor_resistance_calibration_table,
    mordor_resistance_score = mordor_resistance_score,
    mordor_resistance_plot = mordor_resistance_plot
  ),
  file = output_path("mordor_resistance_calibration_outputs.rds")
)

message(
  "MORDOR-style calibration outputs written to: ",
  normalizePath(config$output_dir, mustWork = FALSE)
)
