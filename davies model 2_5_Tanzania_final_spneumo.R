#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(deSolve)
  library(dplyr)
})

# ================================================================
# Core model: Azithromycin MDA, AMR, age structure
# S. pneumoniae carriage-AMR version
# ================================================================
# Key choices in this version:
# - same scaffold as the cleaned E. coli core model
# - stationary demography (births balance all deaths each time step)
# - under-5 targeting uses correct R indexing (1:5 for ages 0-4)
# - routine antibiotic pressure is separate from MDA pressure
# - use_mda actually switches MDA on/off
# - m_contact is taken from the parameter list (not a global object)
# - mda_cov is applied through the intervention vector
# - cumulative resistant incidence is tracked explicitly
# - AMR-attributable deaths are tracked explicitly
# ================================================================

# -----------------------------
# Utilities
# -----------------------------

safe_numeric_column <- function(df, preferred_name = NULL, fallback_index = NULL) {
  if (!is.null(preferred_name) && preferred_name %in% names(df)) {
    return(as.numeric(df[[preferred_name]]))
  }
  if (!is.null(fallback_index) && fallback_index <= ncol(df)) {
    return(as.numeric(df[[fallback_index]]))
  }
  stop("Could not find requested numeric column.")
}

make_ageing_matrix <- function(n_age) {
  dd <- rep(1, n_age)
  ageing <- t(diff(diag(dd), lag = 1) / 365.25)
  ageing <- cbind(ageing, rep(0, n_age))
  ageing
}

make_indices <- function(n_age) {
  list(
    X = 1:n_age,
    S = (1 * n_age + 1):(2 * n_age),
    R = (2 * n_age + 1):(3 * n_age),
    Sr = (3 * n_age + 1):(4 * n_age),
    Rs = (4 * n_age + 1):(5 * n_age),
    D = (5 * n_age + 1):(6 * n_age),
    CumIncR = (6 * n_age + 1):(7 * n_age),
    AMRD = (7 * n_age + 1):(8 * n_age)
  )
}

make_column_names <- function(age_groups) {
  compartments <- c("X", "S", "R", "Sr", "Rs", "D", "CumIncR", "AMRD")
  out <- "time"
  for (comp in compartments) {
    for (age in age_groups) {
      out <- c(out, paste0(comp, "_", age))
    }
  }
  out
}

# -----------------------------
# Data loading
# -----------------------------

load_tanzania_context <- function(data_dir = ".") {
  pop_path <- file.path(data_dir, "Population_Afro_2023_1yearage.csv")
  mort_path <- file.path(
    data_dir, "3.U.1.AFRO_mortality_by_age_group_1yearage.csv"
  )
  contact_path <- file.path(data_dir, "3.U.1.contact_Tanzania_1y.csv")

  pop_df <- read.csv(pop_path, check.names = FALSE)
  mort_df <- read.csv(mort_path, check.names = FALSE)
  contact_df <- read.csv(contact_path, check.names = FALSE)

  # Population: Tanzania, 2023
  pop_tz <- pop_df %>%
    filter(Country == "United Republic of Tanzania", Year %in% c(2023, "2023"))

  if (nrow(pop_tz) == 0) stop("No Tanzania population rows found for 2023.")

  pop_vec <- safe_numeric_column(pop_tz, preferred_name = "Population_age", fallback_index = 5)
  pop_vec <- pop_vec * 1000

  n_age <- length(pop_vec)
  if (n_age < 2) stop("Population vector looks invalid.")

  age_groups <- if (n_age == 101) c(as.character(0:99), "100+") else as.character(seq_len(n_age) - 1)

  # Mortality: convert to per person per day, matching original script logic
  mort_tz <- mort_df %>%
    filter(Country == "United Republic of Tanzania", Year %in% c(2023, "2023"))

  if (nrow(mort_tz) == 0) stop("No Tanzania mortality rows found for 2023.")
  if (nrow(mort_tz) != n_age) stop("Population and mortality age dimensions do not match.")

  mort_raw <- safe_numeric_column(mort_tz, preferred_name = "Percentage", fallback_index = 5)
  mort_per_day <- 1000 * mort_raw / (pop_vec * 365.25)

  # Contact matrix
  contact_mat <- as.matrix(contact_df)
  suppressWarnings(storage.mode(contact_mat) <- "numeric")

  # Some csv exports include an extra first column; drop it if needed
  if (ncol(contact_mat) == n_age + 1) {
    contact_mat <- contact_mat[, -1, drop = FALSE]
  }
  if (nrow(contact_mat) == n_age + 1) {
    contact_mat <- contact_mat[-1, , drop = FALSE]
  }
  if (nrow(contact_mat) != n_age || ncol(contact_mat) != n_age) {
    stop("Contact matrix dimensions do not match age structure.")
  }

  # Match original scaling
  contact_mat <- contact_mat / 25
  colnames(contact_mat) <- age_groups
  rownames(contact_mat) <- age_groups

  list(
    data_dir = data_dir,
    n_age = n_age,
    age_groups = age_groups,
    population_by_age = pop_vec,
    mort = mort_per_day,
    ageing = make_ageing_matrix(n_age),
    m_contact = contact_mat,
    indices = make_indices(n_age)
  )
}

# -----------------------------
# Intervention helpers
# -----------------------------

mda_active <- function(time, mda_starts, duration) {
  any(vapply(mda_starts, function(start) {
    time >= start && time < (start + duration)
  }, logical(1)))
}

make_mda_vector <- function(n_age, mda_cov, targeted_ages = 1:5) {
  out <- rep(0, n_age)
  out[targeted_ages] <- mda_cov
  out
}

annual_mda_starts <- function(n_years) (0:n_years) * 365.25
biannual_mda_starts <- function(n_years) (0:(2 * n_years)) * (365.25 / 2)

# -----------------------------
# Parameters and initial state
# -----------------------------

make_spneumo_parameters <- function(ctx) {
  pars <- list(
    # Pathogen parameters
    # Provisional placeholders for a matched carriage-AMR scaffold;
    # these should be calibrated to pneumococcal carriage/resistance data.
    beta.S = 5.5,
    u.S = 1.5,
    u.R = 1.5,
    u.C = 1.5,
    k = 0.5,
    c = 0.10,

    # MDA parameters
    a = 0.16,
    a.C = 0.16,
    mda_cycle = 365,
    mda_duration = 30,
    mda_cov = 0.6,
    theta = 0,

    # Routine antibiotic pressure
    tau = (25 / 1000) / 7,

    # Country-specific contact matrix
    m_contact = ctx$m_contact,

    # Other
    kappa = 0,
    amrd_rate = 0,

    # Demography and structure
    mort = ctx$mort,
    ageing = ctx$ageing,
    n_age = ctx$n_age,
    use_mda = TRUE,
    mda_start_times = numeric(0),
    targeted_ages = 1:5
  )

  # Convert transmission/clearance parameters from per month to per day,
  # mirroring the E. coli scaffold.
  pars[1:4] <- lapply(pars[1:4], function(x) x * 12 / 365.25)

  pars$azt <- make_mda_vector(ctx$n_age, pars$mda_cov, pars$targeted_ages)
  pars
}

make_no_mda_parameters <- function(pars) {
  out <- pars
  out$use_mda <- FALSE
  out$a <- 0
  out$a.C <- 0
  out$theta <- 0
  out$mda_start_times <- numeric(0)
  out
}

make_initial_state <- function(population_by_age,
                               init_fractions = c(X = 0.95, S = 0.025, R = 0.025, Sr = 0, Rs = 0)) {
  stopifnot(abs(sum(init_fractions) - 1) < 1e-8)

  initX <- init_fractions["X"] * population_by_age
  initS <- init_fractions["S"] * population_by_age
  initR <- init_fractions["R"] * population_by_age
  initSr <- init_fractions["Sr"] * population_by_age
  initRs <- init_fractions["Rs"] * population_by_age
  initD <- 0 * population_by_age
  initCumIncR <- 0 * population_by_age
  initAMRD <- 0 * population_by_age

  c(initX, initS, initR, initSr, initRs, initD, initCumIncR, initAMRD)
}

# -----------------------------
# ODE system (stationary demography)
# -----------------------------

bacteria_odes <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    idx <- make_indices(n_age)

    X <- state[idx$X]
    S <- state[idx$S]
    R <- state[idx$R]
    Sr <- state[idx$Sr]
    Rs <- state[idx$Rs]
    D <- state[idx$D]
    CumIncR <- state[idx$CumIncR]
    AMRD <- state[idx$AMRD]

    N <- X + S + R + Sr + Rs
    denom <- pmax(N, 1e-12)

    S.tot <- S + Sr
    R.tot <- R + Rs

    lamda.S <- beta.S * (m_contact %*% (S.tot / denom))
    lamda.R <- beta.R * (m_contact %*% (R.tot / denom))

    is_mda <- isTRUE(use_mda) && mda_active(t, mda_start_times, mda_duration)

    a_t <- tau + ifelse(is_mda, a * azt, 0)
    a.C_t <- tau + ifelse(is_mda, a.C * azt, 0)

    mort_eff <- mort
    mort_eff[1:5] <- if (is_mda) mort[1:5] * (1 - theta) else mort[1:5]

    # Stationary demography: births balance all deaths at each time step
    death_flow <- mort_eff * N + amrd_rate * (R + Rs)
    births <- rep(0, n_age)
    births[1] <- sum(death_flow)

    dX <- births +
      (u.S + a_t) * S +
      u.R * R +
      u.C * (Sr + Rs) -
      (lamda.S + lamda.R) * X +
      ageing %*% X -
      mort_eff * X

    dS <- lamda.S * X -
      u.S * S -
      k * lamda.R * S -
      a_t * S +
      ageing %*% S -
      mort_eff * S

    dR <- lamda.R * X -
      u.R * R -
      k * lamda.S * R +
      a.C_t * (Sr + Rs) +
      ageing %*% R -
      mort_eff * R -
      amrd_rate * R

    dSr <- k * lamda.R * S - u.C * Sr - a.C_t * Sr + ageing %*% Sr - mort_eff * Sr

    dRs <- k * lamda.S * R - u.C * Rs - a.C_t * Rs + ageing %*% Rs - mort_eff * Rs - amrd_rate * Rs

    dD <- mort_eff * (X + S + R + Sr + Rs) + amrd_rate * (R + Rs)
    dCumIncR <- lamda.R * X + k * lamda.R * S
    dAMRD <- amrd_rate * (R + Rs)

    list(c(dX, dS, dR, dSr, dRs, dD, dCumIncR, dAMRD))
  })
}

solve_model <- function(times, state, parameters) {
  parameters[["beta.R"]] <- parameters[["beta.S"]] * (1 - parameters[["c"]])
  as.data.frame(ode(y = state, times = times, func = bacteria_odes, parms = parameters))
}

# -----------------------------
# Post-processing helpers
# -----------------------------

add_output_names <- function(out_df, age_groups) {
  colnames(out_df) <- make_column_names(age_groups)
  out_df
}

compartment_totals <- function(out_df, n_age) {
  idx <- make_indices(n_age)
  list(
    X = rowSums(out_df[, idx$X + 1, drop = FALSE]),
    S = rowSums(out_df[, idx$S + 1, drop = FALSE]),
    R = rowSums(out_df[, idx$R + 1, drop = FALSE]),
    Sr = rowSums(out_df[, idx$Sr + 1, drop = FALSE]),
    Rs = rowSums(out_df[, idx$Rs + 1, drop = FALSE]),
    D = rowSums(out_df[, idx$D + 1, drop = FALSE]),
    CumIncR = rowSums(out_df[, idx$CumIncR + 1, drop = FALSE]),
    AMRD = rowSums(out_df[, idx$AMRD + 1, drop = FALSE])
  )
}

summarise_output <- function(out_df, ctx) {
  tot <- compartment_totals(out_df, ctx$n_age)

  live_total <- tot$X + tot$S + tot$R + tot$Sr + tot$Rs
  colonised_total <- tot$S + tot$R + tot$Sr + tot$Rs
  resistant_total <- tot$R + tot$Rs
  cum_res_inc <- tot$CumIncR

  data.frame(
    time = out_df$time,
    total_population = live_total,
    carriage_prevalence_pct = 100 * colonised_total / pmax(live_total, 1e-12),
    resistance_among_carriers_pct = 100 * resistant_total / pmax(colonised_total, 1e-12),
    resistant_prevalence_pct = 100 * resistant_total / pmax(live_total, 1e-12),
    cum_resistant_incidence = cum_res_inc,
    resistant_incidence = c(0, diff(cum_res_inc)),
    cumulative_deaths = tot$D,
    daily_deaths = c(0, diff(tot$D)),
    cumulative_amr_deaths = tot$AMRD,
    daily_amr_deaths = c(0, diff(tot$AMRD))
  )
}

extract_under5_deaths <- function(out_df) {
  d_cols <- grep("^D_[0-4]$", names(out_df))
  amrd_cols <- grep("^AMRD_[0-4]$", names(out_df))

  d_cum <- rowSums(out_df[, d_cols, drop = FALSE])
  amrd_cum <- rowSums(out_df[, amrd_cols, drop = FALSE])

  list(
    deaths_total = sum(c(0, diff(d_cum))),
    amr_deaths_total = sum(c(0, diff(amrd_cum)))
  )
}

final_age_distribution <- function(out_df, ctx, organism = "E. coli", scenario = "") {
  idx <- ctx$indices
  final_state <- as.numeric(out_df[nrow(out_df), -1, drop = TRUE])

  X_final <- final_state[idx$X]
  S_final <- final_state[idx$S]
  R_final <- final_state[idx$R]
  Sr_final <- final_state[idx$Sr]
  Rs_final <- final_state[idx$Rs]

  live_totals <- X_final + S_final + R_final + Sr_final + Rs_final
  carriage <- S_final + R_final + Sr_final + Rs_final
  resistant <- R_final + Rs_final

  data.frame(
    organism = organism,
    scenario = scenario,
    age_group = ctx$age_groups,
    population = live_totals,
    carriage_prevalence_pct = 100 * carriage / pmax(live_totals, 1e-12),
    resistance_among_carriers_pct = 100 * resistant / pmax(carriage, 1e-12),
    resistant_prevalence_pct = 100 * resistant / pmax(live_totals, 1e-12)
  )
}

# -----------------------------
# Scenario runner
# -----------------------------

run_equilibrium <- function(ctx, state0, pars_no_mda, years = 70) {
  times <- seq(0, years * 365.25, by = 1)
  out <- solve_model(times, state0, pars_no_mda)
  out <- add_output_names(out, ctx$age_groups)

  eq_state <- as.numeric(out[nrow(out), -1, drop = TRUE])
  eq_summary <- summarise_output(out, ctx)

  list(out = out, state = eq_state, summary = eq_summary)
}

run_single_scenario <- function(ctx, state, pars, years, strategy = c("none", "annual", "biannual"), organism = "S. pneumoniae") {
  strategy <- match.arg(strategy)
  pars_run <- pars

  if (strategy == "none") {
    pars_run <- make_no_mda_parameters(pars_run)
  } else if (strategy == "annual") {
    pars_run$use_mda <- TRUE
    pars_run$mda_start_times <- annual_mda_starts(years)
  } else if (strategy == "biannual") {
    pars_run$use_mda <- TRUE
    pars_run$mda_start_times <- biannual_mda_starts(years)
  }

  times <- seq(0, years * 365.25, by = 1)
  out <- solve_model(times, state, pars_run)
  out <- add_output_names(out, ctx$age_groups)

  scenario_label <- sprintf("%s_%sY", strategy, years)

  list(
    out = out,
    summary = summarise_output(out, ctx),
    under5 = extract_under5_deaths(out),
    age_distribution = final_age_distribution(
      out, ctx,
      organism = organism,
      scenario = scenario_label
    )
  )
}

run_core_scenarios <- function(ctx, equilibrium_state, pars, years = c(1, 5, 10, 20), organism = "S. pneumoniae") {
  results <- list()
  for (yy in years) {
    results[[paste0("NoMDA_", yy, "Y")]] <- run_single_scenario(ctx, equilibrium_state, pars, yy, "none", organism = organism)
    results[[paste0("Annual_", yy, "Y")]] <- run_single_scenario(ctx, equilibrium_state, pars, yy, "annual", organism = organism)
    results[[paste0("Biannual_", yy, "Y")]] <- run_single_scenario(ctx, equilibrium_state, pars, yy, "biannual", organism = organism)
  }
  results
}

# -----------------------------
# Example main run
# -----------------------------

data_dir <- "."

ctx <- load_tanzania_context(data_dir = data_dir)
pars <- make_spneumo_parameters(ctx)
pars_no_mda <- make_no_mda_parameters(pars)

state0 <- make_initial_state(
  ctx$population_by_age,
  init_fractions = c(X = 0.75, S = 0.20, R = 0.05, Sr = 0, Rs = 0)
)

message("Running equilibrium for S. pneumoniae scaffold...")
eq <- run_equilibrium(ctx, state0, pars_no_mda, years = 70)

# Check stability over the last 10 years
total_pop_last <- tail(eq$summary$total_population, 3650)
message(
  "Max |delta total pop| over final 10 years of equilibrium: ",
  signif(max(abs(diff(total_pop_last))), 4)
)

message("Running S. pneumoniae scenarios...")
scenario_results <- run_core_scenarios(ctx, eq$state, pars, years = c(1, 5, 10, 20), organism = "S. pneumoniae")

# Example compact summary table
# For this matched scaffold, mortality outputs remain on the table but
# default to zero until theta/amrd_rate are calibrated or externally set.
scenario_table <- bind_rows(lapply(names(scenario_results), function(nm) {
  ss <- scenario_results[[nm]]$summary
  u5 <- scenario_results[[nm]]$under5
  data.frame(
    scenario = nm,
    final_carriage_prevalence_pct = tail(ss$carriage_prevalence_pct, 1),
    final_resistance_among_carriers_pct = tail(ss$resistance_among_carriers_pct, 1),
    final_resistant_prevalence_pct = tail(ss$resistant_prevalence_pct, 1),
    cumulative_resistant_incidence = tail(ss$cum_resistant_incidence, 1) - ss$cum_resistant_incidence[1],
    under5_total_deaths = u5$deaths_total,
    under5_amr_deaths = u5$amr_deaths_total
  )
}))

print(scenario_table)
