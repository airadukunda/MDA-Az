# ============================================================================
# Rewritten Tanzania MDA-AMR model
# Shared scaffold for:
#   1) Escherichia coli
#   2) Streptococcus pneumoniae
#
# This script is a cleaned and modular rewrite of the original Tanzania model.
# It keeps the same age-structured ODE scaffold so the E. coli and S. pneumoniae
# analyses can sit side-by-side in one report.
#
# -----------------------------------------------------------------------------
# CORRECTIONS MADE RELATIVE TO THE ORIGINAL SCRIPT
# -----------------------------------------------------------------------------
# CORRECTION 1: Fixes age targeting for MDA by using *ages* mapped to R indices,
#               rather than 0-based R indexing such as 0:5.
# CORRECTION 2: Separates routine antibiotic pressure from MDA pressure:
#               a_t   = routine + MDA pulse
#               a.C_t = routine + MDA pulse
#               so background antibiotic pressure is not restricted to targeted ages.
# CORRECTION 3: Uses m_contact from the parameter list rather than the global
#               contact matrix inside the ODE.
# CORRECTION 4: Uses mda_cov from the parameter list rather than a hard-coded 0.8.
# CORRECTION 5: Removes/repurposes unused parameters. r_mda is removed.
#               kappa is now optional within-host emergence under treatment.
# CORRECTION 6: Corrects scenario timing so a "50-year" scenario is actually 50 years.
# CORRECTION 7: Uses output summaries based on the correct column offsets (+1 for time).
# CORRECTION 8: Uses cumulative state variables directly for cumulative outcomes,
#               rather than cumsum() over state prevalence.
# CORRECTION 9: Replaces fragile repeated scenario code with reusable functions.
# CORRECTION 10: Defines observation-style outputs that are closer to trial data:
#                - carriage / isolation prevalence
#                - resistance among carriers / isolates
#                - resistant prevalence in the total population
#
# -----------------------------------------------------------------------------
# S. PNEUMONIAE ADAPTATION MADE HERE
# -----------------------------------------------------------------------------
# S.PNEUMO 1: Keeps the same compartments (X, S, R, Sr, Rs, D, CumR, AMRD),
#             but interprets them as pneumococcal carriage states.
# S.PNEUMO 2: Uses the same MDA scaffold and age-structured contact structure.
# S.PNEUMO 3: Uses an organism-specific parameter preset and initial conditions.
# S.PNEUMO 4: Leaves mortality / AMR-death effects at zero by default for
#             pneumococcus, so the model behaves as a carriage-AMR model unless
#             you later calibrate a disease component.
# S.PNEUMO 5: Keeps outputs parallel to E. coli for a single combined report.
#
# Calibration is intentionally left for a later step.
# ============================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("deSolve", "ggplot2")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Please install required packages before running this script: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

# ---- Data loading -------------------------------------------------------------
load_tanzania_inputs <- function(
  data_dir,
  country = "United Republic of Tanzania",
  population_file = "Population_Afro_2023_1yearage.csv",
  births_file = "3.U.1.Birth_1year_Afro.csv",
  mortality_file = "3.U.1.AFRO_mortality_by_age_group_1yearage.csv",
  contact_file = "3.U.1.contact_Tanzania_1y.csv",
  contact_divisor = 25
) {
  pop <- read.csv(file.path(data_dir, population_file), stringsAsFactors = FALSE)
  births <- read.csv(file.path(data_dir, births_file), stringsAsFactors = FALSE)
  mortality <- read.csv(file.path(data_dir, mortality_file), stringsAsFactors = FALSE)
  contact <- as.matrix(read.csv(file.path(data_dir, contact_file), check.names = FALSE))

  pop <- subset(pop, Year %in% c(2023, "2023") & Country == country)
  births <- subset(births, Year %in% c(2023, "2023") & Country == country)
  mortality <- subset(mortality, Year %in% c(2023, "2023") & Country == country)

  if (nrow(pop) == 0) stop("No population rows found for Tanzania in 2023.")
  if (nrow(births) == 0) stop("No birth rows found for Tanzania in 2023.")
  if (nrow(mortality) == 0) stop("No mortality rows found for Tanzania in 2023.")

  # Population in persons
  pop$Population_age <- pop$Population_age * 1000
  pop$Annual_population <- pop$Annual_population * 1000

  n_age <- nrow(pop)
  age_groups <- if (n_age == 101) c(as.character(0:99), "100+") else as.character(seq_len(n_age) - 1)
  if (n_age != 101) {
    warning("Population table does not have 101 rows. Using row count to define age structure.")
  }

  # Ageing matrix: 1-year age bands, daily rate
  dd <- rep(1, n_age)
  ageing <- t(diff(diag(dd), lag = 1) / 365.25)
  ageing <- cbind(ageing, rep(0, n_age))

  # Births and mortality converted to per-person per-day rates
  births[, 5] <- 1000 * births[, 5] / (pop[, 5] * 365.25)
  mortality[, 5] <- 1000 * mortality[, 5] / (pop[, 5] * 365.25)
  mort <- mortality[, 5]

  if (!all(dim(contact) == c(n_age, n_age))) {
    stop("Contact matrix dimensions do not match the number of age groups.")
  }
  contact <- contact / contact_divisor
  rownames(contact) <- age_groups
  colnames(contact) <- age_groups

  list(
    country = country,
    data_dir = data_dir,
    pop = pop,
    births = births,
    mortality = mortality,
    mort = mort,
    contact = contact,
    ageing = ageing,
    age_groups = age_groups,
    n_age = n_age
  )
}

# ---- Indices and helpers ------------------------------------------------------
make_indices <- function(n_age) {
  list(
    X = 1:n_age,
    S = (1 * n_age + 1):(2 * n_age),
    R = (2 * n_age + 1):(3 * n_age),
    Sr = (3 * n_age + 1):(4 * n_age),
    Rs = (4 * n_age + 1):(5 * n_age),
    D = (5 * n_age + 1):(6 * n_age),
    CumR = (6 * n_age + 1):(7 * n_age),
    AMRD = (7 * n_age + 1):(8 * n_age)
  )
}

get_target_indices <- function(age_groups, target_ages) {
  idx <- match(as.character(target_ages), age_groups)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) {
    warning("No MDA target ages were matched to the age groups.")
  }
  idx
}

mda_active <- function(time, mda_starts, duration) {
  if (length(mda_starts) == 0) {
    return(FALSE)
  }
  any(vapply(mda_starts, function(start) time >= start && time < (start + duration), logical(1)))
}

build_mda_start_times <- function(
  years,
  frequency = c("none", "annual", "biannual"),
  stop_after_years = years
) {
  frequency <- match.arg(frequency)
  if (frequency == "none" || stop_after_years <= 0) {
    return(numeric(0))
  }

  interval_days <- switch(frequency,
    annual = 365.25,
    biannual = 365.25 / 2,
    none = Inf
  )

  max_time <- stop_after_years * 365.25
  seq(0, max_time, by = interval_days)
}

make_initial_state <- function(population_by_age, fractions) {
  required_names <- c("X", "S", "R", "Sr", "Rs")
  if (!all(required_names %in% names(fractions))) {
    stop("fractions must contain names: ", paste(required_names, collapse = ", "))
  }
  if (abs(sum(fractions[required_names]) - 1) > 1e-8) {
    stop("Initial live-state fractions must sum to 1.")
  }

  initX <- fractions[["X"]] * population_by_age
  initS <- fractions[["S"]] * population_by_age
  initR <- fractions[["R"]] * population_by_age
  initSr <- fractions[["Sr"]] * population_by_age
  initRs <- fractions[["Rs"]] * population_by_age
  initD <- 0 * population_by_age
  initCumR <- 0 * population_by_age
  initAMRD <- 0 * population_by_age

  c(initX, initS, initR, initSr, initRs, initD, initCumR, initAMRD)
}

# ---- Shared ODE model ---------------------------------------------------------
# This is the common scaffold for both E. coli and S. pneumoniae.
# The organism-specific adaptation comes from the parameter preset and the
# interpretation of the states, not from changing the ODE structure.

bacteria_odes_shared <- function(t, state, parameters) {
  with(as.list(parameters), {
    idx <- indices

    X <- state[idx$X]
    S <- state[idx$S]
    R <- state[idx$R]
    Sr <- state[idx$Sr]
    Rs <- state[idx$Rs]

    N <- pmax(X + S + R + Sr + Rs, 1e-12)
    S_tot <- S + Sr
    R_tot <- R + Rs

    # CORRECTION 3: use the contact matrix from parameters, not a global object.
    lambda.S <- as.vector(beta.S * (m_contact %*% (S_tot / N)))
    lambda.R <- as.vector(beta.R * (m_contact %*% (R_tot / N)))

    # CORRECTION 1 + 4: correct MDA age targeting and use parameterised coverage.
    target_idx <- get_target_indices(age_groups, mda_target_ages)
    azt <- rep(0, n_age)
    if (length(target_idx) > 0) {
      azt[target_idx] <- mda_cov
    }

    is_mda_on <- as.numeric(mda_active(t, mda_start_times, mda_duration))

    # CORRECTION 2: routine antibiotic pressure is applied to all ages, not only
    # the MDA-targeted ages. MDA adds an age-targeted pulse on top.
    tau <- a.use.eff * a.use
    tau.C <- a.use.eff * a.use
    a_t <- tau + is_mda_on * a * azt
    a.C_t <- tau.C + is_mda_on * a.C * azt

    # CORRECTION 5: kappa is now used as optional within-host emergence / selection
    # from susceptible carriage under treatment pressure.
    emergence <- kappa * a_t * S

    mort_eff <- mort
    if (is_mda_on > 0 && length(target_idx) > 0) {
      mort_eff[target_idx] <- mort[target_idx] * (1 - theta)
    }

    births <- rep(0, n_age)
    births[1] <- sum(mort_eff * N)

    dX <- births + (u.S + a_t) * S + u.R * R + u.C * (Sr + Rs) -
      (lambda.S + lambda.R) * X + ageing %*% X - mort_eff * X

    dS <- lambda.S * X - u.S * S - k * lambda.R * S - a_t * S - emergence +
      ageing %*% S - mort_eff * S

    dR <- lambda.R * X - u.R * R - k * lambda.S * R + a.C_t * (Sr + Rs) + emergence -
      amrd_rate * R + ageing %*% R - mort_eff * R

    dSr <- k * lambda.R * S - u.C * Sr - a.C_t * Sr +
      ageing %*% Sr - mort_eff * Sr

    dRs <- k * lambda.S * R - u.C * Rs - a.C_t * Rs +
      ageing %*% Rs - mort_eff * Rs

    dD <- mort_eff * (X + S + R + Sr + Rs)

    # CORRECTION 8: define cumulative resistant acquisition as a dedicated state.
    dCumR <- lambda.R * X + k * lambda.S * R + a.C_t * (Sr + Rs) + emergence
    dAMRD <- amrd_rate * R

    list(c(dX, dS, dR, dSr, dRs, dD, dCumR, dAMRD))
  })
}

solve_shared_model <- function(times, state, parameters) {
  parameters$beta.R <- parameters$beta.S * (1 - parameters$c)
  as.data.frame(deSolve::ode(y = state, times = times, func = bacteria_odes_shared, parms = parameters))
}

# ---- Summary functions --------------------------------------------------------
extract_compartment_totals <- function(out, indices) {
  list(
    X = rowSums(out[, indices$X + 1, drop = FALSE]),
    S = rowSums(out[, indices$S + 1, drop = FALSE]),
    R = rowSums(out[, indices$R + 1, drop = FALSE]),
    Sr = rowSums(out[, indices$Sr + 1, drop = FALSE]),
    Rs = rowSums(out[, indices$Rs + 1, drop = FALSE]),
    D = rowSums(out[, indices$D + 1, drop = FALSE]),
    CumR = rowSums(out[, indices$CumR + 1, drop = FALSE]),
    AMRD = rowSums(out[, indices$AMRD + 1, drop = FALSE])
  )
}

summarise_solution <- function(out, indices, organism, scenario, horizon_years, stop_after_years) {
  totals <- extract_compartment_totals(out, indices)
  total_pop <- totals$X + totals$S + totals$R + totals$Sr + totals$Rs
  carriage <- totals$S + totals$R + totals$Sr + totals$Rs
  resistant <- totals$R + totals$Rs
  sensitive <- totals$S + totals$Sr

  data.frame(
    organism = organism,
    scenario = scenario,
    horizon_years = horizon_years,
    mda_stop_years = stop_after_years,
    time = out$time,
    total_population = total_pop,
    carriage_prevalence_pct = 100 * carriage / pmax(total_pop, 1e-12),
    resistance_among_carriers_pct = 100 * resistant / pmax(carriage, 1e-12),
    resistant_prevalence_pct = 100 * resistant / pmax(total_pop, 1e-12),
    sensitive_prevalence_pct = 100 * sensitive / pmax(total_pop, 1e-12),
    cumulative_resistant_acquisition = totals$CumR,
    cumulative_deaths = totals$D,
    cumulative_amr_deaths = totals$AMRD
  )
}

final_age_distribution <- function(out, indices, age_groups, organism, scenario) {
  final_state <- as.numeric(out[nrow(out), -1, drop = TRUE]) # drop time column

  X_final <- final_state[indices$X]
  S_final <- final_state[indices$S]
  R_final <- final_state[indices$R]
  Sr_final <- final_state[indices$Sr]
  Rs_final <- final_state[indices$Rs]

  live_totals <- X_final + S_final + R_final + Sr_final + Rs_final
  carriage <- S_final + R_final + Sr_final + Rs_final
  resistant <- R_final + Rs_final

  data.frame(
    organism = organism,
    scenario = scenario,
    age_group = age_groups,
    population = live_totals,
    carriage_prevalence_pct = 100 * carriage / pmax(live_totals, 1e-12),
    resistance_among_carriers_pct = 100 * resistant / pmax(carriage, 1e-12),
    resistant_prevalence_pct = 100 * resistant / pmax(live_totals, 1e-12)
  )
}

# ---- Organism presets ---------------------------------------------------------
make_ecoli_preset <- function(context) {
  list(
    organism = "E. coli",

    # Core transmission / clearance parameters
    beta.S = 5 * 12 / 365.25,
    u.S = 1 * 12 / 365.25,
    u.R = 1 * 12 / 365.25,
    u.C = 1 * 12 / 365.25,

    # Drug-related parameters
    a = 0.16,
    a.C = 0.16,
    a.use = 0.06,
    a.use.eff = 0.05,

    # Competition / resistance parameters
    k = 0.5,
    c = 0.01,
    kappa = 0.00,

    # Demography / contacts
    m_contact = context$contact,
    mort = context$mort,
    ageing = context$ageing,
    age_groups = context$age_groups,
    n_age = context$n_age,
    indices = make_indices(context$n_age),

    # MDA settings
    # CORRECTION 1: target ages are specified as ages, then mapped safely.
    # For children 1-59 months in 1-year age bands, 0:4 is usually the closest.
    mda_target_ages = 0:4,
    mda_cov = 0.60,
    mda_duration = 30,
    mda_start_times = numeric(0),

    # Mortality effect retained to preserve the original scaffold.
    theta = 0.00,
    amrd_rate = 27.3 / (100000 * 365),

    # Initial state fractions
    init_fractions = c(X = 0.95, S = 0.025, R = 0.025, Sr = 0.00, Rs = 0.00)
  )
}

make_spneumo_preset <- function(context) {
  list(
    organism = "S. pneumoniae",

    # S.PNEUMO ADAPTATION 1:
    # Same parameter names and same scaffold, but these are placeholders for a
    # pneumococcal carriage model and should later be calibrated to carriage data.
    beta.S = 4.0 * 12 / 365.25,
    u.S = 1.5 * 12 / 365.25,
    u.R = 1.5 * 12 / 365.25,
    u.C = 1.5 * 12 / 365.25,

    # S.PNEUMO ADAPTATION 2:
    # Same azithromycin-MDA mechanism, reinterpreted as selection on carriage.
    a = 0.16,
    a.C = 0.16,
    a.use = 0.06,
    a.use.eff = 0.05,
    k = 0.5,
    c = 0.02,
    kappa = 0.00,
    m_contact = context$contact,
    mort = context$mort,
    ageing = context$ageing,
    age_groups = context$age_groups,
    n_age = context$n_age,
    indices = make_indices(context$n_age),

    # Same under-5 targeting scaffold for comparability across organisms.
    mda_target_ages = 0:4,
    mda_cov = 0.60,
    mda_duration = 30,
    mda_start_times = numeric(0),

    # S.PNEUMO ADAPTATION 3:
    # Default to a carriage-AMR model: keep mortality effects switched off until
    # a disease layer or empirical mortality mapping is introduced.
    theta = 0.00,
    amrd_rate = 0.00,

    # S.PNEUMO ADAPTATION 4:
    # Higher baseline carriage placeholder than E. coli.
    init_fractions = c(X = 0.80, S = 0.15, R = 0.05, Sr = 0.00, Rs = 0.00)
  )
}

# ---- Scenario runner ----------------------------------------------------------
run_equilibrium <- function(parameters, years = 200) {
  times <- seq(0, years * 365.25, by = 1)
  pop <- parameters$population_by_age
  state0 <- make_initial_state(pop, parameters$init_fractions)

  eq_pars <- parameters
  eq_pars$mda_start_times <- numeric(0)
  eq_pars$theta <- 0

  out_eq <- solve_shared_model(times, state0, eq_pars)
  as.numeric(out_eq[nrow(out_eq), -1])
}

run_one_scenario <- function(parameters, state, years, frequency, stop_after_years = years) {
  pars <- parameters
  pars$mda_start_times <- build_mda_start_times(
    years = years,
    frequency = frequency,
    stop_after_years = stop_after_years
  )

  times <- seq(0, years * 365.25, by = 1)
  solve_shared_model(times, state, pars)
}

run_scenario_set <- function(parameters, horizons = c(1, 5, 10, 50), stop_after_years = NULL) {
  pop <- parameters$population_by_age
  eq_state <- run_equilibrium(parameters, years = 200)

  summaries <- list()
  finals <- list()
  raw <- list()
  counter <- 1

  for (h in horizons) {
    stop_y <- if (is.null(stop_after_years)) h else min(stop_after_years, h)

    scenario_specs <- list(
      list(name = "No MDA", frequency = "none", stop_after = 0),
      list(name = "Annual MDA", frequency = "annual", stop_after = stop_y),
      list(name = "Biannual MDA", frequency = "biannual", stop_after = stop_y)
    )

    for (spec in scenario_specs) {
      out <- run_one_scenario(
        parameters = parameters,
        state = eq_state,
        years = h,
        frequency = spec$frequency,
        stop_after_years = spec$stop_after
      )

      key <- paste(parameters$organism, h, spec$name, sep = " | ")
      raw[[key]] <- out

      summaries[[counter]] <- summarise_solution(
        out = out,
        indices = parameters$indices,
        organism = parameters$organism,
        scenario = spec$name,
        horizon_years = h,
        stop_after_years = spec$stop_after
      )

      finals[[counter]] <- final_age_distribution(
        out = out,
        indices = parameters$indices,
        age_groups = parameters$age_groups,
        organism = parameters$organism,
        scenario = paste0(spec$name, " (", h, "y)")
      )

      counter <- counter + 1
    }
  }

  list(
    equilibrium_state = eq_state,
    raw = raw,
    summary = do.call(rbind, summaries),
    final_age = do.call(rbind, finals)
  )
}

# ---- Plot helper --------------------------------------------------------------
plot_metric_by_scenario <- function(summary_df, organism_name, horizon_years_value, metric = c(
                                      "resistance_among_carriers_pct",
                                      "carriage_prevalence_pct",
                                      "resistant_prevalence_pct"
                                    )) {
  metric <- match.arg(metric)
  dat <- subset(summary_df, organism == organism_name & horizon_years == horizon_years_value)
  if (nrow(dat) == 0) stop("No rows found for the requested organism / horizon.")

  dat$metric_value <- dat[[metric]]
  ggplot2::ggplot(dat, ggplot2::aes(x = time / 365.25, y = metric_value, color = scenario)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = "Time (years)",
      y = metric,
      title = paste0(organism_name, ": ", gsub("_", " ", metric), " over time"),
      color = "Scenario"
    ) +
    ggplot2::theme_classic(base_size = 12)
}

# ---- Main wrapper -------------------------------------------------------------
prepare_organism <- function(context, preset_fn) {
  pars <- preset_fn(context)
  pars$population_by_age <- context$pop[, 5]
  pars
}

run_dual_analysis <- function(data_dir, stop_after_years = NULL) {
  context <- load_tanzania_inputs(data_dir = data_dir)

  ecoli_pars <- prepare_organism(context, make_ecoli_preset)
  spneumo_pars <- prepare_organism(context, make_spneumo_preset)

  ecoli_res <- run_scenario_set(ecoli_pars, horizons = c(1, 5, 10, 50), stop_after_years = stop_after_years)
  spneumo_res <- run_scenario_set(spneumo_pars, horizons = c(1, 5, 10, 50), stop_after_years = stop_after_years)

  list(
    context = context,
    ecoli = ecoli_res,
    spneumo = spneumo_res,
    combined_summary = rbind(ecoli_res$summary, spneumo_res$summary),
    combined_final_age = rbind(ecoli_res$final_age, spneumo_res$final_age)
  )
}

# ---- Example usage ------------------------------------------------------------
# Set this to the folder containing:
#   Population_emro_2023_1yearage.csv
#   3.U.1.Birth_1year_Afro.csv
#   3.U.1.EMRO_mortality_by_age_group_1yearage.csv
#   3.U.1.contact_Tanzania_1y.csv
#
# Example:

results <- run_dual_analysis(data_dir = ".", stop_after_years = 10)
head(results$combined_summary)
tail(subset(
  results$combined_summary,
  organism == "E. coli" &
    scenario == "Annual MDA" &
    horizon_years == 10
)[, c("time", "resistance_among_carriers_pct")])

p1 <- plot_metric_by_scenario(results$combined_summary, "E. coli", 10,
  metric = "resistance_among_carriers_pct"
)
p2 <- plot_metric_by_scenario(results$combined_summary, "S. pneumoniae", 10,
  metric = "resistance_among_carriers_pct"
)
print(p1)
print(p2)

# Write outputs if needed:
write.csv(results$combined_summary, "dual_model_summary.csv", row.names = FALSE)
write.csv(results$combined_final_age, "dual_model_final_age.csv", row.names = FALSE)
