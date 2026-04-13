suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(scales)
})

# ============================================================
# Numeric validation: original vs cleaned E. coli script
# Assumes:
#   1) davies model 2_5_Tanzania_final.R has already been run
#   2) cleaned_core_stationary_mda_model.R is available
# ============================================================

# -----------------------------
# 0) Keep original objects safe
# -----------------------------
orig_objects <- mget(ls(), envir = .GlobalEnv)

# -----------------------------
# 1) Source cleaned script into its own environment
# -----------------------------
clean_env <- new.env(parent = globalenv())
sys.source("~/Downloads/cleaned_core_stationary_mda_model.R", envir = clean_env)

# -----------------------------
# 2) Validation context loader using the SAME AFRO files
#    as the original script
# -----------------------------
clean_env$load_tanzania_context_validation <- function(data_dir = ".") {
    pop_path <- file.path(data_dir, "Population_Afro_2023_1yearage.csv")
    mort_path <- file.path(data_dir, "3.U.1.AFRO_mortality_by_age_group_1yearage.csv")
    contact_path <- file.path(data_dir, "3.U.1.contact_Tanzania_1y.csv")

    pop_df <- read.csv(pop_path, check.names = FALSE)
    mort_df <- read.csv(mort_path, check.names = FALSE)
    contact_df <- read.csv(contact_path, check.names = FALSE)

    pop_tz <- pop_df %>%
        filter(Country == "United Republic of Tanzania", Year %in% c(2023, "2023"))

    pop_vec <- clean_env$safe_numeric_column(pop_tz, preferred_name = "Population_age", fallback_index = 5)
    pop_vec <- pop_vec * 1000

    n_age <- length(pop_vec)
    age_groups <- if (n_age == 101) c(as.character(0:99), "100+") else as.character(seq_len(n_age) - 1)

    mort_tz <- mort_df %>%
        filter(Country == "United Republic of Tanzania", Year %in% c(2023, "2023"))

    mort_raw <- clean_env$safe_numeric_column(mort_tz, preferred_name = "Percentage", fallback_index = 5)
    mort_per_day <- 1000 * mort_raw / (pop_vec * 365.25)

    contact_mat <- as.matrix(contact_df)
    suppressWarnings(storage.mode(contact_mat) <- "numeric")

    if (ncol(contact_mat) == n_age + 1) contact_mat <- contact_mat[, -1, drop = FALSE]
    if (nrow(contact_mat) == n_age + 1) contact_mat <- contact_mat[-1, , drop = FALSE]

    contact_mat <- contact_mat / 25
    colnames(contact_mat) <- age_groups
    rownames(contact_mat) <- age_groups

    list(
        data_dir = data_dir,
        n_age = n_age,
        age_groups = age_groups,
        population_by_age = pop_vec,
        mort = mort_per_day,
        ageing = clean_env$make_ageing_matrix(n_age),
        m_contact = contact_mat,
        indices = clean_env$make_indices(n_age)
    )
}

# -----------------------------
# 3) Run cleaned model with same setup
# -----------------------------
data_dir <- "."

ctx_clean <- clean_env$load_tanzania_context_validation(data_dir = data_dir)
pars_clean <- clean_env$make_ecoli_parameters(ctx_clean)
pars_no_mda_clean <- clean_env$make_no_mda_parameters(pars_clean)

state0_clean <- clean_env$make_initial_state(
    ctx_clean$population_by_age,
    init_fractions = c(X = 0.95, S = 0.025, R = 0.025, Sr = 0, Rs = 0)
)

eq_clean <- clean_env$run_equilibrium(ctx_clean, state0_clean, pars_no_mda_clean, years = 70)

scenario_results_clean <- clean_env$run_core_scenarios(
    ctx_clean,
    equilibrium_state = eq_clean$state,
    pars = pars_clean,
    years = c(1, 5, 10, 20)
)

# -----------------------------
# 4) Helpers for original-script outputs
# -----------------------------
compute_original_summary <- function(out_df, Xindex, Sindex, Rindex, Srindex, Rsindex, CumIncRindex) {
    X_total <- rowSums(out_df[, Xindex + 1, drop = FALSE])
    S_total <- rowSums(out_df[, Sindex + 1, drop = FALSE])
    R_total <- rowSums(out_df[, Rindex + 1, drop = FALSE])
    Sr_total <- rowSums(out_df[, Srindex + 1, drop = FALSE])
    Rs_total <- rowSums(out_df[, Rsindex + 1, drop = FALSE])
    CumIncR_total <- rowSums(out_df[, CumIncRindex + 1, drop = FALSE])

    live_total <- X_total + S_total + R_total + Sr_total + Rs_total
    colonised_total <- S_total + R_total + Sr_total + Rs_total
    resistant_total <- R_total + Rs_total

    data.frame(
        time = out_df[, 1],
        carriage_prevalence_pct = 100 * colonised_total / pmax(live_total, 1e-12),
        resistance_among_carriers_pct = 100 * resistant_total / pmax(colonised_total, 1e-12),
        resistant_prevalence_pct = 100 * resistant_total / pmax(live_total, 1e-12),
        cum_resistant_incidence = CumIncR_total,
        resistant_incidence = c(0, diff(CumIncR_total))
    )
}

extract_original_age_distribution <- function(out_df, age_groups, Xindex, Sindex, Rindex, Srindex, Rsindex) {
    final_state <- as.numeric(out_df[nrow(out_df), -1, drop = TRUE])

    X_final <- final_state[Xindex]
    S_final <- final_state[Sindex]
    R_final <- final_state[Rindex]
    Sr_final <- final_state[Srindex]
    Rs_final <- final_state[Rsindex]

    live_totals <- X_final + S_final + R_final + Sr_final + Rs_final
    carriage <- S_final + R_final + Sr_final + Rs_final
    resistant <- R_final + Rs_final

    data.frame(
        age_group = age_groups,
        carriage_prevalence_pct = 100 * carriage / pmax(live_totals, 1e-12),
        resistance_among_carriers_pct = 100 * resistant / pmax(carriage, 1e-12),
        resistant_prevalence_pct = 100 * resistant / pmax(live_totals, 1e-12)
    )
}

# -----------------------------
# 5) Build matched original summaries
# -----------------------------
age_groups_orig <- orig_objects$age_groups
Xindex_orig <- orig_objects$Xindex
Sindex_orig <- orig_objects$Sindex
Rindex_orig <- orig_objects$Rindex
Srindex_orig <- orig_objects$Srindex
Rsindex_orig <- orig_objects$Rsindex
CumIncRindex_orig <- orig_objects$CumIncRindex

original_runs <- list(
    NoMDA_1Y = orig_objects$out_1_b_Tanzania,
    Annual_1Y = orig_objects$out_1_a_Tanzania,
    Biannual_1Y = orig_objects$out_1_c_Tanzania,
    NoMDA_5Y = orig_objects$out_5_b_Tanzania,
    Annual_5Y = orig_objects$out_5_a_Tanzania,
    Biannual_5Y = orig_objects$out_5_c_Tanzania,
    NoMDA_10Y = orig_objects$out_10_b_Tanzania,
    Annual_10Y = orig_objects$out_10_a_Tanzania,
    Biannual_10Y = orig_objects$out_10_c_Tanzania,
    NoMDA_20Y = orig_objects$out_50_b_Tanzania,
    Annual_20Y = orig_objects$out_50_a_Tanzania,
    Biannual_20Y = orig_objects$out_50_c_Tanzania
)

original_summaries <- lapply(original_runs, function(out_df) {
    compute_original_summary(
        out_df,
        Xindex = Xindex_orig,
        Sindex = Sindex_orig,
        Rindex = Rindex_orig,
        Srindex = Srindex_orig,
        Rsindex = Rsindex_orig,
        CumIncRindex = CumIncRindex_orig
    )
})

original_age_profiles <- lapply(original_runs, function(out_df) {
    extract_original_age_distribution(
        out_df,
        age_groups = age_groups_orig,
        Xindex = Xindex_orig,
        Sindex = Sindex_orig,
        Rindex = Rindex_orig,
        Srindex = Srindex_orig,
        Rsindex = Rsindex_orig
    )
})

# -----------------------------
# 6) Build matched cleaned summaries
# -----------------------------
clean_runs <- list(
    NoMDA_1Y = scenario_results_clean$NoMDA_1Y$out,
    Annual_1Y = scenario_results_clean$Annual_1Y$out,
    Biannual_1Y = scenario_results_clean$Biannual_1Y$out,
    NoMDA_5Y = scenario_results_clean$NoMDA_5Y$out,
    Annual_5Y = scenario_results_clean$Annual_5Y$out,
    Biannual_5Y = scenario_results_clean$Biannual_5Y$out,
    NoMDA_10Y = scenario_results_clean$NoMDA_10Y$out,
    Annual_10Y = scenario_results_clean$Annual_10Y$out,
    Biannual_10Y = scenario_results_clean$Biannual_10Y$out,
    NoMDA_20Y = scenario_results_clean$NoMDA_20Y$out,
    Annual_20Y = scenario_results_clean$Annual_20Y$out,
    Biannual_20Y = scenario_results_clean$Biannual_20Y$out
)

clean_summaries <- lapply(clean_runs, function(out_df) {
    clean_env$summarise_output(out_df, ctx_clean)
})

clean_age_profiles <- lapply(names(clean_runs), function(nm) {
    dd <- clean_env$final_age_distribution(clean_runs[[nm]], ctx_clean, scenario = nm)
    dd$scenario <- nm
    dd
})
names(clean_age_profiles) <- names(clean_runs)

# -----------------------------
# 7) Endpoint comparison table
# -----------------------------
endpoint_comparison <- bind_rows(lapply(names(original_summaries), function(nm) {
    orig <- original_summaries[[nm]]
    cln <- clean_summaries[[nm]]

    data.frame(
        scenario = nm,
        original_final_carriage = tail(orig$carriage_prevalence_pct, 1),
        cleaned_final_carriage = tail(cln$carriage_prevalence_pct, 1),
        diff_final_carriage = tail(cln$carriage_prevalence_pct, 1) - tail(orig$carriage_prevalence_pct, 1),
        original_final_resistance_among_carriers = tail(orig$resistance_among_carriers_pct, 1),
        cleaned_final_resistance_among_carriers = tail(cln$resistance_among_carriers_pct, 1),
        diff_final_resistance_among_carriers = tail(cln$resistance_among_carriers_pct, 1) - tail(orig$resistance_among_carriers_pct, 1),
        original_final_resistant_prevalence = tail(orig$resistant_prevalence_pct, 1),
        cleaned_final_resistant_prevalence = tail(cln$resistant_prevalence_pct, 1),
        diff_final_resistant_prevalence = tail(cln$resistant_prevalence_pct, 1) - tail(orig$resistant_prevalence_pct, 1),
        original_cum_res_inc = tail(orig$cum_resistant_incidence, 1) - orig$cum_resistant_incidence[1],
        cleaned_cum_res_inc = tail(cln$cum_resistant_incidence, 1) - cln$cum_resistant_incidence[1],
        diff_cum_res_inc = (tail(cln$cum_resistant_incidence, 1) - cln$cum_resistant_incidence[1]) -
            (tail(orig$cum_resistant_incidence, 1) - orig$cum_resistant_incidence[1])
    )
}))

print(endpoint_comparison)

write.csv(endpoint_comparison, "validation_endpoint_comparison.csv", row.names = FALSE)

# -----------------------------
# 8) Trajectory comparison table
# -----------------------------
trajectory_comparison <- bind_rows(lapply(names(original_summaries), function(nm) {
    orig <- original_summaries[[nm]]
    cln <- clean_summaries[[nm]]

    merged <- full_join(
        orig %>% transmute(
            time,
            orig_resistance_among_carriers = resistance_among_carriers_pct,
            orig_carriage = carriage_prevalence_pct,
            orig_incidence = resistant_incidence
        ),
        cln %>% transmute(
            time,
            clean_resistance_among_carriers = resistance_among_carriers_pct,
            clean_carriage = carriage_prevalence_pct,
            clean_incidence = resistant_incidence
        ),
        by = "time"
    ) %>%
        mutate(
            scenario = nm,
            diff_resistance_among_carriers = clean_resistance_among_carriers - orig_resistance_among_carriers,
            diff_carriage = clean_carriage - orig_carriage,
            diff_incidence = clean_incidence - orig_incidence
        )

    merged
}))

write.csv(trajectory_comparison, "validation_trajectory_comparison.csv", row.names = FALSE)

# -----------------------------
# 9) Quick diagnostics
# -----------------------------
validation_summary <- trajectory_comparison %>%
    group_by(scenario) %>%
    summarise(
        max_abs_diff_resistance = max(abs(diff_resistance_among_carriers), na.rm = TRUE),
        mean_abs_diff_resistance = mean(abs(diff_resistance_among_carriers), na.rm = TRUE),
        max_abs_diff_carriage = max(abs(diff_carriage), na.rm = TRUE),
        mean_abs_diff_carriage = mean(abs(diff_carriage), na.rm = TRUE),
        max_abs_diff_incidence = max(abs(diff_incidence), na.rm = TRUE),
        mean_abs_diff_incidence = mean(abs(diff_incidence), na.rm = TRUE),
        .groups = "drop"
    )

print(validation_summary)

write.csv(validation_summary, "validation_summary_stats.csv", row.names = FALSE)

# -----------------------------
# 10) Plots of cleaned - original differences
# -----------------------------
plot_diff_resistance <- ggplot(
    trajectory_comparison,
    aes(x = time / 365.25, y = diff_resistance_among_carriers)
) +
    geom_line() +
    facet_wrap(~scenario, scales = "free_x") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(
        x = "Time (years)",
        y = "Cleaned - original resistance among carriers (%)",
        title = "Trajectory difference: resistance among carriers"
    ) +
    theme_minimal()

plot_diff_carriage <- ggplot(
    trajectory_comparison,
    aes(x = time / 365.25, y = diff_carriage)
) +
    geom_line() +
    facet_wrap(~scenario, scales = "free_x") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(
        x = "Time (years)",
        y = "Cleaned - original carriage prevalence (%)",
        title = "Trajectory difference: carriage prevalence"
    ) +
    theme_minimal()

plot_diff_incidence <- ggplot(
    trajectory_comparison,
    aes(x = time / 365.25, y = diff_incidence)
) +
    geom_line() +
    facet_wrap(~scenario, scales = "free_x") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(
        x = "Time (years)",
        y = "Cleaned - original resistant incidence",
        title = "Trajectory difference: resistant incidence"
    ) +
    theme_minimal()

print(plot_diff_resistance)
print(plot_diff_carriage)
print(plot_diff_incidence)

ggsave("validation_diff_resistance.png", plot_diff_resistance, width = 12, height = 8, dpi = 300)
ggsave("validation_diff_carriage.png", plot_diff_carriage, width = 12, height = 8, dpi = 300)
ggsave("validation_diff_incidence.png", plot_diff_incidence, width = 12, height = 8, dpi = 300)

# -----------------------------
# 11) Age-profile comparison
# -----------------------------
age_profile_comparison <- bind_rows(lapply(names(original_age_profiles), function(nm) {
    orig <- original_age_profiles[[nm]] %>%
        transmute(
            age_group,
            scenario = nm,
            source = "Original",
            carriage_prevalence_pct,
            resistance_among_carriers_pct
        )

    cln <- clean_age_profiles[[nm]] %>%
        transmute(
            age_group,
            scenario = nm,
            source = "Cleaned",
            carriage_prevalence_pct,
            resistance_among_carriers_pct
        )

    bind_rows(orig, cln)
}))

p_age_compare_carriage <- ggplot(
    age_profile_comparison,
    aes(x = age_group, y = carriage_prevalence_pct, fill = source)
) +
    geom_col(position = "dodge") +
    facet_wrap(~scenario, scales = "free_y") +
    labs(
        x = "Age group",
        y = "Carriage prevalence (%)",
        title = "Age-profile comparison: original vs cleaned"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

p_age_compare_resistance <- ggplot(
    age_profile_comparison,
    aes(x = age_group, y = resistance_among_carriers_pct, fill = source)
) +
    geom_col(position = "dodge") +
    facet_wrap(~scenario, scales = "free_y") +
    labs(
        x = "Age group",
        y = "Resistance among carriers (%)",
        title = "Age-profile comparison: original vs cleaned"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

print(p_age_compare_carriage)
print(p_age_compare_resistance)

ggsave("validation_age_compare_carriage.png", p_age_compare_carriage, width = 13, height = 8, dpi = 300)
ggsave("validation_age_compare_resistance.png", p_age_compare_resistance, width = 13, height = 8, dpi = 300)
