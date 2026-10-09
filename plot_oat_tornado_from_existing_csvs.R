# =============================================================================
# Plot OAT tornado figures from existing CSV outputs
# =============================================================================
#
# Use this after the calibrated sensitivity analyses have already been run on
# the cluster and the tornado CSV files have been produced. This script does not
# source any model scripts, does not solve ODEs, and does not recalibrate anything.
# It only reads the precomputed tornado CSVs and makes plots.
#
# By default it searches recursively from the current working directory for files
# with names such as:
#   ecoli_oat_tornado_annual_10y_all_ages.csv
#   spneumo_oat_tornado_biannual_10y_all_ages.csv
#
# You can restrict the search by setting OAT_CSV_DIR before running, e.g.
#   OAT_CSV_DIR=/path/to/cluster/outputs Rscript plot_oat_tornado_from_existing_csvs.R
#
# Outputs are written to OAT_PLOT_DIR if set, otherwise to:
#   outputs_combined_report_figures/from_existing_tornado_csvs
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(tibble)
})

input_root <- Sys.getenv("OAT_CSV_DIR", unset = ".")
output_dir <- Sys.getenv(
  "OAT_PLOT_DIR",
  unset = file.path("outputs_combined_report_figures", "from_existing_tornado_csvs")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

find_existing_csv <- function(filename, input_root = ".") {
  # First try exact path in input_root.
  direct <- file.path(input_root, filename)
  if (file.exists(direct)) {
    return(direct)
  }

  # Then search recursively.
  hits <- list.files(
    input_root,
    pattern = paste0("^", filename, "$"),
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(hits) == 0) {
    stop("Could not find required CSV: ", filename, " under ", normalizePath(input_root), call. = FALSE)
  }

  if (length(hits) > 1) {
    message("Multiple matches for ", filename, "; using: ", hits[[1]])
  }

  hits[[1]]
}

parameter_label <- function(x) {
  dplyr::recode(
    x,
    beta.S_start_multiplier = "Transmission",
    beta.S_multiplier = "Transmission",
    beta_multiplier = "Transmission",
    macrolide_use_multiplier = "Background macrolide use",
    background_macrolide_use_multiplier = "Background macrolide use",
    background_antibiotic_multiplier = "Background antibiotic use",
    background_antibiotic_use_multiplier = "Background antibiotic use",
    mda_effect_multiplier = "MDA treatment effect",
    mda_cov = "MDA coverage",
    mda_duration = "MDA duration",
    c = "Resistance fitness cost",
    k = "Co-colonisation efficiency",
    clearance_multiplier = "Carriage clearance",
    u.S_monthly_multiplier = "Sensitive carriage clearance",
    u.R_monthly_multiplier = "Resistant carriage clearance",
    u.C_monthly_multiplier = "Mixed-carriage clearance",
    u.S_multiplier = "Sensitive carriage clearance",
    u.R_multiplier = "Resistant carriage clearance",
    u.C_multiplier = "Mixed-carriage clearance",
    amrd_rate_multiplier = "AMR mortality rate",
    theta = "MDA mortality effect",
    .default = x
  )
}

pathogen_config <- tibble::tribble(
  ~pathogen,  ~pathogen_label,
  "ecoli",    "E. coli",
  "spneumo",  "S. pneumoniae",
  "saureus",  "S. aureus"
)

scenario_config <- tibble::tribble(
  ~scenario_key, ~scenario_label,
  "annual", "Annual MDA",
  "biannual", "Biannual MDA"
)

read_tornado_csv <- function(pathogen, pathogen_label, scenario_key, scenario_label) {
  filename <- paste0(pathogen, "_oat_tornado_", scenario_key, "_10y_all_ages.csv")
  csv_path <- find_existing_csv(filename, input_root)

  dat <- readr::read_csv(csv_path, show_col_types = FALSE)

  required_cols <- c("parameter", "min_effect", "max_effect", "range", "baseline_effect")
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
      scenario_key = scenario_key,
      scenario_label = scenario_label,
      parameter_label = parameter_label(parameter),
      source_file = csv_path,
      .before = 1
    )
}

tornado_all <- tidyr::crossing(pathogen_config, scenario_config) |>
  purrr::pmap_dfr(read_tornado_csv)

tornado_all <- tornado_all |>
  dplyr::filter(
    !(pathogen == "ecoli" & parameter_label %in% c(
      "AMR mortality rate",
      "MDA mortality effect"
    ))
  )

tornado_all <- tornado_all |>
  dplyr::filter(
    !(pathogen == "saureus" & parameter_label %in% c(
      "AMR mortality rate",
      "MDA mortality effect"
    ))
  )

# Drop completely uninformative rows where the range is exactly zero.
# Set DROP_ZERO_RANGE=FALSE if you want to keep these rows.
drop_zero_range <- isTRUE(as.logical(Sys.getenv("DROP_ZERO_RANGE", "TRUE")))
if (drop_zero_range) {
  tornado_all <- tornado_all |>
    dplyr::filter(is.na(range) | range > 0)
}

# Order facets and parameters. The free y-axis means each pathogen only shows
# parameters present in its own CSVs.
tornado_all <- tornado_all |>
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

parameter_order <- tornado_all |>
  dplyr::group_by(parameter_label) |>
  dplyr::summarise(mean_range = mean(abs(range), na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(mean_range) |>
  dplyr::pull(parameter_label)

tornado_all <- tornado_all |>
  dplyr::mutate(parameter_label = factor(parameter_label, levels = parameter_order))

plot_tornado <- function(dat,
                         title,
                         x_label = "Increase in resistance vs matched no-MDA scenario (percentage points)") {
  ggplot2::ggplot(dat, ggplot2::aes(y = parameter_label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = min_effect, xend = max_effect, yend = parameter_label),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(ggplot2::aes(x = min_effect), size = 1.8) +
    ggplot2::geom_point(ggplot2::aes(x = max_effect), size = 1.8) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = baseline_effect),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = "Parameter varied"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", colour = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(1.0, "lines")
    )
}

# Combined 6-panel figure.
combined_plot <- plot_tornado(
  tornado_all,
  title = "One-at-a-time sensitivity analysis on resistance increase after 10 years"
) +
  ggplot2::facet_grid(
    pathogen_label ~ scenario_label,
    scales = "free",
    space = "free_y"
  )

ggplot2::ggsave(
  file.path(output_dir, "combined_oat_tornado_from_existing_csvs.png"),
  combined_plot,
  width = 11,
  height = 8.5,
  dpi = 300
)

readr::write_csv(
  tornado_all,
  file.path(output_dir, "combined_oat_tornado_from_existing_csvs.csv")
)

# Individual pathogen/scenario plots.
for (this_pathogen in unique(tornado_all$pathogen)) {
  for (this_scenario in unique(tornado_all$scenario_key)) {
    dat <- tornado_all |>
      dplyr::filter(pathogen == this_pathogen, scenario_key == this_scenario)

    if (nrow(dat) == 0) next

    pathogen_label_i <- as.character(dat$pathogen_label[[1]])
    scenario_label_i <- as.character(dat$scenario_label[[1]])

    p <- plot_tornado(
      dat,
      title = paste0(pathogen_label_i, " sensitivity: ", scenario_label_i, " at 10 years")
    )

    ggplot2::ggsave(
      file.path(
        output_dir,
        paste0(this_pathogen, "_oat_tornado_", this_scenario, "_10y_all_ages.png")
      ),
      p,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
}

message("Plot-only OAT tornado figures written to: ", normalizePath(output_dir, mustWork = FALSE))
