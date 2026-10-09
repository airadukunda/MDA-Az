# Run one S. aureus recalibrated OAT sensitivity case.
# Intended for SLURM array jobs.

Sys.setenv(RUN_ALL_CASES = "FALSE")

source("saureus_oat_sensitivity_recalibrated_common_clearance.R")

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

if (is.na(task_id) || task_id < 1 || task_id > nrow(oat_specs)) {
  stop("Invalid SLURM_ARRAY_TASK_ID: ", task_id,
       "; expected 1 to ", nrow(oat_specs))
}

case_spec <- oat_specs[task_id, ]

message("Running task ", task_id, " / ", nrow(oat_specs))
message("Case: ", case_spec$case_id)

case_out_dir <- file.path(config$output_dir, "slurm_cases")
dir.create(case_out_dir, recursive = TRUE, showWarnings = FALSE)

case_file <- file.path(
  case_out_dir,
  sprintf("case_%03d_%s.rds", task_id, case_spec$case_id)
)

error_file <- file.path(
  case_out_dir,
  sprintf("case_%03d_%s_ERROR.rds", task_id, case_spec$case_id)
)

tryCatch(
  {
    case_result <- run_saureus_oat_case(case_spec)
    saveRDS(case_result, file = case_file)
    message("Saved case output: ", case_file)
  },
  error = function(e) {
    saveRDS(
      list(
        task_id = task_id,
        case_spec = case_spec,
        error = conditionMessage(e),
        traceback = traceback()
      ),
      file = error_file
    )
    stop(e)
  }
)
