# =============================================================================
# run_all.R — Master Execution Script
# =============================================================================

cat("\n")
cat("================================================================\n")
cat("  Bayesian SAE for Odisha — WorldPop + Sentinel-2 + BYM2\n")
cat("  Ujjwal Kumar Swain | UNFPA India Odisha State Office\n")
cat("================================================================\n\n")

start_time <- Sys.time()

# =============================================================================
# STEP 0: Install Dependencies
# =============================================================================
cat("[Step 0] Checking R dependencies...\n")

# INLA dependencies MUST come first — fmesher is critical
inla_deps <- c("fmesher", "MatrixModels", "numDeriv", "Matrix")
missing_deps <- inla_deps[!sapply(inla_deps, requireNamespace, quietly = TRUE)]
if (length(missing_deps) > 0) {
  cat("Installing INLA dependencies:", paste(missing_deps, collapse = ", "), "\n")
  install.packages(missing_deps, dependencies = TRUE,
                   repos = "https://cloud.r-project.org")
}

# INLA
if (!requireNamespace("INLA", quietly = TRUE)) {
  cat("Installing INLA from r-inla.org...\n")
  install.packages("INLA",
    repos = c(INLA = "https://inla.r-inla-download.org/R/stable",
              CRAN = "https://cloud.r-project.org"),
    dependencies = TRUE)
}

# All other packages
required_packages <- c(
  "terra", "sf", "spdep", "httr", "jsonlite",
  "dplyr", "tidyr", "purrr", "readr", "tibble",
  "ggplot2", "tmap", "patchwork", "viridis", "scales", "forcats",
  "cli", "glue", "tools", "knitr", "kableExtra", "rmarkdown"
)
missing_pkgs <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat("Installing:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs, dependencies = TRUE,
                   repos = "https://cloud.r-project.org")
}

# Optional: rstac + gdalcubes for Sentinel-2 (graceful if unavailable)
for (pkg in c("rstac", "gdalcubes")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch(install.packages(pkg, repos = "https://cloud.r-project.org"),
             error = function(e) cat("Optional:", pkg, "not installed. Fallback will be used.\n"))
  }
}

cat("[Step 0] All dependencies satisfied.\n\n")

# =============================================================================
# Helper: timed step runner
# =============================================================================
run_step <- function(step_name, script_path) {
  cat(sprintf("[%s] Starting...\n", step_name))
  t0 <- Sys.time()
  tryCatch({
    source(script_path, local = new.env(parent = globalenv()))
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    cat(sprintf("[%s] Completed in %.1f min.\n\n", step_name, elapsed))
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("[%s] ERROR: %s\n\n", step_name, conditionMessage(e)))
    return(FALSE)
  })
}

# =============================================================================
# Execute Pipeline
# =============================================================================
status <- list()
status$download   <- run_step("Step 1: Data Download",  "R/01_data_download.R")
status$covariates <- run_step("Step 2: Covariate Prep", "R/02_covariate_prep.R")
status$modelling  <- run_step("Step 3: Model Fitting",  "R/03_model_fitting.R")
status$validation <- run_step("Step 4: Validation",     "R/04_validation.R")
status$visualise  <- run_step("Step 5: Visualisation",  "R/05_visualisation.R")

# Summary
total_time <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
cat("\n================================================================\n")
cat("  PIPELINE SUMMARY\n")
cat("================================================================\n")
for (nm in names(status)) {
  cat(sprintf("  %-20s : %s\n", nm, ifelse(isTRUE(status[[nm]]), "SUCCESS", "FAILED")))
}
cat(sprintf("\n  Total runtime: %.1f minutes\n", total_time))
cat("================================================================\n\n")
