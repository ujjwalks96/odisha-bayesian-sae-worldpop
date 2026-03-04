# =============================================================================
# 03_model_fitting.R
# Bayesian Small Area Estimation using INLA (BYM2 Spatial Prior)
#
# Model:
#   y_i ~ NegBinomial(mu_i, r)
#   log(mu_i) = alpha + X_i * beta + b_i + u_i
#
#   where:
#     y_i  = WorldPop population estimate (aggregated, district i)
#     X_i  = district-level covariates (NDBI, NDVI, NTL, EVI)
#     b_i  = spatially structured random effect (Besag/ICAR)
#     u_i  = unstructured random effect (IID)
#     BYM2 = b_i + u_i combined as per Riebler et al. (2016)
#
# Additionally fits:
#   Model 0: Null (intercept only)
#   Model 1: Fixed effects only (no spatial)
#   Model 2: BYM2 spatial (no covariates)
#   Model 3: BYM2 spatial + covariates [FULL MODEL]
#
# Reference:
#   Riebler et al. (2016) Statistical Methods in Medical Research
#   WorldPop small area methods: Stevens et al. (2015) PLOS ONE
# =============================================================================

source("R/00_worldpop_api.R")

suppressPackageStartupMessages({
  library(INLA)
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cli)
  library(scales)
})

PROC_DIR <- "data/processed"
OUT_DIR  <- "outputs"
FIG_DIR  <- "outputs/figures"
TAB_DIR  <- "outputs/tables"
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# STEP 1: Load Data
# =============================================================================
cli::cli_h1("Step 1: Load Processed Data")

# Boundaries
odisha_dist <- st_read(file.path(PROC_DIR, "odisha_districts.gpkg"), quiet = TRUE)
odisha_dist <- odisha_dist |> arrange(GID_2)   # Consistent ordering

# Covariate table
covs <- readRDS(file.path(PROC_DIR, "district_covariate_table.rds"))

# Merge covariates with spatial data
model_data <- odisha_dist |>
  left_join(covs |> select(-ID), by = c("GID_2" = "district_id")) |>
  arrange(GID_2)

cli::cli_alert_success("Model dataset: {nrow(model_data)} districts")


# =============================================================================
# STEP 2: Build Spatial Neighbourhood Matrix (Queen Contiguity)
# =============================================================================
cli::cli_h1("Step 2: Spatial Neighbourhood Matrix")

# Queen contiguity (shared boundary or vertex)
nb <- poly2nb(model_data, queen = TRUE, snap = 0.01)

# Check connectivity
cli::cli_alert_info("Average neighbours: {round(mean(card(nb)), 2)}")
cli::cli_alert_info("Islands (0 neighbours): {sum(card(nb) == 0)}")

# Handle any islands (snap to nearest)
if (any(card(nb) == 0)) {
  cli::cli_alert_warning("Fixing islands with nearest-neighbour links...")
  nb <- make.sym.nb(nb)
}

# Convert to INLA-compatible adjacency format
nb2INLA(nb, file = file.path(PROC_DIR, "odisha_adj.graph"))
g <- inla.read.graph(file.path(PROC_DIR, "odisha_adj.graph"))
cli::cli_alert_success("Adjacency graph written ({g$n} nodes, {g$ne} edges)")


# =============================================================================
# STEP 3: Prepare Model Data Frame
# =============================================================================
cli::cli_h1("Step 3: Prepare Model Data Frame")

# Response variable: WorldPop population sum per district
# Scale covariates to mean=0, sd=1 for better MCMC mixing
model_df <- model_data |>
  st_drop_geometry() |>
  mutate(
    # Response
    Y = round(pop_sum_worldpop),
    # Spatial index (1-based for INLA)
    sp_idx = row_number(),
    # Scaled covariates
    ndbi_sc  = as.numeric(scale(NDBI,    center = TRUE, scale = TRUE)),
    ndvi_sc  = as.numeric(scale(NDVI,    center = TRUE, scale = TRUE)),
    evi_sc   = as.numeric(scale(EVI,     center = TRUE, scale = TRUE)),
    ntl_sc   = as.numeric(scale(NTL_log, center = TRUE, scale = TRUE)),
    pop_sc   = as.numeric(scale(pop_count, center = TRUE, scale = TRUE)),
    # Log area as offset (population exposure offset)
    log_area = log(st_area(model_data) |> as.numeric())
  ) |>
  # Replace NAs with 0 (for districts with no valid pixel)
  mutate(across(where(is.numeric), ~replace_na(., 0)))

cli::cli_alert_info("Response range: {min(model_df$Y, na.rm=T)} - {max(model_df$Y, na.rm=T)}")
cli::cli_alert_info("N districts: {nrow(model_df)}")

# =============================================================================
# STEP 4: Define INLA Prior Specifications
# =============================================================================
cli::cli_h1("Step 4: INLA Prior Specifications")

# BYM2 hyperpriors - PC priors (Penalised Complexity)
# Riebler et al. (2016): P(sigma > 1) = 0.01; P(phi < 0.5) = 0.5
prior_bym2 <- list(
  prec = list(prior = "pc.prec", param = c(1, 0.01)),
  phi  = list(prior = "pc",      param = c(0.5, 0.5))
)

# Negative binomial overdispersion
prior_nb <- list(hyper = list(
  size = list(prior = "loggamma", param = c(1, 0.01))
))

# Fixed effect priors (weakly informative Normal)
prior_fixed <- list(
  mean = 0,
  prec = 0.001   # wide = sd~31
)


# =============================================================================
# STEP 5: Fit Four Competing Models
# =============================================================================
cli::cli_h1("Step 5: Model Fitting")

# ── Formula definitions ────────────────────────────────────────────────────
f0 <- Y ~ 1 + offset(log_area)

f1 <- Y ~ 1 + ndbi_sc + ndvi_sc + ntl_sc + evi_sc + offset(log_area)

f2 <- Y ~ 1 +
  f(sp_idx,
    model  = "bym2",
    graph  = g,
    scale.model = TRUE,
    hyper  = prior_bym2) +
  offset(log_area)

f3 <- Y ~ 1 +
  ndbi_sc + ndvi_sc + ntl_sc + evi_sc +
  f(sp_idx,
    model  = "bym2",
    graph  = g,
    scale.model = TRUE,
    hyper  = prior_bym2) +
  offset(log_area)


# ── INLA call with sensible defaults ──────────────────────────────────────
inla_fit <- function(formula, data, ...) {
  inla(
    formula,
    family            = "nbinomial",
    data              = data,
    control.family    = prior_nb,
    control.fixed     = prior_fixed,
    control.predictor = list(compute = TRUE, link = 1),
    control.compute   = list(
      dic        = TRUE,
      waic       = TRUE,
      cpo        = TRUE,
      return.marginals = TRUE,
      config     = TRUE      # needed for posterior samples
    ),
    verbose = FALSE,
    ...
  )
}

cli::cli_alert_info("Fitting Model 0: Null model...")
m0 <- inla_fit(f0, model_df)
cli::cli_alert_success("M0 DIC={round(m0$dic$dic,1)}  WAIC={round(m0$waic$waic,1)}")

cli::cli_alert_info("Fitting Model 1: Fixed effects only...")
m1 <- inla_fit(f1, model_df)
cli::cli_alert_success("M1 DIC={round(m1$dic$dic,1)}  WAIC={round(m1$waic$waic,1)}")

cli::cli_alert_info("Fitting Model 2: BYM2 spatial (no covariates)...")
m2 <- inla_fit(f2, model_df)
cli::cli_alert_success("M2 DIC={round(m2$dic$dic,1)}  WAIC={round(m2$waic$waic,1)}")

cli::cli_alert_info("Fitting Model 3: BYM2 + Covariates [FULL MODEL]...")
m3 <- inla_fit(f3, model_df)
cli::cli_alert_success("M3 DIC={round(m3$dic$dic,1)}  WAIC={round(m3$waic$waic,1)}")


# =============================================================================
# STEP 6: Model Comparison Table
# =============================================================================
cli::cli_h1("Step 6: Model Comparison")

model_compare <- tibble(
  Model       = c("M0: Null", "M1: Fixed Effects",
                  "M2: BYM2 Spatial", "M3: BYM2 + Covariates"),
  DIC         = c(m0$dic$dic, m1$dic$dic, m2$dic$dic, m3$dic$dic),
  WAIC        = c(m0$waic$waic, m1$waic$waic, m2$waic$waic, m3$waic$waic),
  p_eff       = c(m0$dic$p.eff, m1$dic$p.eff, m2$dic$p.eff, m3$dic$p.eff),
  LCPO        = c(-mean(log(m0$cpo$cpo), na.rm=T),
                   -mean(log(m1$cpo$cpo), na.rm=T),
                   -mean(log(m2$cpo$cpo), na.rm=T),
                   -mean(log(m3$cpo$cpo), na.rm=T))
) |>
  mutate(across(where(is.numeric), ~round(., 2))) |>
  mutate(delta_DIC  = DIC  - min(DIC),
         delta_WAIC = WAIC - min(WAIC))

print(model_compare)
readr::write_csv(model_compare, file.path(TAB_DIR, "model_comparison.csv"))
cli::cli_alert_success("Best model (lowest DIC): {model_compare$Model[which.min(model_compare$DIC)]}")


# =============================================================================
# STEP 7: Extract Posterior Estimates from Best Model (M3)
# =============================================================================
cli::cli_h1("Step 7: Posterior Estimates from Full Model (M3)")

best_model <- m3

# ── Fixed effects ─────────────────────────────────────────────────────────
fixed_fx <- as.data.frame(best_model$summary.fixed) |>
  tibble::rownames_to_column("Parameter") |>
  rename(
    mean   = mean,
    sd     = sd,
    Q2.5   = `0.025quant`,
    Q50    = `0.5quant`,
    Q97.5  = `0.975quant`,
    mode   = mode
  ) |>
  mutate(across(where(is.numeric), ~round(., 4)))

print(fixed_fx)
readr::write_csv(fixed_fx, file.path(TAB_DIR, "fixed_effects_m3.csv"))

# ── Hyperparameters ────────────────────────────────────────────────────────
hyper_df <- as.data.frame(best_model$summary.hyperpar) |>
  tibble::rownames_to_column("Hyperparameter") |>
  mutate(across(where(is.numeric), ~round(., 4)))

print(hyper_df)
readr::write_csv(hyper_df, file.path(TAB_DIR, "hyperparameters_m3.csv"))


# =============================================================================
# STEP 8: Generate District-Level Posterior Predictions with Uncertainty
# =============================================================================
cli::cli_h1("Step 8: Posterior Predictive Distribution")

# Mean fitted values
model_df$fitted_mean <- best_model$summary.fitted.values$mean
model_df$fitted_sd   <- best_model$summary.fitted.values$sd
model_df$fitted_Q025 <- best_model$summary.fitted.values$`0.025quant`
model_df$fitted_Q975 <- best_model$summary.fitted.values$`0.975quant`

# Coefficient of variation (uncertainty measure)
model_df$fitted_cv <- model_df$fitted_sd / model_df$fitted_mean

# Spatially structured vs unstructured random effects
re_total   <- best_model$summary.random$sp_idx
model_df$re_mean <- re_total$mean[seq_len(nrow(model_df))]
model_df$re_sd   <- re_total$sd[seq_len(nrow(model_df))]

# Rejoin to spatial data
results_sf <- model_data |>
  left_join(model_df |> select(GID_2, fitted_mean, fitted_sd,
                                 fitted_Q025, fitted_Q975,
                                 fitted_cv, re_mean, re_sd,
                                 Y),
            by = "GID_2")

# Save full results
st_write(results_sf,
         file.path(PROC_DIR, "model_results_m3.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)
saveRDS(list(m0=m0, m1=m1, m2=m2, m3=m3),
        file.path(PROC_DIR, "inla_models.rds"))
saveRDS(results_sf, file.path(PROC_DIR, "results_sf.rds"))

cli::cli_alert_success("Model results saved.")
cli::cli_alert_success("Fitted vs Observed correlation: {round(cor(results_sf$fitted_mean, results_sf$Y, use='complete.obs'), 4)}")


# =============================================================================
# STEP 9: Marginal Posterior Plots for Fixed Effects
# =============================================================================
cli::cli_h1("Step 9: Posterior Plots")

# Extract marginals
marginal_names <- c("ndbi_sc", "ndvi_sc", "ntl_sc", "evi_sc")
nice_names     <- c("NDBI (Built-up)", "NDVI (Vegetation)",
                     "NTL (Night Lights)", "EVI (Enhanced Veg.)")

marginal_plots <- map2(marginal_names, nice_names, function(var, label) {
  mg <- inla.smarginal(best_model$marginals.fixed[[var]])
  mg_df <- as.data.frame(mg)

  # HPD interval
  hpd <- inla.hpdmarginal(0.95, best_model$marginals.fixed[[var]])

  ggplot(mg_df, aes(x, y)) +
    geom_ribbon(data = filter(mg_df, x >= hpd[1] & x <= hpd[2]),
                aes(ymin = 0, ymax = y), fill = "#2166AC", alpha = 0.3) +
    geom_line(color = "#2166AC", linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(title = label,
         x     = "Coefficient",
         y     = "Posterior Density") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10))
})

p_marginals <- wrap_plots(marginal_plots, ncol = 2) +
  plot_annotation(
    title    = "Posterior Distributions of Fixed Effects - BYM2 Model",
    subtitle = "Shaded: 95% HPD interval. Red dashed: zero.",
    caption  = "Response: WorldPop district population, 2020"
  )

ggsave(file.path(FIG_DIR, "posterior_fixed_effects.png"),
       p_marginals, width = 10, height = 8, dpi = 150)
cli::cli_alert_success("Posterior plots saved.")
