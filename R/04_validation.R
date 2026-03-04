# =============================================================================
# 04_validation.R - Model Validation and Diagnostics
# FIX: NAME_2 is in spatial object not in res_df directly - join correctly
# =============================================================================

source("R/00_worldpop_api.R")

suppressPackageStartupMessages({
  library(sf)
  library(spdep)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(cli)
  library(tmap)
  library(readr)
})

PROC_DIR <- "data/processed"
FIG_DIR  <- "outputs/figures"
TAB_DIR  <- "outputs/tables"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# Load results
results_sf  <- readRDS(file.path(PROC_DIR, "results_sf.rds"))
inla_models <- readRDS(file.path(PROC_DIR, "inla_models.rds"))
m3          <- inla_models$m3

# Build flat data frame - keep NAME_2 from the sf object
res_df <- results_sf |>
  st_drop_geometry() |>
  mutate(
    residual  = Y - fitted_mean,
    resid_pct = residual / Y * 100,
    abs_resid = abs(residual)
  ) |>
  filter(!is.na(Y) & !is.na(fitted_mean))

# Make sure we have a district name column
if (!"NAME_2" %in% names(res_df)) {
  res_df$NAME_2 <- res_df$district_name
}

# =============================================================================
# SECTION 1: Fitted vs Observed
# =============================================================================
cli::cli_h1("Section 1: Fitted vs Observed")

rmse <- sqrt(mean(res_df$residual^2, na.rm = TRUE))
mae  <- mean(res_df$abs_resid, na.rm = TRUE)
r2   <- cor(res_df$Y, res_df$fitted_mean, use = "complete.obs")^2
bias <- mean(res_df$residual, na.rm = TRUE)

gof_metrics <- tibble(
  Metric = c("RMSE", "MAE", "R-squared", "Mean Bias"),
  Value  = c(rmse, mae, r2, bias)
)
print(gof_metrics)
write_csv(gof_metrics, file.path(TAB_DIR, "goodness_of_fit.csv"))

p_scatter <- ggplot(res_df, aes(x = Y, y = fitted_mean)) +
  geom_abline(slope = 1, intercept = 0, color = "red",
              linetype = "dashed", linewidth = 0.8) +
  geom_errorbar(aes(ymin = fitted_Q025, ymax = fitted_Q975),
                color = "grey70", width = 0, alpha = 0.6) +
  geom_point(aes(color = fitted_cv), size = 3, alpha = 0.85) +
  scale_color_viridis_c(name = "CV\n(uncertainty)", option = "plasma",
                         direction = -1) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = glue::glue("R\u00b2 = {round(r2,3)}  RMSE = {round(rmse/1e3,1)}K  MAE = {round(mae/1e3,1)}K"),
           size = 3.5, fontface = "italic") +
  labs(title    = "WorldPop Observed vs BYM2 Fitted Population \u2014 Odisha Districts",
       subtitle = "Error bars: 95% credible intervals. Dashed: perfect fit.",
       x = "WorldPop Observed Population (2020)",
       y = "Posterior Mean Estimate",
       caption  = "Model: Neg-Bin BYM2 | Covariates: NDBI, NDVI, EVI, NTL") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11))

ggsave(file.path(FIG_DIR, "fitted_vs_observed.png"),
       p_scatter, width = 9, height = 7, dpi = 150)
cli::cli_alert_success("Scatter plot saved.")


# =============================================================================
# SECTION 2: LOO-CV via CPO
# =============================================================================
cli::cli_h1("Section 2: LOO-CV via CPO")

cpo_vals <- m3$cpo$cpo
pit_vals <- m3$cpo$pit
lcpo     <- -sum(log(cpo_vals), na.rm = TRUE)
cli::cli_alert_info("LCPO: {round(lcpo, 2)}")

pit_df <- tibble(PIT = pit_vals[!is.na(pit_vals)])
p_pit  <- ggplot(pit_df, aes(PIT)) +
  geom_histogram(aes(y = after_stat(density)),
                 breaks = seq(0, 1, 0.1),
                 fill = "#4575B4", color = "white") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  labs(title    = "PIT Histogram (LOO Calibration Check)",
       subtitle = "Uniform = well-calibrated model",
       x = "PIT value", y = "Density") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "pit_histogram.png"),
       p_pit, width = 7, height = 5, dpi = 150)
cli::cli_alert_success("PIT histogram saved.")


# =============================================================================
# SECTION 3: Moran's I on Residuals
# =============================================================================
cli::cli_h1("Section 3: Moran's I on Residuals")

odisha_val <- results_sf |>
  left_join(res_df |> select(GID_2, residual, resid_pct) |>
              rename(resid_val = residual, resid_pct_val = resid_pct),
            by = "GID_2") |>
  arrange(GID_2)

nb_val <- poly2nb(odisha_val, queen = TRUE, snap = 0.01)
lw     <- nb2listw(nb_val, style = "W", zero.policy = TRUE)

resids_clean <- odisha_val$resid_val
resids_clean[is.na(resids_clean)] <- 0

moran_test <- moran.test(resids_clean, lw, zero.policy = TRUE)
cli::cli_alert_info("Moran's I: {round(moran_test$estimate[1], 4)}")
cli::cli_alert_info("p-value  : {round(moran_test$p.value, 4)}")
cli::cli_alert_info(
  if (moran_test$p.value > 0.05)
    "No significant spatial autocorrelation - BYM2 adequately captures spatial structure."
  else
    "Residual spatial autocorrelation detected."
)

moran_df <- tibble(
  Statistic      = "Moran's I",
  Estimate       = round(moran_test$estimate[1], 4),
  Expected       = round(moran_test$estimate[2], 4),
  Variance       = round(moran_test$estimate[3], 6),
  p_value        = round(moran_test$p.value, 4),
  Interpretation = ifelse(moran_test$p.value > 0.05,
                          "No significant spatial autocorrelation",
                          "Significant spatial autocorrelation in residuals")
)
write_csv(moran_df, file.path(TAB_DIR, "morans_i_residuals.csv"))


# =============================================================================
# SECTION 4: District Deviation Table
# =============================================================================
cli::cli_h1("Section 4: District Deviation Table")

deviation_table <- res_df |>
  select(any_of(c("NAME_2", "GID_2", "district_name")),
         Y, fitted_mean, fitted_sd, fitted_Q025, fitted_Q975,
         residual, resid_pct, fitted_cv) |>
  arrange(desc(abs(resid_pct))) |>
  mutate(across(c(Y, fitted_mean, fitted_sd), ~round(., 0)),
         across(c(resid_pct, fitted_cv),      ~round(., 2)))

print(head(deviation_table, 10))
write_csv(deviation_table, file.path(TAB_DIR, "district_deviation_table.csv"))
cli::cli_alert_success("Deviation table saved.")


# =============================================================================
# SECTION 5: Predictive Interval Coverage
# =============================================================================
cli::cli_h1("Section 5: 95% CI Coverage")

coverage_df <- res_df |>
  mutate(covered = (Y >= fitted_Q025 & Y <= fitted_Q975)) |>
  summarise(
    coverage_pct = round(mean(covered, na.rm = TRUE) * 100, 1),
    n_districts  = n(),
    n_covered    = sum(covered, na.rm = TRUE)
  )

cli::cli_alert_info("95% CI Coverage: {coverage_df$coverage_pct}% ({coverage_df$n_covered}/{coverage_df$n_districts} districts)")
write_csv(coverage_df, file.path(TAB_DIR, "interval_coverage.csv"))


# =============================================================================
# SECTION 6: Residual Map
# =============================================================================
cli::cli_h1("Section 6: Residual Map")

tmap_mode("plot")
odisha_resid <- results_sf |>
  left_join(res_df |> select(GID_2, resid_pct), by = "GID_2")

tm_resid <- tm_shape(odisha_resid) +
  tm_fill("resid_pct",
          style    = "jenks", n = 7,
          palette  = "RdBu",  midpoint = 0,
          title    = "Residual (%)") +
  tm_borders(col = "white", lwd = 0.5) +
  tm_layout(title          = "Model Residuals (BYM2) - Odisha Districts 2020",
            title.size     = 0.9,
            title.fontface = "bold",
            legend.outside = TRUE,
            frame          = FALSE)

tmap_save(tm_resid, file.path(FIG_DIR, "residual_map.png"),
          width = 8, height = 8, dpi = 150)
cli::cli_alert_success("Residual map saved.")
cli::cli_alert_success("Validation complete.")
