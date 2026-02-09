################################################################################
# Parametric and Semi-Parametric Survival Modeling
# -------------------------------------------------
# Compares Weibull, Log-Normal, and Cox PH models for survival prediction
# with comprehensive model diagnostics and validation.
#
# Dataset: veteran (survival package) — Veterans' Administration Lung Cancer
#   - 137 patients from a randomized trial
#   - Two treatment groups: standard vs. test chemotherapy
#
# Author: Tareq Aldirawi
# Date:   2025
################################################################################

# ── 0. Setup ─────────────────────────────────────────────────────────────────

required_packages <- c("survival", "flexsurv", "ggplot2", "dplyr",
                       "tidyr", "gridExtra", "survminer")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

library(survival)
library(flexsurv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(survminer)

set.seed(2025)

if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")

cat("═══════════════════════════════════════════════════════════════════\n")
cat("  Parametric and Semi-Parametric Survival Modeling\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# ── 1. Load and Explore Data ─────────────────────────────────────────────────

data(veteran, package = "survival")

veteran <- veteran %>%
  mutate(
    trt = factor(trt, labels = c("Standard", "Test")),
    celltype = factor(celltype)
  )

n <- nrow(veteran)

cat("── Dataset Summary ─────────────────────────────────────────────\n")
cat(sprintf("  Total patients:    %d\n", n))
cat(sprintf("  Events (deaths):   %d\n", sum(veteran$status == 1)))
cat(sprintf("  Censored:          %d\n", sum(veteran$status == 0)))
cat(sprintf("  Median survival:   %.0f days\n", median(veteran$time)))
cat(sprintf("  Treatment groups:  %s\n",
            paste(levels(veteran$trt), collapse = " vs. ")))
cat("────────────────────────────────────────────────────────────────\n\n")

# ── 2. Train/Test Split ──────────────────────────────────────────────────────

idx <- sample(1:n)
n_train <- floor(n * 0.7)
train_data <- veteran[idx[1:n_train], ]
test_data  <- veteran[idx[(n_train + 1):n], ]

cat(sprintf("  Training set: %d patients\n", nrow(train_data)))
cat(sprintf("  Test set:     %d patients\n\n", nrow(test_data)))

# ── 3. Fit Models ────────────────────────────────────────────────────────────

cat("── Fitting Models ──────────────────────────────────────────────\n\n")

# 3a. Cox Proportional Hazards (Semi-parametric)
cox_model <- coxph(Surv(time, status) ~ trt + celltype + karno + age,
                   data = train_data)
cat("  Cox PH Model:\n")
print(summary(cox_model)$coefficients[, c(1, 3, 5)])
cat(sprintf("  Concordance: %.3f\n\n", summary(cox_model)$concordance[1]))

# 3b. Weibull AFT Model (Parametric)
weibull_model <- flexsurvreg(Surv(time, status) ~ trt + celltype + karno + age,
                              data = train_data, dist = "weibull")
cat("  Weibull AFT Model:\n")
cat(sprintf("  AIC: %.1f\n", AIC(weibull_model)))
cat(sprintf("  BIC: %.1f\n\n", BIC(weibull_model)))

# 3c. Log-Normal AFT Model (Parametric)
lnorm_model <- flexsurvreg(Surv(time, status) ~ trt + celltype + karno + age,
                            data = train_data, dist = "lognormal")
cat("  Log-Normal AFT Model:\n")
cat(sprintf("  AIC: %.1f\n", AIC(lnorm_model)))
cat(sprintf("  BIC: %.1f\n\n", BIC(lnorm_model)))

# 3d. Log-Logistic AFT Model
llogis_model <- flexsurvreg(Surv(time, status) ~ trt + celltype + karno + age,
                             data = train_data, dist = "llogis")
cat("  Log-Logistic AFT Model:\n")
cat(sprintf("  AIC: %.1f\n", AIC(llogis_model)))
cat(sprintf("  BIC: %.1f\n\n", BIC(llogis_model)))

# ── 4. Model Comparison ─────────────────────────────────────────────────────

# Concordance on test set
cox_risk <- predict(cox_model, newdata = test_data, type = "risk")
cox_conc <- concordance(Surv(time, status) ~ cox_risk, data = test_data)$concordance

# For parametric models, use predicted median survival
weibull_summary <- summary(weibull_model, newdata = test_data, type = "median", tidy = TRUE)
lnorm_summary   <- summary(lnorm_model, newdata = test_data, type = "median", tidy = TRUE)
llogis_summary  <- summary(llogis_model, newdata = test_data, type = "median", tidy = TRUE)

weibull_conc <- concordance(Surv(time, status) ~ I(-weibull_summary$est),
                            data = test_data)$concordance
lnorm_conc   <- concordance(Surv(time, status) ~ I(-lnorm_summary$est),
                            data = test_data)$concordance
llogis_conc  <- concordance(Surv(time, status) ~ I(-llogis_summary$est),
                            data = test_data)$concordance

comparison <- data.frame(
  Model = c("Cox PH", "Weibull AFT", "Log-Normal AFT", "Log-Logistic AFT"),
  Type = c("Semi-parametric", "Parametric", "Parametric", "Parametric"),
  AIC = c(AIC(cox_model), AIC(weibull_model), AIC(lnorm_model), AIC(llogis_model)),
  Concordance = c(cox_conc, weibull_conc, lnorm_conc, llogis_conc)
)

cat("── Model Comparison ────────────────────────────────────────────\n\n")
print(comparison, row.names = FALSE)
cat("\n")

# ── 5. Diagnostics ───────────────────────────────────────────────────────────

cat("Generating diagnostic plots and figures...\n\n")

theme_surv <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(color = "grey40", size = 10),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

# --- 5a. Kaplan-Meier by Treatment ---
km_trt <- survfit(Surv(time, status) ~ trt, data = veteran)
km_df <- data.frame(
  time  = c(km_trt$time),
  surv  = c(km_trt$surv),
  strata = rep(names(km_trt$strata), km_trt$strata)
) %>%
  mutate(strata = gsub("trt=", "", strata))

p1 <- ggplot(km_df, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 0.9) +
  scale_color_manual(values = c("Standard" = "#2C5F8A", "Test" = "#D94F3D")) +
  labs(
    title    = "Kaplan-Meier Survival Curves by Treatment Group",
    subtitle = "Veterans' Administration Lung Cancer Trial",
    x = "Time (days)", y = "Survival Probability",
    color = "Treatment"
  ) +
  theme_surv()

ggsave("figures/01_km_by_treatment.png", p1, width = 8, height = 5, dpi = 300)

# --- 5b. Cox PH Diagnostics: Schoenfeld Residuals ---
schoenfeld_test <- cox.zph(cox_model)

cat("  Schoenfeld Test for PH Assumption:\n")
print(schoenfeld_test)
cat("\n")

png("figures/02_schoenfeld_residuals.png", width = 800, height = 600, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (i in 1:min(4, ncol(schoenfeld_test$y))) {
  plot(schoenfeld_test, var = i, main = colnames(schoenfeld_test$y)[i])
}
dev.off()

# --- 5c. Model Fit Comparison ---
conc_df <- comparison %>%
  select(Model, Concordance) %>%
  mutate(Model = factor(Model, levels = Model[order(Concordance)]))

p3 <- ggplot(conc_df, aes(x = Model, y = Concordance, fill = Model)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = sprintf("%.3f", Concordance)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Cox PH" = "#2C5F8A", "Weibull AFT" = "#E69F00",
                                "Log-Normal AFT" = "#009E73", "Log-Logistic AFT" = "#56B4E9")) +
  labs(
    title    = "Test Set Concordance Index by Model",
    subtitle = "Higher is better · 0.5 = random",
    x = NULL, y = "C-index"
  ) +
  ylim(0, max(conc_df$Concordance) * 1.15) +
  theme_surv() +
  theme(legend.position = "none")

ggsave("figures/03_concordance_comparison.png", p3, width = 8, height = 5, dpi = 300)

# --- 5d. Predicted vs Observed (Calibration) ---
test_data$pred_weibull <- weibull_summary$est
test_data$pred_lnorm   <- lnorm_summary$est

# Bin by predicted quantiles and compare observed median
calibration_plot <- function(data, pred_col, model_name) {
  data$bin <- cut(data[[pred_col]], breaks = 5, labels = FALSE)
  cal <- data %>%
    group_by(bin) %>%
    summarise(
      pred_median = median(.data[[pred_col]]),
      obs_median  = median(time),
      .groups = "drop"
    )
  cal$model <- model_name
  cal
}

cal_weibull <- calibration_plot(test_data, "pred_weibull", "Weibull")
cal_lnorm   <- calibration_plot(test_data, "pred_lnorm", "Log-Normal")
cal_all <- bind_rows(cal_weibull, cal_lnorm)

p4 <- ggplot(cal_all, aes(x = pred_median, y = obs_median, color = model)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Weibull" = "#E69F00", "Log-Normal" = "#009E73")) +
  labs(
    title    = "Calibration: Predicted vs. Observed Median Survival",
    subtitle = "Points near diagonal indicate good calibration",
    x = "Predicted Median Survival (days)",
    y = "Observed Median Survival (days)",
    color = "Model"
  ) +
  theme_surv()

ggsave("figures/04_calibration_plot.png", p4, width = 8, height = 6, dpi = 300)

# ── 6. Save Results ──────────────────────────────────────────────────────────

write.csv(comparison, "results/model_comparison.csv", row.names = FALSE)

cat("\n── Output Files ────────────────────────────────────────────────\n")
cat("  figures/01_km_by_treatment.png\n")
cat("  figures/02_schoenfeld_residuals.png\n")
cat("  figures/03_concordance_comparison.png\n")
cat("  figures/04_calibration_plot.png\n")
cat("  results/model_comparison.csv\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("  Analysis complete.\n")
cat("════════════════════════════════════════════════════════════════\n")
