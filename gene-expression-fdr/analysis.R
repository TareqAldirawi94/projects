################################################################################
# Multiple Hypothesis Testing for Gene Expression Analysis
# -------------------------------------------------------
# Identifies significant biomarkers from gene expression data while
# controlling the false discovery rate (FDR) at 5%.
#
# Compares three correction methods:
#   1. Bonferroni (Family-Wise Error Rate control)
#   2. Benjamini-Hochberg (FDR control)
#   3. Storey's q-value (Adaptive FDR control)
#
# Dataset: Golub et al. (1999) Leukemia microarray dataset
#   - 7,129 genes measured across 72 patients
#   - Two classes: ALL (Acute Lymphoblastic Leukemia) vs AML (Acute Myeloid Leukemia)
#
# Author: [Your Name]
# Date:   2025
################################################################################

# ── 0. Setup ─────────────────────────────────────────────────────────────────

# Install packages if not already installed
required_packages <- c("multtest", "qvalue", "ggplot2", "dplyr", "tidyr",
                       "gridExtra", "ggrepel", "scales")
bioc_packages    <- c("multtest", "qvalue")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
}

# Install CRAN packages
for (pkg in setdiff(required_packages, bioc_packages)) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

library(multtest)
library(qvalue)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggrepel)
library(scales)

# Set seed for reproducibility
set.seed(42)

# Create output directory for figures
if (!dir.exists("figures")) dir.create("figures")
if (!dir.exists("results")) dir.create("results")

# ── 1. Load and Explore the Golub Leukemia Dataset ──────────────────────────

data(golub, package = "multtest")

# golub: matrix of gene expressions (rows = genes, cols = samples)
# golub.cl: class labels (0 = ALL, 1 = AML)
# golub.gnames: gene names

expression_matrix <- golub
class_labels      <- golub.cl
gene_names        <- golub.gnames[, 3]  # Human-readable gene names

n_genes   <- nrow(expression_matrix)
n_samples <- ncol(expression_matrix)
n_all     <- sum(class_labels == 0)
n_aml     <- sum(class_labels == 1)

# ── 2. Conduct Gene-Wise Hypothesis Tests ────────────────────────────────────

# For each gene, perform a two-sample Welch's t-test:
#   H0: No difference in mean expression between ALL and AML
#   H1: Significant difference in mean expression

cat("Running Welch's t-test for each gene...\n")

all_samples <- expression_matrix[, class_labels == 0]
aml_samples <- expression_matrix[, class_labels == 1]

# Compute t-statistics and p-values efficiently
run_ttests <- function(expr_matrix, labels) {
  group0 <- expr_matrix[, labels == 0]
  group1 <- expr_matrix[, labels == 1]

  n0 <- ncol(group0)
  n1 <- ncol(group1)

  mean0 <- rowMeans(group0)
  mean1 <- rowMeans(group1)

  var0 <- apply(group0, 1, var)
  var1 <- apply(group1, 1, var)

  se <- sqrt(var0 / n0 + var1 / n1)
  t_stat <- (mean0 - mean1) / se

  # Welch-Satterthwaite degrees of freedom
  df <- (var0 / n0 + var1 / n1)^2 /
    ((var0 / n0)^2 / (n0 - 1) + (var1 / n1)^2 / (n1 - 1))

  p_values <- 2 * pt(-abs(t_stat), df = df)

  data.frame(
    gene       = gene_names,
    t_stat     = t_stat,
    p_value    = p_values,
    mean_all   = mean0,
    mean_aml   = mean1,
    log2_fc    = mean1 - mean0,  # Already log-scale
    stringsAsFactors = FALSE
  )
}

results <- run_ttests(expression_matrix, class_labels)

# ── 3. Multiple Testing Correction Methods ───────────────────────────────────

alpha <- 0.05  # Significance level / target FDR

cat("── Applying Multiple Testing Corrections (α = 0.05) ──────────\n\n")

# --- 3a. Bonferroni Correction (FWER control) ---
results$p_bonferroni <- p.adjust(results$p_value, method = "bonferroni")
results$sig_bonferroni <- results$p_bonferroni < alpha

n_bonf <- sum(results$sig_bonferroni)
cat(sprintf("  Bonferroni:           %4d significant genes\n", n_bonf))
cat(sprintf("    Adjusted threshold: p < %.2e\n", alpha / n_genes))

# --- 3b. Benjamini-Hochberg (FDR control) ---
results$p_bh <- p.adjust(results$p_value, method = "BH")
results$sig_bh <- results$p_bh < alpha

n_bh <- sum(results$sig_bh)
cat(sprintf("  Benjamini-Hochberg:   %4d significant genes\n", n_bh))

# --- 3c. Storey's q-value (Adaptive FDR control) ---
qobj <- qvalue(results$p_value)
results$q_value <- qobj$qvalues
results$sig_qvalue <- results$q_value < alpha

n_qval <- sum(results$sig_qvalue)
cat(sprintf("  Storey's q-value:     %4d significant genes\n", n_qval))
cat(sprintf("    Estimated π₀:       %.3f\n", qobj$pi0))
cat(sprintf("    (%.1f%% of genes estimated to be truly null)\n\n",
            qobj$pi0 * 100))

# ── 4. Summary Comparison ────────────────────────────────────────────────────

comparison <- data.frame(
  Method = c("No Correction", "Bonferroni (FWER)",
             "Benjamini-Hochberg (FDR)", "Storey's q-value (FDR)"),
  Controls = c("Nothing", "FWER ≤ 0.05", "FDR ≤ 0.05", "FDR ≤ 0.05"),
  Significant = c(sum(results$p_value < alpha), n_bonf, n_bh, n_qval),
  Pct_of_Total = round(c(sum(results$p_value < alpha), n_bonf, n_bh, n_qval)
                        / n_genes * 100, 2)
)

print(comparison, row.names = FALSE)
cat("\n")

# ── 5. Visualizations ────────────────────────────────────────────────────────

cat("Generating publication-quality figures...\n\n")

# Custom theme for all plots
theme_genomics <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(color = "grey40", size = 10),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

# --- 5a. P-value Distribution (Histogram) ---
p1 <- ggplot(results, aes(x = p_value)) +
  geom_histogram(bins = 50, fill = "#2C5F8A", color = "white",
                 linewidth = 0.2, alpha = 0.85) +
  geom_hline(yintercept = n_genes / 50, color = "#D94F3D",
             linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 0.75, y = n_genes / 50 + 15,
           label = "Expected under H₀ (uniform)", color = "#D94F3D",
           size = 3.5) +
  labs(
    title    = "Distribution of Raw P-values",
    subtitle = "Spike near zero indicates true differential expression signals",
    x = "P-value", y = "Frequency"
  ) +
  theme_genomics()

ggsave("figures/01_pvalue_distribution.png", p1, width = 8, height = 5, dpi = 300)

# --- 5b. Adjusted P-value Comparison ---
adj_long <- results %>%
  select(gene, p_value, p_bonferroni, p_bh, q_value) %>%
  pivot_longer(cols = c(p_bonferroni, p_bh, q_value),
               names_to = "method", values_to = "adjusted_p") %>%
  mutate(method = recode(method,
    "p_bonferroni" = "Bonferroni",
    "p_bh"         = "Benjamini-Hochberg",
    "q_value"      = "Storey's q-value"
  ))

p2 <- ggplot(adj_long, aes(x = p_value, y = adjusted_p, color = method)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red", linewidth = 0.6) +
  annotate("text", x = 0.85, y = 0.08, label = "α = 0.05",
           color = "red", size = 3.5) +
  scale_color_manual(values = c("Bonferroni" = "#E69F00",
                                "Benjamini-Hochberg" = "#56B4E9",
                                "Storey's q-value" = "#009E73")) +
  labs(
    title    = "Raw vs. Adjusted P-values by Correction Method",
    subtitle = "Points below α = 0.05 line are deemed significant",
    x = "Raw P-value", y = "Adjusted P-value / q-value",
    color = "Method"
  ) +
  theme_genomics() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

ggsave("figures/02_adjusted_pvalue_comparison.png", p2, width = 9, height = 6, dpi = 300)

# --- 5c. Volcano Plot ---

# Classify significance for coloring
results$significance <- case_when(
  results$sig_qvalue & abs(results$log2_fc) > 0.5 ~ "Significant & |FC| > 0.5",
  results$sig_qvalue ~ "Significant (FDR < 0.05)",
  TRUE ~ "Not significant"
)

# Top genes by q-value for labeling
top_genes <- results %>%
  filter(sig_qvalue, abs(log2_fc) > 0.5) %>%
  arrange(q_value) %>%
  slice_head(n = 15)

p3 <- ggplot(results, aes(x = log2_fc, y = -log10(p_value), color = significance)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(max(results$p_value[results$sig_bh])),
             linetype = "dashed", color = "grey60") +
  geom_text_repel(data = top_genes, aes(label = gene),
                  size = 2.5, max.overlaps = 15, color = "black",
                  segment.color = "grey60", segment.size = 0.3) +
  scale_color_manual(values = c(
    "Not significant"             = "grey70",
    "Significant (FDR < 0.05)"    = "#56B4E9",
    "Significant & |FC| > 0.5"    = "#D94F3D"
  )) +
  labs(
    title    = "Volcano Plot: Differential Expression (ALL vs AML)",
    subtitle = "Top biomarkers labeled by gene name",
    x = "Log₂ Fold Change (AML − ALL)",
    y = "-log₁₀(P-value)",
    color = NULL
  ) +
  theme_genomics()

ggsave("figures/03_volcano_plot.png", p3, width = 10, height = 7, dpi = 300)

# --- 5d. Venn-style Bar Chart: Method Overlap ---
results$method_count <- results$sig_bonferroni + results$sig_bh + results$sig_qvalue

overlap_data <- data.frame(
  Category = c("Bonferroni only",
                "BH only (not Bonferroni)",
                "q-value only (not BH)",
                "All three methods"),
  Count = c(
    sum(results$sig_bonferroni & !results$sig_bh & !results$sig_qvalue),
    sum(results$sig_bh & !results$sig_bonferroni),
    sum(results$sig_qvalue & !results$sig_bh),
    sum(results$sig_bonferroni & results$sig_bh & results$sig_qvalue)
  )
)


# --- 5e. Storey's π₀ estimation plot ---
p5 <- ggplot(data.frame(lambda = qobj$lambda, pi0_lambda = qobj$pi0.lambda),
             aes(x = lambda, y = pi0_lambda)) +
  geom_point(color = "#2C5F8A", size = 2) +
  geom_line(color = "#2C5F8A", alpha = 0.5) +
  geom_hline(yintercept = qobj$pi0, linetype = "dashed", color = "#D94F3D",
             linewidth = 0.8) +
  annotate("text", x = 0.1, y = qobj$pi0 + 0.02,
           label = sprintf("π̂₀ = %.3f", qobj$pi0),
           color = "#D94F3D", size = 4, fontface = "bold") +
  labs(
    title    = "Storey's π₀ Estimation",
    subtitle = "Proportion of truly null hypotheses estimated via spline fit",
    x = "λ (tuning parameter)", y = "π̂₀(λ)"
  ) +
  ylim(0, 1) +
  theme_genomics()

ggsave("figures/04_pi0_estimation.png", p5, width = 8, height = 5, dpi = 300)

# --- 5f. Number of Rejections at Various FDR Thresholds ---
fdr_thresholds <- seq(0.001, 0.20, by = 0.001)

rejection_curves <- data.frame(
  FDR = rep(fdr_thresholds, 3),
  Rejections = c(
    sapply(fdr_thresholds, function(t) sum(results$p_bonferroni < t)),
    sapply(fdr_thresholds, function(t) sum(results$p_bh < t)),
    sapply(fdr_thresholds, function(t) sum(results$q_value < t))
  ),
  Method = rep(c("Bonferroni", "Benjamini-Hochberg", "Storey's q-value"),
               each = length(fdr_thresholds))
)

p6 <- ggplot(rejection_curves, aes(x = FDR, y = Rejections, color = Method)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.055, y = max(rejection_curves$Rejections) * 0.9,
           label = "FDR = 0.05", hjust = 0, color = "grey40", size = 3.5) +
  scale_color_manual(values = c("Bonferroni" = "#E69F00",
                                "Benjamini-Hochberg" = "#56B4E9",
                                "Storey's q-value" = "#009E73")) +
  labs(
    title    = "Number of Discoveries at Various Significance Thresholds",
    subtitle = "Storey's q-value yields the most discoveries at each threshold",
    x = "Significance Threshold", y = "Number of Rejected Hypotheses",
    color = "Method"
  ) +
  theme_genomics()

ggsave("figures/05_rejection_curves.png", p6, width = 9, height = 6, dpi = 300)

# ── 6. Top Biomarker Table ───────────────────────────────────────────────────

top_biomarkers <- results %>%
  filter(sig_qvalue) %>%
  arrange(q_value) %>%
  select(gene, t_stat, log2_fc, p_value, p_bonferroni, p_bh, q_value) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 4))) %>%
  head(30)

cat("── Top 30 Biomarkers (by Storey's q-value) ────────────────────\n\n")
print(top_biomarkers, row.names = FALSE)

# Save full results
write.csv(results, "results/full_results.csv", row.names = FALSE)
write.csv(top_biomarkers, "results/top_biomarkers.csv", row.names = FALSE)
