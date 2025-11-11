# ---- 0. Setup ----

# Define base directory for DDA and DIA comparison (Mac vs Windows)
if (.Platform$OS.type == "unix") {
  Classen_data <- file.path(
    "Filepath_for_macOS",
    "Re-integration-20250802-TTS-including-drugs-including-SIRIUS-PublicationNaming.csv"
  )
} else {
  Classen_data <- file.path(
    "Filepath_for_Windows",
    "Re-integration-20250802-TTS-including-drugs-including-SIRIUS-PublicationNaming.csv"
  )
}


{
  library(tidyverse)
  library(readr)
  library(gtools)
  library(janitor)
  library(lme4)
  library(lmerTest)
  library(performance)
  library(ggeffects)
  library(glue)
  library(stringr)
  library(ggrepel)
  library(scales)
  library(patchwork)  
  library(robustbase)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(performance)
  library(broom.mixed)
  library(emmeans)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
}
setwd(dirname(Classen_data))
# ---- 1. Load data ----
raw_data <- read_csv(Classen_data, skip = 4)
metabolite_col <- "Publication Naming"
sample_start_col <- which(colnames(raw_data) == "20250219_Collabs_R001_F5A_6_MS051_TLC047_MH_20250214_Paul_Classen_Blank_1_")
sample_cols <- colnames(raw_data)[sample_start_col:ncol(raw_data)]

# ---- 2. Sample metadata ----
meta_raw <- read_csv(Classen_data, n_max = 4, col_names = FALSE)
sample_meta <- tibble(
  File_Name = sample_cols,
  Class = as.numeric(meta_raw[1, sample_start_col:ncol(meta_raw)] %>% unlist()),
  File_type = meta_raw[2, sample_start_col:ncol(meta_raw)] %>% unlist() %>% as.character(),
  Injection_order = as.integer(meta_raw[3, sample_start_col:ncol(meta_raw)] %>% unlist()),
  Batch_ID = 1
)

# ---- 3. Reshape and join metadata ----
raw_data_filtered <- raw_data[, c(metabolite_col, sample_cols)]
# Long format
intensity_long <- raw_data_filtered %>%
  pivot_longer(-all_of(metabolite_col), names_to = "File_Name", values_to = "Intensity") %>%
  rename(Metabolite = all_of(metabolite_col)) %>%
  left_join(sample_meta, by = "File_Name") %>%
  mutate(Intensity = as.numeric(Intensity)) %>%
  filter(!is.na(Intensity))  
intensity_long_backup <- intensity_long

# ---- 4. Filter to real samples (Class 4) and parse Patient_ID & Timepoint ----
pad_id <- function(x, width = 2) {
  m <- str_match(x, "^(SLE|HC)(\\d+)$")  # prefix + digits
  out <- ifelse(is.na(m[,1]), x, sprintf("%s%0*d", m[,2], width, as.integer(m[,3])))
  out
}
timepoint_hours <- c(
  "T1" = 1,
  "T2" = 2,
  "T3" = 3,
  "T4" = 4,
  "T5" = 8,
  "T6" = 24
)

intensity_long <- intensity_long %>%
  filter(Class == 4) %>%
  mutate(
    File_Name_clean = str_remove(File_Name, "_$"),
    parts     = str_split(File_Name_clean, "_", simplify = TRUE),
    Patient_ID = parts[, ncol(parts) - 1],
    Timepoint  = parts[, ncol(parts)],
    # map to HC/SLE right here
    Condition  = if_else(str_detect(Patient_ID, "^SLE"), "SLE", "HC"),
    # normalize IDs (SLE03, HC02, etc.)
    Patient_ID = pad_id(Patient_ID),
    # numeric time
    Time_hr    = timepoint_hours[as.character(Timepoint)]
  ) %>%
  # (keep/remove exclusions as needed later in your pipeline)
  select(-File_Name_clean, -parts) %>%
  mutate(
    # nice ordering in legends/facets
    Condition = factor(Condition, levels = c("HC","SLE"))
  )
intensity_T1_all <- intensity_long %>%
  filter(Class == 4) %>%
  mutate(
    File_Name_clean = str_remove(File_Name, "_$"),
    parts     = str_split(File_Name_clean, "_", simplify = TRUE),
    Patient_ID = parts[, ncol(parts) - 1],
    Timepoint  = parts[, ncol(parts)],
    Condition  = if_else(str_detect(Patient_ID, "^SLE"), "SLE", "HC"),
    Patient_ID = pad_id(Patient_ID),
    Time_hr    = timepoint_hours[as.character(Timepoint)],
    Intensity  = as.numeric(Intensity)
  ) %>%
  select(-File_Name_clean, -parts) %>%
  filter(Timepoint == "T1", !is.na(Intensity)) %>%
  mutate(Condition = factor(Condition, levels = c("HC","SLE")))

# ---- 5. QC CVs ----
qc_data <- intensity_long_backup %>%
  filter(File_type == "QC") %>%
  select(Metabolite, File_Name, Intensity)

qc_cv <- qc_data %>%
  group_by(Metabolite) %>%
  summarise(
    Mean = mean(Intensity, na.rm = TRUE),
    SD   = sd(Intensity, na.rm = TRUE),
    CV   = 100 * SD / Mean,
    n_QC = n(),
    .groups = "drop"
  ) %>%
  arrange(CV)

# Export to CSV for documentation
write_csv(qc_cv, "QC_CV_per_metabolite.csv")

# ---- 6. Filter to metabolites with CV < 25% ----
good_mets <- qc_cv %>% filter(CV < 25) %>% pull(Metabolite)
filtered_data <- intensity_long %>% filter(Metabolite %in% good_mets)


# ---- 7. Figure 1 Panel A + B: T1 PCA (no exclusions) + T1 Volcano (no exclusions) ----
# Helpers
eps <- 1e-9
fdr_thr <- 0.05
p_thr   <- 0.05
label_n <- 15
highlight <- c("Arginine")  # optional labels to always keep

# Reduce Matrix to QC filters
use_qc_filter <- TRUE
if (use_qc_filter) {
  intensity_T1_all <- intensity_T1_all %>% filter(Metabolite %in% good_mets)
}

# ---- Panel A: PCA at T1 ----
# Helper for Unicode minus display
unicode_minus_labels <- function(x) {
  s <- format(x, scientific = FALSE, trim = TRUE)
  gsub("-", "\u2212", s, fixed = TRUE)
}

# Wide matrix (samples × metabolites)
pca_input_wide <- intensity_T1_all %>%
  select(File_Name, Patient_ID, Condition, Metabolite, Intensity) %>%
  pivot_wider(names_from = Metabolite, values_from = Intensity) %>%
  drop_na()

pca_mat <- pca_input_wide %>%
  select(-File_Name, -Patient_ID, -Condition) %>%
  scale()

pca_fit <- prcomp(pca_mat, center = TRUE, scale. = TRUE)
var_exp <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)

pca_scores <- as_tibble(pca_fit$x[, 1:2]) %>%
  bind_cols(pca_input_wide %>% select(Patient_ID, Condition))

p_pca <- ggplot(pca_scores, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 2, alpha = 0.9) +
  stat_ellipse(type = "norm", linewidth = 0.5, alpha = 0.6) +
  geom_text_repel(
    data = pca_scores %>% filter(Patient_ID == "SLE02"),
    aes(PC1, PC2, label = "SLE02"),
    inherit.aes = FALSE, size = 3.2, fontface = "bold",
    min.segment.length = 0, box.padding = 0.3, seed = 42
  ) +
  scale_x_continuous(
    labels = unicode_minus_labels
  ) +
  scale_y_continuous(
    labels = unicode_minus_labels
  ) +
  labs(
    title = "(A) PCA at T1 = 1h",
    x = paste0("PC1 (", scales::percent(var_exp[1], accuracy = 0.1), ")"),
    y = paste0("PC2 (", scales::percent(var_exp[2], accuracy = 0.1), ")"),
    color = "Condition"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_pca)



# ---- Panel B: Volcano plot at T1 ----
fdr_thr <- 0.05
fc_thr <- 0.58   # log2 fold-change threshold
label_n <- 10
highlight <- c("Arginine")  # metabolites to always label

# Compute statistics (keep numeric operations with ASCII minus)
t1_stats <- intensity_T1_all %>%
  group_by(Metabolite) %>%
  summarise(
    mean_C   = mean(Intensity[Condition == "HC"],   na.rm = TRUE),
    mean_SLE = mean(Intensity[Condition == "SLE"], na.rm = TRUE),
    log2_fc  = log2((mean_SLE + eps) / (mean_C + eps)),
    p_value  = tryCatch(t.test(Intensity ~ Condition)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    FDR       = p.adjust(p_value, method = "BH"),
    neglog10p = -log10(p_value),
    direction = case_when(
      FDR < fdr_thr & log2_fc >= fc_thr  ~ "Up in SLE (FDR & FC)",
      FDR < fdr_thr & log2_fc <= -fc_thr ~ "Down in SLE (FDR & FC)",
      FDR < fdr_thr                     ~ "Up in SLE (FDR)",
      TRUE                              ~ "NS"
    )
  )

# Determine horizontal FDR threshold line
if (any(t1_stats$FDR < fdr_thr, na.rm = TRUE)) {
  fdr_cutoff_p <- max(t1_stats$p_value[t1_stats$FDR < fdr_thr], na.rm = TRUE)
  fdr_hline_y <- -log10(fdr_cutoff_p)
} else {
  fdr_hline_y <- NA
}

# Select labels: top FDR hits + highlighted metabolites
labels_all <- t1_stats %>%
  filter(FDR < fdr_thr) %>%
  arrange(FDR) %>%
  slice_head(n = label_n) %>%
  bind_rows(t1_stats %>% filter(Metabolite %in% highlight)) %>%
  distinct(Metabolite, .keep_all = TRUE)

# Symmetric x-axis limits
xlim_val <- max(abs(t1_stats$log2_fc), na.rm = TRUE)
if (!is.finite(xlim_val) || xlim_val == 0) xlim_val <- 1

# Colors
cols <- c(
  "Up in SLE (FDR & FC)"   = "#D62728",
  "Down in SLE (FDR & FC)" = "#1F77B4",
  "Up in SLE (FDR)"        = "#FF7F0E",
  "NS"                     = "grey75"
)

# Helper function to render axis ticks with Unicode minus sign
unicode_minus_labels <- function(x) {
  # format numbers to avoid scientific notation; then replace leading hyphen with U+2212
  s <- format(x, scientific = FALSE, trim = TRUE)
  gsub("-", "\u2212", s, fixed = TRUE)
}

# Volcano plot
p_volcano <- ggplot(t1_stats, aes(x = log2_fc, y = neglog10p)) +
  geom_point(aes(color = direction), alpha = 0.9, size = 2) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = fdr_hline_y, linetype = "dotted", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed", linewidth = 0.4) +
  
  # Label repelling: shift based on significance type
  geom_text_repel(
    data = labels_all,
    aes(label = Metabolite),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.25,
    # Move left if FDR-significant but NOT passing FC threshold; otherwise right
    nudge_x = ifelse(
      labels_all$Metabolite == "Phenolethanolamine (RT 5.2)", -3,  # extra left for this one
      ifelse(
        labels_all$FDR < fdr_thr & abs(labels_all$log2_fc) < fc_thr,
        -5,   # left shift (FDR-only)
        5     # right shift (FDR+FC)
      )
    ),
    nudge_y = 0,
    segment.size = 0.3,
    segment.color = "grey40",
    min.segment.length = 0,
    max.overlaps = Inf,
    force = 2,
    seed = 42
  ) +
  
  # Replace tick labels on x-axis with Unicode minus
  scale_x_continuous(
    limits = c(-xlim_val, xlim_val),
    labels = unicode_minus_labels
  ) +
  coord_cartesian(clip = "off", ylim = c(0, max(t1_stats$neglog10p, na.rm = TRUE) * 1.05)) +
  labs(
    title = "(B) Volcano plot at T1 = 1h",
    # use Unicode minus in axis label strings (these are plain strings, not parsed expressions)
    x = "Log\u2082 fold change (SLE / HC)",
    y = "\u2212Log10(p-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = -0.06)
  )

print(p_volcano)




print(p_volcano)
ggsave("Volcano.png", p_volcano, width = 5, height = 4.8, dpi = 300)


# Save combined PCA + Volcano figure
combo <- p_pca + p_volcano + plot_layout(widths = c(1,1))
ggsave("Figure_T1_PCA_and_Volcano.png", combo, width = 10, height = 4.8, dpi = 300)

# Save stats
write_csv(t1_stats, "Volcano_T1_stats_allPatients.csv")



###
###Side by side volcano
###


# ---- Supplementary: Side-by-side volcano (raw vs FDR) ----

library(patchwork)

# Add raw significance column
t1_stats <- t1_stats %>%
  mutate(
    sig_raw = case_when(
      p_value < 0.05 & log2_fc >= fc_thr  ~ "Up in SLE (raw & FC)",
      p_value < 0.05 & log2_fc <= -fc_thr ~ "Down in SLE (raw & FC)",
      p_value < 0.05                     ~ "Up in SLE (raw)",
      TRUE                               ~ "NS"
    )
  )

# Left: raw p-value volcano
p_volcano_raw <- ggplot(t1_stats, aes(x = log2_fc, y = neglog10p)) +
  geom_point(aes(color = sig_raw), alpha = 0.9, size = 2) +
  scale_color_manual(values = c(
    "Up in SLE (raw & FC)" = "#D62728",
    "Down in SLE (raw & FC)" = "#1F77B4",
    "Up in SLE (raw)" = "#FF7F0E",
    "NS" = "grey75"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed", linewidth = 0.4) +
  labs(
    title = "Raw p-value",
    x = expression(Log[2]*" fold change (SLE / HC)"),
    y = expression(-Log[10]*"(p-value)")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

# Right: FDR volcano (existing)
p_volcano_fdr <- ggplot(t1_stats, aes(x = log2_fc, y = neglog10p)) +
  geom_point(aes(color = direction), alpha = 0.9, size = 2) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = fdr_hline_y, linetype = "dotted", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-fc_thr, fc_thr), linetype = "dashed", linewidth = 0.4) +
  labs(
    title = "FDR-adjusted",
    x = expression(Log[2]*" fold change (SLE / HC)"),
    y = expression(-Log[10]*"(p-value)")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

# Combine side by side
p_volcano_side <- p_volcano_raw + p_volcano_fdr + plot_layout(widths = c(1, 1))

# Save as supplementary figure
ggsave("Supplementary_Volcano_Raw_vs_FDR.png", p_volcano_side, width = 10, height = 5, dpi = 300)

# Print to R
print(p_volcano_side)










###
### Power analysis (final fixed version)
###

library(pwr)
library(glue)
library(dplyr)

# --- Get unique sample counts ---
n_HC <- intensity_T1_all %>%
  filter(Condition == "HC") %>%
  distinct(File_Name) %>%
  nrow()

n_SLE <- intensity_T1_all %>%
  filter(Condition == "SLE") %>%
  distinct(File_Name) %>%
  nrow()

cat("Number of unique HC samples:", n_HC, "\n")
cat("Number of unique SLE samples:", n_SLE, "\n")

# --- Power calculations ---
n_eff <- 2 / (1/n_HC + 1/n_SLE)
power_medium <- pwr.t2n.test(n1 = n_HC, n2 = n_SLE, d = 0.5, sig.level = 0.05)
effect_needed <- pwr.t2n.test(n1 = n_HC, n2 = n_SLE, power = 0.8, sig.level = 0.05)

# --- Power curve ---
d_vals <- seq(0.2, 1.2, by = 0.05)
power_vals <- sapply(d_vals, function(d)
  pwr.t2n.test(n1 = n_HC, n2 = n_SLE, d = d, sig.level = 0.05)$power
)
power_vals <- pmin(pmax(as.numeric(power_vals), 0), 1)
max_power <- max(power_vals, na.rm = TRUE)
max_d <- d_vals[which.max(power_vals)]

# --- Plot function with robust label placement ---
plot_power_curve <- function(save = FALSE) {
  if (save) {
    png("Power_Curve_SLE_vs_HC.png", width = 1800, height = 1400, res = 300)
  }
  
  par(mar = c(5, 5, 4, 2) + 0.1)
  plot(
    d_vals, power_vals, type = "b", pch = 19, lwd = 2,
    xlab = "Effect size (Cohen's d)",
    ylab = "Statistical power (1 - β)",
    main = "Power curve for SLE vs HC comparison",
    ylim = c(0, 1.05), xlim = range(d_vals),
    cex.lab = 1.3, cex.main = 1.4
  )
  
  # Horizontal line (80% power)
  abline(h = 0.8, col = "red", lty = 2, lwd = 2)
  
  # Vertical line (effect size for 80% power)
  if (!is.null(effect_needed$d) && is.finite(effect_needed$d)) {
    abline(v = effect_needed$d, col = "blue", lty = 3, lwd = 2)
    text(
      x = effect_needed$d, y = 0.1,
      labels = glue("d = {round(effect_needed$d, 2)} for 80% power"),
      srt = 90, adj = c(0, 0.5), col = "blue", cex = 0.9
    )
  }
  
  # --- Robust "Max power" label positioning (slightly above point) ---
  offset_x <- 0.05 * diff(range(d_vals))
  offset_y <- 0.03  # vertical offset in y-axis units
  
  # Clamp coordinates inside plotting area
  text_x <- pmin(max_d + offset_x, max(d_vals) - offset_x)
  text_y <- pmin(max_power + offset_y, 0.98)  # move label slightly above point
  
  # Choose left vs right placement
  pos <- if (max_d > mean(d_vals)) 2 else 4  # left if on right side, right otherwise
  
  text(
    x = text_x, y = text_y,
    labels = glue("Max power = {round(max_power, 2)} (d = {round(max_d, 2)})"),
    col = "black", cex = 1.0, pos = pos, xpd = TRUE
  )
  
  # Legend
  legend(
    "bottomright",
    legend = c("Power curve", "80% threshold", "d for 80% power"),
    col = c("black", "red", "blue"),
    lty = c(1, 2, 3), pch = c(19, NA, NA),
    bty = "n", pt.cex = 1.2
  )
  
  if (save) {
    dev.off()
    cat("✅ Saved: Power_Curve_SLE_vs_HC.png (300 dpi)\n")
  }
}

# --- Save high-res PNG ---
plot_power_curve(save = TRUE)

# --- View in RStudio ---
plot_power_curve(save = FALSE)








# ---- 8. Outlier justification at T1 (no exclusions), QC CV <25% ----
# 1) Wide matrix per condition
t1_wide <- intensity_T1_all %>%
  select(Sample = File_Name, Patient_ID, Condition, Metabolite, Intensity) %>%
  pivot_wider(names_from = Metabolite, values_from = Intensity) %>%
  drop_na()
# 2) Robust z-score outlier index (within condition)
robust_z_index <- t1_wide %>%
  group_by(Condition) %>%
  group_modify(~{
    # .x does NOT contain the grouping column; it's in .y
    X <- .x %>% select(-Sample, -Patient_ID)   # no -Condition here
    
    med  <- apply(X, 2, median, na.rm = TRUE)
    madv <- apply(X, 2, mad,    na.rm = TRUE); madv[madv == 0] <- NA
    Z <- sweep(sweep(as.matrix(X), 2, med, "-"), 2, (madv * 1.4826), "/")
    
    tibble(
      Sample           = .x$Sample,
      Patient_ID       = .x$Patient_ID,
      outlier_frac_3sd = rowMeans(abs(Z) > 3,   na.rm = TRUE),
      outlier_frac_2_5 = rowMeans(abs(Z) > 2.5, na.rm = TRUE)
    )
  }) %>%
  ungroup()
# 3) Quick plot highlighting SLE2
p_z <- robust_z_index %>%
  mutate(is_SLE02 = Patient_ID == "SLE02") %>%
  ggplot(aes(x = reorder(Patient_ID, outlier_frac_3sd),
             y = outlier_frac_3sd, fill = is_SLE02)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(values = c("FALSE"="grey70","TRUE"="#D62728"), guide = "none") +
  coord_flip() +
  labs(title = "Sample-level outlier index at T1 (robust z within condition)",
       y = "Fraction of metabolites with |z| > 3", x = "Patient ID") +
  theme_minimal(base_size = 11)
ggsave("T1_outlier_index_robustZ.png", p_z, width = 6, height = 4, dpi = 300)
# 4) Robust Mahalanobis distance in PCA space (all samples together)
# 5) Use scaled intensities (all T1 samples) and the first k PCs
X_all <- t1_wide %>% select(-Sample, -Patient_ID, -Condition) %>% scale()
pc <- prcomp(X_all, center = TRUE, scale. = TRUE)
k <- min(5, ncol(pc$x))  # first 5 PCs
PCs <- pc$x[, 1:k, drop = FALSE]
# 6) Robust covariance on PCs
mcd <- covMcd(PCs)
md2 <- mahalanobis(PCs, center = mcd$center, cov = mcd$cov)
# 7) p-values vs chi-square_k
pvals <- pchisq(md2, df = k, lower.tail = FALSE)
md_df <- t1_wide %>% select(Sample, Patient_ID, Condition) %>%
  mutate(Mahalanobis2 = as.numeric(md2),
         p_outlier = as.numeric(pvals))
# 8) Plot with SLE2 labeled
p_md <- md_df %>%
  mutate(is_SLE2 = Patient_ID == "SLE02") %>%
  ggplot(aes(x = Patient_ID, y = Mahalanobis2, color = Condition)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_text_repel(data = subset(md_df, Patient_ID == "SLE02"),
                  aes(label = Patient_ID), seed = 42, size = 3.2, fontface = "bold") +
  geom_hline(yintercept = qchisq(0.975, df = k), linetype = "dashed", alpha = 0.6) +
  labs(title = paste0("Robust Mahalanobis distance at T1 (", k, " PC space)"),
       subtitle = "Dashed = 97.5% χ² threshold",
       y = expression(Mahalanobis^2), x = "Patient ID") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.margin = margin(t = -5, b = 5),  # pull legend closer to subtitle
    plot.title = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # rotate labels
  )
ggsave("T1_outlier_mahalanobis.png", p_md, width = 6.5, height = 4, dpi = 300)
# 9) Save table for SI
write_csv(robust_z_index %>% left_join(md_df, by = c("Sample","Patient_ID")),
          "T1_outlier_scores_robustZ_and_Mahalanobis.csv")


# ---- 9. Transform intensity_long to not include SLE01 and SLE02 -----
intensity_long_all_ID <- intensity_long
intensity_long <- intensity_long %>%
  filter(!(Patient_ID %in% c("SLE01", "SLE02")))

# ---- 10. PCA trajectories -----
# Filter to good metabolites & real samples (intensity_long should already be pre-filtered)
pca_input_long <- intensity_long %>%
  filter(Metabolite %in% good_mets)

# Pivot to wide format
pca_input_wide <- pca_input_long %>%
  select(File_Name, Patient_ID, Time_hr, Condition, Metabolite, Intensity) %>%
  pivot_wider(names_from = Metabolite, values_from = Intensity) %>%
  drop_na()  # remove incomplete entries

# Prepare matrix for PCA
pca_mat <- pca_input_wide %>%
  select(-File_Name, -Patient_ID, -Time_hr, -Condition) %>%
  scale()  # center and scale intensities

# Run PCA
pca_result <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

# Scores + metadata
pca_scores <- as_tibble(pca_result$x) %>%
  bind_cols(pca_input_wide %>% select(Patient_ID, Time_hr, Condition))

# Variance explained
pca_var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Helper for Unicode minus labels
unicode_minus_labels <- function(x) {
  s <- format(x, scientific = FALSE, trim = TRUE)
  gsub("-", "\u2212", s, fixed = TRUE)
}

# PCA trajectories plot (no legend, no title)
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, group = Patient_ID, color = Patient_ID)) +
  geom_path(
    arrow = arrow(type = "closed", length = unit(0.125, "inches")),
    linewidth = 0.75,
    alpha = 0.8
  ) +
  geom_point(size = 2) +
  geom_text_repel(
    data = subset(pca_scores, Time_hr == 1),
    aes(label = Patient_ID),
    size = 3,
    seed = 42,
    box.padding = 0.2,
    point.padding = 0.1,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.curvature = 0.05,
    segment.ncp = 3,
    segment.angle = 20
  ) +
  facet_wrap(~ Condition, ncol = 2) +
  scale_x_continuous(labels = unicode_minus_labels) +
  scale_y_continuous(labels = unicode_minus_labels) +
  labs(
    title = NULL,        # remove title
    x = paste0("PC1 (", round(pca_var_explained[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(pca_var_explained[2] * 100, 1), "%)"),
    color = NULL         # remove legend title
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"  # fully hide legend
  )

ggsave("PCA_Trajectories_By_Patient.png", plot = p, width = 8, height = 5)






# ---- 11. LMEs ----
# Reincluding SLE01
# intensity_long <- intensity_long_all_ID %>%
#   filter(!(Patient_ID %in% c("SLE02")))
# ── PACKAGES ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(performance)
  library(emmeans)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(tibble)
})

# 0) BUILD MODELING DATA
# Exclude only SLE02 (keep SLE01 in LME)
time_levels <- c("T1","T2","T3","T4","T5","T6")

df_fit <- intensity_long %>%
  filter(Metabolite %in% good_mets,
         !(Patient_ID %in% c("SLE02"))) %>%
  mutate(
    Time_f    = factor(Timepoint, levels = time_levels),
    Condition = factor(Condition, levels = c("HC","SLE"))
  )

# ── 1) PER-METABOLITE LME + R² + EMM CONTRASTS ───────────────────────────────
fit_tbl <- df_fit %>%
  group_by(Metabolite) %>%
  group_modify(~{
    dat <- .x
    m <- tryCatch(
      lmer(Intensity ~ Condition * Time_f + (1 | Patient_ID), data = dat, REML = TRUE),
      error = function(e) NULL
    )
    if (is.null(m)) return(tibble())
    
    r2 <- performance::r2(m)
    an <- anova(m, type = 3)
    
    em <- emmeans(m, ~ Condition | Time_f)
    cmp <- pairs(em) %>% as.data.frame()   # SLE - HC at each time
    
    tibble(
      model           = list(m),
      R2_marginal     = r2$R2_marginal,
      R2_conditional  = r2$R2_conditional,
      anova           = list(as.data.frame(an)),
      contrasts       = list(cmp)
    )
  }) %>%
  ungroup()

# ── 2) FLAT TABLES (R², ANOVA, CONTRASTS) + ROBUST P-COL DISCOVERY ───────────
anova_tbl <- fit_tbl %>%
  transmute(Metabolite,
            anova = map(anova, ~ .x %>% tibble::rownames_to_column("term") %>% as_tibble())) %>%
  unnest(anova)

# find the p-value column name robustly and rename to p_value
get_pcol <- function(nm) intersect(c("Pr(>F)","Pr(>Chisq)","p_value","p.value"), nm)[1]
pcol <- get_pcol(names(anova_tbl))
if (is.na(pcol)) stop("Couldn't find a p-value column in ANOVA table.")
anova_tbl <- anova_tbl %>% rename(p_value = !!rlang::sym(pcol))

# detect interaction term label robustly (e.g., "Condition:Time_f")
interaction_term <- unique(anova_tbl$term)[str_detect(unique(anova_tbl$term), "Condition.*Time")]
if (length(interaction_term) == 0) stop("Couldn't detect the Condition×Time interaction term.")
interaction_term <- interaction_term[1]

# contrasts tidy table (HC vs SLE at each time)
contrasts_tbl <- fit_tbl %>%
  select(Metabolite, contrasts) %>%
  unnest(contrasts) %>%
  rename(
    Time     = Time_f,
    estimate_SLE_minus_HC = estimate,
    SE       = SE,
    df       = df,
    p_value  = p.value
  )

# R² table
r2_tbl <- fit_tbl %>% select(Metabolite, R2_marginal, R2_conditional)

# export (optional)
write_csv(r2_tbl,       "LME_TimeFactor_R2.csv")
write_csv(anova_tbl,    "LME_TimeFactor_TypeIII.csv")
write_csv(contrasts_tbl,"LME_TimeFactor_Contrasts_by_Time.csv")

# ── 3) R² SCATTER (half A4, legends side-by-side at bottom) ──────────────────
anv_int <- anova_tbl %>%
  filter(term == interaction_term) %>%
  group_by(term) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  select(Metabolite, p_value, q_value)

r2_join <- r2_tbl %>%
  left_join(anv_int, by = "Metabolite") %>%
  mutate(
    inter_sig = !is.na(q_value) & q_value < 0.05,
    size_var  = if_else(!is.na(q_value) & q_value > 0, -log10(q_value), 0)
  )

# label top by interaction strength
top_labels <- r2_join %>% arrange(q_value) %>% slice_head(n = 15)

# Helper for Unicode minus labels
unicode_minus_labels <- function(x) {
  s <- format(x, scientific = FALSE, trim = TRUE)
  gsub("-", "\u2212", s, fixed = TRUE)
}

# R² scatter plot with Unicode minus
p_r2 <- ggplot(r2_join, aes(R2_marginal, R2_conditional)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(size = size_var, fill = inter_sig),
             shape = 21, color = "grey20", alpha = 0.9, stroke = 0.25) +
  geom_text_repel(
    data = top_labels,
    aes(label = Metabolite),
    size = 3.2,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 42
  ) +
  scale_fill_manual(
    values = c(`FALSE` = "grey75", `TRUE` = "#D62728"),
    name   = "Interaction q < 0.05",
    labels = c("No", "Yes")
  ) +
  scale_size_continuous(
    name = expression("\u2212log"[10]*"(q)"),  # Unicode minus in legend title
    range = c(2.5, 10), breaks = c(1, 2, 3, 4), limits = c(0, NA)
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.02),
                     labels = unicode_minus_labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02),
                     labels = unicode_minus_labels) +
  labs(
    title = NULL,
    x = "Marginal R² (fixed effects: Condition + Time + Interaction)",
    y = "Conditional R² (fixed + random Patient_ID)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.spacing.x = unit(10, "pt"),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    plot.margin = margin(10, 12, 10, 10)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5)),
    size = guide_legend(order = 2)
  )

print(p_r2)
ggsave("LME_R2_scatter_halfA4.pdf", p_r2, width = 8.27, height = 5.85, device = cairo_pdf)
ggsave("LME_R2_scatter_halfA4.png", p_r2, width = 8.27, height = 5.85, dpi = 300)


# ── 4) EMM PLOTS FOR TOP METABOLITES BY MARGINAL R² + 3×4 GRID ───────────────
export_tbl <- r2_tbl %>%
  left_join(anv_int %>% select(Metabolite, q_value), by = "Metabolite")
write_csv(export_tbl, "LME_R2_with_interaction_q.csv")

top12 <- c("Hypoxanthin",
           "Hexose",
           "Choline",
           "Glycerophosphocholine",
           "LPC 18:2 RT7.5",
           "Glutamic acid",
           "Pyroglutamic acid",
           "Glutamine",
           "Arginine",
           "Cystine",
           "Asp-Phe",
           "Phe-Phe")
cond_cols <- c(HC = "#1F77B4", SLE = "#D62728")

scientific_superscript <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) return("")
    if (val == 0) return("0")
    exponent <- floor(log10(abs(val)))
    mantissa <- round(val / (10^exponent), 2)
    parse(text = paste0(mantissa, " %*% 10^", exponent))
  })
}

plot_emm_one <- function(target_met, log_y = FALSE, sci_y = TRUE) {
  dat <- df_fit %>% filter(Metabolite == target_met)
  if (nrow(dat) < 6) return(NULL)
  
  m <- tryCatch(
    lmer(Intensity ~ Condition * Time_f + (1 | Patient_ID), data = dat, REML = TRUE),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  emm_df <- as.data.frame(emmeans(m, ~ Condition | Time_f))
  
  p <- ggplot(dat, aes(x = Time_f, y = Intensity, group = Patient_ID, color = Condition)) +
    geom_point(alpha = 0.35, size = 1.5) +
    geom_line(alpha = 0.25, linewidth = 0.5) +
    geom_point(data = emm_df,
               aes(x = Time_f, y = emmean, color = Condition),
               inherit.aes = FALSE, size = 2.7, shape = 17) +
    geom_line(data = emm_df,
              aes(x = Time_f, y = emmean, group = Condition, color = Condition),
              inherit.aes = FALSE, linewidth = 0.7) +
    geom_errorbar(data = emm_df,
                  aes(x = Time_f, y = emmean, ymin = lower.CL, ymax = upper.CL, color = Condition),
                  inherit.aes = FALSE, width = 0.12, alpha = 0.8) +
    scale_color_manual(values = cond_cols) +
    labs(title = target_met, x = NULL, y = NULL, color = "Condition") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          legend.box.margin = margin(0,0,0,0),
          legend.margin = margin(t = 2, b = 5),
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          plot.margin = margin(4, 6, 4, 6))
  
  if (log_y) {
    p <- p + scale_y_log10(labels = scientific_superscript)
  } else if (sci_y) {
    p <- p + scale_y_continuous(labels = scientific_superscript)
  }
  
  p
}

# Generate all plots
plots <- map(top12, ~ plot_emm_one(.x, log_y = FALSE, sci_y = TRUE)) %>% compact()

# Combine plots with shared legend
combined <- wrap_plots(plots, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

library(cowplot)

# Add title
final <- combined +
  plot_annotation(
    title = "Estimated marginal means over time (selected metabolites)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.margin = margin(5, 5, 5, 30)
    )
  )

# Draw labels dynamically
final_labeled <- ggdraw(final) +
  draw_label(
    "Area",
    x = 0.02, y = 0.5, angle = 90,     # very close to left edge
    vjust = 0.5, hjust = 0.5,
    size = 12
  ) +
  draw_label(
    "Timepoint",
    x = 0.5, y = 0.06,                  # ~8% from bottom → just above legend
    vjust = 0.5, size = 12
  )

print(final_labeled)

ggsave("EMM_top12_grid_with_labels.png", final_labeled,
       width = 8.27, height = 8.75, dpi = 300)


# ── 3x3 grid for top9 metabolites ──
top9 <- c("Hypoxanthine",
          "Hexose",
          "Choline",
          "GPC",
          "Glutamic acid",
          "Pyroglutamic acid",
          "Glutamine",
          "Arginine",
          "Cystine (M+H)")

plots <- purrr::map(top9, ~ plot_emm_one(.x, log_y = FALSE, sci_y = TRUE)) %>% compact()

combined <- wrap_plots(plots, ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

final <- combined +
  plot_annotation(
    title = NULL,
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.margin = margin(5, 5, 5, 30)
    )
  )

# Draw labels dynamically
final_labeled <- ggdraw(final) +
  draw_label(
    "Area",
    x = 0.02, y = 0.5, angle = 90,     # very close to left edge
    vjust = 0.5, hjust = 0.5,
    size = 12
  ) +
  draw_label(
    "Timepoint",
    x = 0.5, y = 0.06,                  # ~8% from bottom → just above legend
    vjust = 0.5, size = 12
  )

print(final_labeled)

ggsave("EMM_top9_grid_with_labels.png", final_labeled,
       width = 8.27, height = 8.75, dpi = 300)







### LME Models for supplementary
# Make sure glue is loaded
library(glue)

plot_emm_one <- function(target_met, log_y = FALSE, sci_y = TRUE) {
  dat <- df_fit %>% dplyr::filter(Metabolite == target_met)
  if (nrow(dat) < 6) return(NULL)
  
  m <- tryCatch(
    lmer(Intensity ~ Condition * Time_f + (1 | Patient_ID), data = dat, REML = TRUE),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  emm_df <- as.data.frame(emmeans(m, ~ Condition | Time_f))
  r2 <- performance::r2(m)
  
  # Coerce to numeric scalars robustly
  r2_m <- tryCatch(as.numeric(r2$R2_marginal[[1]]), silent = TRUE)
  if (!is.finite(r2_m)) r2_m <- tryCatch(as.numeric(r2$R2_marginal), silent = TRUE)
  r2_c <- tryCatch(as.numeric(r2$R2_conditional[[1]]), silent = TRUE)
  if (!is.finite(r2_c)) r2_c <- tryCatch(as.numeric(r2$R2_conditional), silent = TRUE)
  
  # Get q-value from export_tbl if available
  q_val <- NA_real_
  if (exists("export_tbl", inherits = TRUE)) {
    qtmp <- export_tbl %>% dplyr::filter(Metabolite == target_met) %>% dplyr::pull(q_value)
    if (length(qtmp)) q_val <- qtmp[1]
  }
  
  annot <- glue(
    "Marginal R² = {formatC(r2_m, format = 'f', digits = 2)} | ",
    "Conditional R² = {formatC(r2_c, format = 'f', digits = 2)} | ",
    "Interaction q = {ifelse(is.finite(q_val), signif(q_val, 2), 'NA')}"
  )
  
  p <- ggplot(dat, aes(x = Time_f, y = Intensity, group = Patient_ID, color = Condition)) +
    geom_point(alpha = 0.35, size = 1.5) +
    geom_line(alpha = 0.25, linewidth = 0.5) +
    geom_point(data = emm_df, aes(x = Time_f, y = emmean, color = Condition),
               inherit.aes = FALSE, size = 2.7, shape = 17) +
    geom_line(data = emm_df, aes(x = Time_f, y = emmean, group = Condition, color = Condition),
              inherit.aes = FALSE, linewidth = 0.7) +
    geom_errorbar(data = emm_df,
                  aes(x = Time_f, ymin = lower.CL, ymax = upper.CL, color = Condition),
                  inherit.aes = FALSE, width = 0.12, alpha = 0.8) +
    scale_color_manual(values = cond_cols) +
    labs(title = target_met, subtitle = annot, x = NULL, y = NULL, color = "Condition") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  
  if (log_y) {
    p <- p + scale_y_log10()
  } else if (sci_y) {
    p <- p + scale_y_continuous(labels = scales::label_scientific())
  }
  p
}
# Build ALL plots (choose your metabolite list: good_mets or unique(df_fit$Metabolite))
met_list <- sort(unique(df_fit$Metabolite))
plots_all <- purrr::map(met_list, ~ plot_emm_one(.x, log_y = FALSE, sci_y = TRUE)) %>% purrr::compact()
# Save as multi-page PDF (one metabolite per page)
pdf("Supplementary_Figure_EMM_all_metabolites.pdf", width = 7, height = 5)
for (p in plots_all) print(p)
dev.off()


# ---- 12. Boxplot Function: grouped x-axis (Condition_Time) ----
plot_metabolite_boxplot_grouped <- function(met) {
  time_order <- c(1, 2, 3, 4, 8, 24)
  condition_order <- c("HC", "SLE")
  
  # Create ordered combinations: Control first, then SLE
  group_levels <- unlist(lapply(condition_order, function(cond) {
    paste(cond, time_order, sep = "_")
  }))
  
  intensity_long %>%
    filter(Metabolite == met, !is.na(Condition)) %>%
    mutate(
      Time_hr = as.numeric(Time_hr),
      Group = paste(Condition, Time_hr, sep = "_"),
      Group = factor(Group, levels = group_levels)
    ) %>%
    ggplot(aes(x = Group, y = Intensity, fill = Condition)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.8) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
    scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
    labs(
      title = paste("Boxplot by Condition and Time (hours) -", met),
      x = "Condition & Time (h)", y = "Raw Intensity"
    ) +
    theme_minimal()
}

# Create CV-based output folders
dir.create("metabolite_boxplots_by_condition_CV_low", showWarnings = FALSE)
dir.create("metabolite_boxplots_by_condition_CV_high", showWarnings = FALSE)

# Loop through all metabolites and sort boxplots into folders by CV
walk(unique(intensity_long$Metabolite), ~{
  met <- .x
  safe_name <- make_clean_names(met) %>% str_trunc(100, side = "right", ellipsis = "")
  
  # Lookup CV from previously calculated qc_cv table
  cv_val <- qc_cv %>% filter(Metabolite == met) %>% pull(CV)
  
  # Decide folder based on CV
  folder <- if (!is.na(cv_val) && cv_val < 25) {
    "metabolite_boxplots_by_condition_CV_low"
  } else {
    "metabolite_boxplots_by_condition_CV_high"
  }
  
  # Save boxplot to the appropriate folder
  ggsave(
    filename = file.path(folder, paste0(safe_name, ".pdf")),
    plot = plot_metabolite_boxplot_grouped(met),
    width = 10, height = 5
  )
})

# ---- 13. Line plots ---------------
plot_metabolite_lines_faceted <- function(met) {
  intensity_long %>%
    filter(Metabolite == met, !is.na(Condition)) %>%
    mutate(
      Time_hr = as.numeric(Time_hr)  # ensures x-axis reflects true time distances
    ) %>%
    ggplot(aes(x = Time_hr, y = Intensity, group = Patient_ID, color = Patient_ID)) +
    geom_line(alpha = 0.7, linewidth = 0.8) +
    geom_point(size = 1) +
    facet_wrap(~ Condition, ncol = 2, scales = "fixed") +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 8, 24)) +
    labs(
      title = paste("Time Courses by Condition -", met),
      x = "Time (hours)", y = "Raw Intensity"
    ) +
    theme_minimal()
}

# Create output folders for lineplots only
dir.create("metabolite_lineplots_by_condition_CV_low", showWarnings = FALSE)
dir.create("metabolite_lineplots_by_condition_CV_high", showWarnings = FALSE)

# Loop over metabolites to save lineplots into CV-binned folders
walk(unique(intensity_long$Metabolite), ~{
  met <- .x
  safe_name <- make_clean_names(met) %>% str_trunc(100, side = "right", ellipsis = "")
  
  # Look up CV value
  cv_val <- qc_cv %>% filter(Metabolite == met) %>% pull(CV)
  
  # Decide folder based on CV
  suffix <- if (!is.na(cv_val) && cv_val <= 25) "CV_low" else "CV_high"
  output_folder <- paste0("metabolite_lineplots_by_condition_", suffix)
  
  # Save plot
  ggsave(
    filename = file.path(output_folder, paste0(safe_name, ".pdf")),
    plot = plot_metabolite_lines_faceted(met),
    width = 10, height = 5
  )
})

# ---- 14.0 Drug ----
# ---- 14.1 Drug Specific Effects ----
suppressPackageStartupMessages({
  library(tidyverse); library(janitor); library(stringr); library(emmeans)
  library(lme4); library(lmerTest); library(broom.mixed); library(performance)
  library(glue); library(ggrepel); library(ggplot2)
})

# ---- 0) Inputs & prechecks -----------------------------------------------
# Assumes earlier in your script you created:
#   - intensity_long  (possibly filtered)
#   - intensity_long_all_ID  (unfiltered backup; if not, we'll fall back)
#   - good_mets (QC CV < 25% list)
# If the backup does not exist, we'll use intensity_long.

# Optional: list features that are *drugs themselves* to de-emphasize in summaries
drug_feature_names <- c(
  "Prednisolone","Hydroxychloroquine","Mycophenolic acid","Mycophenolic acid RT6.2",
  "MPA","MPA_glucuronide","Paraxanthine (caffeine metabolite)" # adapt to your naming
)

# Where to write outputs
dir.create("drug_effects", showWarnings = FALSE)
dir.create("drug_effects/plots", showWarnings = FALSE)

# ---- 1) Read + harmonize drug flags (NO grouping; keep drugs as-is) -------
drug_path <- "drug_info.csv"

pad_id2 <- function(x) {
  x <- toupper(gsub("\\s+", "", as.character(x)))
  m <- stringr::str_match(x, "^(SLE|HC)(\\d+)$")
  ifelse(is.na(m[,1]), x, sprintf("%s%02d", m[,2], as.integer(m[,3])))
}
binify <- function(v) {
  if (is.numeric(v)) return(as.integer(replace_na(v, 0) > 0))
  x <- tolower(trimws(as.character(v)))
  y <- dplyr::case_when(
    x %in% c("1","y","yes","true","t","x","+","on","present") ~ 1L,
    x %in% c("0","n","no","false","f","-","off","absent","")  ~ 0L,
    TRUE ~ NA_integer_
  )
  replace_na(y, 0L)
}

drug_raw <- readr::read_csv(drug_path, show_col_types = FALSE) %>% janitor::clean_names()
pid_col <- intersect(names(drug_raw), c("patient_id","patient","id","pid","patientid","patient_id_"))
stopifnot("drug_info.csv must have an ID column like patient_id / patient / id / pid" = length(pid_col) > 0)

drug_meta <- drug_raw %>%
  rename(patient_id = all_of(pid_col[1])) %>%
  mutate(Patient_ID = pad_id2(patient_id)) %>%
  select(-patient_id)

# Convert all non-ID columns to binary 0/1 flags (as recorded)
flag_cols <- setdiff(names(drug_meta), "Patient_ID")
drug_meta <- drug_meta %>% mutate(across(all_of(flag_cols), binify))

# Keep ONLY the as-recorded drug flags (no artificial groupings)
drug_flags <- setdiff(names(drug_meta), "Patient_ID")

# Choose the metabolomics base: prefer unfiltered backup if available
int_base <- if (exists("intensity_long_all_ID")) intensity_long_all_ID else intensity_long

# Join with long data (avoid duplicates if re-running)
already <- intersect(names(int_base), drug_flags)
if (length(already)) message("Note: dropping pre-existing drug columns: ", paste(already, collapse=", "))
int_drug <- int_base %>%
  select(-any_of(drug_flags)) %>%
  left_join(drug_meta, by = "Patient_ID") %>%
  mutate(across(all_of(drug_flags), ~ tidyr::replace_na(., 0L)))

# ---- 2) SLE-only baseline (T1) split per drug ------------------------------
# For each drug: within SLE at T1, drug+ vs drug-.
# Returns log2FC, Welch t p, Wilcoxon p, Hedges' g; BH-correct p_t per drug.
min_per_group <- 3

hedges_g <- function(x, y) {
  nx <- sum(!is.na(x)); ny <- sum(!is.na(y))
  sx <- sd(x, na.rm = TRUE); sy <- sd(y, na.rm = TRUE)
  if (nx < 2 || ny < 2 || (sx == 0 && sy == 0)) return(NA_real_)
  sp <- sqrt(((nx-1)*sx^2 + (ny-1)*sy^2) / (nx + ny - 2))
  g <- (mean(x, na.rm=TRUE) - mean(y, na.rm=TRUE)) / sp
  J <- 1 - (3 / (4*(nx+ny) - 9))  # small-sample correction
  g * J
}

t1_sle <- int_drug %>% filter(Condition == "SLE", Timepoint == "T1")

t1_results_all <- purrr::map_dfr(drug_flags, function(dv) {
  dat <- t1_sle %>% select(Metabolite, Patient_ID, Intensity, !!sym(dv))
  out <- dat %>%
    group_by(Metabolite) %>%
    group_modify(~{
      g <- .x[[dv]]
      n1 <- sum(g == 1, na.rm = TRUE); n0 <- sum(g == 0, na.rm = TRUE)
      if (n1 < min_per_group || n0 < min_per_group) {
        tibble(n_pos = n1, n_neg = n0,
               log2FC = NA_real_, p_t = NA_real_, p_wilcox = NA_real_, hedges_g = NA_real_)
      } else {
        x <- .x$Intensity[g == 1]; y <- .x$Intensity[g == 0]
        tibble(
          n_pos = n1, n_neg = n0,
          log2FC  = log2(mean(x, na.rm = TRUE) / mean(y, na.rm = TRUE)),
          p_t     = tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_),
          p_wilcox= tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_),
          hedges_g= hedges_g(x, y)
        )
      }
    }) %>%
    ungroup() %>%
    mutate(drug = dv)
  out
})

t1_results_all <- t1_results_all %>%
  group_by(drug) %>%
  mutate(q_t = p.adjust(p_t, method = "BH")) %>%
  ungroup()

readr::write_csv(t1_results_all, file.path("drug_effects", "T1_SLE_split_all_drugs.csv"))

# Per-drug patient counts table (useful for SI)
per_drug_counts <- int_drug %>%
  filter(Condition == "SLE") %>%
  distinct(Patient_ID, across(all_of(drug_flags))) %>%
  pivot_longer(cols = all_of(drug_flags), names_to = "drug", values_to = "flag") %>%
  group_by(drug) %>%
  summarise(n_pos = sum(flag == 1, na.rm = TRUE),
            n_neg = sum(flag == 0, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n_pos))
readr::write_csv(per_drug_counts, file.path("drug_effects", "drug_group_sizes_SLE.csv"))

# ---- 3) SLE-only LME: Time × Drug interaction per metabolite ---------------
time_levels <- c("T1","T2","T3","T4","T5","T6")

df_fit_sle <- int_drug %>%
  filter(Condition == "SLE", Metabolite %in% good_mets) %>%
  mutate(Time_f = factor(Timepoint, levels = time_levels))

time_by_drug_all <- purrr::map_dfr(drug_flags, function(dv) {
  dat <- df_fit_sle %>% filter(.data[[dv]] %in% c(0,1))
  # require minimal representation across patients
  min_pos <- dat %>%
    distinct(Patient_ID, .keep_all = TRUE) %>%
    summarise(n_pos = sum(.data[[dv]] == 1, na.rm = TRUE),
              n_neg = sum(.data[[dv]] == 0, na.rm = TRUE)) %>% as.list()
  if (min_pos$n_pos < min_per_group || min_pos$n_neg < min_per_group) return(tibble())
  res <- dat %>%
    group_by(Metabolite) %>%
    group_modify(~{
      fml <- as.formula(paste0("Intensity ~ Time_f * ", dv, " + (1|Patient_ID)"))
      m <- tryCatch(lmer(fml, data = .x, REML = TRUE), error = function(e) NULL)
      if (is.null(m)) return(tibble())
      an <- tryCatch(anova(m, type = 3), error = function(e) NULL)
      if (is.null(an)) return(tibble())
      an_tidy <- broom.mixed::tidy(an) %>% janitor::clean_names()
      p_int <- an_tidy %>% filter(str_detect(term, "Time_f:")) %>% pull(p_value) %>% dplyr::first()
      r2 <- performance::r2(m)
      tibble(
        p_interaction = p_int,
        R2_marginal   = r2$R2_marginal,
        R2_conditional= r2$R2_conditional
      )
    }) %>% ungroup() %>%
    mutate(drug = dv)
  res
})

time_by_drug_all <- time_by_drug_all %>%
  group_by(drug) %>%
  mutate(q_interaction = p.adjust(p_interaction, method = "BH")) %>%
  ungroup()

readr::write_csv(time_by_drug_all, file.path("drug_effects", "LME_TimeXDrug_SLE_all.csv"))

# ---- 4) Plots: per-drug quick volcano (T1 split) ---------------------------
plot_drug_t1_volcano <- function(dv, fdr = 0.1) {
  tab <- t1_results_all %>% filter(drug == dv, !is.na(p_t))
  if (!nrow(tab)) return(NULL)
  ggplot(tab, aes(x = log2FC, y = -log10(p_t))) +
    geom_point(aes(color = q_t < fdr), size = 2) +
    scale_color_manual(values = c(`TRUE` = "#E1812C", `FALSE` = "grey70"), guide = "none") + # orange for sig
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      title = paste0("T1 SLE split by ", dv),
      x = "log2(drug+ / drug-)",
      y = expression(-log[10]~italic(p)[Welch])
    ) +
    theme_minimal(base_size = 12)
}
for (dv in unique(t1_results_all$drug)) {
  p <- plot_drug_t1_volcano(dv)
  if (!is.null(p)) ggsave(file.path("drug_effects/plots", paste0("T1_volcano_", dv, ".pdf")),
                          p, width = 6.5, height = 5)
}
drug_labels <- c(
  "cell_cept_myfortic" = "CellCept/Myfortic",
  "prednisolone"       = "Prednisolone",
  "hydroxychloroquine" = "Hydroxychloroquine"
  # add others as needed
)
print(p)

# ---- 5) Plots: EMM curves for top hits per drug (interaction) --------------
# Neutral palette to avoid HC/SLE colors: grey vs orange
drug_palette <- c("No" = "grey40", "Yes" = "#E1812C")

plot_emm_sle_by_drug <- function(target_met, dv) {
  dat <- df_fit_sle %>%
    filter(Metabolite == target_met) %>%
    filter(.data[[dv]] %in% c(0,1)) %>%
    mutate(Drug = factor(.data[[dv]], levels = c(0,1), labels = c("No","Yes")))
  if (nrow(dat) < 8) return(NULL)
  m <- tryCatch(lmer(as.formula(paste0("Intensity ~ Time_f * Drug + (1|Patient_ID)")),
                     data = dat, REML = TRUE), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  emm_df <- as.data.frame(emmeans(m, ~ Drug | Time_f))
  ggplot(dat, aes(x = Time_f, y = Intensity, group = Patient_ID, color = Drug)) +
    geom_point(alpha = 0.35, size = 1.5) +
    geom_line(alpha = 0.25, linewidth = 0.5) +
    geom_point(data = emm_df, aes(x = Time_f, y = emmean, color = Drug),
               inherit.aes = FALSE, size = 2.7, shape = 17) +
    geom_line(data = emm_df, aes(x = Time_f, y = emmean, group = Drug, color = Drug),
              inherit.aes = FALSE, linewidth = 0.8) +
    geom_errorbar(data = emm_df,
                  aes(x = Time_f, y = emmean, ymin = lower.CL, ymax = upper.CL, color = Drug),
                  inherit.aes = FALSE, width = 0.1, alpha = 0.8) +
    scale_color_manual(values = drug_palette) +
    labs(
      title = glue("{target_met} — EMMs by {drug_labels[[dv]] %||% dv} (SLE only)"),
      subtitle = subtitle_txt,
      ...
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

top_k <- 6
#for (dv in unique(time_by_drug_all$drug)) {
#  tab <- time_by_drug_all %>%
#    filter(drug == dv, !Metabolite %in% drug_feature_names) %>%
#    arrange(q_interaction)
#  if (!nrow(tab)) next
#  mets <- head(tab$Metabolite, top_k)
#  for (met in mets) {
#    p <- plot_emm_sle_by_drug(met, dv)
#    if (!is.null(p)) {
#      safe <- janitor::make_clean_names(glue("{dv}_{met}")) %>% stringr::str_trunc(120, ellipsis = "")
#      ggsave(file.path("drug_effects/plots", paste0("EMM_", safe, ".pdf")), p, width = 7.2, height = 5.2)
#    }
#  }
#}
# gives error: Error in labs(title = glue("{target_met} — EMMs by {drug_labels[[dv]] %||% dv} (SLE only)"),  : '...' used in an incorrect context



#####################################################################################################
# ==== Supplementary multi-page PDF: EMM curves for top hits per drug ====
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(emmeans)
  library(lme4)
  library(glue)
  library(grid)
})

# Unicode superscript mapping
superscript <- function(x) {
  x <- as.character(x)
  chars <- strsplit(x, "")[[1]]
  paste0(sapply(chars, function(ch) {
    switch(ch,
           "0" = "\u2070",
           "1" = "\u00B9",
           "2" = "\u00B2",
           "3" = "\u00B3",
           "4" = "\u2074",
           "5" = "\u2075",
           "6" = "\u2076",
           "7" = "\u2077",
           "8" = "\u2078",
           "9" = "\u2079",
           "-" = "\u207B",
           ch)
  }), collapse = "")
}

scientific_10x <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      "NA"                # will appear as text
    } else if (val == 0) {
      "0"                 # zero as valid expression
    } else {
      exponent <- floor(log10(abs(val)))
      mantissa <- round(val / 10^exponent, 2)
      paste0(mantissa, " %*% 10^", exponent)
    }
  }) |> parse(text = _)
}


# Helper: nice drug label
nice_drug_label <- function(dv) {
  drug_labels <- c(
    "cell_cept_myfortic" = "CellCept/Myfortic",
    "prednisolone"       = "Prednisolone",
    "hydroxychloroquine" = "Hydroxychloroquine"
  )
  if (dv %in% names(drug_labels)) drug_labels[[dv]] else dv
}

plot_emm_sle_by_drug_annotated <- function(target_met, dv) {
  dat <- df_fit_sle %>%
    filter(Metabolite == target_met) %>%
    filter(.data[[dv]] %in% c(0,1)) %>%
    mutate(Drug = factor(.data[[dv]], levels = c(0,1), labels = c("No","Yes")))
  
  if (nrow(dat) < 8) return(NULL)
  
  # Fit linear mixed model
  m <- tryCatch(
    lmer(Intensity ~ Time_f * Drug + (1 | Patient_ID), data = dat, REML = TRUE),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  # Estimated marginal means
  emm_df <- as.data.frame(emmeans(m, ~ Drug | Time_f))
  
  # R² and q-value (still calculated but not shown in plot)
  r2 <- performance::r2(m)
  r2_m <- suppressWarnings(if (!is.null(r2$R2_marginal)) as.numeric(r2$R2_marginal[[1]]) else NA_real_)
  r2_c <- suppressWarnings(if (!is.null(r2$R2_conditional)) as.numeric(r2$R2_conditional[[1]]) else NA_real_)
  
  q_val <- time_by_drug_all %>%
    filter(drug == dv, Metabolite == target_met) %>%
    pull(q_interaction)
  q_val <- if (length(q_val)) q_val[1] else NA_real_
  
  drug_label <- nice_drug_label(dv)
  
  # Plot — no title, no subtitle
  ggplot(dat, aes(x = Time_f, y = Intensity, group = Patient_ID, color = Drug)) +
    geom_point(alpha = 0.35, size = 1.5) +
    geom_line(alpha = 0.25, linewidth = 0.5) +
    geom_point(data = emm_df, aes(x = Time_f, y = emmean, color = Drug),
               inherit.aes = FALSE, size = 2.7, shape = 17) +
    geom_line(data = emm_df, aes(x = Time_f, y = emmean, group = Drug, color = Drug),
              inherit.aes = FALSE, linewidth = 0.8) +
    geom_errorbar(data = emm_df,
                  aes(x = Time_f, y = emmean, ymin = lower.CL, ymax = upper.CL, color = Drug),
                  inherit.aes = FALSE, width = 0.10, alpha = 0.8) +
    scale_color_manual(values = c("No" = "grey40", "Yes" = "#E1812C")) +
    scale_y_continuous(labels = scientific_10x) +
    labs(
      title = NULL,        # remove main title
      subtitle = NULL,     # remove subtitle
      x = NULL, y = "Intensity", color = "Drug"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}



# Open the PDF and capture the device number
outfile <- "Supplementary_Figure_Sx_EMM_DrugTime_ALL_perDrug.pdf"
if (capabilities("cairo")) {
  grDevices::cairo_pdf(outfile, width = 7.2, height = 5.2, onefile = TRUE)
} else {
  grDevices::pdf(outfile, width = 7.2, height = 5.2, useDingbats = FALSE, onefile = TRUE)
}
pdf_dev <- grDevices::dev.cur()  # <-- remember which device is the PDF

page_ct <- 0L; fail_ct <- 0L
for (dv in unique(time_by_drug_all$drug)) {
  tab <- time_by_drug_all %>%
    dplyr::filter(drug == dv, !Metabolite %in% drug_feature_names) %>%
    dplyr::arrange(is.na(q_interaction), q_interaction)
  if (!nrow(tab)) next
  
  for (met in tab$Metabolite) {
    p <- try(plot_emm_sle_by_drug_annotated(met, dv), silent = TRUE)
    if (inherits(p, "try-error") || is.null(p)) { fail_ct <- fail_ct + 1L; next }
    grDevices::dev.set(pdf_dev)   # <-- make sure the PDF is active
    print(p)                      # each print() adds a new page
    page_ct <- page_ct + 1L
  }
}

# Close THAT device (don’t rely on the "current" device)
grDevices::dev.off(pdf_dev)
message(sprintf("Wrote %d pages to %s (skipped %d)", page_ct, normalizePath(outfile), fail_ct))

table(time_by_drug_all$drug)
dv <- unique(time_by_drug_all$drug)[1]
tab <- subset(time_by_drug_all, drug == dv & !(Metabolite %in% drug_feature_names))
sum(vapply(tab$Metabolite, function(m) !is.null(plot_emm_sle_by_drug_annotated(m, dv)), logical(1)))
print(plot_emm_sle_by_drug_annotated(tab$Metabolite[1], dv))
print(plot_emm_sle_by_drug_annotated(tab$Metabolite[2], dv))

#######
## AMP Cell cept vis
# Pick metabolite and drug
target_met <- "AMP"
dv <- "cell_cept_myfortic"

# Generate the plot
p_amp <- plot_emm_sle_by_drug_annotated(target_met, dv)

# Save as PNG for main manuscript (same dimensions as supplement PDFs)
if (!is.null(p_amp)) {
  ggsave("Figure_AMP_CellCept-Myfortic.png",
         plot   = p_amp,
         width  = 7.2, 
         height = 5.2,
         dpi    = 300)
}

