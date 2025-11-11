# ---- Load packages ----
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# ---- Define base directory ----
if (.Platform$OS.type == "unix") {
  Classen_data <- file.path(
    "Filepath_for_macOS",
    "mass_error_distribution.csv"
  )
} else {
  Classen_data <- file.path(
    "Filepath_for_Windows",
    "mass_error_distribution.csv"
  )
}

# ---- Import CSV ----
mass_data <- read_csv(Classen_data)

# ---- Clean and convert Dalton Mass Error safely ----
mass_data <- mass_data %>%
  mutate(
    `Dalton Mass Error` = str_trim(`Dalton Mass Error`),
    # Replace Unicode minus (U+2212) with ASCII minus
    `Dalton Mass Error` = str_replace_all(`Dalton Mass Error`, "\u2212", "-"),
    # Replace variants like "×10^" or "x10^" with scientific notation 'e'
    `Dalton Mass Error` = str_replace_all(`Dalton Mass Error`, "[×xX]\\s*10\\^", "e"),
    # Remove spaces and commas
    `Dalton Mass Error` = str_replace_all(`Dalton Mass Error`, "[ ,]", ""),
    # Convert to numeric (after all formatting normalized)
    `Dalton Mass Error` = as.numeric(`Dalton Mass Error`)
  )

# ---- Verify conversion ----
summary(mass_data$`Dalton Mass Error`)
range(mass_data$`Dalton Mass Error`, na.rm = TRUE)

# ---- Compute mass-error statistics ----
# ---- Compute mass-error statistics (optional for annotation) ----
mass_summary <- mass_data %>%
  summarise(
    mean_error = mean(`Dalton Mass Error`, na.rm = TRUE),
    median_error = median(`Dalton Mass Error`, na.rm = TRUE),
    sd_error = sd(`Dalton Mass Error`, na.rm = TRUE),
    n = sum(!is.na(`Dalton Mass Error`))
  )

# ---- Refined publication-quality plot ----
p <- ggplot(mass_data, aes(x = `Dalton Mass Error`)) +
  geom_histogram(
    bins = 40,
    fill = "#4C72B0",        # subtle blue
    color = "black",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "#DD1C77",
    size = 1
  ) +
  labs(
    title = "Mass-Error Distribution of Annotated Features",
    x = "Mass Error (Da)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )
print(p)
# ---- Save as high-resolution PNG ----
ggsave(
  filename = "mass_error_distribution.png",
  plot = p,
  width = 6,      # width in inches
  height = 4,     # height in inches
  dpi = 300       # publication-quality resolution
)
