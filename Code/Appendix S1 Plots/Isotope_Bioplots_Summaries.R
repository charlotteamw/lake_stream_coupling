
# =============================================================================
# SUPPLEMENTARY ISOTOPE ANALYSES AND BIPLOTS
#
# Common and Golden Shiner Comparisons
# Lake–Stream Coupling Manuscript
#
# Analyses:
#   1. Raw isotope summary statistics
#   2. Species comparisons: delta13C and delta15N
#   3. Model diagnostics and estimated marginal means
#
# Plots:
#   - Raw delta13C and delta15N distributions
#   - Shiner species isotope biplot (mean ± SD)
#   - Shiner species biplot by habitat and season
#   - Combined shiner individuals + baseline means ± SD
#   - All consumer species means ± SD + baselines
#
# Data: Data/iso_metadata.csv
#
# Run from the lake_creek_coupling-main project root.
# =============================================================================


# =============================================================================
# 1. SETUP
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(ggpubr)
library(ggrepel)

file_path <- getwd()

data_path <- file.path(
  file_path,
  "Data",
  "iso_metadata.csv"
)

output_dir <- file.path(
  file_path,
  "Figures",
  "Supplementary_Shiner_Isotopes"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

options(
  contrasts = c("contr.sum", "contr.poly")
)

month_levels <- c(
  "may",
  "august",
  "october"
)

month_labels <- c(
  "Spring",
  "Summer",
  "Fall"
)

# Lake–stream colours
lake_stream_colors <- c(
  "Stream" = "#F06C57",
  "Lake" = "#4A90B8"
)

# Species colours
species_colors <- c(
  "common shiner" = "slategray3",
  "golden shiner" = "goldenrod2"
)

species_labels <- c(
  "common shiner" = "Common shiner",
  "golden shiner" = "Golden shiner"
)

custom_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 12)
  )


# =============================================================================
# 2. LOAD DATA
# =============================================================================

if (!file.exists(data_path)) {
  stop(
    "Cannot find Data/iso_metadata.csv. ",
    "Run this script from the project root."
  )
}

isotope_data <- read.csv(
  data_path,
  header = TRUE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "sample_ID",
  "organism",
  "tissue",
  "month",
  "location",
  "year",
  "d13C",
  "d15N"
)



isotope_data <- isotope_data %>%
  mutate(
    organism = str_to_lower(str_trim(organism)),
    tissue = str_to_lower(str_trim(tissue)),
    location = str_to_lower(str_trim(location)),
    location = if_else(location == "creek", "stream", location),
    month = str_to_lower(str_trim(month)),
    d13C = as.numeric(as.character(d13C)),
    d15N = as.numeric(as.character(d15N))
  )



# =============================================================================
# 3. SAMPLING AREA ASSIGNMENT
# =============================================================================

# Area assignments from the original isotope analysis.
# Sampling identifiers are retained as recorded in the dataset.

assign_area <- function(sample_id) {
  
  sample_id <- trimws(as.character(sample_id))
  
  sample_prefix <- toupper(
    stringr::str_extract(
      sample_id,
      "^[A-Za-z]+[0-9]*"
    )
  )
  
  sample_prefix <- case_when(
    
    str_detect(
      toupper(sample_id),
      "^(WU|WD)"
    ) |
      sample_prefix %in% c(
        "DN2", "UP1", "UP3"
      ) ~ "CR5",
    
    str_detect(
      toupper(sample_id),
      "^(FU|FD)|^SMB1[[:space:]]+FU(?:[-_]|$)"
    ) ~ "CR1",
    
    sample_prefix == "FN" ~ "CR1",
    
    TRUE ~ sample_prefix
  )
  
  factor(
    case_when(
      
      is.na(sample_id) ~ NA_character_,
      
      sample_prefix %in% c(
        "CR8",
        "CR10",
        paste0("CC", 18:24)
      ) ~ "CR_Area3",
      
      sample_prefix %in% c(
        "FYKE",
        "CR1",
        "CR2",
        "CR3",
        "CC",
        "CC2"
      ) ~ "CR_Area1",
      
      sample_prefix %in% c(
        "CR5",
        "EF4",
        "WEIR",
        "UP",
        "UP2",
        "CR4",
        "CR6",
        "CR7",
        "CR9",
        paste0("CC", 5:14),
        "CC16",
        "CC25"
      ) ~ "CR_Area2",
      
      sample_prefix %in% c(
        "LK4",
        "LK10",
        "LK5",
        "LK6",
        "LK3",
        "OP3",
        "OP4",
        "OP5",
        "OP6",
        "OP10",
        "OP11",
        "OP12"
      ) ~ "LK_Area1",
      
      sample_prefix %in% c(
        "OP1",
        "OP2",
        "LK",
        "LK1",
        "LK2"
      ) ~ "LK_Area2",
      
      sample_prefix %in% c(
        "OP9",
        "OP8",
        "LK9",
        "LK8",
        "OP17",
        "LK7",
        "OP7"
      ) ~ "LK_Area3",
      
      TRUE ~ NA_character_
    )
  )
}


# =============================================================================
# 4. PREPARE SHINER DATA
# =============================================================================

shiners_raw <- isotope_data %>%
  filter(
    organism %in% c(
      "common shiner",
      "golden shiner"
    ),
    tissue == "liver"
  ) %>%
  mutate(
    area = assign_area(sample_ID),
    
    organism = factor(
      organism,
      levels = c(
        "common shiner",
        "golden shiner"
      )
    ),
    
    location = factor(
      location,
      levels = c(
        "stream",
        "lake"
      )
    ),
    
    month = factor(
      month,
      levels = month_levels,
      labels = month_labels
    ),
    
    year = factor(year)
  )


# Check for duplicate sample identifiers
duplicate_samples <- shiners_raw %>%
  count(sample_ID) %>%
  filter(n > 1)

if (nrow(duplicate_samples) > 0) {
  warning("Duplicate sample IDs detected.")
  print(duplicate_samples)
}


# Check unresolved area assignments
missing_areas <- shiners_raw %>%
  filter(is.na(area)) %>%
  distinct(
    sample_ID,
    organism,
    location
  )

if (nrow(missing_areas) > 0) {
  
  print(missing_areas)
  
  stop(
    "Some shiner samples have unresolved areas. ",
    "Verify area assignments before modelling."
  )
}


# Sample sizes
sample_sizes <- shiners_raw %>%
  count(
    organism,
    month,
    location,
    year,
    name = "n"
  )

print(sample_sizes, n = Inf)

write.csv(
  sample_sizes,
  file.path(
    output_dir,
    "shiner_species_sample_sizes.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 5. RAW ISOTOPE SUMMARY STATISTICS
# =============================================================================

species_summary <- shiners_raw %>%
  group_by(
    organism,
    location,
    month
  ) %>%
  summarise(
    n_d13C = sum(!is.na(d13C)),
    mean_d13C = mean(d13C, na.rm = TRUE),
    sd_d13C = sd(d13C, na.rm = TRUE),
    
    n_d15N = sum(!is.na(d15N)),
    mean_d15N = mean(d15N, na.rm = TRUE),
    sd_d15N = sd(d15N, na.rm = TRUE),
    
    .groups = "drop"
  )

print(species_summary, n = Inf)

write.csv(
  species_summary,
  file.path(
    output_dir,
    "shiner_species_isotope_summary.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 6. SPECIES COMPARISON MODELS
# =============================================================================

# Models test whether raw isotope values differ between species,
# accounting for season, habitat, year and sampling area.
#
# The species effect represents an adjusted average difference.
# Interactions involving species are not included in these models.

model_data <- shiners_raw %>%
  filter(
    !is.na(d13C),
    !is.na(d15N),
    !is.na(organism),
    !is.na(month),
    !is.na(location),
    !is.na(year),
    !is.na(area)
  ) %>%
  droplevels()

if (n_distinct(model_data$organism) != 2) {
  stop("Both shiner species must be represented.")
}

if (n_distinct(model_data$area) < 2) {
  stop("At least two sampling areas are required.")
}

cat("\nModel sample sizes:\n")

print(
  model_data %>%
    count(
      organism,
      month,
      location,
      year
    ),
  n = Inf
)


# -----------------------------------------------------------------------------
# 6A. DELTA13C
# -----------------------------------------------------------------------------

model_d13C <- glmmTMB(
  d13C ~ organism + month + location + year + (1 | area),
  data = model_data,
  family = gaussian()
)

cat("\n========== DELTA13C ==========\n")

print(summary(model_d13C))

anova_d13C <- as.data.frame(
  car::Anova(
    model_d13C,
    type = 2
  )
) %>%
  tibble::rownames_to_column("Predictor")

print(anova_d13C)

emm_d13C <- emmeans(
  model_d13C,
  ~ organism
)

contrast_d13C <- as.data.frame(
  summary(
    pairs(emm_d13C),
    infer = c(TRUE, TRUE)
  )
)

print(contrast_d13C)

write.csv(
  anova_d13C,
  file.path(
    output_dir,
    "model_d13C_anova.csv"
  ),
  row.names = FALSE
)

write.csv(
  contrast_d13C,
  file.path(
    output_dir,
    "model_d13C_species_contrast.csv"
  ),
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# 6B. DELTA15N
# -----------------------------------------------------------------------------

model_d15N <- glmmTMB(
  d15N ~ organism + month + location + year + (1 | area),
  data = model_data,
  family = gaussian()
)

cat("\n========== DELTA15N ==========\n")

print(summary(model_d15N))

anova_d15N <- as.data.frame(
  car::Anova(
    model_d15N,
    type = 2
  )
) %>%
  tibble::rownames_to_column("Predictor")

print(anova_d15N)

emm_d15N <- emmeans(
  model_d15N,
  ~ organism
)

contrast_d15N <- as.data.frame(
  summary(
    pairs(emm_d15N),
    infer = c(TRUE, TRUE)
  )
)

print(contrast_d15N)

write.csv(
  anova_d15N,
  file.path(
    output_dir,
    "model_d15N_anova.csv"
  ),
  row.names = FALSE
)

write.csv(
  contrast_d15N,
  file.path(
    output_dir,
    "model_d15N_species_contrast.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 7. MODEL DIAGNOSTICS
# =============================================================================

set.seed(123)

sim_d13C <- DHARMa::simulateResiduals(
  model_d13C,
  n = 1000
)

sim_d15N <- DHARMa::simulateResiduals(
  model_d15N,
  n = 1000
)

pdf(
  file.path(
    output_dir,
    "model_d13C_diagnostics.pdf"
  ),
  width = 9,
  height = 7
)

plot(sim_d13C)

dev.off()


pdf(
  file.path(
    output_dir,
    "model_d15N_diagnostics.pdf"
  ),
  width = 9,
  height = 7
)

plot(sim_d15N)

dev.off()


cat("\nDelta13C diagnostics:\n")

print(DHARMa::testUniformity(sim_d13C))
print(DHARMa::testDispersion(sim_d13C))
print(DHARMa::testOutliers(sim_d13C))

cat("\nDelta15N diagnostics:\n")

print(DHARMa::testUniformity(sim_d15N))
print(DHARMa::testDispersion(sim_d15N))
print(DHARMa::testOutliers(sim_d15N))


# =============================================================================
# 8. RAW ISOTOPE DISTRIBUTIONS
# =============================================================================

plot_d13C <- ggplot(
  model_data,
  aes(
    x = organism,
    y = d13C,
    fill = organism
  )
) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.35,
    size = 1.5
  ) +
  scale_fill_manual(
    values = species_colors
  ) +
  scale_x_discrete(
    labels = species_labels
  ) +
  labs(
    x = NULL,
    y = expression(delta^13*C~("\u2030"))
  ) +
  custom_theme +
  theme(
    legend.position = "none"
  )


plot_d15N <- ggplot(
  model_data,
  aes(
    x = organism,
    y = d15N,
    fill = organism
  )
) +
  geom_boxplot(
    width = 0.5,
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.35,
    size = 1.5
  ) +
  scale_fill_manual(
    values = species_colors
  ) +
  scale_x_discrete(
    labels = species_labels
  ) +
  labs(
    x = NULL,
    y = expression(delta^15*N~("\u2030"))
  ) +
  custom_theme +
  theme(
    legend.position = "none"
  )


figure_species_isotopes <- ggpubr::ggarrange(
  plot_d13C,
  plot_d15N,
  ncol = 2,
  labels = c("A", "B")
)

print(figure_species_isotopes)

ggsave(
  file.path(
    output_dir,
    "raw_shiner_isotopes.png"
  ),
  plot = figure_species_isotopes,
  width = 9,
  height = 4.5,
  dpi = 600
)


# =============================================================================
# 9. SHINER SPECIES BIPLOT
# Individual observations + species means ± SD
# =============================================================================

species_biplot_summary <- model_data %>%
  group_by(organism) %>%
  summarise(
    d13C_mean = mean(d13C),
    d13C_sd = sd(d13C),
    d15N_mean = mean(d15N),
    d15N_sd = sd(d15N),
    .groups = "drop"
  )


plot_species_biplot <- ggplot() +
  
  # Individual isotope observations
  geom_point(
    data = model_data,
    aes(
      x = d13C,
      y = d15N,
      color = organism
    ),
    alpha = 0.35,
    size = 1.8
  ) +
  
  # Vertical SD
  geom_errorbar(
    data = species_biplot_summary,
    aes(
      x = d13C_mean,
      ymin = d15N_mean - d15N_sd,
      ymax = d15N_mean + d15N_sd,
      color = organism
    ),
    width = 0,
    linewidth = 0.9
  ) +
  
  # Horizontal SD
  geom_segment(
    data = species_biplot_summary,
    aes(
      x = d13C_mean - d13C_sd,
      xend = d13C_mean + d13C_sd,
      y = d15N_mean,
      yend = d15N_mean,
      color = organism
    ),
    linewidth = 0.9
  ) +
  
  # Species mean isotope values
  geom_point(
    data = species_biplot_summary,
    aes(
      x = d13C_mean,
      y = d15N_mean,
      fill = organism
    ),
    shape = 21,
    size = 4,
    color = "black",
    stroke = 0.6
  ) +
  
  scale_color_manual(
    values = species_colors,
    labels = species_labels,
    name = "Species"
  ) +
  
  scale_fill_manual(
    values = species_colors,
    guide = "none"
  ) +
  
  labs(
    x = expression(delta^13*C~("\u2030")),
    y = expression(delta^15*N~("\u2030"))
  ) +
  
  custom_theme

print(plot_species_biplot)

ggsave(
  file.path(
    output_dir,
    "shiner_species_biplot.png"
  ),
  plot = plot_species_biplot,
  width = 6,
  height = 5,
  dpi = 600
)


# =============================================================================
# 10. SHINER SPECIES BIPLOT BY HABITAT AND SEASON
# Individual observations + species means ± SD
# =============================================================================

# Calculate species-specific isotope means and SDs
# within each habitat and season.

habitat_species_summary <- model_data %>%
  group_by(
    location,
    organism
  ) %>%
  summarise(
    n = n(),
    d13C_mean = mean(d13C),
    d13C_sd = sd(d13C),
    d15N_mean = mean(d15N),
    d15N_sd = sd(d15N),
    .groups = "drop"
  )

print(
  habitat_species_summary,
  n = Inf
)

write.csv(
  habitat_species_summary,
  file.path(
    output_dir,
    "shiner_isotope_means_by_habitat_species.csv"
  ),
  row.names = FALSE
)


plot_habitat_sp_biplot <- ggplot() +
  
  # Individual isotope observations
  geom_point(
    data = model_data,
    aes(
      x = d13C,
      y = d15N,
      color = organism
    ),
    alpha = 0.35,
    size = 1.7
  ) +
  
  # Vertical SD: delta15N
  geom_errorbar(
    data = habitat_species_summary,
    aes(
      x = d13C_mean,
      ymin = d15N_mean - d15N_sd,
      ymax = d15N_mean + d15N_sd,
      color = organism
    ),
    width = 0,
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  
  # Horizontal SD: delta13C
  geom_segment(
    data = habitat_species_summary,
    aes(
      x = d13C_mean - d13C_sd,
      xend = d13C_mean + d13C_sd,
      y = d15N_mean,
      yend = d15N_mean,
      color = organism
    ),
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  
  # Species mean isotope values
  geom_point(
    data = habitat_species_summary,
    aes(
      x = d13C_mean,
      y = d15N_mean,
      fill = organism
    ),
    shape = 21,
    color = "black",
    stroke = 0.6,
    size = 3.5
  ) +
  
  # Habitat panels, 
  facet_grid(
    cols  = vars(location),
    labeller = labeller(
      location = c(
        stream = "Stream",
        lake = "Lake"
      )
    )
  ) +
  
  scale_color_manual(
    values = species_colors,
    labels = species_labels,
    name = "Species"
  ) +
  
  scale_fill_manual(
    values = species_colors,
    guide = "none"
  ) +
  
  labs(
    x = expression(delta^13*C~("\u2030")),
    y = expression(delta^15*N~("\u2030"))
  ) +
  
  custom_theme

print(plot_habitat_sp_biplot)

ggsave(
  file.path(
    output_dir,
    "shiner_species_biplot_by_habitat.png"
  ),
  plot = plot_habitat_sp_biplot,
  width = 10,
  height = 7,
  dpi = 600
)


# =============================================================================
# 11. SHARED DATA FOR GENERAL ISOTOPE BIPLOTS
# =============================================================================

# Baselines: mayfly and mussel, means ± SD
# Shiners: individual liver isotope values
# Consumers: species-specific means ± SD within each habitat


# Baseline means ± SD
baseline_means_sd <- isotope_data %>%
  filter(
    organism %in% c(
      "mayfly",
      "mussel"
    )
  ) %>%
  mutate(
    organism = str_to_title(organism),
    location = str_to_title(location)
  ) %>%
  group_by(
    location,
    organism
  ) %>%
  summarise(
    d13C_mean = mean(d13C, na.rm = TRUE),
    d13C_sd = sd(d13C, na.rm = TRUE),
    d15N_mean = mean(d15N, na.rm = TRUE),
    d15N_sd = sd(d15N, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = organism
  )


# Combined shiner individual isotope values
shiner_plot_df <- isotope_data %>%
  filter(
    organism %in% c(
      "golden shiner",
      "common shiner"
    ),
    tissue == "liver"
  ) %>%
  mutate(
    location = str_to_title(location)
  )


# Original consumer species abbreviations
abbr <- c(
  "smallmouth bass" = "SMB",
  "golden shiner" = "GS",
  "creek chub" = "CC",
  "common shiner" = "CS",
  "pumpkinseed" = "PS",
  "yellow perch" = "YP",
  "brook trout" = "BT",
  "pearl dace" = "PD",
  "blacknose shiner" = "BNS",
  "northern redbelly dace" = "NRD",
  "bluntnose minnow" = "BNM"
)


# Consumer species means ± SD
species_means_sd <- isotope_data %>%
  filter(
    organism %in% names(abbr)
  ) %>%
  mutate(
    location = str_to_title(location)
  ) %>%
  group_by(
    location,
    organism
  ) %>%
  summarise(
    d13C_mean = mean(d13C, na.rm = TRUE),
    d13C_sd = sd(d13C, na.rm = TRUE),
    d15N_mean = mean(d15N, na.rm = TRUE),
    d15N_sd = sd(d15N, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label_abbr = unname(abbr[organism])
  )


# Save descriptive biplot summaries
write.csv(
  baseline_means_sd,
  file.path(
    output_dir,
    "baseline_means_sd.csv"
  ),
  row.names = FALSE
)

write.csv(
  species_means_sd,
  file.path(
    output_dir,
    "all_consumer_species_means_sd.csv"
  ),
  row.names = FALSE
)


# =============================================================================
# 12. COMBINED SHINER INDIVIDUALS + BASELINE MEANS ± SD
# =============================================================================

plot_biplot_shiners <- ggplot() +
  
  # Individual common and golden shiners combined
  geom_point(
    data = shiner_plot_df,
    aes(
      x = d13C,
      y = d15N,
      color = location
    ),
    alpha = 0.6,
    size = 2
  ) +
  
  # Baseline vertical error bars
  geom_errorbar(
    data = baseline_means_sd,
    aes(
      x = d13C_mean,
      ymin = d15N_mean - d15N_sd,
      ymax = d15N_mean + d15N_sd,
      color = location
    ),
    width = 0,
    linewidth = 0.6
  ) +
  
  # Baseline horizontal error bars
  geom_segment(
    data = baseline_means_sd,
    aes(
      x = d13C_mean - d13C_sd,
      xend = d13C_mean + d13C_sd,
      y = d15N_mean,
      yend = d15N_mean,
      color = location
    ),
    linewidth = 0.6
  ) +
  
  # Baseline mean points
  geom_point(
    data = baseline_means_sd,
    aes(
      x = d13C_mean,
      y = d15N_mean,
      color = location,
      shape = organism
    ),
    size = 2.8,
    stroke = 1.1,
    fill = "white"
  ) +
  
  scale_color_manual(
    values = lake_stream_colors,
    name = "Location"
  ) +
  
  scale_shape_manual(
    values = c(
      "Mayfly" = 21,
      "Mussel" = 24
    ),
    name = "Baseline"
  ) +
  
  labs(
    x = expression(delta^13*C~("\u2030")),
    y = expression(delta^15*N~("\u2030"))
  ) +
  
  custom_theme

print(plot_biplot_shiners)

ggsave(
  file.path(
    output_dir,
    "combined_shiners_with_baselines.png"
  ),
  plot = plot_biplot_shiners,
  width = 6,
  height = 5,
  dpi = 600
)


# =============================================================================
# 13. ALL CONSUMER SPECIES MEANS ± SD + BASELINE MEANS ± SD
# =============================================================================

plot_biplot_allspecies <- ggplot() +
  
  # Baseline vertical error bars
  geom_errorbar(
    data = baseline_means_sd,
    aes(
      x = d13C_mean,
      ymin = d15N_mean - d15N_sd,
      ymax = d15N_mean + d15N_sd,
      color = location
    ),
    width = 0,
    linewidth = 0.4,
    alpha = 0.6
  ) +
  
  # Baseline horizontal error bars
  geom_segment(
    data = baseline_means_sd,
    aes(
      x = d13C_mean - d13C_sd,
      xend = d13C_mean + d13C_sd,
      y = d15N_mean,
      yend = d15N_mean,
      color = location
    ),
    linewidth = 0.4,
    alpha = 0.6
  ) +
  
  # Baseline mean points
  geom_point(
    data = baseline_means_sd,
    aes(
      x = d13C_mean,
      y = d15N_mean,
      color = location
    ),
    shape = 21,
    size = 2.4,
    stroke = 0.9,
    fill = "white"
  ) +
  
  # Baseline labels
  ggrepel::geom_label_repel(
    data = baseline_means_sd,
    aes(
      x = d13C_mean,
      y = d15N_mean,
      label = label,
      color = location
    ),
    size = 3.0,
    label.size = 0,
    fill = scales::alpha("white", 0.7),
    box.padding = 0.3,
    point.padding = 0.2,
    seed = 123,
    show.legend = FALSE
  ) +
  
  # Consumer vertical error bars
  geom_errorbar(
    data = species_means_sd,
    aes(
      x = d13C_mean,
      ymin = d15N_mean - d15N_sd,
      ymax = d15N_mean + d15N_sd,
      color = location
    ),
    width = 0,
    linewidth = 0.6,
    alpha = 0.75
  ) +
  
  # Consumer horizontal error bars
  geom_segment(
    data = species_means_sd,
    aes(
      x = d13C_mean - d13C_sd,
      xend = d13C_mean + d13C_sd,
      y = d15N_mean,
      yend = d15N_mean,
      color = location
    ),
    linewidth = 0.6,
    alpha = 0.75
  ) +
  
  # Consumer mean points
  geom_point(
    data = species_means_sd,
    aes(
      x = d13C_mean,
      y = d15N_mean,
      color = location
    ),
    size = 3.2,
    stroke = 0.7
  ) +
  
  # Consumer labels
  ggrepel::geom_label_repel(
    data = species_means_sd,
    aes(
      x = d13C_mean,
      y = d15N_mean,
      label = label_abbr,
      color = location
    ),
    size = 3.2,
    label.size = 0,
    segment.size = 0.3,
    fill = scales::alpha("white", 0.65),
    box.padding = 0.4,
    point.padding = 0.35,
    max.overlaps = 200,
    seed = 42,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = lake_stream_colors,
    name = "Location"
  ) +
  
  labs(
    x = expression(delta^13*C~("\u2030")),
    y = expression(delta^15*N~("\u2030"))
  ) +
  
  custom_theme

print(plot_biplot_allspecies)

ggsave(
  file.path(
    output_dir,
    "all_species_with_baselines.png"
  ),
  plot = plot_biplot_allspecies,
  width = 7,
  height = 5.5,
  dpi = 600
)


# =============================================================================
# 14. SAVE MODEL RESULTS
# =============================================================================

sink(
  file.path(
    output_dir,
    "species_model_results.txt"
  )
)

cat(
  "RAW ISOTOPE COMPARISONS\n",
  "COMMON VS GOLDEN SHINERS\n\n"
)

cat("\nDELTA13C MODEL\n")
print(summary(model_d13C))
print(anova_d13C)
print(contrast_d13C)

cat("\nDELTA15N MODEL\n")
print(summary(model_d15N))
print(anova_d15N)
print(contrast_d15N)

sink()



