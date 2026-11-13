# Script for Figure 3 — Seasonal Comparisons of Stable Isotopes & Telemetry variables using glms, for manuscript "Habitat coupling by small fish across lake-stream habitats: a missing link in freshwater food webs?"

# Version: 2026-10-25
# Author(s): Charlotte Ward & Reilly O'Connor

# load Pkgs
library(tidyverse)
library(lubridate)
library(purrr)
library(glatos)
library(ggpubr)
library(cowplot)
library(car)
library(emmeans)
library(glmmTMB)
library(DHARMa)

#Set file paths n ze working directory
file_path <- getwd()

source(file.path(file_path, "Code/0 - Functions.R"))

options(contrasts = c("contr.sum", "contr.poly"))

ensure_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

#Assign ggplot themes, colours and labels...
custom_theme <- theme_bw() +
  theme(
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "black"),
    legend.position = "right",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    axis.title      = element_text(size = 14),
    axis.text.x     = element_text(size = 14),
    axis.text.y     = element_text(size = 12)
  )

fill_scheme <- c("creek" = "#F06C57", "lake" = "#4A90B8")
season_levels <- c("Spring", "Summer", "Fall")
month_levels  <- c("may", "august", "october")
month_labels  <- c("Spring", "Summer", "Fall")

legend_override <- guide_legend(
  override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
)

make_bottom_legend <- function(plot_obj) {
  get_legend(
    plot_obj + theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 12)))
}


# load Data
isotope_data <- read.csv(file.path(file_path, "Data/iso_metadata.csv"), header = TRUE)

shiner_species <- c("golden shiner", "common shiner")

isotope_data <- isotope_data %>%
  mutate(
    sample_prefix = toupper(sub("-.*$", "", trimws(sample_ID))),
    area = assign_area(sample_ID)
  )

print(isotope_data %>% 
        filter(is.na(area)) %>%
        count(organism, sample_prefix, location))

excluded_baselines <- isotope_data %>%
  filter(organism %in% c("mayfly", "mussel"), is.na(area))

isotope_data <- isotope_data %>%
  filter(!(organism %in% c("mayfly", "mussel") & is.na(area)))

shiners_liver <- calculate_prop_tp(isotope_data, shiner_species, "liver") %>%
  distinct() %>%
  mutate(location = factor(location, levels = c("creek", "lake")),
         month = factor(month, levels = month_levels, labels = month_labels),
         habitat_coupling = (0.5 - abs(lake_carbon - 0.5)) / 2,
         lake_carbon_adj = adjust_beta_response(lake_carbon))

shiners_liver %>%
  mutate(sample_prefix = toupper(sub("-.*$", "", trimws(sample_ID)))) %>%
  count(sample_prefix, area, location, name = "n_records") %>%
  arrange(area, sample_prefix, location) 

shiners_liver %>%
  select(sample_ID, organism, month, location, area) %>%
  arrange(area, sample_ID) %>%
  as_tibble() %>%
  print(n = Inf)

#Check number of sampling records for the assigned areas
shiners_liver %>%
  count(area, location, name = "n_records") %>%
  arrange(area, location)

##### Isotope Models #####
##### Baseline C13, N15 #####
mayfly_stats <- isotope_data %>%
  filter(organism == "mayfly") %>%
  group_by(location, month) %>%
  summarise(n = n(),
            mean_d13C = round(mean(d13C, na.rm = TRUE), 2),
            sd_d13C = round(sd(d13C, na.rm = TRUE), 2),
            mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
            sd_d15N = round(sd(d15N, na.rm = TRUE), 2),
            .groups = "drop")

mayfly_stats

mussel_stats <- isotope_data %>%
  filter(organism == "mussel") %>%
  group_by(location, month) %>%
  summarise(n = n(),
            mean_d13C = round(mean(d13C, na.rm = TRUE), 2),
            sd_d13C = round(sd(d13C, na.rm = TRUE), 2),
            mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
            sd_d15N = round(sd(d15N, na.rm = TRUE), 2),
            .groups = "drop")

mussel_stats

baseline_data <- isotope_data %>%
  filter(organism %in% c("mayfly", "mussel")) %>%
  mutate(location = factor(location, levels = c("creek", "lake")),
         month = factor(month, levels = month_levels, labels = month_labels))

count(baseline_data, organism, area, location, month)

baseline_mf <- filter(baseline_data, organism == "mayfly")
baseline_muss <- filter(baseline_data, organism == "mussel")

#Glm for mayfly carbon 13 with interaction
glm_d13C_mf <- glmmTMB(d13C ~ month*location + (1 | area),
                       data = baseline_mf, family = gaussian())

Anova(glm_d13C_mf, type = 3)
#close to significant...

#Glm for mayfly carbon 13 without interaction
glm_d13C_mf <- glmmTMB(d13C ~ month + location + (1|area),
                       data = baseline_mf, family = gaussian())

Anova(glm_d13C_mf, type = 2)
summary(glm_d13C_mf)
performance::r2_nakagawa(glm_d13C_mf, tolerance = 1e-12)


d13C_contrasts_mf <- as.data.frame(pairs(emmeans(glm_d13C_mf, ~ month + location, type = "response"), adjust = "tukey"))
d13C_contrasts_mf

#glm for mayfly nitrogen 15 with interaction
glm_d15N_mf <- glmmTMB(d15N ~ month * location + (1|area),
                       data = baseline_mf, family = gaussian())

Anova(glm_d15N_mf, type = 3)
performance::r2_nakagawa(glm_d15N_mf, tolerance = 1e-12)


d15N_contrasts_mf <- as.data.frame(pairs(emmeans(glm_d15N_mf, ~ month * location, type = "response"), adjust = "tukey"))
d15N_contrasts_mf

#glm for mussel carbon 13 with interaction
glm_d13C_muss <- glmmTMB(d13C ~ month*location + (1 | area),
                         data = baseline_muss, family = gaussian())

Anova(glm_d13C_muss, type = 3)

#glm for mussel carbon 13 without interaction
glm_d13C_muss <- glmmTMB(d13C ~ month + location + (1 | area),
                         data = baseline_muss, family = gaussian())

Anova(glm_d13C_muss, type = 2)
performance::r2_nakagawa(glm_d13C_muss, tolerance = 1e-12)

d13C_contrasts_muss <- as.data.frame(pairs(emmeans(glm_d13C_muss, ~ month * location), adjust = "tukey"))
d13C_contrasts_muss

#glm for mussel nitrogen 15 with interaction
glm_d15N_muss <- glmmTMB(d15N ~ month*location + (1 | area),
                         data = baseline_muss, family = gaussian())

Anova(glm_d15N_muss, type = 3)

#glm for mussel nitrogen 15 without interaction
glm_d15N_muss <- glmmTMB(d15N ~ month + location + (1 | area),
                         data = baseline_muss, family = gaussian())

Anova(glm_d15N_muss, type = 2)
performance::r2_nakagawa(glm_d15N_muss, tolerance = 1e-12)

d15N_contrasts_muss <- as.data.frame(pairs(emmeans(glm_d15N_muss, ~ month * location), adjust = "tukey"))

#Check model outputs and evaluate residual diagnostics
summary(glm_d13C_mf)
d13C_contrasts_mf
sim_d13C_mf <- simulateResiduals(glm_d13C_mf, n = 1000, seed = 123)
plot(sim_d13C_mf)

summary(glm_d15N_mf)
d15N_contrasts_mf
sim_d15N_mf <- simulateResiduals(glm_d15N_mf, n = 1000, seed = 123)
plot(sim_d15N_mf)

summary(glm_d13C_muss)
d13C_contrasts_muss
sim_d13C_muss <- simulateResiduals(glm_d13C_muss, n = 1000, seed = 123)
plot(sim_d13C_muss)

summary(glm_d15N_muss)
print(d15N_contrasts_muss)
sim_d15N_muss <- simulateResiduals(glm_d15N_muss, n = 1000, seed = 123)
plot(sim_d15N_muss)

##### Shiner Isotope Models #####
shiner_iso_stats<- shiners_liver %>%
  group_by(organism, location, month) %>%
  summarise(
    n         = n(),
    mean_d13C = round(mean(d13C, na.rm = TRUE), 2),
    sd_d13C   = round(sd(d13C, na.rm = TRUE), 2),
    mean_d15N = round(mean(d15N, na.rm = TRUE), 2),
    sd_d15N   = round(sd(d15N, na.rm = TRUE), 2),
    .groups = "drop"
  )

common_shiners <- shiners_liver %>% filter(organism == "common shiner")
golden_shiners <- shiners_liver %>% filter(organism == "golden shiner")

#glm carbon 13 common shiners
glm_d13C_common <- glmmTMB(d13C ~ month * location + (1 | area),
                           data = common_shiners, family = gaussian())

Anova(glm_d13C_common, type = 3)
performance::r2_nakagawa(glm_d13C_common, tolerance = 1e-12)

d13C_contrasts_common <- as.data.frame(pairs(emmeans(glm_d13C_common, ~ month * location), adjust = "tukey"))
d13C_contrasts_common
summary(glm_d13C_common)

sim_d13C_common <- simulateResiduals(glm_d13C_common, n = 1000, seed = 123)
plot(sim_d13C_common)

#glm nitrogen 15 common shiners
glm_d15N_common <- glmmTMB(d15N ~ month * location + (1 | area),
                           data = common_shiners, family = gaussian())

Anova(glm_d15N_common, type = 3)
performance::r2_nakagawa(glm_d15N_common, tolerance = 1e-12)

d15N_contrasts_common <- as.data.frame(pairs(emmeans(glm_d15N_common, ~ month * location), adjust = "tukey"))
d15N_contrasts_common
summary(glm_d15N_common)

sim_d15N_common <- simulateResiduals(glm_d15N_common, n = 1000, seed = 123)
plot(sim_d15N_common)

#glm carbon 13 golden shiners
glm_d13C_golden <- glmmTMB(d13C ~ month * location + (1 | area),
                           data = golden_shiners, family = gaussian())

Anova(glm_d13C_golden, type = 3)
performance::r2_nakagawa(glm_d13C_golden, tolerance = 1e-12)


d13C_contrasts_golden <- as.data.frame(pairs(emmeans(glm_d13C_golden, ~ month * location), adjust = "tukey"))
d13C_contrasts_golden
summary(glm_d13C_golden)

sim_d13C_golden <- simulateResiduals(glm_d13C_golden, n = 1000, seed = 123)
plot(sim_d13C_golden)

#glm nitrogen 15 golden shiners - interaction not significant
glm_d15N_golden <- glmmTMB(d15N ~ month * location + (1 | area),
                           data = golden_shiners, family = gaussian())
Anova(glm_d15N_golden, type = 3)

glm_d15N_golden <- glmmTMB(d15N ~ month + location + (1 | area),
                           data = golden_shiners, family = gaussian())

Anova(glm_d15N_golden, type = 2)
performance::r2_nakagawa(glm_d15N_golden, tolerance = 1e-12)

d15N_contrasts_golden <- as.data.frame(pairs(emmeans(glm_d15N_golden, ~ month * location), adjust = "tukey"))
d15N_contrasts_golden
summary(glm_d15N_golden)

sim_d15N_golden <- simulateResiduals(glm_d15N_golden, n = 1000, seed = 123)
plot(sim_d15N_golden)

##### Shiner Prop Lake Carbon Model #####
#Shiners pooled together...
#use glmTMB for lake_carbon as proportion (0, 1)
glmb_carbon_shnr <- glmmTMB(lake_carbon ~ month * location + (1| area),
                            data = shiners_liver, family = beta_family(link = "logit"))

car::Anova(glmb_carbon_shnr, type = 3)
performance::r2_nakagawa(glmb_carbon_shnr, tolerance = 1e-12)

emmeans(glmb_carbon_shnr, ~ month * location, type = "response")
pairs(emmeans(glmb_carbon_shnr, ~ month | location, type = "response"), adjust = "tukey")
pairs(emmeans(glmb_carbon_shnr, ~ location | month, type = "response"), adjust = "tukey")

prop_carbon_contrasts <- pairs(emmeans(glmb_carbon_shnr, ~ month * location, type = "response")) %>%
  as.data.frame()

prop_carbon_contrasts

sim_carb <- simulateResiduals(glmb_carbon_shnr)
plot(sim_carb)

testOutliers(sim_carb, type = "bootstrap", nBoot = 1000)

##### Shiners Trophic Position Model #####
#Shiners pooled together...
#Interaction not significant...
glmb_tp_shnr <- glmmTMB(trophic_position ~ month * location + (1|area),
                        data = shiners_liver, family = gaussian())

car::Anova(glmb_tp_shnr, type = 3)

glmb_tp_shnr <- glmmTMB(trophic_position ~ month + location + (1|area),
                        data = shiners_liver, family = gaussian())

car::Anova(glmb_tp_shnr, type = 2)
performance::r2_nakagawa(glmb_tp_shnr, tolerance = 1e-15)

tp_contrasts_season <- pairs(emmeans(glmb_tp_shnr, ~ month, type = "response")) %>%
  as.data.frame()
tp_contrasts_season

tp_contrasts_habitat <- pairs(emmeans(glmb_tp_shnr, ~ location, type = "response")) %>%
  as.data.frame()
tp_contrasts_habitat

sim_tp <- simulateResiduals(glmb_tp_shnr)
plot(sim_tp) 

##### Shiner Summary Stats, Delta Calculations #####
mixing_summary <- shiners_liver %>%
  group_by(month, location) %>%
  summarise(
    n       = n(),
    mean_lc = round(mean(lake_carbon,      na.rm = TRUE), 3),
    se_lc   = round(standard_error(lake_carbon), 3),
    sd_lc   = round(sd(lake_carbon,        na.rm = TRUE), 3),
    mean_tp = round(mean(trophic_position, na.rm = TRUE), 3),
    se_tp   = round(standard_error(trophic_position), 3),
    sd_tp   = round(sd(trophic_position,   na.rm = TRUE), 3),
    .groups = "drop"
  )

mixing_summary


for (loc in c("creek", "lake")) {
  d <- mixing_summary %>% filter(location == loc)
  spr <- d %>% filter(month == "Spring") %>% pull(mean_lc)
  sumr <- d %>% filter(month == "Summer") %>% pull(mean_lc)
  fal <- d %>% filter(month == "Fall")   %>% pull(mean_lc)
  cat("\n--- Delta: season-to-season changes within location (lake carbon) ---\n")
  cat(sprintf("\n%s:\n", toupper(loc)))
  cat(sprintf("  Spring -> Summer: Δ = %+.3f\n", sumr - spr))
  cat(sprintf("  Summer -> Fall:   Δ = %+.3f\n", fal - sumr))
  cat(sprintf("  Spring -> Fall:   Δ = %+.3f\n", fal - spr))
}

for (seas in c("Spring", "Summer", "Fall")) {
  lk <- mixing_summary %>% filter(month == seas, location == "lake")  %>% pull(mean_lc)
  cr <- mixing_summary %>% filter(month == seas, location == "creek") %>% pull(mean_lc)
  cat("\n--- Delta: lake vs creek per season (lake carbon) ---\n")
  cat(sprintf("  %s - lake minus creek: Δ = %+.3f\n", seas, lk - cr))
}

for (loc in c("creek", "lake")) {
  d <- mixing_summary %>% filter(location == loc)
  spr <- d %>% filter(month == "Spring") %>% pull(mean_tp)
  sumr <- d %>% filter(month == "Summer") %>% pull(mean_tp)
  fal <- d %>% filter(month == "Fall")   %>% pull(mean_tp)
  cat("\n--- Delta: season-to-season changes within location (trophic position) ---\n")
  cat(sprintf("\n%s:\n", toupper(loc)))
  cat(sprintf("  Spring -> Summer: Δ = %+.3f\n", sumr - spr))
  cat(sprintf("  Summer -> Fall:   Δ = %+.3f\n", fal - sumr))
  cat(sprintf("  Spring -> Fall:   Δ = %+.3f\n", fal - spr))
}

for (seas in c("Spring", "Summer", "Fall")) {
  lk <- mixing_summary %>% filter(month == seas, location == "lake")  %>% pull(mean_tp)
  cr <- mixing_summary %>% filter(month == seas, location == "creek") %>% pull(mean_tp)
  cat("\n--- Delta: lake vs creek per season (trophic position) ---\n")
  cat(sprintf("  %s - lake minus creek: Δ = %+.3f\n", seas, lk - cr))
}


##### Telemetry #####
# load data
dets <- read_csv(
  file.path(file_path, "Data/dets_clean.csv"),
  show_col_types = FALSE
)

#Remove fish assumed to have experienced post-release mortality
transmitters_to_remove <- c("34905", "32999")

#Remove receivers in pre-defined transition habitat
recs_to_remove <- c("R018", "R017", "R016")

dets <- dets %>%
  filter(
    !transmitter_id %in% transmitters_to_remove,
    !rec_ID %in% recs_to_remove
  ) %>%
  mutate(
    tl = ensure_numeric(tl)
  )

dets <- dets %>%
  mutate(
    step_time = case_when(
      step_number == 1 ~ step1_dur,
      step_number == 2 ~ step2_dur,
      step_number == 3 ~ step3_dur
    ),
    time_delay = case_when(
      step_number == 1 ~ 0,
      step_number == 2 ~ step1_dur,
      step_number == 3 ~ step1_dur + step2_dur
    ),
    release_date = as.Date(release_date),
    tag_on_date = as.POSIXct(
      release_date + days(time_delay), tz = "UTC"
    ),
    tag_off_date = as.POSIXct(
      release_date + days(time_delay) + days(step_time),
      tz = "UTC"
    ),
    detection_date = as.Date(detection_timestamp_utc),
    detection_timestamp = as.POSIXct(
      detection_timestamp_utc, tz = "UTC"
    )
  )

#Important to note that:
#transmitter_id identifies the fish and
#id_time identifies a deployment step

metadata_records <- dets %>%
  transmute(
    id_time = ifelse(
      id_time == "1576142_3", "1576142_2", id_time),
    transmitter_id = na_if(
      trimws(as.character(transmitter_id)), ""),
    release_site = na_if(
      trimws(as.character(release_site)), ""),
    recorded_release_location =
      tolower(trimws(as.character(release_location))),
    tl = ensure_numeric(tl)) %>%
  distinct()

metadata_conflicts <- metadata_records %>%
  group_by(id_time) %>%
  summarise(n_fish = n_distinct(transmitter_id, na.rm = TRUE),
            n_sites = n_distinct(release_site, na.rm = TRUE),
            n_release_habitats = n_distinct(recorded_release_location, na.rm = TRUE),
            n_lengths = n_distinct(tl, na.rm = TRUE),
            .groups = "drop") %>%
  filter(n_fish > 1 | n_sites > 1 | n_release_habitats > 1 |n_lengths > 1)

first_known <- function(x) {
  if (all(is.na(x))) return(x[NA_integer_][1])
  x[which(!is.na(x))[1]]
}

fish_metadata_lookup <- metadata_records %>%
  group_by(id_time) %>%
  summarise(across(c(transmitter_id,
                     release_site,
                     recorded_release_location,
                     tl),
                   first_known),
            .groups = "drop")

detections_filtered <- dets %>%
  false_detections(tf = 3600, show_plot = FALSE) %>%
  filter(passed_filter == 1)

window_check <- detections_filtered %>%
  group_by(id_time) %>%
  summarise(first_det = min(detection_timestamp, na.rm = TRUE),
            last_det = max(detection_timestamp, na.rm = TRUE),
            tag_on = first(tag_on_date),
            tag_off = first(tag_off_date),
            step_dur = first(step_time),
            .groups = "drop")

invalid_id_times <- window_check %>%
  filter(is.na(tag_on) | is.na(tag_off) | tag_off <= tag_on | last_det < tag_on | first_det > tag_off)

if (nrow(invalid_id_times) > 0) {
  message(
    "\n--- Dropping ",
    nrow(invalid_id_times),
    " id_time(s) with invalid/mismatched windows ---")
  
  invalid_id_times %>%
    dplyr::select(
      id_time, tag_on, tag_off,
      first_det, last_det, step_dur
    ) %>%
    print(n = Inf)
}

detections_filtered <- detections_filtered %>%
  filter(id_time %in% setdiff(window_check$id_time, invalid_id_times$id_time))

id_times <- unique(detections_filtered$id_time)

time_series_list <- purrr::map(id_times, ~ create_time_series(detections_filtered, .x))

names(time_series_list) <- id_times

time_series_list <- time_series_list[!sapply(time_series_list, is.null)]

time_series_all <- bind_rows(time_series_list) %>%
  mutate(Month_Year = format(date_time, "%Y-%m"),
         id_time = ifelse(id_time == "1576142_3", "1576142_2", id_time)) %>%
  filter(!(id_time == "1576142_2" & Month_Year < "2024-04")) %>%
  group_by(id_time) %>%
  filter(date_time >= min(date_time, na.rm = TRUE)) %>%
  ungroup()

#Assign a season and a calendar year to all of the observations
time_series_all <- time_series_all %>%
  mutate(
    Season = case_when(
      (month(date_time) == 3 & day(date_time) >= 21) |
        month(date_time) %in% 4:5 |
        (month(date_time) == 6 & day(date_time) <= 20) ~ "Spring",
      
      (month(date_time) == 6 & day(date_time) >= 21) |
        month(date_time) %in% 7:8 |
        (month(date_time) == 9 & day(date_time) <= 20) ~ "Summer",
      
      (month(date_time) == 9 & day(date_time) >= 21) |
        month(date_time) %in% 10:11 |
        (month(date_time) == 12 & day(date_time) <= 20) ~ "Fall",
      
      (month(date_time) == 12 & day(date_time) >= 21) |
        month(date_time) %in% 1:2 |
        (month(date_time) == 3 & day(date_time) <= 20) ~ "Winter",
      
      TRUE ~ NA_character_),
    year = lubridate::year(date_time))

time_series_all %>%
  count(year, Season)

##### Telemetry Data Cleaning and Calculating Residency & Recurrence #####
reference_habitat_lookup <- time_series_all %>%
  group_by(id_time) %>%
  summarise(target_habitat = first(location[!is.na(location)],
                                   default = NA_character_),
            .groups = "drop")

#Summarize Residency separately for each year and season...
residency_seasonal <- time_series_all %>%
  left_join(reference_habitat_lookup,by = "id_time") %>%
  group_by(id_time, year, Season, target_habitat) %>%
  summarise(total_minutes = sum(!is.na(location)),
            habitat_minutes = sum(location == target_habitat, na.rm = TRUE),
            habitat_residency = ifelse(total_minutes > 0, 
                                       habitat_minutes / total_minutes,
                                       NA_real_),
            .groups = "drop") %>%
  rename(reference_habitat = target_habitat)

#Identify return events across each individual fish's full time series
#Recurrence events are then assigned to the year and season in which they occur

recurrence_seasonal_list <- list()

for (fish_id in unique(time_series_all$id_time)) {
  
  fish_data <- time_series_all %>%
    filter(id_time == fish_id) %>%
    arrange(date_time)
  
  fish_data$return_event <- FALSE
  fish_data$recurrence_interval <- NA_real_
  
  target_habitat <- fish_data %>%
    filter(!is.na(location)) %>%
    dplyr::slice(1) %>%
    pull(location)
  
  if (
    length(target_habitat) == 0 ||
    is.na(target_habitat)
  ) {
    recurrence_seasonal_list[[fish_id]] <- fish_data
    next
  }
  
  exited_target <- FALSE
  
  last_exit_time <- as.POSIXct(
    NA,
    origin = "1970-01-01",
    tz = "UTC"
  )
  
  for (i in seq_len(nrow(fish_data))) {
    
    this_loc <- fish_data$location[i]
    this_time <- fish_data$date_time[i]
    
    if (is.na(this_loc)) next
    
    if (
      this_loc != target_habitat &&
      !exited_target
    ) {
      exited_target <- TRUE
      last_exit_time <- this_time
    }
    
    if (
      this_loc == target_habitat &&
      exited_target
    ) {
      
      interval <- as.numeric(
        difftime(this_time, last_exit_time, units = "mins"))
      
      fish_data$return_event[i] <- TRUE
      
      fish_data$recurrence_interval[i] <- ifelse(
        interval > 0,
        interval,
        NA_real_)
      
      exited_target <- FALSE
      
      last_exit_time <- as.POSIXct(
        NA,
        origin = "1970-01-01",
        tz = "UTC")
    }
  }
  
  recurrence_seasonal_list[[fish_id]] <- fish_data
}

#Summarize return events by individual fish, year, and season

recurrence_seasonal <- bind_rows(
  recurrence_seasonal_list) %>% 
  left_join(reference_habitat_lookup, by = "id_time") %>%
  group_by(id_time, year, Season, target_habitat) %>%
  summarise(total_returns = sum(return_event, na.rm = TRUE),
            mean_recurrence_interval = mean(recurrence_interval, na.rm = TRUE), .groups = "drop") %>%
  rename(reference_habitat = target_habitat)


#Put together final data set for anaylsis
final_results <- recurrence_seasonal %>%
  left_join(residency_seasonal,by = c("id_time", "year", "Season", "reference_habitat")) %>%
  left_join(fish_metadata_lookup, by = "id_time") %>%
  distinct(id_time, year, Season,
           .keep_all = TRUE) %>%
  filter(!is.na(Season), Season != "Winter", !is.na(habitat_minutes),
         !is.na(total_minutes), total_minutes > 0, is.finite(habitat_residency)) %>%
  mutate(year = factor(year),
         tl = ensure_numeric(tl),
         release_site = factor(release_site),
         transmitter_id = factor(transmitter_id),
         Season = factor(
           Season,
           levels = season_levels),
         reference_habitat = factor(reference_habitat, levels = c("creek", "lake")),
         season_loc = paste(Season, reference_habitat, sep = "_"),
         habitat_residency_adj =
           (habitat_residency * (total_minutes - 1) + 0.5) / total_minutes,
         release_area = assign_area(release_site))

count(final_results, release_site, release_area, sort = TRUE)

final_results %>%
  count(year, Season, reference_habitat)

count(final_results, release_site, release_area, sort = TRUE)

telemetry_summary <- final_results %>%
  group_by(Season, reference_habitat) %>%
  summarise(n = n(),
            mean_habitat_residency = round(mean(habitat_residency, na.rm = TRUE), 3),
            se_habitat_residency = round(standard_error(habitat_residency), 3),
            mean_recurrence = round(mean(total_returns, na.rm = TRUE), 1),
            se_recurrence = round(standard_error(total_returns), 1),
            .groups = "drop")

print(telemetry_summary)

names(time_series_all)
names(recurrence_seasonal)
names(residency_seasonal)

##### Telemetry Residency & Occurrence Models #####
glmb_residency <- glmmTMB(
  habitat_residency_adj ~ Season * reference_habitat + (1 | release_area) + (1 | transmitter_id),
  data = final_results, family = beta_family(link = "logit"))

Anova(glmb_residency, type = 3)

glmb_residency <- glmmTMB(
  habitat_residency_adj ~ Season + reference_habitat + (1 | release_area) + (1 | transmitter_id),
  data = final_results, family = beta_family(link = "logit"))

Anova(glmb_residency, type = 2)
performance::r2_nakagawa(glmb_residency, tolerance = 1e-15)


sim_residency <- simulateResiduals(glmb_residency, n = 1000, seed = 123)
plot(sim_residency)

#Observation time is measured in days with an assigned habitat (1,440 minutes).
final_results <- final_results %>%
  mutate(observation_days = total_minutes / 1440)

print(final_results %>% group_by(Season, reference_habitat) %>%
        summarise(n = n(), min_days = min(observation_days),
                  median_days = median(observation_days),
                  max_days = max(observation_days), .groups = "drop"))

glmb_recurrence <- glmmTMB(total_returns ~ Season * reference_habitat 
                           + offset(log(observation_days)) + (1 | release_area) + (1 | transmitter_id),
                           data = final_results, family = nbinom2(link = "log"))

Anova(glmb_recurrence, type = 3)

glmb_recurrence <- glmmTMB(total_returns ~ Season + reference_habitat 
                           + offset(log(observation_days)) + (1 | release_area) + (1 | transmitter_id),
                           data = final_results, family = nbinom2(link = "log"))

Anova(glmb_recurrence, type = 2)
performance::r2_nakagawa(glmb_recurrence, tolerance = 1e-12)

sim_recurrence <- simulateResiduals(glmb_recurrence, n = 1000, seed = 123)
plot(sim_recurrence)

#Telemetry comparisons: additive models, averaged over the other fixed factor.
#Tests use the link scale; reported effects are odds ratios (residency) or rate ratios (recurrence)
#P-values are Tukey-adjusted within each set.
residency_season_contrasts <- as.data.frame(summary(pairs(emmeans(glmb_residency, ~ Season), adjust = "tukey"), type = "response", infer = c(FALSE, TRUE), adjust = "tukey"))

print(residency_season_contrasts)

residency_habitat_contrasts <- as.data.frame(summary(pairs(emmeans(glmb_residency, ~ reference_habitat), adjust = "tukey"), type = "response", infer = c(FALSE, TRUE), adjust = "tukey"))

print(residency_habitat_contrasts)

recurrence_season_contrasts <- as.data.frame(summary(
  pairs(emmeans(glmb_recurrence, ~ Season, offset = 0), adjust = "tukey"),
  type = "response", infer = c(FALSE, TRUE), adjust = "tukey"
))
print(recurrence_season_contrasts)

recurrence_habitat_contrasts <- as.data.frame(summary(
  pairs(emmeans(glmb_recurrence, ~ reference_habitat, offset = 0), adjust = "tukey"),
  type = "response", infer = c(FALSE, TRUE), adjust = "tukey"
))

print(recurrence_habitat_contrasts)

##### Estimate Marginal Means from Isotope and Telemetry Top Models #####
emm_carbon <- emmeans(glmb_carbon_shnr, ~ month * location, type = "response") %>%
  as.data.frame() %>%
  rename(est = response) %>%
  mutate(lo = pmax(0, est - SE), hi = pmin(1, est + SE))

emm_tp <- emmeans(glmb_tp_shnr, ~ month + location, type = "response") %>%
  as.data.frame() %>%
  rename(est = emmean) %>%
  mutate(lo = est - SE, hi = est + SE)

emm_residency <- emmeans(glmb_residency, ~ Season + reference_habitat, type = "response") %>%
  as.data.frame() %>%
  rename(est = response) %>%
  mutate(lo = pmax(0, est - SE), hi = pmin(1, est + SE))

emm_recurrence <- emmeans(
  glmb_recurrence, ~ Season + reference_habitat, type = "response", offset = 0) %>%
  as.data.frame() %>%
  rename(est = response) %>%
  mutate(lo = est - SE, hi = est + SE)


##### Figure 3 #####
#3A - Lake Derived Carbon
plot_shiner_carbon <- ggplot(
  emm_carbon,
  aes(x = month, y = est, group = location)
) +
  geom_line(alpha = 0.6) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7,
    # position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = location),
    size = 2.8, shape = 21, color = "#323332", stroke = 0.6,
    # position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name = "Location",
    values = fill_scheme,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  labs(y = "Proportion Lake-derived Carbon", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0, 1)) +
  custom_theme +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
  ))

#3B - Trophic position
plot_shiner_tp <- ggplot(
  emm_tp,
  aes(x = month, y = est, group = location)
) +
  geom_line(alpha = 0.6) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7,
    # position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = location),
    size = 2.8, shape = 21, color = "#323332", stroke = 0.6,
    # position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name = "Location",
    values = fill_scheme,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  labs(y = "Trophic Position", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  custom_theme +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
  ))

#3C - Recurrence
plot_recurrence <- ggplot(
  emm_recurrence,
  aes(x = Season, y = est, group = reference_habitat)
) +
  geom_line(alpha = 0.6) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7,
    # position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = reference_habitat),
    size = 2.8, shape = 21, color = "#323332", stroke = 0.6,
    # position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name = "Location",
    values = fill_scheme,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(y = "Habitat Recurrence (returns/day)", x = NULL) +
  custom_theme +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
  ))

#3D - Habitat residency
plot_residency <- ggplot(
  emm_residency,
  aes(x = Season, y = est, group = reference_habitat)
) +
  geom_line(alpha = 0.6) +
  geom_errorbar(
    aes(ymin = lo, ymax = hi),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7,
    # position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = reference_habitat),
    size = 2.8, shape = 21, color = "#323332", stroke = 0.6,
    # position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name = "Location",
    values = fill_scheme,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  labs(y = "Habitat Residency", x = NULL) +
  scale_y_continuous(breaks = scales::pretty_breaks(), limits = c(0, 1)) +
  custom_theme +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 3, color = "#323332", stroke = 0.6)
  ))


left_iso <- ggarrange(
  plot_shiner_carbon + theme(legend.position = "none"),
  plot_shiner_tp     + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v"
)

right_tel <- ggarrange(
  plot_recurrence + theme(legend.position = "none"),
  plot_residency  + theme(legend.position = "none"),
  ncol = 1, nrow = 2, align = "v"
)

legend_shared <- get_legend(
  plot_shiner_carbon +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 12)
    )
)

top_row <- ggarrange(left_iso, right_tel, ncol = 2)

figure_3 <- ggarrange(top_row, legend_shared, ncol = 1, heights = c(10, 1))
figure_3

ggsave(file.path(file_path, "Figures/raw_r_output_figs/figure_3.jpg"), plot = figure_3, width = 8, height = 7, dpi = 900)


##### Supplementary Analysis & Sensitivity Checks for Species & Year #####
#Ensure year is treated as a categorical variable
shiners_liver <- shiners_liver %>%
  mutate(year = factor(year))

#Check sample sizes across years, species, seasons and habitats
shiners_liver %>% 
  count(year, organism, month, location)


#Models for each response variable with organism and year added
#Lake carbon glm
carbon_species <- glmmTMB(lake_carbon ~ month * location + organism + year + (1 | area),
                          data = shiners_liver, family = beta_family(link = "logit"))

Anova(carbon_species, type = 3)
performance::r2_nakagawa(carbon_species, tolerance = 1e-12)

#Year and Organisms not significant (p > 0.05)
summary(carbon_species)

#TP glm
tp_species <- glmmTMB(trophic_position ~ month + location + organism + year + (1 | area),
                      data = shiners_liver, family = gaussian())

Anova(tp_species, type = 2)
performance::r2_nakagawa(tp_species, tolerance = 1e-15)

#Year and organism not significant (p > 0.05)
summary(tp_species)

#Residency glm (only testing year - all individuals are golden shiners...)
#Ensure year is treated as a categorical variable
final_results <- final_results %>%
  mutate(year = factor(year)) 

glmb_residency_year <- glmmTMB(habitat_residency_adj ~ Season + reference_habitat + year + (1 | release_area) + (1 | transmitter_id),
                               data = final_results, family = beta_family(link = "logit"))

Anova(glmb_residency_year, type = 2)
performance::r2_nakagawa(glmb_residency_year, tolerance = 1e-20)

#Year not significant (p > 0.05)
summary(glmb_residency_year)

#Recurrence glm
glmb_recurrence_year <- glmmTMB(total_returns ~ Season + reference_habitat + year + 
                                  offset(log(observation_days)) + (1 | release_area) + (1 | transmitter_id),
                                data = final_results, family = nbinom2(link = "log"))

Anova(glmb_recurrence_year, type = 2)
performance::r2_nakagawa(glmb_recurrence_year, tolerance = 1e-15)

#Year not significant (p > 0.05)
summary(glmb_residency_year)


