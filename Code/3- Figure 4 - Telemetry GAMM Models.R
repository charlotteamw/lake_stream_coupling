# Script for Figure 4 — Generalized Additive Models for Shiner Telemetry, for manuscript "Habitat coupling by small fish across lake-stream habitats: a missing link in freshwater food webs?"

# Version: 2026-10-25
# Code Co-Authors: Charlotte Ward & Reilly O'Connor

# load Pkgs
library(glatos)
library(tidyverse)
library(lubridate)
library(mgcv)
library(gratia)
#unhash and install below if not already installed
# pak::pak(mgcvUtils)
library(mgcvUtils)
library(nlme)

file_path <- getwd()

source(file.path(file_path, "Code/0 - Functions.R"))

# load data
dets <- read_csv("Data/dets_clean.csv")

##### Filter Data, assign areas, calculate response variables #####
transmitters_to_remove <- c("34905", "32999")
recs_to_remove <- c("R018", "R017", "R016")

dets_1 <- dets %>%
  filter(!transmitter_id %in% transmitters_to_remove) %>%
  filter(!rec_ID %in% recs_to_remove)

dets_2 <- dets_1 %>%
  mutate(step_time = case_when(
    step_number == 1 ~ step1_dur,
    step_number == 2 ~ step2_dur,
    step_number == 3 ~ step3_dur),
    time_delay = case_when(
      step_number == 1 ~ 0,
      step_number == 2 ~ step1_dur,
      step_number == 3 ~ step1_dur + step2_dur),
    release_date = as.Date(release_date),
    tag_on_date = release_date + days(time_delay),
    tag_off_date = tag_on_date + days(step_time),
    detection_date = as.Date(detection_timestamp_utc, tz = "UTC"),
    detection_timestamp = as.POSIXct(detection_timestamp_utc, tz = "UTC"),
    tag_on_date = as.POSIXct(tag_on_date, tz = "UTC"),
    tag_off_date = as.POSIXct(tag_off_date, tz = "UTC"))

detections_filtered <- dets_2 %>%
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
  
  message("\nDropping ", nrow(invalid_id_times), " id_time(s) with invalid/mismatched windows.")
  
  invalid_id_times %>%
    select(id_time, tag_on, tag_off,
           first_det, last_det, step_dur) %>%
    print(n = Inf)
}

detections_filtered <- detections_filtered %>%
  filter(id_time %in% (window_check %>% 
                         filter(!id_time %in% invalid_id_times$id_time) %>%
                         pull(id_time)))%>%
  mutate(release_area = assign_area(release_site))



#Double check the site-to-area assignments...
detections_filtered %>%
  distinct(release_site, release_area) %>%
  arrange(release_area, release_site) %>%
  print(n = Inf)

#Calculate daily movement for multiple locations...
daily_movements <- detections_filtered %>%
  group_by(transmitter_id, detection_date, release_location, release_area, tl) %>%
  summarise(unique_locations = n_distinct(location, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(multiple_locations = ifelse(unique_locations > 1, 1, 0),
         date = as.Date(detection_date),
         day_of_year = yday(date))

#Check active tags and make final data frame for analysis...
#Preserve the original earliest-on/latest-off active-tag calculation and use date keys on both sides of the join
active_tag_sequence <- detections_filtered %>%
  filter(!is.na(tag_on_date) & !is.na(tag_off_date)) %>%
  group_by(transmitter_id) %>%
  summarise(tag_on = as.Date(min(tag_on_date), tz = "UTC"),
            tag_off = as.Date(max(tag_off_date), tz = "UTC"),
            .groups = "drop") %>%
  mutate(dates = map2(tag_on, tag_off, ~ seq(.x, .y, by = "day"))) %>%
  unnest(dates) %>%
  group_by(dates) %>%
  summarise(active_tags = n(), .groups = "drop") %>%
  rename(date = dates)

daily_movements_active <- daily_movements %>%
  left_join(active_tag_sequence, by = "date") %>%
  filter(active_tags > 2)

required_variables <- c("multiple_locations", "day_of_year","tl",
                        "release_location", "release_area", "transmitter_id", "date")

daily_movements_final <- daily_movements_active %>%
  filter(if_all(all_of(required_variables), ~ !is.na(.x))) %>%
  mutate(transmitter_id = factor(transmitter_id),
         release_area = factor(release_area),
         time_day = as.integer(date - min(date))) %>%
  arrange(release_area, transmitter_id, date) %>%
  droplevels()

str(daily_movements_final$detection_date)

cat("Model observations:", nrow(daily_movements_final), "\n")
cat("Fish:", nlevels(daily_movements_final$transmitter_id), "\n")
cat("Release areas:", nlevels(daily_movements_final$release_area), "\n")


##### Fit GAMMs (using bam() for testing results sensitivity to AR1 addition) and Model Selection #####
bam_data <- daily_movements_final %>%
  ungroup() %>%
  arrange(transmitter_id, date) %>%
  group_by(transmitter_id) %>%
  mutate(ar_start = row_number() == 1L | as.integer(date - lag(date)) != 1L) %>%
  ungroup()

m_bam_mod_no_re <- mgcv::bam(
  multiple_locations ~
    s(day_of_year, bs = "cc", k = 20) 
  + s(tl, bs = 'tp', k = 20) 
  # + s(release_area, bs = "re")
  # + s(transmitter_id, bs = "re")
  ,
  family = binomial(link = "logit"),
  data = bam_data,
  knots = list(day_of_year = c(0.5, 366.5)),
  method = "ML",
  discrete = F,
  na.action = na.fail)

m_bam_mod_rel_area <- mgcv::bam(
  multiple_locations ~
    s(day_of_year, bs = "cc", k = 20) 
  + s(tl, bs = 'tp', k = 20) 
  + s(release_area, bs = "re")
  # + s(transmitter_id, bs = "re")
  ,
  family = binomial(link = "logit"),
  data = bam_data,
  knots = list(day_of_year = c(0.5, 366.5)),
  method = "ML",
  discrete = F,
  na.action = na.fail)

m_bam_mod_trans_id <- mgcv::bam(
  multiple_locations ~
    s(day_of_year, bs = "cc", k = 20) 
  + s(tl, bs = 'tp', k = 20) 
  # + s(release_area, bs = "re")
  + s(transmitter_id, bs = "re")
  ,
  family = binomial(link = "logit"),
  data = bam_data,
  knots = list(day_of_year = c(0.5, 366.5)),
  method = "ML",
  discrete = F,
  na.action = na.fail)

m_bam_mod_trans_id_rel_area <- mgcv::bam(
  multiple_locations ~
    s(day_of_year, bs = "cc", k = 20) 
  + s(tl, bs = 'tp', k = 20) 
  + s(release_area, bs = "re")
  + s(transmitter_id, bs = "re")
  ,
  family = binomial(link = "logit"),
  data = bam_data,
  knots = list(day_of_year = c(0.5, 366.5)),
  method = "ML",
  discrete = F,
  na.action = na.fail)

gam.check(m_bam_mod_no_re)
gam.check(m_bam_mod_rel_area)
gam.check(m_bam_mod_trans_id)
gam.check(m_bam_mod_trans_id_rel_area)

appraise(m_bam_mod_no_re)
appraise(m_bam_mod_rel_area)
appraise(m_bam_mod_trans_id)
appraise(m_bam_mod_trans_id_rel_area)

draw(m_bam_mod_no_re, residuals = T)
draw(m_bam_mod_rel_area, residuals = T)
draw(m_bam_mod_trans_id, residuals = T)
draw(m_bam_mod_trans_id_rel_area, residuals = T)

summary(m_bam_mod_no_re)
summary(m_bam_mod_rel_area)
summary(m_bam_mod_trans_id)
summary(m_bam_mod_trans_id_rel_area)

bbmle::AICtab(m_bam_mod_no_re, m_bam_mod_rel_area, m_bam_mod_trans_id, m_bam_mod_trans_id_rel_area)

### Retain only transmitter ID random effect; release area is not improving model fit and explains little variation in the data

# full model, with linear term for total length (tl) instead of smooth term
m_bam_mod_trans_id_tl <- mgcv::bam(
  multiple_locations ~
    s(day_of_year, bs = "cc", k = 20) 
  + tl
  # + s(release_area, bs = "re")
  + s(transmitter_id, bs = "re")
  ,
  family = binomial(link = "logit"),
  data = bam_data,
  knots = list(day_of_year = c(0.5, 366.5)),
  method = "ML",
  discrete = FALSE,
  na.action = na.fail)

# no total length 
m_bam_mod_trans_id_no_tl <- mgcv::bam(
  multiple_locations ~
    s(day_of_year, bs = "cc", k = 20) 
  # + tl
  # + s(release_area, bs = "re")
  + s(transmitter_id, bs = "re")
  ,
  family = binomial(link = "logit"),
  data = bam_data,
  knots = list(day_of_year = c(0.5, 366.5)),
  method = "ML",
  discrete = FALSE,
  na.action = na.fail)

summary(m_bam_mod_trans_id)
summary(m_bam_mod_trans_id_tl)
summary(m_bam_mod_trans_id_no_tl)

bbmle::AICtab(m_bam_mod_trans_id, m_bam_mod_trans_id_tl, m_bam_mod_trans_id_no_tl)

candidate_models <- list(
  m_bam_mod_no_re,
  m_bam_mod_rel_area,
  m_bam_mod_trans_id,
  m_bam_mod_trans_id_rel_area,
  m_bam_mod_trans_id_tl,
  m_bam_mod_trans_id_no_tl)

aic <- vapply(candidate_models, AIC, numeric(1))
delta_aic <- aic - min(aic)
relative_likelihood <- exp(-delta_aic / 2)

model_comparison <- data.frame(
  Model = vapply(
    candidate_models, function(m) deparse1(formula(m)[[3]]), character(1)
  ),
  n = vapply(candidate_models, nobs, numeric(1)),
  EDF = round(vapply(candidate_models, function(m) sum(m$edf), numeric(1)), 2),
  LogLik = round(vapply(
    candidate_models, function(m) as.numeric(logLik(m)), numeric(1)
  ), 2),
  AIC = round(aic, 2),
  Delta_AIC = round(delta_aic, 2),
  Akaike_weight = signif(relative_likelihood / sum(relative_likelihood), 3),
  Deviance_explained_pct = round(vapply(
    candidate_models, function(m) 100 * summary(m)$dev.expl, numeric(1)
  ), 2)
)

model_comparison <- model_comparison[order(aic), ]
print(model_comparison, row.names = FALSE)

write.csv(model_comparison, "Figures/raw_r_output_figs/GAMM_model_comparison.csv", row.names = FALSE)

##### Test Influence of Autocorrelation #####
#refit with freml and discrete = T...
m_bam_mod_trans_id_no_tl <- mgcv::bam(
  multiple_locations ~
    s(day_of_year, bs = "cc", k = 20) 
  # + s(tl, bs = 'tp', k = 20) 
  # + s(release_area, bs = "re")
  + s(transmitter_id, bs = "re")
  ,
  family = binomial(link = "logit"),
  data = bam_data,
  knots = list(day_of_year = c(0.5, 366.5)),
  method = "fREML",
  discrete = TRUE,
  na.action = na.fail)

new_data <- tibble(
  day_of_year = 1:366,
  transmitter_id = factor(
    levels(bam_data$transmitter_id)[1],
    levels = levels(bam_data$transmitter_id)
  )
)


pred <- predict(
  m_bam_mod_trans_id_no_tl,
  newdata = new_data,
  type = "link",
  se.fit = TRUE,
  exclude = "s(transmitter_id)"
)

new_data <- new_data %>%
  mutate(
    predicted_prob = plogis(as.numeric(pred$fit)),
    lower_CI = plogis(as.numeric(pred$fit - 1.96 * pred$se.fit)),
    upper_CI = plogis(as.numeric(pred$fit + 1.96 * pred$se.fit)),
    doy_rotated = (day_of_year - 79) %% 366 + 1
  ) %>%
  arrange(doy_rotated)

axis_breaks <- seq(1, 361, by = 30)

gg_gam <- ggplot(new_data, aes(doy_rotated, predicted_prob)) +
  geom_ribbon(
    aes(ymin = lower_CI, ymax = upper_CI),
    alpha = 0.4, fill = "gray50"
  ) +
  geom_line(linewidth = 1.1, colour = "black") +
  geom_vline(
    xintercept = c(94, 186, 277),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_x_continuous(
    name = "Julian Day",
    breaks = axis_breaks,
    labels = (axis_breaks + 77) %% 366 + 1
  ) +
  labs(y = "Prob. of Lake-Stream Detection") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    plot.margin = grid::unit(c(0.3, 0.3, 1.2, 0.3), "cm"),
    axis.title.x = element_text(vjust = -1.5)
  ) +
  coord_cartesian(ylim = c(0, 0.8))

gg_gam

#Use PACF to generate estimate for rho...
pacf(resid(m_bam_mod_trans_id_no_tl)) 

#it is around 0.25, so account for auto-correlation with rho = 0.25
m_bam_mod_trans_id_no_tl_ar1 <- mgcv::bam(
  multiple_locations ~
    s(day_of_year, bs = "cc", k = 20) 
  # + s(tl, bs = 'tp', k = 5) 
  # + s(release_area, bs = "re")
  + s(transmitter_id, bs = "re")
  ,
  family = binomial(link = "logit"),
  data = bam_data,
  rho = 0.25,
  AR.start = bam_data$ar_start,
  knots = list(day_of_year = c(0.5, 366.5)),
  method = "fREML",
  discrete = TRUE,
  na.action = na.fail)

#See how model predictions change...

pred <- predict(
  m_bam_mod_trans_id_no_tl_ar1,
  newdata = new_data,
  type = "link",
  se.fit = TRUE,
  exclude = "s(transmitter_id)"
)

new_data <- new_data %>%
  mutate(
    predicted_prob = plogis(as.numeric(pred$fit)),
    lower_CI = plogis(as.numeric(pred$fit - 1.96 * pred$se.fit)),
    upper_CI = plogis(as.numeric(pred$fit + 1.96 * pred$se.fit)),
    doy_rotated = (day_of_year - 79) %% 366 + 1
  ) %>%
  arrange(doy_rotated)

axis_breaks <- seq(1, 361, by = 30)

gg_gam_ar1 <- ggplot(new_data, aes(doy_rotated, predicted_prob)) +
  geom_ribbon(
    aes(ymin = lower_CI, ymax = upper_CI),
    alpha = 0.4, fill = "gray50"
  ) +
  geom_line(linewidth = 1.1, colour = "black") +
  geom_vline(
    xintercept = c(94, 186, 277),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_x_continuous(
    name = "Julian Day",
    breaks = axis_breaks,
    labels = (axis_breaks + 77) %% 366 + 1
  ) +
  labs(y = "Prob. of Lake-Stream Detection") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    plot.margin = grid::unit(c(0.3, 0.3, 1.2, 0.3), "cm"),
    axis.title.x = element_text(vjust = -1.5)
  ) +
  coord_cartesian(ylim = c(0, 0.8))

gg_gam
gg_gam_ar1

summary(m_bam_mod_trans_id_no_tl_ar1)


ggsave(file.path(file_path, "Figures/raw_r_output_figs/figure_4.png"), gg_gam_ar1,
       width = 8.5, height = 4.5, dpi = 600, bg = "transparent")



