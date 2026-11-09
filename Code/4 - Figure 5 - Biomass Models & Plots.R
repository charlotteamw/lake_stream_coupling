#Script for Figure 5 - Minnow Trapping data and models, for manuscript "Habitat coupling by small fish across lake-stream habitats: a missing link in freshwater food webs?"

#Author(s): Charlotte Ward & Reilly O'Connor
# Version: 2026-10-25

# Load Pkgs
library(tidyverse) 
library(patchwork)
library(grid)
library(car)
library(easystats)
library(cowplot)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(emmeans)

file_path <- getwd()

source(file.path(file_path, "Code/0 - Functions.R"))

#Set theme
habitat_colors <- c("creek" = "#F06C57", "lake" = "#4A90B8")
base_theme <- theme_bw() + theme(
  panel.grid = element_blank(), panel.border = element_blank(),
  axis.line = element_line(color = "black"), legend.position = "none",
  axis.title = element_text(size = 12), axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 10))

# Load Data
mt_data <- read.csv(file.path(file_path, "Data/minnowtrap_data.csv"), header = TRUE)

#Filter for shiners and clean
shiner_mt_data <- mt_data %>%
  filter(Species %in% c("Golden Shiner", "Common Shiner")) %>%
  mutate(wt = Wt,
         tl = TL,
         Location = tolower(Location))

all_sites  <- c(paste0("CR", 1:10), paste0("LK", 1:10))
all_months <- c("May", "August", "October")

site_counts <- shiner_mt_data %>%
  mutate(Site = str_to_upper(Site)) %>%
  group_by(Site, Month, Location) %>%
  summarise(wt = if (anyNA(wt)) NA_real_ else sum(wt),
            .groups = "drop")

site_counts_complete <- site_counts %>%
  complete(Site  = all_sites, Month = all_months, fill  = list(wt = 0), explicit = FALSE) %>%
  mutate(Location = case_when(
    str_starts(Site, "CR") ~ "creek",
    str_starts(Site, "LK") ~ "lake"),
    Site = factor(Site),
    Location = factor(Location, levels = c("creek", "lake")),
    Month = factor(Month, levels = c("May", "August", "October"))) %>%
  arrange(Location, Site, Month) %>%
  #Add a single positive biomass observation for august, shiners likely just below detectability in lake, 3 represents half the minimum biomass of shiners caught in traps... 
  mutate(wt = case_when(Site == "LK6" & Month == "August" ~ 3, TRUE ~ wt))

cpue_data <- site_counts_complete %>% 
  filter(!is.na(wt)) %>% 
  droplevels()

count(cpue_data, Month, Location)

##### Shiner CPUE glm models #####
glm_cpue <- glmmTMB(wt ~ Month*Location + (1 | Site), 
                    data = cpue_data, family = tweedie(link = "log"),
                    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

Anova(glm_cpue, type = 3)
summary(glm_cpue)

#Check conditional and marginal R2 - random effect influence...
performance::r2_nakagawa(glm_cpue)

#Check residuals with Dharma
sim_cpue <- simulateResiduals(glm_cpue, n = 1000, seed = 27)
plot(sim_cpue)

#emmeans marginal means
emm_cpue <- emmeans(glm_cpue, ~ Month + Location, type = "response")
cpue_contrasts <- as.data.frame(pairs(emm_cpue, adjust = "tukey"))
print(cpue_contrasts)


# These points are model estimates, not observed catches. The additive model
# can predict positive summer lake biomass even though its observed catches are zero.
df_emm_cpue <- emmeans(glm_cpue, ~ Month:Location, type = "response") %>%
  as.data.frame() %>%
  mutate(
    Month    = factor(Month, levels = c("May", "August", "October"),
                      labels = c("Spring", "Summer", "Fall")),
    Location = factor(Location, levels = c("creek", "lake"))
  ) %>%
  mutate(log_response = log10(response),
         log_lcl = log10(asymp.LCL),
         log_ucl = log10(asymp.UCL))


plot_cpue <- ggplot(
  df_emm_cpue, aes(x = Month, y = log_response, group = Location)
) +
  geom_errorbar(
    aes(ymin = log_lcl, ymax = log_ucl),
    width = 0.15, linewidth = 0.6,
    color = "#323332", alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = Location),
    size = 3, shape = 21,
    color = "#323332", stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name = "Location",
    values = habitat_colors,
    breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50", linewidth = 0.5
  ) +
  labs(x = NULL, y = "Log Biomass per Trap (g)") +
  coord_cartesian(clip = "off") +
  base_theme

plot_cpue

##### Load Directional Trap Data - Estimate Total Directional Flux (Upstream) #####
trap_data <- read.csv(file.path(file_path, "Data/directional_traps_all.csv"))

shiner_data <- trap_data %>%
  filter(species %in% c("golden shiner", "common shiner"))

shiner_data_fac <- shiner_data %>%
  mutate(season_order = factor(season, levels = c("spring", "summer", "fall")),
         type = factor(type, levels = c("fyke", "weir")),
         direction = factor(direction, levels = c("upstream", "downstream")),
         wt = as.numeric(wt))

shiner_data_bs <- shiner_data_fac %>%
  filter(number == 1) %>%
  filter(!is.na(tl), !is.na(wt), wt > 0)

ggplot(shiner_data_bs, aes(x = tl, y = wt, color = species)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~species) +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw(base_size = 14)

common_bs <- shiner_data_bs %>%
  filter(species == 'common shiner')

golden_bs <- shiner_data_bs %>%
  filter(species == 'golden shiner')

#make length weight relationships for common and golden shiners...
lm_common_len_wt <- lm(log10(wt) ~ log10(tl), data = common_bs)
summary(lm_common_len_wt)

lm_goolden_len_wt <- lm(log10(wt) ~ log10(tl), data = golden_bs)
summary(lm_goolden_len_wt)

shiner_data_fac_na <- shiner_data_fac %>%
  filter(is.na(wt), !is.na(tl)) %>%
  mutate(wt = case_when(
    species == "common shiner" ~ 10^predict(lm_common_len_wt, newdata = .),
    species == "golden shiner" ~ 10^predict(lm_goolden_len_wt, newdata = .),
      TRUE ~ wt)) %>%
  mutate(wt = wt*as.numeric(number))

shiner_data_fac_na_fixed <- shiner_data_fac %>%
  filter(!is.na(wt)) %>%
  bind_rows(shiner_data_fac_na)

#Add groups for sampling days and events...
trap_clean <- shiner_data_fac_na_fixed %>%
  mutate(sampling_day = interaction(year, month, day, drop = TRUE),
         sampling_event = interaction(type, sampling_day, drop = TRUE))

#Add back sampling days where no shiners were recorded, but other species were - assign true zeros...
keys <- c("year", "season", "month", "day", "type", "direction")

#All represented samples, including those containing only other species
observed_samples <- trap_data %>%
  distinct(across(all_of(keys)))

#Sampling days where shiners were recorded...
shiner_samples <- shiner_data %>%
  distinct(across(all_of(keys))) %>%
  mutate(shiners_recorded = TRUE)

#Sum measured and estimated shiner weights
shiner_totals <- shiner_data_fac_na_fixed %>%
  group_by(across(all_of(keys))) %>%
  summarise(total_daily_biomass = sum(wt),
            .groups = "drop")

flux_data <- observed_samples %>%
  left_join(shiner_samples, by = keys) %>%
  left_join(shiner_totals, by = keys) %>%
  mutate(total_daily_biomass = if_else(
    is.na(shiners_recorded), 0, total_daily_biomass)) %>%
  # filter(!(year == "2024" & month == "april" &
  #            day == 24 & type == "fyke")) %>%
  group_by(year, season, month, day, direction, type) %>%
  summarise(total_daily_biomass = sum(total_daily_biomass),
            .groups = "drop") %>%
  mutate(season = factor(season, levels = c("spring", "summer", "fall")),
         direction = factor(direction, levels = c("upstream", "downstream")),
         year = factor(year),
         sampling_day = interaction(year, month, day, drop = TRUE),
         sampling_event = interaction(sampling_day, type, drop = TRUE)) %>%
  droplevels()

#Model
glm_flux <- glmmTMB(total_daily_biomass ~ season * direction + (1 | sampling_event), 
                    family = tweedie(link = "log"), data = flux_data)

Anova(glm_flux, type = 3)
summary(glm_flux)
r2_nakagawa(glm_flux)

sim_flux <- simulateResiduals(glm_flux, n = 1000, seed = 27)
plot(sim_flux)

emm_flux <- emmeans(glm_flux, ~ season * direction, type = "response")

df_emm_flux <- as.data.frame(emm_flux) %>%
  mutate(season = factor(season, levels = c("spring", "summer", "fall"),
                         labels = c("Spring", "Summer", "Fall")),
         direction = factor(direction, levels = c("upstream", "downstream")))

direction_colors <- c("upstream" = "grey75", "downstream" = "grey25")

plot_flux <- ggplot(
  df_emm_flux, aes(x = season, y = log10(response), group = direction)
) +
  geom_errorbar(
    aes(ymin = log10(asymp.LCL), ymax = log10(asymp.UCL)),
    width = 0.15, linewidth = 0.6,
    color = "#323332", alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = direction),
    size = 3, shape = 21,
    color = "#323332", stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name = "Direction",
    values = direction_colors,
    breaks = c("upstream", "downstream"),
    labels = c("Upstream", "Downstream")
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5),
    linetype = "dashed",
    color = "grey50", linewidth = 0.5
  ) +
  labs(x = NULL, y = "Log Daily Biomass (g)") +
  coord_cartesian(clip = "off") +
  base_theme

plot_flux

figure_5_combined <- ggarrange(
  plot_flux, plot_cpue,
  ncol = 1,
  heights = c(1, 1),
  align = "v"
)

figure_5_combined

ggsave(
  file.path(file_path, "Figures/raw_r_output_figs/figure_5.jpg"),
  figure_5_combined,
  width = 6, height = 5, dpi = 900
)


##### Supplementary Analysis - Check Influence of Year & Species #####
#CPUE Data
cpue_species_data <- shiner_mt_data %>%
  mutate(Site = str_to_upper(Site)) %>%
  group_by(Site, Month, Species) %>%
  summarise(wt = if (anyNA(wt)) NA_real_ else sum(wt),
            .groups = "drop") %>%
  complete(Site = all_sites, Month = all_months,
           Species = c("Common Shiner", "Golden Shiner"),
           fill = list(wt = 0), explicit = FALSE) %>%
  mutate(Location = factor(
    if_else(str_starts(Site, "CR"), "creek", "lake"), levels = c("creek", "lake")),
    Month = factor(Month, levels = all_months),
    Species = factor(Species),
    Site = factor(Site),
    sample_id = interaction(Site, Month, drop = TRUE)) %>%
  filter(!is.na(wt)) %>%
  mutate(wt = if_else(Site == "LK6" & Month == "August", 1.5, wt))

glm_cpue_sens <- glmmTMB(wt ~ Month * Location + Species + (1 | Site),
                         family = tweedie(link = "log"), data = cpue_species_data)

Anova(glm_cpue_sens, type = 3)
r2_nakagawa(glm_cpue_sens)

#Look at trap data
flux_species_samples <- crossing(observed_samples, species = c("common shiner", "golden shiner"))

species_keys <- c(keys, "species")

flux_species_recorded <- shiner_data %>%
  distinct(across(all_of(species_keys))) %>%
  mutate(species_recorded = TRUE)

flux_species_totals <- shiner_data_fac_na_fixed %>%
  group_by(across(all_of(species_keys))) %>%
  summarise(total_daily_biomass = sum(wt), .groups = "drop")

flux_species_data <- flux_species_samples %>%
  left_join(flux_species_recorded, by = species_keys) %>%
  left_join(flux_species_totals, by = species_keys) %>%
  mutate(total_daily_biomass = if_else(
    is.na(species_recorded), 0, total_daily_biomass)) %>%
  group_by(year, season, month, day, type, direction, species) %>%
  summarise(total_daily_biomass = sum(total_daily_biomass),
            .groups = "drop") %>%
  mutate(season = factor(season, levels = c("spring", "summer", "fall")),
         direction = factor(direction, levels = c("upstream", "downstream")),
         species = factor(species),
         year = factor(year),
         sampling_day = interaction(year, month, day, drop = TRUE),
         sampling_event = interaction(sampling_day, type, drop = TRUE))

glm_flux_sens <- glmmTMB(total_daily_biomass ~ season * direction 
                         + species + year + (1 | sampling_event),
                         data = flux_species_data, family = tweedie(link = "log"))

Anova(glm_flux_sens, type = 3)
r2_nakagawa(glm_flux_sens)


