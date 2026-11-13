# Script for Figure 6 — Seasonal Comparisons of Fish Body Size for manuscript "Habitat coupling by small fish across lake-stream habitats: a missing link in freshwater food webs?"

# Version: 2026-10-25
#Author(s): Charlotte Ward & Reilly O'Connor

# load Pkgs
library(tidyverse)
library(glmmTMB)
library(car)
library(DHARMa)
library(emmeans)
library(ggpubr)
library(performance)


file_path <- getwd()

source(file.path(file_path, "Code/0 - Functions.R"))

# load Data
fish_file <- file.path(file_path, "Data/fish_catches_all_merged.csv")

df_tl <- read.csv(fish_file, header = TRUE, stringsAsFactors = FALSE) 

df_tl_num <- df_tl %>%
  mutate(
  across(c(fl, tl, wt), ~ suppressWarnings(as.numeric(trimws(.x))))) %>%
  filter(!is.na(tl), tl > 0) %>%
  mutate(number = (as.numeric(as.character(number))),
         species = str_to_lower(trimws(species)),
         location = str_to_lower(trimws(location)),
         season = str_to_title(trimws(season)),
         area = factor(area))

#Set theme
habitat_colors <- c("creek" = "#F06C57", "lake" = "#4A90B8")
base_theme <- theme_bw() + theme(
  panel.grid = element_blank(), panel.border = element_blank(),
  axis.line = element_line(color = "black"), legend.position = "none",
  axis.title = element_text(size = 12), axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 10))


##### Shiner body length model #####
df_size_data <- df_tl_num %>%
  filter(species %in% c("common shiner", "golden shiner", "smallmouth bass"), 
         is.finite(tl), tl > 0, number == 1, !is.na(area), !site %in% paste0("EF", 1:5), 
         location %in% c("creek", "lake"), year %in% 2022:2024,
         season %in% c("Spring", "Summer", "Fall"), 
         !is.na(id_conflict), !id_conflict, !is.na(location_conflict), !location_conflict) %>%
  mutate(season = factor(season, levels = c("Spring", "Summer", "Fall")),
         location = factor(location, levels = c("creek", "lake")))

df_shiner_size <- df_size_data %>%
  filter(species %in% c("common shiner", "golden shiner")) %>% 
  droplevels()

#Body Length glm model
glm_size_shiner <- glmmTMB(tl ~ season * location + (1 | area),
                           data = df_shiner_size, family = gaussian())

Anova(glm_size_shiner, type = 3)
summary(glm_size_shiner)
r2_nakagawa(glm_size_shiner)

sim_size_shiner <- simulateResiduals(glm_size_shiner)
plot(sim_size_shiner)

emm_size_shiner <- emmeans(glm_size_shiner, ~ season * location)
shiner_size_contrasts <- as.data.frame(summary(pairs(emm_size_shiner, adjust = "tukey"), infer = c(FALSE, TRUE)))

shiner_size_contrasts

df_emm_size_shiner <- as.data.frame(emm_size_shiner)

plot_size_shiner <- ggplot(
  df_emm_size_shiner, aes(x = season, y = emmean, group = location)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.15, linewidth = 0.6, 
                color = "#323332", alpha = 0.7, position = position_dodge(width = 0.4)) +
  geom_point(aes(fill = location), size = 3, shape = 21, 
             color = "#323332", stroke = 0.6, position = position_dodge(width = 0.4)) +
  scale_fill_manual(name = "Location", values = habitat_colors, 
                    breaks = c("creek", "lake"), labels = c("Creek", "Lake")) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  labs(x = NULL, y = "Shiner Total Length (mm)") +
  coord_cartesian(clip = "off") + 
  base_theme

plot_size_shiner

##### Smallmouth Bass Model #####
df_smb_size <- df_size_data %>%
  filter(species == "smallmouth bass") %>%
  droplevels()

glm_size_smb <- glmmTMB(tl ~ season*location + (1 | area),
                        data = df_smb_size, family = Gamma(link = "log"),
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

Anova(glm_size_smb, type = 3)

glm_size_smb <- glmmTMB(tl ~ season + location + (1 | area),
                        data = df_smb_size, family = Gamma(link = "log"),
                        control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

Anova(glm_size_smb, type = 2)
r2_nakagawa(glm_size_smb, tolerance = 1e-12)

summary(glm_size_smb)

sim_size_smb <- simulateResiduals(glm_size_smb)
plot(sim_size_smb)

emm_size_smb_season <- emmeans(glm_size_smb, ~ season, type = "response")
smb_size_season_contrasts <- as.data.frame(summary(
  pairs(emm_size_smb_season, adjust = "tukey"), infer = c(FALSE, TRUE), type = "response"
))
smb_size_season_contrasts

emm_size_smb_habitat <- emmeans(glm_size_smb, ~ location, type = "response")
smb_size_habitat_contrasts <- as.data.frame(summary(
  pairs(emm_size_smb_habitat, adjust = "tukey"), infer = c(FALSE, TRUE), type = "response"
))

smb_size_habitat_contrasts

df_emm_size_smb <- emmeans(glm_size_smb, ~ season + location, type = "response") %>%
  as.data.frame()


plot_size_smb <- ggplot(
  df_emm_size_smb, aes(x = season, y = response, group = location)
) +
  geom_errorbar(
    aes(ymin = response - SE, ymax = response + SE),
    width = 0.15, linewidth = 0.6, color = "#323332", alpha = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    aes(fill = location), size = 3, shape = 21, color = "#323332", stroke = 0.6,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(
    name = "Location", values = habitat_colors, breaks = c("creek", "lake"),
    labels = c("Creek", "Lake")
  ) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  labs(x = NULL, y = "SMB Total Length (mm)") +
  coord_cartesian(clip = "off") + base_theme

plot_size_smb

figure_6_combined <- ggarrange(
  plot_size_shiner, plot_size_smb, ncol = 1, heights = c(1, 1), align = "v"
)

figure_6_combined

ggsave(file.path(file_path, "Figures/raw_r_output_figs/figure_6.jpg"), figure_6_combined,
       width = 6, height = 5, dpi = 900)


##### Supplementary Sensitivity Analyses #####
#Check influence of species and year on model
shiner_size_species_year <- glmmTMB(tl ~ season * location + species + year + (1 | area),
                                    data = df_shiner_size, family = gaussian(),
                                    control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

Anova(shiner_size_species_year, type = 3)
r2_nakagawa(smb_size_year)

summary(shiner_size_species_year)

summary(pairs(emmeans(shiner_size_species_year, ~ location | season)), by = NULL, adjust = "holm")

#Check year in SMB model
df_smb_size <- df_smb_size %>%
  mutate(year = factor(year)) %>%
  droplevels()

smb_size_year <- glmmTMB(tl ~ season + location + year + (1 | area),
                         data = df_smb_size, family = Gamma(link = "log"),
  control = glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000)))

#No significant influence of year
Anova(smb_size_year, type = 2)
r2_nakagawa(smb_size_year)

emmeans(smb_size_year,~ year,type = "response")

#Tukey-adjusted pairwise comparisons among years
smb_size_year_contrasts <- as.data.frame(pairs(emmeans(smb_size_year, ~ year), adjust = "tukey"))
smb_size_year_contrasts



