
##Read in fyke data from 2023 
fyke_data<- read.csv("cleaned_fyke_data_2023.csv", header=T)


head(fyke_data)

library(tidyverse)
library(lubridate)
library(stringr)


fyke_summ_data <- fyke_data %>%
  group_by(Month, Direction, Species) %>%
  summarize(total_biomass = sum(Wt, na.rm = TRUE)) %>%
  ungroup()

# Calculate average daily biomass per month
daily_summ_data <- fyke_summ_data %>%
  group_by(Month, Direction, Species) %>%
  mutate(daily_avg_biomass = mean(total_biomass))

# Filter the data to include only "Creek Chub"
creek_chub_data <- fyke_summ_data %>%
    filter(Species == c("Creek Chub"))

month_order <- factor(creek_chub_data$Month, levels = c("May", "August", "October"))

##CC 
ggplot(data = creek_chub_data, aes(x = month_order, y = total_biomass, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  theme_minimal() +  # Using a minimal theme as a base
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10))
  )+
  scale_x_discrete(labels = c("May", "August", "October")) +
  scale_fill_manual(values = c('azure3', 'azure4')) +
  labs(x = "Month", y = "Total Biomass (g)")+
  ylim(0,3000)


# Filter the data to include only "Common Shiner"
common_shiner_data <- fyke_summ_data %>%
    filter(Species == "Common Shiner")

month_order <- factor(common_shiner_data$Month, levels = c("May", "August", "October"))


##CS 
ggplot(data = common_shiner_data, aes(x = month_order, y = total_biomass, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  theme_minimal() +  # Using a minimal theme as a base
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10))
  )+
  scale_x_discrete(labels = c("May", "August", "October")) +
  scale_fill_manual(values = c('azure3', 'azure4')) +
  labs(x = "Month", y = "Total Biomass (g)")+
  ylim(0,3000)

  

# Filter the data to include only "Golden Shiner"
golden_shiner_data <- fyke_summ_data %>%
  filter(Species == "Golden Shiner")

# Define the order of months
month_order <- factor(golden_shiner_data$Month, levels = c("May", "August", "October"))

##GS
ggplot(data = golden_shiner_data, aes(x = month_order, y = total_biomass, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  theme_minimal() +  # Using a minimal theme as a base
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10))
  )+
  scale_x_discrete(labels = c("May", "August", "October")) +
  scale_fill_manual(values = c('azure3', 'azure4')) +
  labs(x = "Month", y = "Total Biomass (g)")+
  ylim(0,3000)


