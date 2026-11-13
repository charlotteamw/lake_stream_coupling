
##Read in MT data from 2023 
mt_data<- read.csv("/Users/charlotteward/Documents/Field Study Files/Opeongo_Costello_2023/minnowtrap_data_2023.csv", header=T)


head(mt_data)


library(tidyverse)
library(lubridate)
library(stringr)

summ_data <- mt_data %>%
  group_by(Month, Location, Species) %>%
  summarize(total_count = sum(Number, na.rm = TRUE), total_biomass = sum(Wt, na.rm = TRUE))         

head(summ_data)


# Filter the data to include only "creek chub"
  creek_chub_data <- summ_data %>%
    filter(Species == "Creek Chub")

# Define the order of months
  month_order <- factor(creek_chub_data$Month, levels = c("May", "August", "October"))
# Create the plot
ggplot(data = creek_chub_data, aes(x = month_order, y = total_biomass, fill = Location)) +
    geom_bar(stat = "identity" , position = position_dodge2(preserve = "single")) +
    theme_minimal() +
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
    scale_fill_manual(values = c('cadetblue', 'sandybrown')) +
    labs(x = "Month", y = "Total Biomass (g)")+
    ylim(0,2500)

# Filter the data to include only "Golden Shiner"
cs_mt_data <- summ_data %>%
  filter(Species == "Common Shiner")

# Define the order of months
month_order <- factor(cs_mt_data$Month, levels = c("May", "August", "October"))
# Create the plot
ggplot(data = cs_mt_data, aes(x = month_order, y = total_biomass, fill = Location)) +
  geom_bar(stat = "identity" , position = position_dodge2(preserve = "single")) +
  theme_minimal() +
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
  scale_fill_manual(values = c('cadetblue', 'sandybrown')) +
  labs(x = "Month", y = "Total Biomass (g)")+
  ylim(0,1500)

 
  # Filter the data to include only "Golden Shiner"
  gs_mt_data <- summ_data %>%
    filter(Species == "Golden Shiner")
  
  # Define the order of months
  month_order <- factor(gs_mt_data$Month, levels = c("May", "August", "October"))
  # Create the plot
  ggplot(data = gs_mt_data, aes(x = month_order, y = total_biomass, fill = Location)) +
    geom_bar(stat = "identity" , position = position_dodge2(preserve = "single")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "grey"),
      axis.text = element_text(size = 16),
      axis.title.x = element_text(size = 18, margin = margin(t = 10)),
      axis.title.y = element_text(size = 18, margin = margin(r = 10)),
      axis.title.y.right = element_text(size = 18, margin = margin(l = 10))
    )+
    scale_x_discrete(labels = c("May", "August", "October")) +
    scale_fill_manual(values = c('cadetblue', 'sandybrown')) +
    labs(x = "Month", y = "Total Biomass (g)")+
    ylim(0,500)
  

  # Filter the data to include only "Yellow Perch"
  yp_mt_data <- summ_data %>%
    filter(Species == "Yellow Perch")
  
  # Define the order of months
  month_order <- factor(yp_mt_data$Month, levels = c("May", "August", "October"))
  # Create the plot
  ggplot(data = yp_mt_data, aes(x = month_order, y = total_biomass, fill = Location)) +
    geom_bar(stat = "identity" , position = position_dodge2(preserve = "single")) +
    theme_minimal() +
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
    scale_fill_manual(values = c('cadetblue', 'sandybrown')) +
    labs(x = "Month", y = "Total Biomass (g)")+
    ylim(0,1500)
  

