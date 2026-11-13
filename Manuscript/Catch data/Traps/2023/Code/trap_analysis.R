
##Read in weir data from 2023 
trap_data<- read.csv("/Users/charlotteward/Documents/Field Study Files/Opeongo_Costello_2023/trap_data_2023.csv", header=T)
head(trap_data)

library(tidyverse)
library(lubridate)
library(stringr)


# List of species to include
selected_species <- c("Golden Shiner", "Common Shiner", "Northern Pearl Dace", "Creek Chub", "Yellow Perch")

# Filter the data for the selected species
filtered_data <- trap_data %>%
  filter(Species %in% selected_species) %>%
  mutate(
    tl = as.numeric(tl),
    month = Month,
    month_order = factor(month, levels = c("May", "August", "October"))
  )


# Create the plot
my_plot <- ggplot(filtered_data, aes(x = tl, color = direction)) +
  geom_density(linewidth = 0.9) +
  facet_wrap(~ Species + month_order, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Downstream" = "azure3", "Upstream" = "azure4")) +
  labs(x = "Total Length (mm)", y = "Density", color = "Direction") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18, margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 18, margin = margin(l = 10)),
    panel.spacing = unit(2, "lines"),
    strip.text = element_text(size = 12, face = "bold")
  )

my_plot
# Save the plot
ggsave("selected_species_plot.png", my_plot, width = 12, height = 8, dpi = 300)


## summarize data for biomass plots 
trap_summ_data <- trap_data %>%
  mutate(
         wt = as.numeric(wt)) %>%  # Convert wt to numeric
  group_by(Month, Species, direction) %>%
  summarize(
    total_biomass = sum(wt, na.rm = TRUE),  # Calculate total biomass for the month
    average_mass = mean(wt, na.rm = TRUE),# Calculate average daily biomass
    average_fl = mean(fl, na.rm = TRUE), 
    average_tl = mean(tl, na.rm = TRUE)
  ) %>%
  ungroup()

##### Golden Shiner Data

# Filter the data to include only "Golden Shiner"
golden_shiner_data <- trap_summ_data %>%
  filter(Species == "Golden Shiner")

# Define the order of months
month_order <- factor(golden_shiner_data$Month, levels = c("May", "August", "October"))

# Create the plot for total monthly biomass
ggplot(data = golden_shiner_data, aes(x = month_order, y = total_biomass, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  theme_minimal() +  # Using a minimal theme as a base
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
  scale_fill_manual(values = c('azure3', 'azure4')) +
  labs(x = "Month", y = "Total Biomass (g)")+
  ylim(0, 6000)


##### Yellow Perch Data

# Filter the data to include only "Golden Shiner"
yellow_perch_data <- trap_summ_data %>%
  filter(Species == "Yellow Perch")

# Define the order of months
month_order <- factor(yellow_perch_data$Month, levels = c("May", "August", "October"))

# Create the plot for total monthly biomass
ggplot(data = yellow_perch_data, aes(x = month_order, y = total_biomass, fill = direction)) +
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
  labs(x = "Month", y = "Total Biomass (g)")


