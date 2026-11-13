
##Read in weir data from 2023 
weir_data<- read.csv("weir_data.csv", header=T)
head(weir_data)

library(tidyverse)
library(lubridate)
library(stringr)

weir_data$date <- ymd(weir_data$date) # Convert to Date format if not already

weir_data <- weir_data %>%
  mutate(date = as.character(date),  # Convert to character if it's in Date format
         date = substring(date, 3),  # Remove the first two digits of the year
         date = dmy(date))  # Convert to Date format with dd-mm-yy structure


weir_data <- weir_data %>%
  filter(date != as.Date("2023-05-19"))

weir_data <- weir_data %>%
  mutate(
    tl = as.numeric(tl),  # Convert tl to numeric if it's not already
    month = month(date),
    month = case_when(
      month == 9 ~ "October",  # Treat September as October
      TRUE ~ month(date, label = TRUE, abbr = FALSE)  # Extract month as full name
    ),
    month_order = factor(month, levels = c("May", "August", "October"))
  )

# List of species to include
selected_species <- c("Golden Shiner", "Common Shiner", "Northern Pearl Dace", "Creek Chub")

# Filter the data for the selected species
filtered_data <- weir_data %>%
  filter(species %in% selected_species) %>%
  mutate(
    tl = as.numeric(tl),
    month = month(date),
    month = case_when(
      month == 9 ~ "October",
      TRUE ~ month(date, label = TRUE, abbr = FALSE)
    ),
    month_order = factor(month, levels = c("May", "August", "October"))
  )

# Create the plot
my_plot <- ggplot(filtered_data, aes(x = tl, color = direction)) +
  geom_density(linewidth = 0.9) +
  facet_wrap(~ species + month_order, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Downstream" = "azure3", "Upstream" = "azure4")) +
  labs(x = "Total Length (mm)", y = "Density", color = "Direction") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)),
    panel.spacing = unit(2, "lines"),
    strip.text = element_text(size = 12, face = "bold")
  )

# Save the plot
ggsave("selected_species_plot.png", my_plot, width = 12, height = 8, dpi = 300)


## summarize data for biomass plots 
weir_summ_data <- weir_data %>%
  mutate(month = month(date), 
         day = day(date),
         wt = as.numeric(wt)) %>%  # Convert wt to numeric
  mutate(month = case_when(
    month == 9 ~ "October",  # Treat September as October
    TRUE ~ month.name[month]  # Convert other months to full names
  )) %>%
  group_by(month, species, direction) %>%
  summarize(
    total_biomass = sum(wt, na.rm = TRUE),  # Calculate total biomass for the month
    average_mass = mean(wt, na.rm = TRUE),# Calculate average daily biomass
    average_fl = mean(fl, na.rm = TRUE), 
    average_tl = mean(tl, na.rm = TRUE)
  ) %>%
  ungroup()

##### Golden Shiner Data

# Filter the data to include only "Golden Shiner"
golden_shiner_data <- weir_summ_data %>%
  filter(species == "Golden Shiner")

# Define the order of months
month_order <- factor(golden_shiner_data$month, levels = c("May", "August", "October"))

# Create the plot for total monthly biomass
ggplot(data = golden_shiner_data, aes(x = month_order, y = total_biomass, fill = direction)) +
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
  ylim(0, 5000)

##size distributions

weir_data <- weir_data %>%
  mutate(tl = as.numeric(tl),  # Convert tl to numeric if it's not already
         month = month(date),
         month = case_when(
           month == 9 ~ "October",  # Treat September as October
           TRUE ~ month(date, label = TRUE, abbr = FALSE)  # Extract month as full name
         ))

weir_data$month_order <- factor(weir_data$month, levels = c("May", "August", "October"))


ggplot(weir_data, aes(x = tl, color = direction)) +
  geom_density(linewidth = 0.9) +
  facet_wrap(~ month_order) +  
  scale_color_manual(values = c("Downstream" = "azure3", "Upstream" = "azure4")) +
  labs(
    x = "Total Length (mm)",
    y = "Frequency",
    color = "direction") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey"),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)),
    panel.spacing = unit(2, "lines"), 
    strip.text = element_text(size = 12, face = "bold")
  )


##### Common Shiner Data

# Filter the data to include only "Golden Shiner"
common_shiner_data <- weir_summ_data %>%
  filter(species == "Common Shiner")

# Define the order of months
month_order <- factor(common_shiner_data$month, levels = c("May", "August", "October"))

# Create the plot for total monthly biomass
ggplot(data = common_shiner_data, aes(x = month_order, y = total_biomass, fill = direction)) +
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
  ylim(0, 8500)

##### Common Shiner Data

# Filter the data to include only "Golden Shiner"
creek_chub_data <- weir_summ_data %>%
  filter(species == "Creek Chub")

# Define the order of months
month_order <- factor(creek_chub_data$month, levels = c("May", "August", "October"))

# Create the plot for total monthly biomass
ggplot(data = creek_chub_data, aes(x = month_order, y = total_biomass, fill = direction)) +
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
  ylim(0, 8500)
  

##### Yellow Perch Data

# Filter the data to include only "Golden Shiner"
yellow_perch_data <- weir_summ_data %>%
  filter(species == "Yellow Perch")

# Define the order of months
month_order <- factor(yellow_perch_data$month, levels = c("May", "August", "October"))

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

