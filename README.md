# lake_stream_coupling

This repository contains the data processing, statistical analysis, and figuregeneration code for the manuscript:

**Habitat coupling by small fish across lake–stream habitats: a missing link in freshwater food webs?**

The analyses evaluate seasonal habitat coupling by small-bodied fishes, focusing on golden shiner (*Notemigonus crysoleucas*) and common shiner (*Luxilus cornutus*) in a connected lake–stream ecosystem. The scripts combine stable isotope data, acoustic telemetry detections, and seasonal/directional trapping data to quantify trophic, behavioural, and biomass-flux dimensions of habitat coupling.

Note that in scripts, stream is coded as "creek" 

## Repository structure

```text
lake_stream_coupling/
├── Data/                  # Input data used by the analysis scripts
├── Code/                  # R scripts for data processing, models, and figures
│   ├── 0 - Functions.R
│   ├── 1 - Figure 2 - Stream Map & Detections.R
│   ├── 2 - Figure 3 - Seasonal Isotopes & Telemetry.R
│   ├── 3 - Figure 4 - Telemetry GAMM Models.R
│   ├── 4 - Figure 5 -  Biomass Models & Plots.R
│   ├── 5 - Figure 6 - Fish Size Models & Plots.R
│   └── 6 - Appendix S1 Plots.R 
├── lake_stream_coupling.Rproj
└── README.md
```

## How to run the analyses

Open `lake_stream_coupling.Rproj` in RStudio and run scripts from the project root. Each script sets:

```r
file_path <- getwd()
```

and reads files using paths such as `Data/...` and `Code/...`. Running from the project root is therefore important for the relative paths to work correctly.

Recommended script order:

1. `Code/0 - Functions.R`
2. `Code/1 - Figure 2 - Stream Map & Detections.R`
3. `Code/2 - Figure 3 - Seasonal Isotopes & Telemetry Analysis.R`
4. `Code/3 - Figure 4 - GAM Models.R`
5. `Code/4 - Figure 5 -  Biomass Flux.R`
6. `Code/5 - Supplement.R`


## R packages used

The scripts use the following R packages:

```r
library(tidyverse)
library(lubridate)
library(purrr)
library(glatos)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
library(sf)
library(terra)
library(IsoriX)
library(igraph)
library(vegan)
library(ggspatial)
library(mgcv)
library(gratia)
library(car)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library(bbmle)
library(easystats)
library(circlize)
library(grid)
```

The `glatos` package is used to remove likely false acoustic telemetry detections. If needed, it can be installed from the Ocean Tracking Network R-universe repository:

```r
install.packages(
  "glatos",
  repos = c("https://ocean-tracking-network.r-universe.dev", "https://cloud.r-project.org")
)
```

## Script descriptions

### `0 - Functions.R`

This script contains helper functions shared across the repository. It is sourced by the main analysis scripts and is not intended to generate figures by itself.

Key functions include:

- `standard_error()`: calculates standard error while ignoring missing values.
- `print_emm_pairs_iso()`: prints estimated marginal means and Tukey-adjusted pairwise comparisons for isotope models across month and location.
- `print_emm_pairs_tel()`: prints estimated marginal means and Tukey-adjusted pairwise comparisons for telemetry models across season and release location.
- `calculate_prop_tp()`: calculates the proportion of lake-derived carbon and trophic position for selected fish species and tissue types using stream and lake mayfly isotope baselines.
- `create_time_series()`: creates a one-minute telemetry time series for each tagged fish, joining detections to the deployment period and carrying the last known location forward.
- `enriched_time_series()`: similar to `create_time_series()`, but used in supplementary telemetry-processing plots.



---

### `1 - Figure 2 - Stream Map & Detections.R`

This script generates the study-system map and receiver-detection visualization used for **Figure 2**.

Main steps:

1. Loads lake and stream spatial layers from `Data/Lakes_Streams/`.
2. Loads receiver metadata from `Data/receivers.csv`.
3. Loads cleaned telemetry detections from `Data/detections_clean_alldata.csv`.
4. Removes specified transmitter IDs from the detection data.
5. Calculates detection frequency by receiver.
6. Builds a movement network based on sequential changes in receiver detections for each fish.
7. Loads isotope baseline-coordinate data from `Data/baselines_coords_d13C.csv`.
8. Uses `IsoriX` to model and map spatial variation in baseline δ13C across the lake–stream study area.
9. Overlays receiver locations, scaled by detection frequency, on the waterbody/isoscape map.



---

### `2 - Figure 3 - Seasonal Isotopes & Telemetry Analysis.R`

This script performs the main seasonal isotope and telemetry analyses and generates **Figure 3**.

Main isotope steps:

1. Loads isotope metadata from `Data/iso_metadata.csv`.
2. Filters to golden shiner and common shiner liver samples.
3. Uses `calculate_prop_tp()` to estimate:
   - proportion lake-derived carbon, and
   - trophic position.
4. Summarizes mayfly and mussel baseline isotope values by location and season.
5. Fits baseline isotope models for δ13C and δ15N.
6. Fits shiner isotope models for raw δ13C and δ15N.
7. Fits a beta-family model for proportion lake-derived carbon.
8. Fits a Gaussian model for trophic position.
9. Calculates estimated marginal means, pairwise contrasts, and seasonal/location deltas.
10. Uses DHARMa residual diagnostics for model checking.

Main telemetry steps:

1. Loads cleaned acoustic telemetry detections from `Data/detections_clean_alldata.csv`.
2. Removes selected transmitter IDs and receivers from the analysis.
3. Converts tag-step timing, release dates, detection dates, and tag-on/tag-off windows to date-time objects.
4. Applies a false-detection filter using `glatos::false_detections()`.
5. Removes fish/tag records with invalid or mismatched detection windows.
6. Builds one-minute location time series for each tagged fish.
7. Assigns observations to seasons.
8. Calculates seasonal habitat residency as the proportion of detected minutes spent in the fish’s reference habitat (where it was first detected).
9. Calculates seasonal habitat recurrence as the number of return events to the reference habitat after an excursion to the alternative habitat.
10. Fits telemetry models for:
   - habitat residency, using a beta-family model; and
   - habitat recurrence, using a negative-binomial model.

Figure output:

The script combines four panels into `Figures/figure_3.jpg`:

- proportion lake-derived carbon,
- trophic position,
- habitat recurrence, and
- habitat residency.

---

### `3 - Figure 4 - GAM Models.R`

This script fits generalized additive models to telemetry detections and generates **Figure 4**.

Main steps:

1. Loads cleaned telemetry detections from `Data/detections_clean_alldata.csv`.
2. Removes specified transmitters and receivers.
3. Converts detection timestamps, release dates, and tag windows to date-time objects.
4. Applies the `glatos` false-detection filter.
5. Removes invalid tag windows.
6. Builds one-minute location time series for each tagged fish.
7. Calculates a binary response indicating whether a fish was detected in the stream and lake in the same day.
8. Fits candidate GAM/GAMM models using cyclic smooths over day-of-year.
9. Compares models with and without random effects and AR(1) autocorrelation structure.
10. Selects the final model based on model fit. 
11. Generates predicted seasonal probabilities of stream detection with confidence intervals.
12. Identifies peak and trough periods in predicted lake-stream use probability.
13. Runs residual diagnostics, including checks for autocorrelation.

Figure output:

The script saves `Figures/figure_4.jpg`, showing the predicted probability of lake–stream detection/use through the annual cycle.


---

### `4 - Figure 5 -  Biomass Flux.R`

This script analyzes seasonal minnow-trapping data and directional trap data and generates **Figure 5**.

Main CPUE steps:

1. Loads minnow trap data from `Data/minnowtrap_data.csv`.
2. Filters to golden shiner and common shiner.
3. Aggregates catch and biomass by month, day, and location.
4. Calculates log biomass per trap as a CPUE response.
5. Fits a `glmmTMB` model testing month, location, and their interaction.
6. Uses Type II or III tests, estimated marginal means, and pairwise contrasts to assess seasonal and habitat differences.
7. Uses DHARMa residual diagnostics for model checking.

Main directional-flux steps:

1. Loads directional trap data from `Data/directional_traps_all.csv`.
2. Filters to golden shiner and common shiner.
3. Aggregates daily biomass by year, season, day, and movement direction.
4. Calculates log biomass flux.
5. Fits model testing season, direction, and their interaction, while including year.
6. Uses Type II or III tests, estimated marginal means, and pairwise contrasts to assess seasonal and directional differences.

Figure output:

The script saves `Figures/figure_5.jpg`, combining:

- seasonal biomass per trap by lake versus creek location, and
- directional biomass flux

---

### `5 - Figure 6 -  Fish Size Models & PLots.R`

This script compares fish total length (mm) between creek and lake habitats across spring, summer, and fall and generates Figure 6.

Data preparation

1. Loads merged fish catch records from `Data/fish_catches_all_merged.csv`.
2. Retains individually recorded common shiner, golden shiner, and smallmouth bass with valid total lengths from 2022 to 2024. 

Main steps:

1. Filters to individually recorded common shiner and golden shiner and smallmouth bass total lengths.
2. Fits model testing season, habitat, and their interaction, with sampling area as a random intercept.
3. Uses Type II or III test, estimated marginal means, and pairwise contrasts to assess seasonal and habitat differences.
4. Uses DHARMa residual diagnostics for model checking.


Figure output

The script saves `Figures/figure_6.jpg`, combining:
- estimated mean total length (mm) by season and habitat for pooled shiners, and
- estimated mean total length (mm) by season and habitat for smallmouth bass.

