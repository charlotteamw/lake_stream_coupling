# =============================================================================
# Combine catch records and processed fish (2022-2024)
# Updated: 15 October 2026
#
# INPUTS:
#   Data/Catch data/Traps
#   Data/Catch data/processed_all.csv
#
# OUTPUT:
#   Data/fish_catches_all_merged_recoded.csv
#
# CHECKS:
#   Data/Merge checks recoded
#
# One row per specimen code.
# Processed measurements take priority.
# Catch records supply missing information.
# Catches without specimen codes remain separate records.
# Bulk catches remain bulk catches.
#
# Includes all 21 confirmed A-coded smallmouth bass.
# =============================================================================


# =============================================================================
# 1. PACKAGES AND FOLDERS
# =============================================================================

library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)

project_dir <- getwd()

source(file.path(project_dir, "Code/0 - Functions.R"))

data_dir <- file.path(project_dir, "Data", "Catch data")
trap_dir <- file.path(data_dir, "Traps")

processed_file <- file.path(data_dir, "processed_all.csv")

output_file <- file.path(
  project_dir,
  "Data",
  "fish_catches_all_merged_recoded.csv"
)

checks_dir <- file.path(project_dir, "Data", "Merge checks recoded")


# =============================================================================
# 2. CATCH FILES
# =============================================================================

# Only these catch files are imported.
# Excel files are read from their first sheet.
# Other files are excluded to avoid importing catches twice.

catch_files <- c(
  "2023/Data/minnowtrap_data_2023.csv",
  "2023/Data/trap_data_2023.csv",
  "2024/Data/Fyke_2024.xlsx",
  "2024/Data/Weir_2024.xlsx",
  "2022/August Digitized Data Sheets/August Weir Data_2022.xlsx",
  "2022/August Digitized Data Sheets/August Fyke Net Data_2022.xlsx",
  "2022/August Digitized Data Sheets/August Efishing Data_2022.xlsx",
  "2022/September Digitized Data Sheets/September Minnow Trap Data_2022.xlsx",
  "2022/September Digitized Data Sheets/September Weir_2022.xlsx",
  "2022/May Digitized Data Sheets/May Electrofishing Data_2022.xlsx",
  "2022/May Digitized Data Sheets/May Fyke Data_2022.xlsx",
  "2022/May Digitized Data Sheets/May Minnow Trap Data_2022.xlsx",
  "2022/May Digitized Data Sheets/May Weir Data_2022.xlsx"
)

catch_paths <- file.path(trap_dir, catch_files)

missing_files <- c(
  catch_paths,
  processed_file
)[
  !file.exists(c(catch_paths, processed_file))
]

if (length(missing_files)) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n")
  )
}


# =============================================================================
# 3. CONFIRMED A-CODED SMALLMOUTH BASS
# =============================================================================

# The A prefix does not identify the capture site.
# Confirmed capture sites are assigned individually.
#
# The two FYKE fish were caught in the creek.

bass_sites <- tribble(
  ~id,        ~site, ~catch_method, ~month,   ~capture_site,
  
  "A-22-01",  "LK5",  "angling", "May",    "Across from Op store",
  "A-22-02",  "LK10", "angling", "May",    "Around corner from Harkness",
  "A-22-03",  "LK10", "angling", "May",    "Around corner from Harkness",
  "A-22-04",  "FYKE", "fyke",    "May",    "US",
  "A-22-05",  "LK9",  "angling", "May",    "EF3",
  "A-22-06",  "LK9",  "angling", "May",    "EF3",
  "A-22-07",  "FYKE", "fyke",    "May",    "US",
  "A-22-08",  "LK2",  "angling", "May",    "EF1",
  "A-22-09",  "LK2",  "angling", "May",    "EF1",
  "A-22-10",  "LK1",  "angling", "May",    "Harkness dock",
  
  "A-A22-01", "LK",   "angling", "August", "Sproule",
  "A-A22-02", "LK",   "angling", "August", "Sproule",
  "A-A22-03", "LK",   "angling", "August", "Sproule",
  "A-A22-04", "LK",   "angling", "August", "Sproule",
  "A-A22-05", "LK",   "angling", "August", "Sproule",
  "A-A22-06", "LK",   "angling", "August", "Sproule",
  "A-A22-07", "LK",   "angling", "August", "Sproule",
  "A-A22-08", "LK5",  "angling", "August", "Sproule; caught OP15",
  "A-A22-09", "LK",   "angling", "August", "Sproule",
  "A-A22-10", "LK",   "angling", "August", "Sproule",
  "A-A22-11", "LK",   "angling", "August", "Sproule"
)

# Preserve confirmed dates.
# Missing or conflicting dates are not inferred.

bass_sites <- bass_sites %>%
  mutate(
    date = case_when(
      id == "A-22-01" ~ "2022-05-15",
      id == "A-22-04" ~ "2022-05-26",
      
      id %in% c(
        "A-22-05",
        "A-22-06",
        "A-22-07"
      ) ~ "2022-05-24",
      
      id %in% c(
        "A-22-08",
        "A-22-09",
        "A-22-10"
      ) ~ "2022-05-18",
      
      id %in% sprintf(
        "A-A22-%02d",
        1:5
      ) ~ "2022-08-29",
      
      TRUE ~ NA_character_
    )
  )

stopifnot(
  nrow(bass_sites) == 21,
  !anyDuplicated(bass_sites$id)
)


# =============================================================================
# 4. HELPER FUNCTIONS
# =============================================================================

clean_text <- function(x) {
  
  x <- str_squish(as.character(x))
  
  x[toupper(x) %in% c(
    "",
    "NA",
    "N/A",
    "NULL"
  )] <- NA_character_
  
  x
}


clean_id <- function(x) {
  
  str_replace_all(
    toupper(clean_text(x)),
    "_",
    "-"
  )
}


as_number <- function(x) {
  
  suppressWarnings(
    as.numeric(clean_text(x))
  )
}


first_known <- function(x) {
  
  if (all(is.na(x))) {
    return(x[1])
  }
  
  x[which(!is.na(x))[1]]
}


# Read whichever equivalent column exists in a source file.

get_column <- function(data, names_to_try) {
  
  value <- rep(
    NA_character_,
    nrow(data)
  )
  
  for (name in intersect(
    names_to_try,
    names(data)
  )) {
    
    value <- dplyr::coalesce(
      value,
      clean_text(data[[name]])
    )
  }
  
  value
}


# Convert month abbreviations and numbers to full month names.

get_month <- function(x) {
  
  month_number <- match(
    tolower(substr(x, 1, 3)),
    tolower(month.abb)
  )
  
  numeric_month <- as_number(x)
  
  use_number <- is.na(month_number) &
    !is.na(numeric_month) &
    numeric_month %in% 1:12
  
  month_number[use_number] <- numeric_month[use_number]
  
  month.name[month_number]
}


# Parse dates from CSV and Excel files.

get_date <- function(x, month_first = FALSE) {
  
  result <- as.Date(
    rep(NA_character_, length(x))
  )
  
  for (i in seq_along(x)) {
    
    if (is.na(x[i])) next
    
    if (grepl(
      "^[0-9]{5}(\\.[0-9]+)?$",
      x[i]
    )) {
      
      result[i] <- as.Date(
        as.numeric(x[i]),
        origin = "1899-12-30"
      )
      
    } else {
      
      formats <- if (grepl(
        "^[0-9]{4}-",
        x[i]
      )) {
        
        "%Y-%m-%d"
        
      } else if (grepl(
        "/[0-9]{4}$",
        x[i]
      )) {
        
        if (month_first) {
          c("%m/%d/%Y", "%d/%m/%Y")
        } else {
          c("%d/%m/%Y", "%m/%d/%Y")
        }
        
      } else {
        
        c("%d/%m/%y", "%d-%m-%Y")
      }
      
      for (format in formats) {
        
        parsed <- as.Date(
          x[i],
          format = format
        )
        
        if (!is.na(parsed)) {
          
          result[i] <- parsed
          break
        }
      }
    }
  }
  
  result
}


# =============================================================================
# 6. READ AND STANDARDIZE CATCH FILES
# =============================================================================

read_catches <- function(path) {
  
  filename <- basename(path)
  
  is_csv <- grepl(
    "\\.csv$",
    path,
    ignore.case = TRUE
  )
  
  data <- if (is_csv) {
    
    read_csv(
      path,
      col_types = cols(
        .default = col_character()
      ),
      show_col_types = FALSE
    )
    
  } else {
    
    read_excel(
      path,
      sheet = 1,
      col_types = "text",
      .name_repair = "unique_quiet"
    )
  }
  
  names(data) <- make.unique(
    str_replace_all(
      str_replace_all(
        tolower(trimws(names(data))),
        "[^a-z0-9]+",
        "_"
      ),
      "^_|_$",
      ""
    )
  )
  
  source_row <- seq_len(nrow(data)) + 1L
  
  species <- str_to_lower(
    get_column(
      data,
      c(
        "species",
        "species_name",
        "sp_name"
      )
    )
  )
  
  keep <- !is.na(species)
  
  data <- data[keep, ]
  species <- species[keep]
  source_row <- source_row[keep]
  
  
  # ---------------------------------------------------------------------------
  # STANDARDIZE SPECIES NAMES
  #
  # FIXED:
  # Replaces recode(species, !!!species_names).
  #
  # Uses named-vector matching to avoid the !!! evaluation error.
  # ---------------------------------------------------------------------------
  
  species_names <- c(
    "c shiner" = "common shiner",
    "common shin" = "common shiner",
    "g shiner" = "golden shiner",
    "smb" = "smallmouth bass",
    "c chub" = "creek chub",
    "y perch" = "yellow perch",
    "perch" = "yellow perch",
    "p seed" = "pumpkinseed",
    "pumpkin seed" = "pumpkinseed",
    "pumpkin" = "pumpkinseed",
    "w sucker" = "white sucker",
    "bnm" = "bluntnose minnow",
    "blunt nose" = "bluntnose minnow",
    "bluntnose" = "bluntnose minnow",
    "black nose shiner" = "blacknose shiner",
    "fall fish" = "fallfish"
  )
  
  species <- dplyr::coalesce(
    unname(species_names[species]),
    species
  )
  
  
  # ---------------------------------------------------------------------------
  # YEAR AND MONTH
  # ---------------------------------------------------------------------------
  
  year <- as.integer(
    str_extract(
      filename,
      "20[0-9]{2}"
    )
  )
  
  month <- coalesce(
    get_month(
      get_column(data, "month")
    ),
    str_extract(
      filename,
      paste(month.name, collapse = "|")
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # CATCH METHOD
  # ---------------------------------------------------------------------------
  
  default_method <- case_when(
    
    str_detect(
      tolower(filename),
      "minnow"
    ) ~ "minnow trap",
    
    str_detect(
      tolower(filename),
      "fyke"
    ) ~ "fyke",
    
    str_detect(
      tolower(filename),
      "weir"
    ) ~ "weir",
    
    str_detect(
      tolower(filename),
      "efishing|electrofishing"
    ) ~ "electrofishing",
    
    TRUE ~ NA_character_
  )
  
  method <- coalesce(
    str_to_lower(
      get_column(
        data,
        c(
          "trap",
          "type",
          "catch_method"
        )
      )
    ),
    default_method
  )
  
  
  # ---------------------------------------------------------------------------
  # SPECIMEN ID
  # ---------------------------------------------------------------------------
  
  minnow_code <- get_column(
    data,
    "minnow_trap"
  )
  
  id <- coalesce(
    get_column(
      data,
      c(
        "id",
        "id_code",
        "sample_id"
      )
    ),
    if_else(
      str_detect(minnow_code, "-"),
      minnow_code,
      NA_character_
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # SITE
  # ---------------------------------------------------------------------------
  
  site <- coalesce(
    
    toupper(
      get_column(
        data,
        c(
          "site",
          "site_code"
        )
      )
    ),
    
    toupper(
      str_extract(
        coalesce(id, minnow_code),
        "^[^-]+"
      )
    ),
    
    case_when(
      method == "fyke" ~ "FYKE",
      method == "weir" ~ "WEIR",
      TRUE ~ NA_character_
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # HABITAT
  # ---------------------------------------------------------------------------
  
  location <- coalesce(
    
    str_to_lower(
      get_column(
        data,
        c(
          "location",
          "habitat"
        )
      )
    ),
    
    case_when(
      
      str_detect(
        site,
        "^(CR|CC|UP|DN|WU|WD|FU|FD|FN)"
      ) ~ "creek",
      
      str_detect(
        site,
        "^(LK|OP)"
      ) ~ "lake",
      
      method %in% c(
        "fyke",
        "weir"
      ) ~ "creek",
      
      TRUE ~ NA_character_
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # DATE
  # ---------------------------------------------------------------------------
  
  date_raw <- get_column(
    data,
    "date"
  )
  
  date_source <- get_date(
    date_raw,
    month_first = filename == "August Weir Data_2022.xlsx"
  )
  
  date_source <- tibble(
    date_source = date_source
  ) %>%
    fill(date_source) %>%
    pull(date_source)
  
  month <- coalesce(
    month,
    month.name[
      as.integer(format(date_source, "%m"))
    ]
  )
  
  day <- as_number(
    get_column(data, "day")
  )
  
  has_day <- !is.na(day) &
    !is.na(month) &
    !is.na(year)
  
  date <- date_source
  
  date[has_day] <- as.Date(
    sprintf(
      "%04d-%02d-%02d",
      year,
      match(month[has_day], month.name),
      day[has_day]
    )
  )
  
  date_conflict <- !is.na(date_source) & (
    
    as.integer(
      format(date_source, "%Y")
    ) != year |
      
      (
        !is.na(month) &
          month.name[
            as.integer(
              format(date_source, "%m")
            )
          ] != month
      ) |
      
      (
        !is.na(date) &
          date_source != date
      )
  )
  
  date_conflict[is.na(date_conflict)] <- FALSE
  
  date[
    date_conflict & !has_day
  ] <- as.Date(NA)
  
  message(
    filename,
    ": ",
    nrow(data),
    " records"
  )
  
  
  # ---------------------------------------------------------------------------
  # RETURN STANDARDIZED DATA
  # ---------------------------------------------------------------------------
  
  tibble(
    
    id = clean_id(id),
    
    species = species,
    
    location = location,
    
    site = site,
    
    catch_method = method,
    
    number = as_number(
      get_column(data, "number")
    ),
    
    fl = as_number(
      get_column(
        data,
        c("fl", "fl_mm")
      )
    ),
    
    tl = as_number(
      get_column(
        data,
        c("tl", "tl_mm")
      )
    ),
    
    wt = as_number(
      get_column(
        data,
        c("wt", "wt_g")
      )
    ),
    
    month = month,
    
    year = year,
    
    date = as.character(date),
    
    date_raw = date_raw,
    
    date_conflict = date_conflict,
    
    notes = get_column(
      data,
      "notes"
    ),
    
    source_file = filename,
    
    source_row = source_row,
    
    data_source = "catch",
    
    priority = 2L
  )
}


# =============================================================================
# 7. IMPORT ALL CATCH RECORDS
# =============================================================================

catches <- bind_rows(
  lapply(
    catch_paths,
    read_catches
  )
) %>%
  mutate(
    across(
      c(fl, tl, wt),
      as.character
    )
  )

message(
  "Total catch records imported: ",
  nrow(catches)
)


# =============================================================================
# 8. READ ALL PROCESSED FISH
# =============================================================================

processed_raw <- read_csv(
  processed_file,
  col_types = cols(
    .default = col_character()
  ),
  na = c(
    "",
    "NA",
    "N/A",
    "NULL"
  ),
  name_repair = "unique_quiet",
  show_col_types = FALSE
)

# The source file contains two notes columns.
# Preserve both where available.

notes_columns <- grep(
  "^notes($|\\.\\.\\.)",
  names(processed_raw),
  value = TRUE
)

get_processed_notes <- function(index) {
  
  if (length(notes_columns) < index) {
    return(
      rep(
        NA_character_,
        nrow(processed_raw)
      )
    )
  }
  
  clean_text(
    processed_raw[[notes_columns[index]]]
  )
}


processed <- processed_raw %>%
  transmute(
    
    id = clean_id(sample_ID),
    
    species = str_to_lower(
      clean_text(species)
    ),
    
    site = str_extract(
      id,
      "^[A-Z]+[0-9]*"
    ),
    
    month_code = str_match(
      id,
      "^[^-]+-([MASO])[0-9]{2}(?:-|$)"
    )[, 2],
    
    year_code = str_match(
      id,
      "^[^-]+-(?:[MASO])?([0-9]{2})(?:-|$)"
    )[, 2],
    
    month = unname(
      c(
        M = "May",
        A = "August",
        S = "September",
        O = "October"
      )[month_code]
    ),
    
    year = 2000L +
      as.integer(year_code),
    
    number = 1,
    
    # Preserve text entries such as "No tail".
    fl = clean_text(fl),
    tl = clean_text(tl),
    wt = clean_text(wt),
    
    notes = get_processed_notes(1),
    
    notes_secondary = get_processed_notes(2),
    
    catch_method = case_when(
      
      str_detect(
        id,
        "^(FU|FD|FN)-"
      ) ~ "fyke",
      
      str_detect(
        id,
        "^(WU|WD)-"
      ) ~ "weir",
      
      TRUE ~ NA_character_
    ),
    
    source_file = "processed_all.csv",
    
    source_row = row_number() + 1L,
    
    data_source = "processed",
    
    priority = 1L
    
  ) %>%
  select(
    -month_code,
    -year_code
  )


# =============================================================================
# 9. APPLY CONFIRMED BASS METADATA
# =============================================================================

is_a_bass <- !is.na(processed$species) &
  processed$species == "smallmouth bass" &
  !is.na(processed$id) &
  str_starts(processed$id, "A")

unmatched_bass <- processed$id[
  is_a_bass &
    !processed$id %in% bass_sites$id
]

if (length(unmatched_bass) > 0) {
  
  stop(
    "A-coded bass missing from confirmed site table:\n",
    paste(
      unique(unmatched_bass),
      collapse = "\n"
    )
  )
}


# Match confirmed metadata directly to processed records.

bass_match <- match(
  processed$id,
  bass_sites$id
)

matched_bass <- !is.na(bass_match) &
  processed$species == "smallmouth bass"

processed$site[matched_bass] <-
  bass_sites$site[
    bass_match[matched_bass]
  ]

processed$catch_method[matched_bass] <-
  bass_sites$catch_method[
    bass_match[matched_bass]
  ]

processed$month[matched_bass] <-
  bass_sites$month[
    bass_match[matched_bass]
  ]

processed$year[matched_bass] <- 2022L

processed$date <- rep(
  NA_character_,
  nrow(processed)
)

processed$date[matched_bass] <-
  bass_sites$date[
    bass_match[matched_bass]
  ]


# =============================================================================
# 10. COMBINE PROCESSED FISH AND CATCH RECORDS
# =============================================================================

# Processed measurements have priority.
# Catch records supply missing information.
#
# Records without IDs receive unique row keys.
# This prevents unrelated bulk catches from being combined.

all_records <- bind_rows(
  processed,
  catches
) %>%
  mutate(
    
    record_key = if_else(
      is.na(id),
      paste0(
        "ROW_",
        row_number()
      ),
      paste0(
        "ID_",
        id
      )
    )
    
  ) %>%
  arrange(
    priority,
    source_file,
    source_row
  )


# =============================================================================
# 11. CHECK FOR CONFLICTING RECORDS
# =============================================================================

id_checks <- all_records %>%
  group_by(record_key) %>%
  filter(
    priority == min(priority)
  ) %>%
  summarise(
    
    id_conflict = any(
      vapply(
        
        pick(
          species,
          fl,
          tl,
          wt,
          site,
          month,
          year
        ),
        
        function(x) {
          
          n_distinct(
            x,
            na.rm = TRUE
          ) > 1
          
        },
        
        logical(1)
      )
    ),
    
    .groups = "drop"
  )


# =============================================================================
# 12. CREATE FINAL MERGED DATASET
# =============================================================================

fish_data <- all_records %>%
  group_by(record_key) %>%
  summarise(
    across(
      everything(),
      first_known
    ),
    .groups = "drop"
  ) %>%
  left_join(
    id_checks,
    by = "record_key"
  ) %>%
  left_join(
    select(
      bass_sites,
      id,
      capture_site
    ),
    by = "id"
  ) %>%
  mutate(
    
    area = as.character(assign_area(if_else(site %in% "SMB1", id, site))),
    
    mapped_habitat = case_when(
      
      str_starts(
        area,
        "CR_"
      ) ~ "creek",
      
      str_starts(
        area,
        "LK_"
      ) ~ "lake",
      
      TRUE ~ NA_character_
    ),
    
    # Confirmed fyke bass belong to the creek.
    
    location = if_else(
      id %in% c(
        "A-22-04",
        "A-22-07"
      ),
      "creek",
      location
    ),
    
    location_conflict =
      !is.na(location) &
      !is.na(mapped_habitat) &
      location != mapped_habitat,
    
    location = coalesce(
      location,
      mapped_habitat
    ),
    
    season = case_when(
      
      month %in% month.name[3:5] ~ "Spring",
      
      month %in% month.name[6:8] ~ "Summer",
      
      month %in% month.name[9:11] ~ "Fall",
      
      month %in% month.name[
        c(12, 1, 2)
      ] ~ "Winter",
      
      TRUE ~ NA_character_
    ),
    
    # Preserve unresolved measurement.
    # A-22-03 has recorded TL = 114.
    # TL = 414 is unverified.
    
    length_unresolved =
      id %in% "A-22-03",
    
    proposed_tl = if_else(
      length_unresolved,
      414,
      NA_real_
    )
    
  ) %>%
  select(
    
    id,
    species,
    location,
    site,
    area,
    catch_method,
    number,
    fl,
    tl,
    wt,
    month,
    season,
    year,
    date,
    notes,
    notes_secondary,
    capture_site,
    id_conflict,
    location_conflict,
    length_unresolved,
    proposed_tl,
    date_raw,
    date_conflict,
    source_file,
    source_row,
    data_source
  )


# =============================================================================
# 13. VALIDATE MERGED DATA
# =============================================================================

bass_check <- fish_data %>%
  filter(
    species == "smallmouth bass",
    id %in% bass_sites$id
  )

# Check for duplicated specimen codes.

duplicate_ids <- fish_data$id[
  !is.na(fish_data$id) &
    duplicated(fish_data$id)
]

if (length(duplicate_ids) > 0) {
  
  stop(
    "Duplicate specimen codes remain:\n",
    paste(
      unique(duplicate_ids),
      collapse = "\n"
    )
  )
}

# Confirm all processed fish with IDs are represented.

stopifnot(
  
  all(
    processed$id[
      !is.na(processed$id)
    ] %in% fish_data$id
  ),
  
  all(
    bass_sites$id %in% bass_check$id
  ),
  
  all(
    !is.na(bass_check$area)
  ),
  
  nrow(bass_check) == 21
)


# =============================================================================
# 14. SUMMARY COUNTS
# =============================================================================

summary_counts <- tibble(
  
  check = c(
    
    "Catch records read",
    
    "Processed records read",
    
    "Merged rows saved",
    
    "Repeated-code rows combined",
    
    "Duplicate codes remaining",
    
    "A-coded bass included"
  ),
  
  count = c(
    
    nrow(catches),
    
    nrow(processed),
    
    nrow(fish_data),
    
    nrow(all_records) -
      nrow(fish_data),
    
    length(duplicate_ids),
    
    nrow(bass_check)
  )
)


# =============================================================================
# 15. RECORDS REQUIRING REVIEW
# =============================================================================

review <- fish_data %>%
  filter(
    
    coalesce(id_conflict, FALSE) |
      
      coalesce(location_conflict, FALSE) |
      
      coalesce(length_unresolved, FALSE) |
      
      is.na(area) |
      
      is.na(location) |
      
      is.na(season)
  )


# Preserve duplicate source records for inspection.

duplicates <- all_records %>%
  group_by(record_key) %>%
  filter(
    n() > 1
  ) %>%
  ungroup()


# =============================================================================
# 16. SAVE OUTPUTS
# =============================================================================

dir.create(
  checks_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  fish_data,
  output_file,
  na = ""
)

write_csv(
  review,
  file.path(
    checks_dir,
    "records_to_review.csv"
  ),
  na = ""
)

write_csv(
  duplicates,
  file.path(
    checks_dir,
    "duplicate_source_records.csv"
  ),
  na = ""
)

write_csv(
  summary_counts,
  file.path(
    checks_dir,
    "merge_counts.csv"
  )
)


# =============================================================================
# 17. FINAL CHECKS
# =============================================================================

print(summary_counts)

message(
  "Merged CSV saved to: ",
  output_file
)

message(
  "Review files saved to: ",
  checks_dir
)

message(
  "Confirmed A-coded bass included: ",
  nrow(bass_check)
)

message(
  "Number of records requiring review: ",
  nrow(review)
)


# =============================================================================
# 18. OPTIONAL: INSPECT A-CODED BASS
# =============================================================================

bass_check %>%
  arrange(id) %>%
  select(
    id,
    species,
    location,
    site,
    area,
    catch_method,
    fl,
    tl,
    wt,
    month,
    season,
    year,
    date
  ) %>%
  print(n = Inf)


# =============================================================================
# END
# =============================================================================