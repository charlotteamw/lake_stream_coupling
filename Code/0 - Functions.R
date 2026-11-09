# Functions used in this repository

# Version: 2026-03-14

#Calculate Standard Error
standard_error <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

#Function to make sure proportional data using beta-family responses are not exactly 0s or 1s
adjust_beta_response <- function(x) {
  n_obs <- sum(!is.na(x))
  if (n_obs <= 1) return(x)
  (x * (n_obs - 1) + 0.5) / n_obs
}

#Generate emmeans pairwise results for isotopes
print_emm_pairs_iso <- function(model, response_name) {
  cat("\nEstimated marginal means:\n")
  print(emmeans(model, ~ month * location, type = "response"))
  
  cat("\nPairwise comparisons for month within location:\n")
  print(pairs(
    emmeans(model, ~ month | location, type = "response"),
    adjust = "tukey"
  ))
  
  cat("\nPairwise comparisons for location within month:\n")
  print(pairs(
    emmeans(model, ~ location | month, type = "response"),
    adjust = "tukey"
  ))
}

#Generate emmeans pairwise results for Telemetry
print_emm_pairs_tel <- function(model) {
  cat("\nEstimated marginal means:\n")
  print(emmeans(model, ~ Season * release_location, type = "response"))
  
  cat("\nPairwise comparisons for Season within release location:\n")
  print(pairs(
    emmeans(model, ~ Season | release_location, type = "response"),
    adjust = "tukey"
  ))
  
  cat("\nPairwise comparisons for release location within Season:\n")
  print(pairs(
    emmeans(model, ~ release_location | Season, type = "response"),
    adjust = "tukey"
  ))
}


#Calculate lake-derived carbon and trophic position (uses d13C)
calculate_prop_tp <- function(iso_data, chosen_species, chosen_tissue) {
  months  <- c("may", "august", "october")
  results <- list()
  
  cr_baseline <- iso_data %>% filter(organism == "mayfly" | organism == "mussel", location == "creek")
  lk_baseline <- iso_data %>% filter(organism == "mayfly" | organism == "mussel", location == "lake")
  
  cr_mean_dC <- mean(cr_baseline$d13C, na.rm = TRUE)
  cr_mean_dN <- mean(cr_baseline$d15N,      na.rm = TRUE)
  lk_mean_dC <- mean(lk_baseline$d13C, na.rm = TRUE)
  lk_mean_dN <- mean(lk_baseline$d15N,      na.rm = TRUE)
  
  for (month_i in months) {
    df <- iso_data %>%
      filter(
        organism %in% chosen_species,
        tissue == chosen_tissue,
        month  == month_i
      )
    
    lake_carbon <- (df$d13C - cr_mean_dC) / (lk_mean_dC - cr_mean_dC)
    lake_carbon <- pmin(pmax(lake_carbon, 0.001), 0.999)
    
    trophic_position <- 2 + (
      (1-lake_carbon)       * ((df$d15N - cr_mean_dN) / 3.4) +
      lake_carbon * ((df$d15N - lk_mean_dN) / 3.4)
    )
    
    df$lake_carbon      <- lake_carbon
    df$trophic_position <- trophic_position
    results[[month_i]]  <- df
  }
  
  bind_rows(results)
}


create_time_series <- function(data, id_time_value) {
  fish_data <- data %>% dplyr::filter(id_time == id_time_value)
  
  if (nrow(fish_data) == 0 ||
      is.na(fish_data$tag_on_date[1]) ||
      is.na(fish_data$tag_off_date[1])) {
    warning(paste("Skipping id_time:", id_time_value, "due to missing or invalid dates"))
    return(NULL)
  }
  
  start_time <- fish_data$tag_on_date[1] + hours(12)
  end_time   <- fish_data$tag_off_date[1] + hours(12)
  
  if (!is.finite(start_time) || !is.finite(end_time) || start_time > end_time) {
    warning(paste("Skipping id_time:", id_time_value, "due to invalid time range"))
    return(NULL)
  }
  
  time_series <- tibble(
    id_time          = id_time_value,
    date_time        = seq(from = start_time, to = end_time, by = "1 min"),
    release_location = fish_data$release_location[1]
  )
  
  detection_data <- fish_data %>%
    dplyr::mutate(detection_minute = floor_date(detection_timestamp, "minute")) %>%
    dplyr::select(detection_minute, location) %>%
    dplyr::distinct(detection_minute, .keep_all = TRUE)
  
  time_series %>%
    left_join(detection_data, by = c("date_time" = "detection_minute")) %>%
    tidyr::fill(location, .direction = "down")
}


enriched_time_series <- function(data, id_time) {
  fish_data <- data %>% filter(id_time == !!id_time)
  
  if (nrow(fish_data) == 0 || is.na(fish_data$tag_on_date[1]) || is.na(fish_data$tag_off_date[1])) {
    warning(paste("Skipping id_time:", id_time, "due to missing or invalid dates"))
    return(NULL)
  }
  
  start_time <- fish_data$tag_on_date[1] + hours(12)
  end_time   <- fish_data$tag_off_date[1] + hours(12)
  
  if (!is.finite(start_time) || !is.finite(end_time) || start_time > end_time) {
    warning(paste("Skipping id_time:", id_time, "due to invalid time range"))
    return(NULL)
  }
  
  release_loc <- fish_data$release_location[1]
  
  time_series <- tibble(
    id_time          = id_time,
    date_time        = seq(from = start_time, to = end_time, by = "1 min"),
    release_location = release_loc
  )
  
  detection_data <- fish_data %>%
    dplyr::mutate(detection_minute = floor_date(detection_timestamp, "minute")) %>%
    dplyr::select(detection_minute, location) %>%
    dplyr::distinct(detection_minute, .keep_all = TRUE)
  
  enriched_series <- time_series %>%
    left_join(detection_data, by = c("date_time" = "detection_minute")) %>%
    fill(location, .direction = "down")
  
  return(enriched_series)
}


# =============================================================================
# MODEL RESULTS AND SUPPLEMENTARY TABLE EXPORT
# =============================================================================

# Collect results from a fitted glmmTMB model.
#
# Arguments:
# model   = fitted model
# data    = dataframe used to fit the model
# specs   = emmeans formula, e.g. ~ season * location
# species = optional species column for sample-size reporting

model_entry <- function(model, data, specs, species = NULL) {
  
  # ---------------------------------------------------------------------------
  # Model information
  # ---------------------------------------------------------------------------
  
  model_info <- data.frame(
    formula = paste(
      deparse(formula(model)),
      collapse = " "
    ),
    family = family(model)$family,
    link = family(model)$link,
    n = nobs(model),
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = as.numeric(logLik(model)),
    convergence = model$fit$convergence,
    positive_definite_Hessian = isTRUE(model$sdr$pdHess)
  )
  
  # ---------------------------------------------------------------------------
  # Fixed-effect estimates
  # ---------------------------------------------------------------------------
  
  fixed_effects <- as.data.frame(
    summary(model)$coefficients$cond
  )
  
  fixed_effects$term <- rownames(fixed_effects)
  rownames(fixed_effects) <- NULL
  
  fixed_effects <- fixed_effects %>%
    dplyr::relocate(term)
  
  # ---------------------------------------------------------------------------
  # Random-effect variances
  # ---------------------------------------------------------------------------
  
  random_effects <- as.data.frame(
    VarCorr(model)
  )
  
  # ---------------------------------------------------------------------------
  # Wald tests
  #
  # Type III for models containing an interaction.
  # Type II for additive models.
  # ---------------------------------------------------------------------------
  
  model_terms <- attr(
    terms(model),
    "term.labels"
  )
  
  has_interaction <- any(
    grepl(":", model_terms, fixed = TRUE)
  )
  
  anova_type <- if (has_interaction) 3 else 2
  
  anova_results <- as.data.frame(
    car::Anova(
      model,
      type = anova_type
    )
  )
  
  anova_results$term <- rownames(anova_results)
  
  rownames(anova_results) <- NULL
  
  anova_results <- anova_results %>%
    dplyr::relocate(term)
  
  # ---------------------------------------------------------------------------
  # Estimated marginal means
  # ---------------------------------------------------------------------------
  
  emm <- emmeans::emmeans(
    model,
    specs
  )
  
  emm_results <- as.data.frame(
    summary(
      emm,
      type = "response",
      infer = c(TRUE, FALSE)
    )
  )
  
  # ---------------------------------------------------------------------------
  # Pairwise comparisons
  # ---------------------------------------------------------------------------
  
  pairwise_results <- as.data.frame(
    summary(
      pairs(
        emm,
        adjust = "tukey"
      ),
      type = "response",
      infer = c(FALSE, TRUE)
    )
  )
  
  # ---------------------------------------------------------------------------
  # Sample sizes
  # ---------------------------------------------------------------------------
  
  if (!is.null(species) && species %in% names(data)) {
    
    sample_sizes <- data %>%
      dplyr::count(
        .data[[species]],
        name = "n"
      )
    
  } else {
    
    sample_sizes <- data.frame(
      n = nrow(data)
    )
  }
  
  # ---------------------------------------------------------------------------
  # Return results
  # ---------------------------------------------------------------------------
  
  list(
    model_info = model_info,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    anova = anova_results,
    anova_type = anova_type,
    emmeans = emm_results,
    pairwise = pairwise_results,
    sample_sizes = sample_sizes
  )
}


# =============================================================================
# SAVE MODEL RESULTS
# =============================================================================

save_analysis <- function(entries, analysis_name, file_path) {
  
  # Create output directory.
  
  output_dir <- file.path(
    file_path,
    "Results",
    analysis_name
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Export each model separately.
  
  for (model_name in names(entries)) {
    
    entry <- entries[[model_name]]
    
    # Safe filename.
    safe_name <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      model_name
    )
    
    # Model information
    readr::write_csv(
      entry$model_info,
      file.path(
        output_dir,
        paste0(safe_name, "_model_info.csv")
      )
    )
    
    # Fixed effects
    readr::write_csv(
      entry$fixed_effects,
      file.path(
        output_dir,
        paste0(safe_name, "_fixed_effects.csv")
      )
    )
    
    # Random effects
    readr::write_csv(
      entry$random_effects,
      file.path(
        output_dir,
        paste0(safe_name, "_random_effects.csv")
      )
    )
    
    # Wald tests
    readr::write_csv(
      entry$anova,
      file.path(
        output_dir,
        paste0(safe_name, "_anova.csv")
      )
    )
    
    # Estimated marginal means
    readr::write_csv(
      entry$emmeans,
      file.path(
        output_dir,
        paste0(safe_name, "_emmeans.csv")
      )
    )
    
    # Pairwise comparisons
    readr::write_csv(
      entry$pairwise,
      file.path(
        output_dir,
        paste0(safe_name, "_pairwise.csv")
      )
    )
    
    # Sample sizes
    readr::write_csv(
      entry$sample_sizes,
      file.path(
        output_dir,
        paste0(safe_name, "_sample_sizes.csv")
      )
    )
  }
  
  message(
    "Analysis results saved to: ",
    output_dir
  )
  
  invisible(output_dir)
}
