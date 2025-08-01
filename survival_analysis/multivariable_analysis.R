# ==============================================================================
# Combined Multivariable Analysis and Forest Plot
# "An AI based Pathology Model to Predict Docetaxel Benefit in Prostate Cancer"
# ==============================================================================

# Load required libraries
packages <- c("survival", "forestplot", "dplyr", "ggplot2", "tidyr")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# MULTIVARIABLE ANALYSIS FUNCTIONS
# ==============================================================================

prepare_chaarted_data <- function(filepath, endpoint = "OS") {
  """Prepare CHAARTED data with proper factor coding"""
  data <- read.csv(filepath) %>%
    mutate(
      risk_group = factor(risk_group, levels = c("Low Risk", "High Risk")),
      treatment = factor(treatment, levels = c("ADT", "ADT+CT")),
      dz_extent = factor(dz_extent, levels = c("Low", "High")),
      visc_dz = case_when(visc_dz == 0 ~ 0, visc_dz == 1 ~ 1, TRUE ~ NA_real_),
      visc_dz = factor(visc_dz, levels = c(0, 1), labels = c("No", "Yes")),
      metachronous = case_when(
        prior_local_tx == 0 ~ 1, 
        prior_local_tx %in% c(1, 2) ~ 0, 
        TRUE ~ NA_real_
      ),
      metachronous = factor(metachronous, levels = c(0, 1), 
                            labels = c("Synchronous", "Metachronous"))
    )
  
  if(endpoint == "OS") {
    data <- data %>% rename(time = os, event = dead)
  } else {
    data <- data %>% rename(time = TT_CRPC, event = CRPC)
  }
  
  return(data)
}

prepare_rtog_data <- function(filepath) {
  """Prepare RTOG data with T-stage binarization"""
  data <- read.csv(filepath) %>%
    mutate(
      risk_group = factor(risk_group),
      chemo = factor(chemo),
      t_stage_bin = case_when(
        t_stage < 8 ~ 0,
        t_stage >= 8 & t_stage <= 11 ~ 1,
        TRUE ~ NA_real_
      ),
      t_stage_bin = factor(t_stage_bin, levels = c(0, 1), labels = c("<T3", "≥T3"))
    )
  return(data)
}

perform_multivariable_analysis <- function(data, trial_name, endpoint = "OS") {
  """Perform multivariable Cox regression with treatment-APIC interaction"""
  
  cat("\n=== MULTIVARIABLE ANALYSIS:", trial_name, "===\n")
  
  if(trial_name == "CHAARTED") {
    # CHAARTED model
    cox_model <- coxph(
      Surv(time, event) ~ treatment * risk_group + dz_extent + 
        metachronous + gleason + baseline_PSA,
      data = data
    )
  } else {
    # RTOG model  
    cox_model <- coxph(
      Surv(time, event) ~ chemo * risk_group + gleason + psa + t_stage_bin,
      data = data
    )
  }
  
  summary_cox <- summary(cox_model)
  
  result <- data.frame(
    Trial = trial_name,
    Variable = rownames(summary_cox$coefficients),
    HR = summary_cox$coefficients[, "exp(coef)"],
    Lower_CI = summary_cox$conf.int[, "lower .95"],
    Upper_CI = summary_cox$conf.int[, "upper .95"],
    P_Value = summary_cox$coefficients[, "Pr(>|z|)"]
  )
  
  cat("Results:\n")
  print(result)
  
  return(result)
}

# ==============================================================================
# FOREST PLOT FUNCTIONS
# ==============================================================================

prepare_forest_data <- function(df, dataset_name) {
  """Prepare data for forest plot analysis"""
  df$risk_group <- ifelse(df$risk_group == 0, "APIC-Negative", "APIC-Positive")
  df$trial <- dataset_name
  df$chemo <- ifelse(df$chemo == 0, "SOC", "Docetaxel")
  df$chemo <- factor(df$chemo, levels = c("SOC", "Docetaxel"))
  df$risk_group <- factor(df$risk_group)
  df$Docetaxel <- ifelse(df$chemo == "Docetaxel", 1, 0)
  return(df)
}

analyze_for_forest <- function(df, dataset_name) {
  """Analyze dataset for forest plot components"""
  
  risk_groups <- levels(df$risk_group)
  
  # Fit Cox models for each risk group
  plot_data_list <- lapply(risk_groups, function(group) {
    group_data <- filter(df, risk_group == group)
    cox_model <- coxph(Surv(time, event) ~ Docetaxel, data = group_data)
    summary_cox <- summary(cox_model)
    
    data.frame(
      Trial = dataset_name,
      Category = group,
      HR = summary_cox$coef[,"exp(coef)"],
      Lower_CI = summary_cox$conf.int[,"lower .95"],
      Upper_CI = summary_cox$conf.int[,"upper .95"],
      P_Value = summary_cox$coef[,"Pr(>|z|)"]
    )
  })
  
  plot_data <- bind_rows(plot_data_list)
  
  # Add event counts
  incidence_data <- df %>%
    group_by(risk_group, chemo) %>%
    summarise(Events = sum(event), N = n(), .groups = 'drop') %>%
    mutate(Incidence_N = paste0(Events, "/", N)) %>%
    select(risk_group, chemo, Incidence_N) %>%
    pivot_wider(names_from = chemo, values_from = Incidence_N)
  
  plot_data <- left_join(plot_data, incidence_data, by = c("Category" = "risk_group"))
  
  # Calculate survival benefits
  time_point <- ifelse(dataset_name == "RTOG 0521", 10, 60)
  
  surv_data <- lapply(risk_groups, function(group) {
    group_data <- filter(df, risk_group == group)
    surv_fit <- survfit(Surv(time, event) ~ Docetaxel, data = group_data)
    summary_surv <- summary(surv_fit, times = time_point, extend = TRUE)
    
    data.frame(
      risk_group = group,
      Docetaxel = gsub("Docetaxel=", "", summary_surv$strata),
      surv_prob = summary_surv$surv,
      surv_se = summary_surv$std.err
    )
  }) %>% bind_rows()
  
  surv_data$Docetaxel <- ifelse(surv_data$Docetaxel == "0", "No", "Yes")
  
  surv_wide <- surv_data %>%
    select(risk_group, Docetaxel, surv_prob, surv_se) %>%
    pivot_wider(names_from = Docetaxel, values_from = c(surv_prob, surv_se)) %>%
    mutate(
      Abs_Benefit = (`surv_prob_Yes` - `surv_prob_No`) * 100,
      SE_Diff = sqrt(`surv_se_Yes`^2 + `surv_se_No`^2) * 100,
      Abs_Benefit_Lower = Abs_Benefit - 1.96 * SE_Diff,
      Abs_Benefit_Upper = Abs_Benefit + 1.96 * SE_Diff,
      Abs_Benefit_CI = sprintf("%.1f%% (%.1f%% - %.1f%%)", 
                               Abs_Benefit, Abs_Benefit_Lower, Abs_Benefit_Upper)
    )
  
  plot_data <- left_join(plot_data, surv_wide %>% select(risk_group, Abs_Benefit_CI), 
                         by = c("Category" = "risk_group"))
  
  # Calculate interaction p-value
  interaction_model <- coxph(Surv(time, event) ~ Docetaxel * risk_group, data = df)
  summary_interaction <- summary(interaction_model)
  interaction_term <- rownames(summary_interaction$coefficients)[
    grepl("Docetaxel.*risk_group", rownames(summary_interaction$coefficients))
  ]
  p_interaction <- summary_interaction$coefficients[interaction_term, "Pr(>|z|)"]
  
  plot_data$Interaction_P_Value <- ""
  apic_positive_idx <- which(plot_data$Category == "APIC-Positive")
  if(length(apic_positive_idx) > 0) {
    plot_data$Interaction_P_Value[apic_positive_idx] <- sprintf("%.3f", p_interaction)
  }
  
  return(plot_data)
}

create_forest_plot <- function(chaarted_data, rtog_data) {
  """Create combined forest plot"""
  
  # Analyze both datasets
  plot_data_chaarted <- analyze_for_forest(chaarted_data, "CHAARTED")
  plot_data_rtog <- analyze_for_forest(rtog_data, "RTOG 0521")
  
  # Combine and order results
  plot_data <- bind_rows(plot_data_chaarted, plot_data_rtog) %>%
    mutate(Trial_Category = paste(Trial, Category)) %>%
    arrange(Trial, Category)
  
  plot_data$Trial_Category <- factor(plot_data$Trial_Category,
                                     levels = c("CHAARTED APIC-Negative", 
                                                "CHAARTED APIC-Positive",
                                                "RTOG 0521 APIC-Negative", 
                                                "RTOG 0521 APIC-Positive"))
  
  # Format p-values
  format_p_values <- function(p_vals) {
    ifelse(p_vals < 0.01, "< 0.01**",
           ifelse(p_vals < 0.05, "< 0.05*", sprintf("%.3f", p_vals)))
  }
  
  formatted_p_values <- format_p_values(plot_data$P_Value)
  
  numeric_interaction_p <- suppressWarnings(as.numeric(plot_data$Interaction_P_Value))
  formatted_interaction_p <- ifelse(
    is.na(numeric_interaction_p), "", format_p_values(numeric_interaction_p)
  )
  
  # Create trial names (show only once per trial)
  trial_names <- as.character(plot_data$Trial)
  for(trial in unique(trial_names)) {
    first_occurrence <- which(trial_names == trial)[1]
    trial_names[-first_occurrence][trial_names[-first_occurrence] == trial] <- ""
  }
  
  # Create table for forest plot
  tabletext <- cbind(
    c("Clinical Trial", trial_names),
    c("Group", sub(".* (APIC-.*)", "\\1", plot_data$Trial_Category)),
    c("Docetaxel\nEvents/N", plot_data$Docetaxel),
    c("SOC\nEvents/N", plot_data$SOC),
    c("HR (95% CI)", sprintf("%.2f (%.2f - %.2f)", 
                             plot_data$HR, plot_data$Lower_CI, plot_data$Upper_CI)),
    c("P value", formatted_p_values),
    c("Interaction\nP value", formatted_interaction_p),
    c("Benefit of Docetaxel (95% CI)\u2020", plot_data$Abs_Benefit_CI)
  )
  
  # Generate forest plot
  forestplot(
    labeltext = tabletext,
    mean = c(NA, plot_data$HR),
    lower = c(NA, plot_data$Lower_CI),
    upper = c(NA, plot_data$Upper_CI),
    zero = 1,
    vertices = TRUE,
    boxsize = 0.5,
    col = fpColors(box = "#000000", line = "#000000"),
    xlab = "Favors Docetaxel           Favors SOC",
    title = "Effect of Docetaxel on Overall Survival",
    xticks = seq(0, 2.5, by = 0.5),
    graph.pos = ncol(tabletext),
    graphwidth = unit(6, "cm")
  )
}

# ==============================================================================
# MAIN ANALYSIS PIPELINE
# ==============================================================================

run_complete_analysis <- function(chaarted_os_path, chaarted_crpc_path, rtog_path,
                                  chaarted_forest_path, rtog_forest_path) {
  """Complete analysis pipeline for APIC study"""
  
  cat("=== APIC COMPLETE ANALYSIS PIPELINE ===\n")
  
  # Multivariable Analyses
  cat("\n1. MULTIVARIABLE ANALYSES\n")
  
  # CHAARTED
  chaarted_os_data <- prepare_chaarted_data(chaarted_os_path, "OS")
  chaarted_crpc_data <- prepare_chaarted_data(chaarted_crpc_path, "CRPC")
  
  chaarted_os_results <- perform_multivariable_analysis(chaarted_os_data, "CHAARTED", "OS")
  chaarted_crpc_results <- perform_multivariable_analysis(chaarted_crpc_data, "CHAARTED", "CRPC")
  
  # RTOG
  rtog_data <- prepare_rtog_data(rtog_path)
  rtog_results <- perform_multivariable_analysis(rtog_data, "RTOG 0521", "OS")
  
  # Forest Plot
  cat("\n2. FOREST PLOT GENERATION\n")
  
  chaarted_forest_data <- prepare_forest_data(read.csv(chaarted_forest_path), "CHAARTED")
  rtog_forest_data <- prepare_forest_data(read.csv(rtog_forest_path), "RTOG 0521")
  
  create_forest_plot(chaarted_forest_data, rtog_forest_data)
  
  # Return results
  return(list(
    chaarted_os = chaarted_os_results,
    chaarted_crpc = chaarted_crpc_results,
    rtog_os = rtog_results
  ))
}