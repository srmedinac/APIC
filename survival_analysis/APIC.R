# ==============================================================================
# APIC Development and Validation
# Part of: "An AI based Pathology Model to Predict Docetaxel Benefit in Prostate Cancer"
# ==============================================================================

# Load required libraries
packages <- c("survival", "survminer", "dplyr", "caret", "ggplot2", "gridExtra")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Define feature sets for each trial
CHAARTED_FEATURES <- c(
  'Area.Energy.var', 'Area.InvDiffMom.Skewness', 'MinorAxisLength.Energy.Prcnt90',
  'Area.DiffAvg.Prcnt10', 'X341', 'X51'
)

RTOG_FEATURES <- c(
  'Area.Entropy.IqntlRange', 'Area.InvVar.Median', 'Area.ClusterShade.Kurtosis', 
  'X51', 'X341', 'X336', 'X0'
)

# Function to normalize features using min-max scaling
normalize_minmax <- function(x) {
  """
  Apply min-max normalization to feature vector
  
  Args:
    x: Numeric vector to normalize
  
  Returns:
    Normalized vector with values between 0 and 1
  """
  if (max(x) == min(x)) return(x)  # Handle case where all values are identical
  return((x - min(x)) / (max(x) - min(x)))
}

# Function to prepare and split data
prepare_data <- function(filepath, features, control_arm_filter, treatment_arm_filter) {
  """
  Prepare data for model development and validation
  
  Args:
    filepath: Path to CSV file with clinical and feature data
    features: Vector of feature names to use for modeling
    control_arm_filter: Filter expression for control arm (e.g., ASSIGNED_TX_ARM == "B")
    treatment_arm_filter: Filter expression for treatment arm (e.g., ASSIGNED_TX_ARM == "A")
  
  Returns:
    List containing training, validation, and treatment arm datasets
  """
  
  data <- read.csv(filepath)
  
  # Filter control arm data
  control_data <- data %>%
    filter(has_features == TRUE, !!enquo(control_arm_filter)) %>%
    select(group_uid, TT_CRPC, CRPC, os, dead, all_of(features)) %>%
    na.omit()
  
  # Filter treatment arm data
  treatment_data <- data %>%
    filter(has_features == TRUE, !!enquo(treatment_arm_filter)) %>%
    select(group_uid, TT_CRPC, CRPC, os, dead, all_of(features)) %>%
    na.omit()
  
  # Create train/validation split from control arm (50/50 split)
  set.seed(1)
  train_index <- createDataPartition(control_data$CRPC, p = 0.5, list = FALSE)
  
  train_data <- control_data[train_index, ]
  validation_data <- control_data[-train_index, ]
  
  # Normalize features using training data statistics
  for (feat in features) {
    feat_min <- min(train_data[[feat]])
    feat_max <- max(train_data[[feat]])
    
    # Normalize training data
    train_data[[feat]] <- normalize_minmax(train_data[[feat]])
    
    # Apply same transformation to validation and treatment data
    validation_data[[feat]] <- (validation_data[[feat]] - feat_min) / (feat_max - feat_min)
    treatment_data[[feat]] <- (treatment_data[[feat]] - feat_min) / (feat_max - feat_min)
  }
  
  return(list(
    train = train_data,
    validation = validation_data,
    treatment = treatment_data,
    features = features
  ))
}

# Function to develop APIC model
develop_apic_model <- function(train_data, features, endpoint = "CRPC") {
  """
  Develop APIC model using Cox regression on training data
  
  Args:
    train_data: Training dataset with features and outcomes
    features: Vector of feature names to include in model
    endpoint: Endpoint type ("CRPC" or "OS")
  
  Returns:
    List containing fitted Cox model and risk threshold
  """
  
  # Create Cox model formula
  if (endpoint == "CRPC") {
    cox_formula <- as.formula(paste("Surv(TT_CRPC, CRPC) ~", paste(features, collapse = " + ")))
  } else {
    cox_formula <- as.formula(paste("Surv(os, dead) ~", paste(features, collapse = " + ")))
  }
  
  # Fit Cox model
  cox_model <- coxph(cox_formula, data = train_data)
  
  # Calculate risk scores on training data
  train_risks <- predict(cox_model, type = "risk", newdata = train_data)
  
  # Set threshold at 33rd percentile for CHAARTED, 50th percentile for RTOG
  threshold <- quantile(train_risks, 0.33)  # Adjust as needed per trial
  
  cat("Model Summary:\n")
  print(summary(cox_model))
  cat("\nRisk threshold (33rd percentile):", threshold, "\n")
  
  return(list(
    model = cox_model,
    threshold = threshold
  ))
}

# Function to validate APIC model
validate_apic_model <- function(model_object, validation_data, treatment_data, endpoint = "CRPC") {
  """
  Validate APIC model on independent datasets
  
  Args:
    model_object: List containing Cox model and threshold from develop_apic_model()
    validation_data: Validation dataset (control arm)
    treatment_data: Treatment arm dataset
    endpoint: Endpoint type ("CRPC" or "OS")
  
  Returns:
    List containing validation results and processed datasets
  """
  
  cox_model <- model_object$model
  threshold <- model_object$threshold
  
  # Calculate risk scores for validation and treatment datasets
  validation_risks <- predict(cox_model, type = "risk", newdata = validation_data)
  treatment_risks <- predict(cox_model, type = "risk", newdata = treatment_data)
  
  # Assign risk groups based on threshold
  validation_data$risk_score <- validation_risks
  validation_data$risk_group <- factor(
    ifelse(validation_risks > threshold, "High Risk", "Low Risk"),
    levels = c("Low Risk", "High Risk")
  )
  
  treatment_data$risk_score <- treatment_risks
  treatment_data$risk_group <- factor(
    ifelse(treatment_risks > threshold, "High Risk", "Low Risk"),
    levels = c("Low Risk", "High Risk")
  )
  
  # Validate prognostic value in control arm
  if (endpoint == "CRPC") {
    cox_validation <- coxph(Surv(TT_CRPC, CRPC) ~ risk_group, data = validation_data)
  } else {
    cox_validation <- coxph(Surv(os, dead) ~ risk_group, data = validation_data)
  }
  
  cat("\nValidation Results (Control Arm):\n")
  print(summary(cox_validation))
  
  # Test predictive value (treatment interaction)
  combined_data <- bind_rows(
    validation_data %>% mutate(treatment = "Control"),
    treatment_data %>% mutate(treatment = "Treatment")
  )
  
  if (endpoint == "CRPC") {
    cox_interaction <- coxph(Surv(TT_CRPC, CRPC) ~ treatment * risk_group, data = combined_data)
  } else {
    cox_interaction <- coxph(Surv(os, dead) ~ treatment * risk_group, data = combined_data)
  }
  
  cat("\nPredictive Value (Treatment Interaction):\n")
  print(summary(cox_interaction))
  
  return(list(
    validation_data = validation_data,
    treatment_data = treatment_data,
    combined_data = combined_data,
    cox_validation = cox_validation,
    cox_interaction = cox_interaction
  ))
}

# Function to create survival plots
create_survival_plots <- function(data, title, endpoint = "CRPC") {
  """
  Create Kaplan-Meier survival plots
  
  Args:
    data: Dataset with risk_group, time, and event columns
    title: Plot title
    endpoint: Endpoint type ("CRPC" or "OS")
  
  Returns:
    ggsurvplot object
  """
  
  if (endpoint == "CRPC") {
    fit <- survfit(Surv(TT_CRPC, CRPC) ~ risk_group, data = data)
    ylab <- "CRPC-Free Survival Probability"
  } else {
    fit <- survfit(Surv(os, dead) ~ risk_group, data = data)
    ylab <- "Overall Survival Probability"
  }
  
  # Calculate statistics for annotation
  cox_temp <- coxph(Surv(TT_CRPC, CRPC) ~ risk_group, data = data)
  cox_summary <- summary(cox_temp)
  
  hr <- cox_summary$conf.int[1, "exp(coef)"]
  hr_lower <- cox_summary$conf.int[1, "lower .95"]
  hr_upper <- cox_summary$conf.int[1, "upper .95"]
  hr_text <- sprintf("HR: %.2f (95%% CI: %.2f-%.2f)", hr, hr_lower, hr_upper)
  
  c_index <- cox_summary$concordance[1]
  c_index_text <- sprintf("C-index: %.3f", c_index)
  
  # Create plot
  p <- ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    title = title,
    conf.int.style = "ribbon",
    xlab = "Time (Months)",
    ylab = ylab,
    legend.title = "APIC Group",
    legend.labs = c("APIC-Negative", "APIC-Positive"),
    palette = c("#00BFC4", "#F8766D"),
    surv.median.line = "hv",
    risk.table.height = 0.2,
    font.main = c(18, "bold"),
    font.x = c(16),
    font.y = c(16),
    font.legend = c(16),
    font.tickslab = c(16),
    risk.table.fontsize = 6,
    pval.size = 6
  )
  
  # Add annotations
  p$plot <- p$plot +
    annotate("text", x = 0.25, y = 0.15, label = c_index_text, hjust = 0, size = 6) +
    annotate("text", x = 0.25, y = 0.1, label = hr_text, hjust = 0, size = 6)
  
  return(p)
}

# Function to create treatment comparison plots
create_treatment_comparison <- function(control_data, treatment_data, risk_level, title, endpoint = "CRPC") {
  """
  Create survival plots comparing treatment arms within specific risk group
  
  Args:
    control_data: Control arm dataset
    treatment_data: Treatment arm dataset
    risk_level: Risk group to analyze ("Low Risk" or "High Risk")
    title: Plot title
    endpoint: Endpoint type ("CRPC" or "OS")
  
  Returns:
    ggsurvplot object
  """
  
  # Combine data for specified risk group
  combined_data <- bind_rows(
    control_data %>% mutate(treatment = "Control"),
    treatment_data %>% mutate(treatment = "Treatment")
  ) %>%
    filter(risk_group == risk_level)
  
  if (endpoint == "CRPC") {
    fit <- survfit(Surv(TT_CRPC, CRPC) ~ treatment, data = combined_data)
    ylab <- "CRPC-Free Survival Probability"
  } else {
    fit <- survfit(Surv(os, dead) ~ treatment, data = combined_data)
    ylab <- "Overall Survival Probability"
  }
  
  # Calculate treatment effect statistics
  cox_temp <- coxph(Surv(TT_CRPC, CRPC) ~ treatment, data = combined_data)
  cox_summary <- summary(cox_temp)
  
  hr <- cox_summary$conf.int[1, "exp(coef)"]
  hr_lower <- cox_summary$conf.int[1, "lower .95"]
  hr_upper <- cox_summary$conf.int[1, "upper .95"]
  hr_text <- sprintf("HR: %.2f (95%% CI: %.2f-%.2f)", hr, hr_lower, hr_upper)
  
  # Create plot
  p <- ggsurvplot(
    fit,
    data = combined_data,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    title = title,
    xlab = "Time (Months)",
    ylab = ylab,
    legend.title = "Treatment",
    legend.labs = c("Control", "Treatment"),
    palette = c("#ff8439", "#03c03c"),
    surv.median.line = "hv",
    risk.table.height = 0.2,
    font.main = c(18, "bold"),
    font.x = c(16),
    font.y = c(16),
    font.legend = c(16),
    font.tickslab = c(16),
    risk.table.fontsize = 6,
    pval.size = 6
  )
  
  # Add HR annotation
  p$plot <- p$plot +
    annotate("text", x = 0.25, y = 0.1, label = hr_text, hjust = 0, size = 6)
  
  return(p)
}

# Main analysis pipeline
run_apic_analysis <- function(chaarted_filepath = NULL, rtog_filepath = NULL) {
  """
  Complete APIC analysis pipeline for both trials
  
  Args:
    chaarted_filepath: Path to CHAARTED dataset CSV
    rtog_filepath: Path to RTOG dataset CSV
  
  Returns:
    List containing all results and plots
  """
  
  results <- list()
  
  # CHAARTED Analysis
  if (!is.null(chaarted_filepath)) {
    cat("=== CHAARTED ANALYSIS ===\n")
    
    # Prepare data
    chaarted_data <- prepare_data(
      chaarted_filepath, 
      CHAARTED_FEATURES,
      ASSIGNED_TX_ARM == "B",  # Control arm
      ASSIGNED_TX_ARM == "A"   # Treatment arm
    )
    
    # Develop model
    chaarted_model <- develop_apic_model(
      chaarted_data$train, 
      CHAARTED_FEATURES, 
      "CRPC"
    )
    
    # Validate model
    chaarted_validation <- validate_apic_model(
      chaarted_model,
      chaarted_data$validation,
      chaarted_data$treatment,
      "CRPC"
    )
    
    # Create plots
    plots_chaarted <- list(
      training = create_survival_plots(chaarted_data$train, "CHAARTED Training", "CRPC"),
      validation = create_survival_plots(chaarted_validation$validation_data, "CHAARTED Validation", "CRPC"),
      treatment_high_risk = create_treatment_comparison(
        chaarted_validation$validation_data,
        chaarted_validation$treatment_data,
        "High Risk",
        "CHAARTED: APIC-Positive",
        "CRPC"
      ),
      treatment_low_risk = create_treatment_comparison(
        chaarted_validation$validation_data,
        chaarted_validation$treatment_data,
        "Low Risk",
        "CHAARTED: APIC-Negative",
        "CRPC"
      )
    )
    
    results$chaarted <- list(
      data = chaarted_data,
      model = chaarted_model,
      validation = chaarted_validation,
      plots = plots_chaarted
    )
  }
  
  # RTOG Analysis
  if (!is.null(rtog_filepath)) {
    cat("\n=== RTOG 0521 ANALYSIS ===\n")
    
    # Prepare data (adjust filters for RTOG data structure)
    rtog_data <- prepare_data(
      rtog_filepath,
      RTOG_FEATURES,
      RX == 1,  # Control arm
      RX == 2   # Treatment arm
    )
    
    # Convert years to months for RTOG data
    rtog_data$train$os <- rtog_data$train$os * 12
    rtog_data$validation$os <- rtog_data$validation$os * 12
    rtog_data$treatment$os <- rtog_data$treatment$os * 12
    
    # Develop model
    rtog_model <- develop_apic_model(
      rtog_data$train,
      RTOG_FEATURES,
      "OS"
    )
    
    # Use median threshold for RTOG
    rtog_model$threshold <- quantile(predict(rtog_model$model, type = "risk", newdata = rtog_data$train), 0.5)
    
    # Validate model
    rtog_validation <- validate_apic_model(
      rtog_model,
      rtog_data$validation,
      rtog_data$treatment,
      "OS"
    )
    
    # Create plots
    plots_rtog <- list(
      training = create_survival_plots(rtog_data$train, "RTOG Training", "OS"),
      validation = create_survival_plots(rtog_validation$validation_data, "RTOG Validation", "OS"),
      treatment_high_risk = create_treatment_comparison(
        rtog_validation$validation_data,
        rtog_validation$treatment_data,
        "High Risk",
        "RTOG: APIC-Positive",
        "OS"
      ),
      treatment_low_risk = create_treatment_comparison(
        rtog_validation$validation_data,
        rtog_validation$treatment_data,
        "Low Risk",
        "RTOG: APIC-Negative",
        "OS"
      )
    )
    
    results$rtog <- list(
      data = rtog_data,
      model = rtog_model,
      validation = rtog_validation,
      plots = plots_rtog
    )
  }
  
  return(results)
}