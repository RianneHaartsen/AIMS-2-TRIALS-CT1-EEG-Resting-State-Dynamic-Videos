# AIMS-CT1 dynamic videos & resting state: Statistical analyses

# This scripts includes the statistical analysis for the AIMS-CT1 study for the 
# EEG resting state and dynamic videos tasks. See the paper for more details. 
# This script contains contributions from Rianne Haartsen, Teresa Del Bianco, and Emily J.H. Jones. 

# Sections 1 through 3 aim to test modulations of arbaclofen on change in EEG:
# The following EEG metrics are tested:
# 1) 1/f (narrow range, 3-28Hz)
# a) 1/f offset
# b) 1/f slope
# 2) Spectral power
# a1) Low frequencies in absolute power
# a2) Follow up in canonical frequency bands
# a3) Follow up in relative and log power
# b1) High frequencies in absolute power
# b2) Follow up in canonical frequency bands
# b3) Follow up in relative and log power
# 3) Functional connectivity
# a) Theta frequency band
# b) Alpha frequency band
# 4) Clinical prediction
# a) Alpha frequency band
# b) High frequency bands
# Section 4 aims to test if baseline EEG can predict change in clinical scores. 
# Section 5 extracts the characteristics of the included sample for reporting in the main paper. 
# Section 6 contains explotary analyses on mental health meaures. 

# created by Rianne Haartsen, PhD. on 12th of September, 2025
# updated figures on 9th of March, 2026
# additional analyses added spring 2026


#### load R packages ####
packs <- c("readr", "dplyr", "tidyr", "reshape2", "lmerTest", "kableExtra", "ggplot2", "tibble", "emmeans", "readxl","ggeffects","flextable", "officer", "grid", "gridExtra", "car")
lapply(packs, require, character.only = TRUE)
setwd("xxx/")

## Handy functions for later

# Function for plots to check assumptions of lmer (EEG change)
check_lmer_assumptions <- function(model, data, title = "Model Diagnostics") {
  
  # Add residuals and fitted values
  data$fitted <- fitted(model)
  data$residuals <- residuals(model)
  data$sqrt_abs_resid <- sqrt(abs(data$residuals))
  
  # Create plots
  # 1) Residuals vs fitted values
  p1 <- ggplot(data, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "cornflowerblue", linetype = "dashed") +
    geom_smooth(se = FALSE, color = "cyan3") +
    labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") +
    theme_minimal()
  # 2) Q-Q plot
  p2 <- ggplot(data, aes(sample = residuals)) +
    stat_qq(alpha = 0.5) +
    stat_qq_line(color = "cornflowerblue") +
    labs(title = "Q-Q Plot", x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_minimal()
  # 3) Scale-Location plot
  p3 <- ggplot(data, aes(x = fitted, y = sqrt_abs_resid)) +
    geom_point(alpha = 0.5) +
    geom_smooth(se = FALSE, color = "cyan3") +
    labs(title = "Scale-Location", x = "Fitted", y = "√|Residuals|") +
    theme_minimal()
  # 4) Histogram of residuals
  p4 <- ggplot(data, aes(x = residuals)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    labs(title = "Histogram of Residuals", x = "Residuals", y = "Count") +
    theme_minimal()
  
  # Combine plots
  grid.arrange(p1, p2, p3, p4, ncol = 2,
               top = textGrob(title, 
                              gp = gpar(fontsize = 16, fontface = "bold")))
}

# Function for plots to check assumptions of lm (clinical change)
check_lm_clinchange_assumptions <- function(model, data, title = "Model Diagnostics") {
  
  # Add residuals and fitted values
  data$fitted <- fitted(model)
  data$residuals <- residuals(model)
  data$sqrt_abs_resid <- sqrt(abs(data$residuals))
  # calculate influence statistics for Residuals vs Leverage plot
  data$leverage <- hatvalues(model)
  data$cooks_d <- cooks.distance(model)
  data$std_residuals <- rstandard(model)
  
  # Create plots
  # 1) Residuals vs fitted values
  p1 <- ggplot(data, aes(x = fitted, y = residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "cornflowerblue", linetype = "dashed") +
    geom_smooth(se = FALSE, color = "cyan3") +
    labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") +
    theme_minimal()
  # 2) Q-Q plot
  p2 <- ggplot(data, aes(sample = residuals)) +
    stat_qq(alpha = 0.5) +
    stat_qq_line(color = "cornflowerblue") +
    labs(title = "Q-Q Plot", x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_minimal()
  # 3) Scale-Location plot
  p3 <- ggplot(data, aes(x = fitted, y = sqrt_abs_resid)) +
    geom_point(alpha = 0.5) +
    geom_smooth(se = FALSE, color = "cyan3") +
    labs(title = "Scale-Location", x = "Fitted", y = "√|Residuals|") +
    theme_minimal()
  # 4) Residuals vs Leverage plot
  # Create Cook's distance contour lines
  cook_levels <- c(0.5, 1.0)
  p <- length(coef(model))  # number of parameters
  n <- nrow(data)
  # Calculate limits dynamically based on your data
  y_range <- range(data$std_residuals, na.rm = TRUE)
  y_limits <- c(floor(y_range[1])-1, ceiling(y_range[2])+1)
  p4 <- ggplot(data, aes(x = leverage, y = std_residuals)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "cornflowerblue", linetype = "dashed") +
    geom_smooth(se = FALSE, color = "cyan3", linewidth = 0.5) +
    # Add Cook's distance contour lines
    geom_function(fun = function(x) sqrt(cook_levels[1] * p * (1 - x) / x), 
                  color = "blue4", linetype = "dotted", linewidth = 0.5) +
    geom_function(fun = function(x) -sqrt(cook_levels[1] * p * (1 - x) / x), 
                  color = "blue4", linetype = "dotted", linewidth = 0.5) +
    # Add Cook's distance legend
    annotate("text", x = Inf, y = Inf, 
             label = "Cook's distance", 
             hjust = 1.1, vjust = 1.5, 
             size = 3, color = "blue4", fontface = "italic") +
    annotate("text", x = Inf, y = Inf, 
             label = paste0(cook_levels[1]), 
             hjust = 1.1, vjust = 3, 
             size = 3, color = "blue4") +
    coord_cartesian(ylim = y_limits) + 
    labs(title = "Residuals vs Leverage", 
         x = "Leverage", 
         y = "Standardized Residuals") +
    theme_minimal()
  
  # Combine plots
  grid.arrange(p1, p2, p3, p4, ncol = 2,
               top = textGrob(title, 
                              gp = gpar(fontsize = 16, fontface = "bold")))
}

# Table for reporting lmer results
lmer_results_table <- function(m_int) {
  #calculate CIs
  ci_int <- confint(m_int, method = "boot", nsim = 1000)
  tidy_ci <- ci_int %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "Term") %>%
    filter(!grepl("^\\.sig", Term)) %>%
    mutate_if(is.numeric, function(x)
      round(x, 2))
  
  # Values for reporting in paper
  summary_output <- summary(m_int)
  # Extract fixed-effects coefficients
  tidy_summary <- as.data.frame(coef(summary_output))
  # Rename columns for clarity
  colnames(tidy_summary) <- c("Estimate", "SE", "df","t-value", "p-value")
  # Convert row names (terms) into a column
  tidy_summary <- tidy_summary %>% tibble::rownames_to_column(var = "Term")
  # round values off
  tidy_summary[, c(2, 3, 4, 5)] <- round(tidy_summary[, c(2, 3, 4, 5)],2)
  tidy_summary[, c(4)] <- round(tidy_summary[, c(4)],0)
  tidy_summary[, c(6)] <- round(tidy_summary[, c(6)],3)
  
  # Rename columns for clarity
  colnames(tidy_ci) <- c("Term", "2.5% CI", "97.5% CI")
  tidy_summary2 <- tidy_summary[, 1:2] %>%
    left_join(tidy_ci, by = "Term")
  tidy_summary2 <- tidy_summary2 %>%
    left_join(tidy_summary[, c(1,3:6)], by = "Term")
  tidy_summary2
  
  # Create a nice formatted table
  nice_table <- flextable(tidy_summary2)
}

# Table for reporting lmer results for functional connectivity, using scientific notation
lmer_FCresults_table <- function(m_int) {
  #calculate CIs
  ci_int <- confint(m_int, method = "boot", nsim = 1000)
  tidy_ci <- ci_int %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "Term") 
  
  # Values for reporting in paper
  summary_output <- summary(m_int)
  # Extract fixed-effects coefficients
  tidy_summary <- as.data.frame(coef(summary_output))
  # Rename columns for clarity
  colnames(tidy_summary) <- c("Estimate", "SE", "df","t-value", "p-value")
  # Convert row names (terms) into a column
  tidy_summary <- tidy_summary %>% tibble::rownames_to_column(var = "Term")
  # round values off
  tidy_summary[, c(4)] <- round(tidy_summary[, c(4)],0)
  tidy_summary[, c(5)] <- round(tidy_summary[, c(5)],2)
  tidy_summary[, c(6)] <- round(tidy_summary[, c(6)],3)
  
  # Rename columns for clarity
  colnames(tidy_ci) <- c("Term", "2.5% CI", "97.5% CI")
  tidy_summary2 <- tidy_summary[, 1:2] %>%
    left_join(tidy_ci, by = "Term")
  tidy_summary2 <- tidy_summary2 %>%
    left_join(tidy_summary[, c(1,3:6)], by = "Term")

  # Create a nice formatted table using scientific notation
  nice_table <- flextable(tidy_summary2 %>%
    mutate(across(c(2, 3, 4, 5), ~ format(., scientific = TRUE)))) 
}


# VIF table 
lmer_VIF_table <- function(m_int) {
  vif_values <- car::vif(m_int)
  vif_table <- vif_values %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "Term") %>%
    mutate_if(is.numeric, function(x)
      round(x, 2))
  vif_table
  VIFtable <- flextable(vif_table)
}

# Table for reporting lm results
lm_results_table <- function(m_int) {
  #calculate CIs
  ci_int <- confint(m_int, method = "boot", nsim = 1000)
  tidy_ci <- ci_int %>% as.data.frame() %>%
    tibble::rownames_to_column(var = "Term") %>%
    filter(!grepl("^\\.sig", Term)) %>%
    mutate_if(is.numeric, function(x)
      round(x, 2))
  
  # Values for reporting in paper
  summary_output <- summary(m_int)
  # Extract fixed-effects coefficients
  tidy_summary <- as.data.frame(coef(summary_output))
  # Rename columns for clarity
  colnames(tidy_summary) <- c("Estimate", "SE","t-value", "p-value")
  # Convert row names (terms) into a column
  tidy_summary <- tidy_summary %>% tibble::rownames_to_column(var = "Term")
  # round values off
  tidy_summary[, c(2, 3, 4)] <- round(tidy_summary[, c(2, 3, 4)],2)
  tidy_summary[, c(5)] <- round(tidy_summary[, c(5)],3)
  
  # Rename columns for clarity
  colnames(tidy_ci) <- c("Term", "2.5% CI", "97.5% CI")
  tidy_summary2 <- tidy_summary[, 1:2] %>%
    left_join(tidy_ci, by = "Term")
  tidy_summary2 <- tidy_summary2 %>%
    left_join(tidy_summary[, c(1,3:5)], by = "Term")
  tidy_summary2
  
  # Create a nice formatted table
  nice_table <- flextable(tidy_summary2)
  
  # Extract model statistics
  summary_output <- summary(m_int)
  r_sq <- round(summary_output$r.squared, 3)
  adj_r_sq <- round(summary_output$adj.r.squared, 3)
  f_stat <- summary_output$fstatistic
  p_val <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
  
  # Create flextable with model stats in footer
  nice_table <- flextable(tidy_summary2) %>%
    add_footer_lines(values = paste0(
      "Model: R² = ", r_sq, 
      ", Adj. R² = ", adj_r_sq,
      ", F(", f_stat[2], ", ", f_stat[3], ") = ", round(f_stat[1], 2),
      ", p ", ifelse(p_val < 0.001, "< .001", paste0("= ", round(p_val, 3)))
    ))
  
  return(nice_table)
}

#### 1) 1/f (narrow range, 3-28Hz) ####
#### 1.a) 1/f offset ####
# read in data
d_eeg <- read.csv(
  "RCT1_1oF_longformat_CxRxFB.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_1oF_longformat_CxRxFB.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Ranges == "Narrow") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Offset = ifelse(is.nan(Offset), NA, Offset),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0) %>%
  filter(Rsq >= 0.98) 

# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Offset) %>%
  dcast(Participant_ID + Condition + Regions ~ Sessions, value.var = "Offset")
# calculate average baseline on test sessions
Offset_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Offset, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(Offset_test_avg_val = Offset_test_avg) %>% 
  mutate(Offset_Bs_centred = test - Offset_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(Offset_diff = retest - test) %>%
  ungroup() 


#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")

#select variables for model
fixed_cols <- c(
  "test",
  "retest",
  "Offset_diff",
  "treatment_group_code",
  "Condition",
  "Regions",
  "Offset_Bs_centred",
  "Age_patient_v0"
)
data_oof_o <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]

# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, Offset_test_avg, fixed_cols)


# Interaction model:
m_int <- lmerTest::lmer(
  Offset_diff ~
    (treatment_group_code * Offset_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_oof_o,
  control = lmerControl(optimizer = "bobyqa")
)

summary(m_int)
# create table for reporting
oof_offset_table <- lmer_results_table(m_int)

# Check the sample size
ppts <- (rownames(coef(m_int)$Participant_ID))
save(ppts, file = "IDs_incl_OOFfeatures.RData")
ppts_incl <- data_oof_o %>%
  filter(Participant_ID %in% ppts) %>%
  distinct(Participant_ID, treatment_group_code) %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl, ppts)


# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(Offset_diff ~ (Age_patient_v0) + (Condition) + (1 + Regions | Participant_ID), 
                     data = data_oof_o, control = lmerControl(optimizer = "bobyqa"))
# Now calculate partial residuals
data_oof_o$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_oof_o %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_oof_o$Offset_Bs_centred, na.rm = TRUE), max(data_oof_o$Offset_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "Offset_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(Offset_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(Offset_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (Offset_Bs_centred.trend * Offset_Bs_centred))
# For individual points: average within participant
participant_averages <- data_oof_o %>%
  group_by(Participant_ID) %>%
  reframe(
    Offset_Bs_centred_avg = mean(Offset_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = Offset_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = Offset_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_OoF_offset_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_oof_o, 
             aes(x = Offset_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = Offset_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_OoF_offset_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int)

# 2) Check direction of effect of Age
# Get predictions across age range
age_predictions <- ggpredict(m_int, terms = "Age_patient_v0")
# Plot
plot(age_predictions) +
  labs(title = "Effect of Patient Age on Offset Change",
       x = "Patient Age (years)",
       y = "Predicted Offset Difference") +
  theme_minimal()
rm(age_predictions)

# 2) Check direction of effect of baseline
bl_predictions <- ggpredict(m_int, terms = "Offset_Bs_centred")
# Plot
plot(bl_predictions) +
  labs(title = "Effect of baseline Offset on Offset Change",
       x = "EEG at baseline",
       y = "Predicted Offset Difference") +
  theme_minimal()
rm(bl_predictions)



# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_oof_o, title = "Model assumptions for 1/f offset")
ggsave("OoF_offset_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") # high-quality resolution
# VIF table to check collinearity
oof_offset_VIFtable <- lmer_VIF_table(m_int)

# Clean up variables
rm(nb_vals, slopes, slopes_df, lines_df, plot_ass)
rm(m_plot, m_int, data_oof_o)
rm(emm, emm_df, participant_averages)


#### 1.b) 1/f slope ####
# read in data
d_eeg <- read.csv(
  "RCT1_1oF_longformat_CxRxFB.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_1oF_longformat_CxRxFB.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Ranges == "Narrow") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Offset = ifelse(is.nan(Slope), NA, Slope),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0) %>%
  filter(Rsq >= 0.98) 

# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Slope) %>%
  dcast(Participant_ID + Condition + Regions ~ Sessions, value.var = "Slope")
# calculate average baseline on test sessions
Slope_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Slope, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(Slope_test_avg_val = Slope_test_avg) %>% 
  mutate(Slope_Bs_centred = test - Slope_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(Slope_diff = retest - test) %>%
  ungroup() 

#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")

#select variables for model
fixed_cols <- c(
  "test",
  "retest",
  "Slope_diff",
  "treatment_group_code",
  "Condition",
  "Regions",
  "Slope_Bs_centred",
  "Age_patient_v0"
)
data_oof_s <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]

# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_filt, d_eeg_filt, Slope_test_avg, fixed_cols)

# Run statistical model
m_int <- lmerTest::lmer(
  Slope_diff ~
    (treatment_group_code * Slope_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_oof_s,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
oof_slope_table <- lmer_results_table(m_int)



# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(Slope_diff ~ (Age_patient_v0) + (Condition) + (1 + Regions | Participant_ID), 
                     data = data_oof_s, control = lmerControl(optimizer = "bobyqa"))
# Now calculate partial residuals
data_oof_s$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_oof_s %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_oof_s$Slope_Bs_centred, na.rm = TRUE), max(data_oof_s$Slope_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "Slope_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(Slope_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(Slope_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (Slope_Bs_centred.trend * Slope_Bs_centred))
# For individual points: average within participant
participant_averages <- data_oof_s %>%
  group_by(Participant_ID) %>%
  reframe(
    Slope_Bs_centred_avg = mean(Slope_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = Slope_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = Slope_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_OoF_slope_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_oof_s, 
             aes(x = Slope_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = Slope_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_OoF_slope_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int)

# 2) Check direction of effect of Age
# Get predictions across age range
age_predictions <- ggpredict(m_int, terms = "Age_patient_v0")
# Plot
plot(age_predictions) +
  labs(title = "Effect of Patient Age on Slope Change",
       x = "Patient Age (years)",
       y = "Predicted Slope Difference") +
  theme_minimal()
rm(age_predictions)

# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_oof_s, title = "Model assumptions for 1/f slope")
ggsave("OoF_slope_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
oof_slope_VIFtable <- lmer_VIF_table(m_int)

# clear up
rm(nb_vals, slopes, slopes_df, lines_df, plot_ass)
rm(m_plot, m_int, data_oof_s)
rm(emm, emm_df, participant_averages)


#### Save 1/f tables into Word document ####
doc <- read_docx() %>%
  body_add_par("Table 1: 1/f offset - Model output", style = "heading 1") %>%
  body_add_flextable(oof_offset_table) %>%
  body_add_par("") %>%
  body_add_par("Table 2: 1/f offset - VIF table", style = "heading 1") %>%
  body_add_flextable(oof_offset_VIFtable) %>%
  body_add_par("") %>%
  body_add_par("Table 3: 1/f slope - Model output", style = "heading 1") %>%
  body_add_flextable(oof_slope_table) %>%
  body_add_par("") %>%
  body_add_par("Table 4: 1/f slope - VIF table", style = "heading 1") %>%
  body_add_flextable(oof_slope_VIFtable) 

print(doc, target = "OoF_model_tables.docx")

# clean up
rm(oof_offset_table, oof_offset_VIFtable, oof_slope_table, oof_slope_VIFtable, doc)



####  2) Spectral power #### 
#### 2.a.1) Low frequencies in log power ####
d_eeg <- read.csv(
  "RCT1_Power_LvsHfreqs_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_LvsHfreqs_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Low") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Frequencies = as.factor(Frequencies),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(lp_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "lp_diff", "treatment_group_code", "Condition", "Regions", "LPow_Bs_centred", "Age_patient_v0"
)
data_lf_lpow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, LPow_test_avg, fixed_cols)

# Run interaction model and check the interaction term
m_int <- lmerTest::lmer(
  lp_diff ~
    (treatment_group_code * LPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_lf_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
lf_lpow_table <- lmer_results_table(m_int)

# Check the sample size and save IDs for clinical predictions
ppts <- (rownames(coef(m_int)$Participant_ID))
save(ppts, file = "IDs_incl_LFreq_logpow.RData")
ppts_incl <- data_lf_lpow %>%
  filter(Participant_ID %in% ppts) %>%
  distinct(Participant_ID, treatment_group_code) %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl, ppts)

# Additional check for group differences in baseline EEG
m_bsEEGdiffs <- lmerTest::lmer(
  LPow_Bs_centred ~
    treatment_group_code + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_lf_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_bsEEGdiffs)
# create table for reporting
bsEEG_lf_lpow_table <- lmer_results_table(m_bsEEGdiffs)

# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(lp_diff ~ (Age_patient_v0) + (Condition) + (1 + Regions | Participant_ID), 
                     data = data_lf_lpow, control = lmerControl(optimizer = "bobyqa"))
# Now calculate partial residuals
data_lf_lpow$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_lf_lpow %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_lf_lpow$LPow_Bs_centred, na.rm = TRUE), max(data_lf_lpow$LPow_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "LPow_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(LPow_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(LPow_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (LPow_Bs_centred.trend * LPow_Bs_centred))
# For individual points: average within participant
participant_averages <- data_lf_lpow %>%
  group_by(Participant_ID) %>%
  reframe(
    LPow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = LPow_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = LPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_LF_logpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_lf_lpow, 
             aes(x = LPow_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = LPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_LF_logpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int)

# 2) Check direction of effect of Age
age_predictions <- ggpredict(m_int, terms = "Age_patient_v0")
# Plot
plot(age_predictions) +
  labs(title = "Effect of Patient Age on Lf log power Change",
       x = "Patient Age (years)",
       y = "Predicted power Difference") +
  theme_minimal()
rm(age_predictions)

# 3) Check direction of effect of baseline
bl_predictions <- ggpredict(m_int, terms = "LPow_Bs_centred")
# Plot
plot(bl_predictions) +
  labs(title = "Effect of baseline EEG on EEG Change",
       x = "EEG at baseline",
       y = "Predicted Offset Difference") +
  theme_minimal()
rm(bl_predictions)


# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_lf_lpow, title = "Model assumptions for low frequency log power")
ggsave("LF_logpow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
lf_lpow_VIFtable <- lmer_VIF_table(m_int)

rm(nb_vals, lines_df, plot, slopes, slopes_d, plot_ass)
rm(plot, data_lf_lpow, m_int)
rm(emm, emm_df, slopes_df, participant_averages)


#### 2.a.2) Follow up in canonical frequency bands ####
## For log power (better distributions)
#### 2.a.2.1) Delta frequency band ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Delta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(dp_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest","dp_diff", "treatment_group_code", "Condition", "Regions", "LPow_Bs_centred", "Age_patient_v0"
)
data_delta_lpow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, LPow_test_avg, fixed_cols)

# Run statistical model
m_int <- lmerTest::lmer(
  dp_diff ~
    (treatment_group_code * LPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_delta_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
delta_lpow_table <- lmer_results_table(m_int)

# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_delta_lpow, title = "Model assumptions for delta log power")
ggsave("FU_Delta_logpow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
delta_lpow_VIFtable <- lmer_VIF_table(m_int)

rm(plot_ass)
rm(m_int, data_delta_lpow)



#### 2.a.2.2) Theta frequency band ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Theta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(tp_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest","tp_diff", "treatment_group_code", "Condition", "Regions", "LPow_Bs_centred", "Age_patient_v0"
)
data_theta_lpow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, LPow_test_avg, fixed_cols)
# Run interaction model and check the interaction term
m_int <- lmerTest::lmer(
  tp_diff ~
    (treatment_group_code * LPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_theta_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
theta_lpow_table <- lmer_results_table(m_int)

# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_theta_lpow, title = "Model assumptions for theta log power")
ggsave("FU_Theta_logpow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
theta_lpow_VIFtable <- lmer_VIF_table(m_int)

rm(plot_ass)
rm(m_int, data_theta_lpow)


#### 2.a.2.3) Alpha frequency band ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Alpha") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(ap_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest","ap_diff", "treatment_group_code", "Condition", "Regions", "LPow_Bs_centred", "Age_patient_v0"
)
data_alpha_lpow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_filt, d_eeg_filt, LPow_test_avg, fixed_cols)

# Run interaction model and check the interaction term
m_int <- lmerTest::lmer(
  ap_diff ~
    (treatment_group_code * LPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_alpha_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
alpha_lpow_table <- lmer_results_table(m_int)


# Additional check for group differences in baseline EEG
m_bsEEGdiffs <- lmerTest::lmer(
  LPow_Bs_centred ~
    treatment_group_code + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_alpha_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_bsEEGdiffs)
# create table for reporting
bsEEG_alpha_lpow_table <- lmer_results_table(m_bsEEGdiffs)


# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(ap_diff ~ (Age_patient_v0) + (Condition) + (1 + Regions | Participant_ID), 
                     data = data_alpha_lpow, control = lmerControl(optimizer = "bobyqa"))
# Now calculate partial residuals
data_alpha_lpow$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_alpha_lpow %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_alpha_lpow$LPow_Bs_centred, na.rm = TRUE), max(data_alpha_lpow$LPow_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "LPow_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(LPow_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(LPow_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (LPow_Bs_centred.trend * LPow_Bs_centred))
# For individual points: average within participant
participant_averages <- data_alpha_lpow %>%
  group_by(Participant_ID) %>%
  reframe(
    LPow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = LPow_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = LPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_alpha_logpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_alpha_lpow, 
             aes(x = LPow_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = LPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_alpha_logpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)

# 2) Check direction of effect of Age
age_predictions <- ggpredict(m_int, terms = "Age_patient_v0")
# Plot
plot(age_predictions) +
  labs(title = "Effect of Patient Age on alpha power change",
       x = "Patient Age (years)",
       y = "Predicted power Difference") +
  theme_minimal()
rm(age_predictions)


# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_alpha_lpow, title = "Model assumptions for alpha log power")
ggsave("FU_Alpha_logpow_Check_model.jpg", plot_ass,  width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
alpha_lpow_VIFtable <- lmer_VIF_table(m_int)

rm(plot_ass)
rm(m_int, data_alpha_lpow)



#### 2.a.3) Follow up for alpha frequency band in other versions of power ####
#### 2.a.3.1) Alpha absolute power ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Alpha") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Abs_Power = ifelse(is.nan(Abs_Power), NA, Abs_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Abs_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Abs_Power")
# calculate average baseline on test sessions
APow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Abs_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(APow_test_avg_val = APow_test_avg) %>% 
  mutate(APow_Bs_centred = test - APow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(ap_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest","ap_diff", "treatment_group_code", "Condition", "Regions", "APow_Bs_centred", "Age_patient_v0"
)
data_alpha_apow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_filt, d_eeg_filt, APow_test_avg, fixed_cols)

# Run interaction model and check the interaction term
m_int <- lmerTest::lmer(
  ap_diff ~
    (treatment_group_code * APow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_alpha_apow,
  control = lmerControl(optimizer = "bobyqa")
)
# above model fails to converge, so using (1 | Region ) to simplify the model
m_int <- lmerTest::lmer(
  ap_diff ~
    (treatment_group_code * APow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 | Regions) + (1 | Participant_ID),
  data = data_alpha_apow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
alpha_apow_table <- lmer_results_table(m_int)

# Additional check for group differences in baseline EEG
m_bsEEGdiffs <- lmerTest::lmer(
  APow_Bs_centred ~
    treatment_group_code + (Age_patient_v0) +
    (Condition) + (1 | Regions) + (1 | Participant_ID),
  data = data_alpha_apow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_bsEEGdiffs)
# create table for reporting
bsEEG_alpha_apow_table <- lmer_results_table(m_bsEEGdiffs)




# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(ap_diff ~ (Age_patient_v0) + (Condition) + (1 | Regions) + (1 | Participant_ID), 
                     data = data_alpha_apow, control = lmerControl(optimizer = "bobyqa")
)
# Now calculate partial residuals
data_alpha_apow$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_alpha_apow %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_alpha_apow$APow_Bs_centred, na.rm = TRUE), max(data_alpha_apow$APow_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "APow_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(APow_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(APow_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (APow_Bs_centred.trend * APow_Bs_centred))
# For individual points: average within participant
participant_averages <- data_alpha_apow %>%
  group_by(Participant_ID) %>%
  reframe(
    APow_Bs_centred_avg = mean(APow_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# find min and max for the axes
vals_X <- participant_averages %>%
  summarize(APow_Bs_centred_avg.min = min(APow_Bs_centred_avg, na.rm = TRUE),
            APow_Bs_centred_avg.max = max(APow_Bs_centred_avg, na.rm = TRUE))
vals_X
vals_Y <- participant_averages %>%
  summarize(partial_residuals_avg.min = min(partial_residuals_avg, na.rm = TRUE),
            partial_residuals_avg.max = max(partial_residuals_avg, na.rm = TRUE))
vals_Y
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = APow_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = APow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  xlim(-2.5, 6) +
  ylim(-6, 4) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10b_alpha_abspow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(vals_X, vals_Y)

# Plot with all individual points
ggplot() +
  geom_point(data = data_alpha_apow, 
             aes(x = APow_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = APow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_alpha_abspow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)

# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_alpha_apow, title = "Model assumptions for alpha absolute power")
ggsave("FU_Alpha_abspow_Check_model.jpg", plot_ass,  width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
alpha_apow_VIFtable <- lmer_VIF_table(m_int)

rm(nb_vals, lines_df, plot, slopes, slopes_df, plot_ass)
rm(plot, m_int, data_alpha_apow)



#### 2.a.3.2) Alpha relative power ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Alpha") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Rel_Power = ifelse(is.nan(Rel_Power), NA, Rel_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Frequencies = as.factor(Frequencies),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Rel_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Rel_Power")
# calculate average baseline on test sessions
RPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Rel_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(RPow_test_avg_val = RPow_test_avg) %>% 
  mutate(RPow_Bs_centred = test - RPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(rp_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "rp_diff", "treatment_group_code", "Condition", "Regions", "RPow_Bs_centred", "Age_patient_v0"
)
data_alpha_rpow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, RPow_test_avg, fixed_cols)

# Run interaction model and check the interaction term
m_int <- lmerTest::lmer(
  rp_diff ~
    (treatment_group_code * RPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_alpha_rpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
alpha_rpow_table <- lmer_results_table(m_int)

# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_alpha_rpow, title = "Model assumptions for alpha relative power")
ggsave("FU_Alpha_relpow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
alpha_rpow_VIFtable <- lmer_VIF_table(m_int)

# Clear up
rm(nb_vals, lines_df, plot, slopes, slopes_df, plot_ass)
rm(plot, data_alpha_rpow, m_int)


#### 2.a.3.3) Alpha aperiodic-adjusted power ####
d_eeg <- read.csv(
  "RCT1_ApAdjPower_CanBands_longformat_CxRxFB.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_ApAdjPower_CanBands_longformat_CxRxFB.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Alpha") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Aa_Power = ifelse(is.nan(Nr_Aapow), NA, Nr_Aapow),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Aa_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Aa_Power")
# calculate average baseline on test sessions
AaPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Aa_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(AaPow_test_avg_val = AaPow_test_avg) %>% 
  mutate(AaPow_Bs_centred = test - AaPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(aap_diff = retest - test) %>%
  ungroup() 
# get ppts included in change EEG analyses
load("IDs_incl_OOFfeatures.RData")
d_eeg_filt_ids <- d_filt %>%
  filter(Participant_ID %in% ppts)
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_eeg_filt_ids, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest","aap_diff", "treatment_group_code", "Condition", "Regions", "AaPow_Bs_centred", "Age_patient_v0"
)
data_alpha_aapow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_filt, d_eeg_filt, d_eeg_filt_ids, AaPow_test_avg, fixed_cols)

# Run interaction model and check the interaction term
m_int <- lmerTest::lmer(
  aap_diff ~
    (treatment_group_code * AaPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_alpha_aapow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
alpha_aapow_table <- lmer_results_table(m_int)

# Additional check for group differences in baseline EEG
m_bsEEGdiffs <- lmerTest::lmer(
  AaPow_Bs_centred ~
    (treatment_group_code) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_alpha_aapow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_bsEEGdiffs)
# create table for reporting
bsEEG_alpha_aapow_table <- lmer_results_table(m_bsEEGdiffs)


# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(aap_diff ~ (Age_patient_v0) + (Condition) + (1 | Regions) + (1 | Participant_ID), 
                     data = data_alpha_aapow, control = lmerControl(optimizer = "bobyqa")
)
# Now calculate partial residuals
data_alpha_aapow$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_alpha_aapow %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_alpha_aapow$AaPow_Bs_centred, na.rm = TRUE), max(data_alpha_aapow$AaPow_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "AaPow_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(AaPow_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(AaPow_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (AaPow_Bs_centred.trend * AaPow_Bs_centred))
# For individual points: average within participant
participant_averages <- data_alpha_aapow %>%
  group_by(Participant_ID) %>%
  reframe(
    AaPow_Bs_centred_avg = mean(AaPow_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = AaPow_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = AaPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_alpha_ApAdjpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_alpha_aapow, 
             aes(x = AaPow_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = AaPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_alpha_ApAdjpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)


# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_alpha_aapow, title = "Model assumptions for alpha aperiodic-adjusted power")
ggsave("FU_Alpha_ApAdjpow_Check_model.jpg", plot_ass,  width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
alpha_aapow_VIFtable <- lmer_VIF_table(m_int)

rm(nb_vals, lines_df, plot, slopes, slopes_df, plot_ass)
rm(plot, m_int, data_alpha_aapow)


#### Save tables into Word document ####
doc <- read_docx() %>%
  body_add_par("Table 1.1: Low frequency log power - model output", style = "heading 1") %>%
  body_add_flextable(lf_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 1.2: Low frequency log power - VIF table", style = "heading 2") %>%
  body_add_flextable(lf_lpow_VIFtable) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 2.1: Delta log power - model output", style = "heading 1") %>%
  body_add_flextable(delta_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 2.2: Delta log power - VIF table", style = "heading 2") %>%
  body_add_flextable(delta_lpow_VIFtable) %>%
  body_add_par("") %>%
  body_add_par("Table 3.1: Theta log power - model output", style = "heading 1") %>%
  body_add_flextable(theta_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 3.2: Theta log power - VIF table", style = "heading 2") %>%
  body_add_flextable(theta_lpow_VIFtable) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 4.1: Alpha log power - model output", style = "heading 1") %>%
  body_add_flextable(alpha_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 4.2: Alpha log power - VIF table", style = "heading 2") %>%
  body_add_flextable(alpha_lpow_VIFtable) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 5.1: Alpha absolute power - model output", style = "heading 1") %>%
  body_add_flextable(alpha_apow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 5.2: Alpha absolute power - VIF table", style = "heading 2") %>%
  body_add_flextable(alpha_apow_VIFtable) %>%
  body_add_par("") %>%
  body_add_par("Table 6.1: Alpha relative power - model output", style = "heading 1") %>%
  body_add_flextable(alpha_rpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 6.2: Alpha relative power - VIF table", style = "heading 2") %>%
  body_add_flextable(alpha_rpow_VIFtable) %>%
  body_add_par("") %>%
  body_add_par("Table 7.1: Alpha aperiodic-adjusted power - model output", style = "heading 1") %>%
  body_add_flextable(alpha_aapow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 7.2: Alpha aperiodic-adjusted power - VIF table", style = "heading 2") %>%
  body_add_flextable(alpha_aapow_VIFtable) 

print(doc, target = "LowFrequency_power_tables.docx")

# Clear up
rm(lf_lpow_table, lf_lpow_VIFtable)
rm(delta_lpow_table, delta_lpow_VIFtable, theta_lpow_table, theta_lpow_VIFtable, alpha_lpow_table, alpha_lpow_VIFtable)
rm(alpha_apow_table, alpha_apow_VIFtable, alpha_rpow_table, alpha_rpow_VIFtable, alpha_aapow_table, alpha_aapow_VIFtable)
rm(doc)




#### 2.b.1) High frequencies in log power ####
d_eeg <- read.csv(
  "RCT1_Power_LvsHfreqs_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_LvsHfreqs_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "High") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Frequencies = as.factor(Frequencies),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(hf_lp_diff = retest - test) %>%
  ungroup() 

#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "hf_lp_diff", "treatment_group_code", "Condition", "Regions", "LPow_Bs_centred", "Age_patient_v0"
)
data_hf_lpow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, LPow_test_avg, fixed_cols)

# Run statistical model
m_int <- lmerTest::lmer(
  hf_lp_diff ~
    (treatment_group_code * LPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_hf_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
hf_lpow_table <- lmer_results_table(m_int)

# Additional check for group differences in baseline EEG
m_bsEEGdiffs <- lmerTest::lmer(
  LPow_Bs_centred ~
    (treatment_group_code) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_hf_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_bsEEGdiffs)
# create table for reporting
bsEEG_hf_lpow_table <- lmer_results_table(m_bsEEGdiffs)


# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(hf_lp_diff ~ (Age_patient_v0) + (Condition) + (1 + Regions | Participant_ID), 
                     data = data_hf_lpow, control = lmerControl(optimizer = "bobyqa"))
# Now calculate partial residuals
data_hf_lpow$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_hf_lpow %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_hf_lpow$LPow_Bs_centred, na.rm = TRUE), max(data_hf_lpow$LPow_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "LPow_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(LPow_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(LPow_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (LPow_Bs_centred.trend * LPow_Bs_centred))
# For individual points: average within participant
participant_averages <- data_hf_lpow %>%
  group_by(Participant_ID) %>%
  reframe(
    LPow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = LPow_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = LPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_HF_logpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_hf_lpow, 
             aes(x = LPow_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = LPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_HF_logpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df, emm, emm_df, participant_averages)

# 2) Check direction of effect of Age
age_predictions <- ggpredict(m_int, terms = "Age_patient_v0")
# Plot
plot(age_predictions) +
  labs(title = "Effect of Patient Age on Lf abs power Change",
       x = "Patient Age (years)",
       y = "Predicted power Difference") +
  theme_minimal()
rm(age_predictions)


# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_hf_lpow, title = "Model assumptions for high frequency log power")
ggsave("FU_HighF_logpow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
hf_lpow_VIFtable <- lmer_VIF_table(m_int)

rm(nb_vals, slopes, slopes_df, lines_df, plot, plot_ass)
rm(plot, m_int, data_hf_lpow)


#### 2.b.2) Follow up in canonical frequency bands ####
#### 2.b.2,1) Beta frequency band ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Beta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(bp_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "bp_diff", "treatment_group_code", "Condition", "Regions", "LPow_Bs_centred", "Age_patient_v0"
)
data_beta_lpow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, LPow_test_avg, fixed_cols)

# Run statistical model
m_int <- lmerTest::lmer(
  bp_diff ~
    (treatment_group_code * LPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_beta_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
beta_lpow_table <- lmer_results_table(m_int)

# Additional check for group differences in baseline EEG
m_bsEEGdiffs <- lmerTest::lmer(
  LPow_Bs_centred ~
    (treatment_group_code) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_beta_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_bsEEGdiffs)
# create table for reporting
bsEEG_beta_lpow_table <- lmer_results_table(m_bsEEGdiffs)


# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(bp_diff ~ (Age_patient_v0) + (Condition) + (1 + Regions | Participant_ID), 
                     data = data_beta_lpow, control = lmerControl(optimizer = "bobyqa"))
# Now calculate partial residuals
data_beta_lpow$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_beta_lpow %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_beta_lpow$LPow_Bs_centred, na.rm = TRUE), max(data_beta_lpow$LPow_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "LPow_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(LPow_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(LPow_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (LPow_Bs_centred.trend * LPow_Bs_centred))
# For individual points: average within participant
participant_averages <- data_beta_lpow %>%
  group_by(Participant_ID) %>%
  reframe(
    LPow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = LPow_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = LPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_beta_logpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_beta_lpow, 
             aes(x = LPow_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = LPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_beta_logpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)


# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_beta_lpow, title = "Model assumptions for beta log power")
ggsave("FU_Beta_logpow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
beta_lpow_VIFtable <- lmer_VIF_table(m_int)

rm(plot_ass, m_int, data_beta_lpow)


#### 2.b.2.2) Lower gamma frequency band ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Lower Gamma") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(lgp_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "lgp_diff", "treatment_group_code", "Condition", "Regions", "LPow_Bs_centred", "Age_patient_v0"
)
data_lgamma_lpow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, LPow_test_avg, fixed_cols)

# Run statistical model
m_int <- lmerTest::lmer(
  lgp_diff ~
    (treatment_group_code * LPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 | Regions) + (1 | Participant_ID),
  data = data_lgamma_lpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
lgamma_lpow_table <- lmer_results_table(m_int)

# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_lgamma_lpow, title = "Model assumptions for lower gamma log power")
ggsave("FU_LGamma_logpow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
lgamma_lpow_VIFtable <- lmer_VIF_table(m_int)

# Clear up
rm(plot_ass, m_int, data_lgamma_lpow)



#### 2.b.3) Follow up for beta frequency band in other versions of power ####
#### 2.b.3.1) Beta absolute power ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Beta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Abs_Power = ifelse(is.nan(Abs_Power), NA, Abs_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Abs_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Abs_Power")
# calculate average baseline on test sessions
APow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Abs_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(APow_test_avg_val = APow_test_avg) %>% 
  mutate(APow_Bs_centred = test - APow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(bp_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "bp_diff", "treatment_group_code", "Condition", "Regions", "APow_Bs_centred", "Age_patient_v0"
)
data_beta_apow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, APow_test_avg, fixed_cols)

# Run statistical model
m_int <- lmerTest::lmer(
  bp_diff ~
    (treatment_group_code * APow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_beta_apow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
beta_apow_table <- lmer_results_table(m_int)


# Additional check for group differences in baseline EEG
m_bsEEGdiffs <- lmerTest::lmer(
  APow_Bs_centred ~
    (treatment_group_code) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_beta_apow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_bsEEGdiffs)
# create table for reporting
bsEEG_beta_apow_table <- lmer_results_table(m_bsEEGdiffs)


# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(bp_diff ~ (Age_patient_v0) + (Condition) + (1 + Regions | Participant_ID), 
                     data = data_beta_apow, control = lmerControl(optimizer = "bobyqa"))
# Now calculate partial residuals
data_beta_apow$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_beta_apow %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_beta_apow$APow_Bs_centred, na.rm = TRUE), max(data_beta_apow$APow_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "APow_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(APow_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(APow_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (APow_Bs_centred.trend * APow_Bs_centred))
# For individual points: average within participant
participant_averages <- data_beta_apow %>%
  group_by(Participant_ID) %>%
  reframe(
    APow_Bs_centred_avg = mean(APow_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# find min and max for the axes
vals_X <- participant_averages %>%
  summarize(APow_Bs_centred_avg.min = min(APow_Bs_centred_avg, na.rm = TRUE),
            APow_Bs_centred_avg.max = max(APow_Bs_centred_avg, na.rm = TRUE))
vals_X
vals_Y <- participant_averages %>%
  summarize(partial_residuals_avg.min = min(partial_residuals_avg, na.rm = TRUE),
            partial_residuals_avg.max = max(partial_residuals_avg, na.rm = TRUE))
vals_Y
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = APow_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = APow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  xlim(-0.5, 1.5) +
  ylim(-1, 0.5) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10b_beta_abspow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(vals_X, vals_Y)

# Plot with all individual points
ggplot() +
  geom_point(data = data_beta_apow, 
             aes(x = APow_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = APow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_beta_abspow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)


# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_beta_apow, title = "Model assumptions for beta absolute power")
ggsave("FU_Beta_abspow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
beta_apow_VIFtable <- lmer_VIF_table(m_int)

# Clear up
rm(nb_vals, slopes, slopes_df, lines_df, plot, plot_ass)
rm(m_int, data_beta_apow)


#### 2.b.3.2) Beta relative power ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Beta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Rel_Power = ifelse(is.nan(Rel_Power), NA, Rel_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Rel_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Rel_Power")
# calculate average baseline on test sessions
RPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Rel_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(RPow_test_avg_val = RPow_test_avg) %>% 
  mutate(RPow_Bs_centred = test - RPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(bp_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "bp_diff", "treatment_group_code", "Condition", "Regions", "RPow_Bs_centred", "Age_patient_v0"
)
data_beta_rpow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, RPow_test_avg, fixed_cols)

# Run statistical model
m_int <- lmerTest::lmer(
  bp_diff ~
    (treatment_group_code * RPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_beta_rpow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
beta_rpow_table <- lmer_results_table(m_int)


# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_beta_rpow, title = "Model assumptions for beta relative power")
ggsave("FU_Beta_relpow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
beta_rpow_VIFtable <- lmer_VIF_table(m_int)

# Clear up
rm(nb_vals, slopes, slopes_df, lines_df, plot, plot_ass)
rm(m_int, data_beta_rpow)




#### 2.b.3.3) Beta aperiodic-adjusted power ####
d_eeg <- read.csv(
  "RCT1_ApAdjPower_CanBands_longformat_CxRxFB.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_ApAdjPower_CanBands_longformat_CxRxFB.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Beta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Aa_Power = ifelse(is.nan(Nr_Aapow), NA, Nr_Aapow),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Aa_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Aa_Power")
# calculate average baseline on test sessions
AaPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Aa_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(AaPow_test_avg_val = AaPow_test_avg) %>% 
  mutate(AaPow_Bs_centred = test - AaPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(aap_diff = retest - test) %>%
  ungroup() 
# get ppts included in change EEG analyses
load("IDs_incl_OOFfeatures.RData")
d_eeg_filt_ids <- d_filt %>%
  filter(Participant_ID %in% ppts)
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_eeg_filt_ids, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "aap_diff", "treatment_group_code", "Condition", "Regions", "AaPow_Bs_centred", "Age_patient_v0"
)
data_beta_aapow <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_eeg_filt_ids, d_filt, AaPow_test_avg, fixed_cols)

# Run statistical model
m_int <- lmerTest::lmer(
  aap_diff ~
    (treatment_group_code * AaPow_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_beta_aapow,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
beta_aapow_table <- lmer_results_table(m_int)




# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(aap_diff ~ (Age_patient_v0) + (Condition) + (1 + Regions | Participant_ID), 
                     data = data_beta_aapow, control = lmerControl(optimizer = "bobyqa"))
# Now calculate partial residuals
data_beta_aapow$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_beta_aapow %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_beta_aapow$AaPow_Bs_centred, na.rm = TRUE), max(data_beta_aapow$AaPow_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "AaPow_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(AaPow_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(AaPow_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (AaPow_Bs_centred.trend * AaPow_Bs_centred))
# For individual points: average within participant
participant_averages <- data_beta_aapow %>%
  group_by(Participant_ID) %>%
  reframe(
    AaPow_Bs_centred_avg = mean(AaPow_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = AaPow_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = AaPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_beta_ApAdjpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_beta_aapow, 
             aes(x = AaPow_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = AaPow_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_beta_ApAdjpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)









# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model on complete data
mplot_no_int <- lmer(aap_diff ~ (Age_patient_v0) + (Condition) + (1 + Regions | Participant_ID), 
                     data = data_beta_aapow, control = lmerControl(optimizer = "bobyqa")
)
# Now calculate partial residuals
data_beta_aapow$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_beta_aapow %>%
  mutate(case_match(treatment_group_code, 
                    "1" ~ "Arbaclofen",
                    "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_beta_aapow$AaPow_Bs_centred, na.rm = TRUE), max(data_beta_aapow$AaPow_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "AaPow_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Expand over mean-centred test values
lines_df <- slopes_df %>%
  tidyr::expand_grid(AaPow_Bs_centred = nb_vals) %>%
  mutate(predicted = AaPow_Bs_centred.trend * AaPow_Bs_centred)
lines_df %>%
  mutate(case_match(treatment_group_code, 
                    "1" ~ "Arbaclofen",
                    "2" ~ "Placebo"))
# Plot the datasets together
ggplot() +
  geom_point(data = data_beta_aapow, aes(x = AaPow_Bs_centred, y = partial_residuals, color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, aes(x = AaPow_Bs_centred, y = predicted, color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 1) +
  labs(title = "Beta frequency aperiodic-adjusted power",
       x = "Baseline EEG",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V9_Beta_ApAdjpow_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)


# Check the model assumptions
# Plots
plot_ass <- check_lmer_assumptions(m_int, data_beta_aapow, title = "Model assumptions for beta aperiodic-adjusted power")
ggsave("FU_Beta_ApAdjpow_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
beta_aapow_VIFtable <- lmer_VIF_table(m_int)

# Clear up
rm(nb_vals, slopes, slopes_df, lines_df, plot, plot_ass)
rm(m_int, data_beta_aapow)





#### Save tables into Word document ####
doc <- read_docx() %>%
  body_add_par("Table 1.1: High frequency log power - model output", style = "heading 1") %>%
  body_add_flextable(hf_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 1.2: High frequency log power - VIF table", style = "heading 2") %>%
  body_add_flextable(hf_lpow_VIFtable) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 2.1: Beta log power - model output", style = "heading 1") %>%
  body_add_flextable(beta_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 2.2: Beta log power - VIF table", style = "heading 2") %>%
  body_add_flextable(beta_lpow_VIFtable) %>%
  body_add_par("") %>%
  body_add_par("Table 3.1: Lower gamma log power - model output", style = "heading 1") %>%
  body_add_flextable(lgamma_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 3.2: Lower gamma log power - VIF table", style = "heading 2") %>%
  body_add_flextable(lgamma_lpow_VIFtable) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 4.1: Beta absolute power - model output", style = "heading 1") %>%
  body_add_flextable(beta_apow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 4.2: Beta absolute power - VIF table", style = "heading 2") %>%
  body_add_flextable(beta_apow_VIFtable) %>%
  body_add_par("") %>%
  body_add_par("Table 5.1: Beta relative power - model output", style = "heading 1") %>%
  body_add_flextable(beta_rpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 5.2: Beta relative power - VIF table", style = "heading 2") %>%
  body_add_flextable(beta_rpow_VIFtable) %>%
  body_add_par("") %>%
  body_add_par("Table 6.1: Beta aperiodic-adjusted power - model output", style = "heading 1") %>%
  body_add_flextable(beta_aapow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 6.2: Beta aperiodic-adjusted power - VIF table", style = "heading 2") %>%
  body_add_flextable(beta_aapow_VIFtable) 

print(doc, target = "HighFrequency_power_tables.docx")

# Clear up
rm(hf_lpow_table, hf_lpow_VIFtable)
rm(beta_lpow_table, beta_lpow_VIFtable, lgamma_lpow_table, lgamma_lpow_VIFtable)
rm(beta_apow_table, beta_apow_VIFtable, beta_rpow_table, beta_rpow_VIFtable, beta_aapow_table, beta_aapow_VIFtable)

rm(doc)

# additional checks of baseline EEG differences between groups
doc <- read_docx() %>%
  body_add_par("Table 1.1: Low frequency log power - model output", style = "heading 1") %>%
  body_add_flextable(bsEEG_lf_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 1.2: Alpha log power - model output", style = "heading 1") %>%
  body_add_flextable(bsEEG_alpha_lpow_table) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 1.3: Alpha absolute power - model output", style = "heading 1") %>%
  body_add_flextable(bsEEG_alpha_apow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 1.4: Alpha aperiodic-adjusted power - model output", style = "heading 1") %>%
  body_add_flextable(bsEEG_alpha_aapow_table) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 2.1: High frequency log power - model output", style = "heading 1") %>%
  body_add_flextable(bsEEG_hf_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 2.2: Beta log power - model output", style = "heading 1") %>%
  body_add_flextable(bsEEG_beta_lpow_table) %>%
  body_add_par("") %>%
  body_add_par("Table 2.3: Beta absolute power - model output", style = "heading 1") %>%
  body_add_flextable(bsEEG_beta_apow_table) 

print(doc, target = "BaselineEEG_check_tables.docx")

# Clear up
rm(hf_lpow_table, hf_lpow_VIFtable)
rm(beta_lpow_table, beta_lpow_VIFtable, lgamma_lpow_table, lgamma_lpow_VIFtable)
rm(beta_apow_table, beta_apow_VIFtable, beta_rpow_table, beta_rpow_VIFtable, beta_aapow_table, beta_aapow_VIFtable)

rm(doc)



#### 3) Functional connectivity ####
#### 3a) Theta band functional connectivity ####
d_eeg <- read.csv(
  "RCT1_FC_longformat_CxRxFB.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_FC_longformat_CxRxFB.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Theta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    FC = ifelse(is.nan(FC), NA, FC),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Frequencies = as.factor(Frequencies),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 90)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, FC) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "FC")
# calculate average baseline on test sessions
FC_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(FC, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(FC_test_avg_val = FC_test_avg) %>% 
  mutate(FC_Bs_centred = test - FC_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(fc_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "fc_diff", "treatment_group_code", "Condition", "Regions", "FC_Bs_centred", "Age_patient_v0"
)
data_th_fc <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]

# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, FC_test_avg, fixed_cols)


# Run statistical model
m_int <- lmerTest::lmer(
  fc_diff ~
    (treatment_group_code * FC_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_th_fc,
  control = lmerControl(optimizer = "bobyqa")
)
summary(m_int)
# create table for reporting
fc_theta_table <- lmer_FCresults_table(m_int)

# Check the sample size
ppts <- (rownames(coef(m_int)$Participant_ID))
ppts_incl <- data_th_fc %>%
  filter(Participant_ID %in% ppts) %>%
  distinct(Participant_ID, treatment_group_code) %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl, ppts)


# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(fc_diff ~ (Age_patient_v0) +
                       (Condition) + (1 + Regions | Participant_ID), 
                     data = data_th_fc, control = lmerControl(optimizer = "bobyqa")
)
# Now calculate partial residuals
data_th_fc$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_th_fc %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_th_fc$FC_Bs_centred, na.rm = TRUE), max(data_th_fc$FC_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "FC_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(FC_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(FC_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (FC_Bs_centred.trend * FC_Bs_centred))
# For individual points: average within participant
participant_averages <- data_th_fc %>%
  group_by(Participant_ID) %>%
  reframe(
    FC_Bs_centred_avg = mean(FC_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = FC_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = FC_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_theta_FC_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_th_fc, 
             aes(x = FC_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = FC_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_theta_FC_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)

# Check the model assumptions
plot_ass <- check_lmer_assumptions(m_int, data_th_fc, title = "Model assumptions for theta connectivity")
ggsave("Theta_FC_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
fc_theta_VIFtable <- lmer_VIF_table(m_int)

# Clear up
rm(nb_vals, slopes, slopes_df, lines_df, plot, plot_ass)
rm(m_int)




# On Yeo-Johnson transformed data
# transform FC values using the Yeo-Johnson transform
# Estimate the optimal lambda for test data
pt <- powerTransform(data_th_fc$test, family = "yjPower")
# Apply the transformation to test data using that lambda
data_th_fc$test_yj <- yjPower(data_th_fc$test, pt$lambda)
# Apply the transformation to retest data using same lambda
data_th_fc$retest_yj <- yjPower(data_th_fc$retest, pt$lambda)
# calculate average baseline on test sessions
FC_test_avg_yj <- data_th_fc %>% 
  summarise(val = mean(test_yj, na.rm = TRUE))
# centre the baseline around the average
d_eeg_1 <- data_th_fc %>%
  mutate(FC_test_avg_val = FC_test_avg_yj) %>% 
  mutate(FC_Bs_centred_yj = test_yj - FC_test_avg_val$val)
# calculate change scores
data_th_fc_yj <-d_eeg_1 %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(fc_diff_yj = retest_yj - test_yj) %>%
  ungroup() 


par(mfrow = c(4, 2))
hist(data_th_fc$test, main = "Original test")
hist(data_th_fc_yj$test_yj, main = "Transformed test")
hist(data_th_fc$retest, main = "Original retest")
hist(data_th_fc_yj$retest_yj, main = "Transformed retest")
hist(data_th_fc$fc_diff, main = "Original change")
hist(data_th_fc_yj$fc_diff_yj, main = "Transformed change")
hist(data_th_fc$FC_Bs_centred, main = "Original centred baseline")
hist(data_th_fc_yj$FC_Bs_centred_yj, main = "Transformed centred baseline")

# clear up
rm(pt, FC_test_avg_yj, d_eeg_1)

# Run statistical model
m_int_yj <- lmerTest::lmer(
  fc_diff_yj ~
    (treatment_group_code * FC_Bs_centred_yj) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_th_fc_yj,
  control = lmerControl(optimizer = "bobyqa")
)
# model failed to converge, so removing nested effect
m_int_yj <- lmerTest::lmer(
  fc_diff_yj ~
    (treatment_group_code * FC_Bs_centred_yj) + (Age_patient_v0) +
    (Condition) + (1 | Regions) + (1 | Participant_ID),
  data = data_th_fc_yj,
  control = lmerControl(optimizer = "bobyqa")
)
# model failed to converge, so removing nested effect
m_int_yj <- lmerTest::lmer(
  fc_diff_yj ~
    (treatment_group_code * FC_Bs_centred_yj) + (Age_patient_v0) +
    (Condition) + (Regions) + (1 | Participant_ID),
  data = data_th_fc_yj,
  control = lmerControl(optimizer = "bobyqa")
)

summary(m_int_yj)
# create table for reporting
fc_theta_table_yj <- lmer_FCresults_table(m_int_yj)

# Check the sample size
ppts <- (rownames(coef(m_int_yj)$Participant_ID))
ppts_incl <- data_th_fc %>%
  filter(Participant_ID %in% ppts) %>%
  distinct(Participant_ID, treatment_group_code) %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl, ppts)

# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(fc_diff_yj ~ (Age_patient_v0) + (Condition) + (Regions) + (1 | Participant_ID), 
                     data = data_th_fc_yj, control = lmerControl(optimizer = "bobyqa")
)
# Now calculate partial residuals
data_th_fc_yj$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_th_fc_yj %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_th_fc_yj$FC_Bs_centred_yj, na.rm = TRUE), max(data_th_fc_yj$FC_Bs_centred_yj, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int_yj, specs = ~ treatment_group_code, var = "FC_Bs_centred_yj")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int_yj, specs = ~ treatment_group_code, 
               at = list(FC_Bs_centred_yj = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(FC_Bs_centred_yj = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (FC_Bs_centred_yj.trend * FC_Bs_centred_yj))
# For individual points: average within participant
participant_averages <- data_th_fc_yj %>%
  group_by(Participant_ID) %>%
  reframe(
    FC_Bs_centred_yj_avg = mean(FC_Bs_centred_yj, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = FC_Bs_centred_yj_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = FC_Bs_centred_yj, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_theta_FCyj_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_th_fc_yj, 
             aes(x = FC_Bs_centred_yj, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = FC_Bs_centred_yj, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_theta_FCyj_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)


# Check the model assumptions
plot_ass <- check_lmer_assumptions(m_int_yj, data_th_fc_yj, title = "Model assumptions for theta connectivity (YJt)")
ggsave("Theta_FCyj_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
fc_theta_VIFtable_yj <- lmer_VIF_table(m_int_yj)

# Clear up
rm(nb_vals, slopes, slopes_df, lines_df, plot, plot_ass)
rm(data_th_fc_yj, m_int_yj)
rm(data_th_fc)







#### 3b) Alpha band functional connectivity ####
d_eeg <- read.csv(
  "RCT1_FC_longformat_CxRxFB.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_FC_longformat_CxRxFB.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Alpha") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    FC = ifelse(is.nan(FC), NA, FC),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Frequencies = as.factor(Frequencies),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 90)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, FC) %>%
  dcast(Participant_ID + Condition + Regions ~ Sessions, value.var = "FC")
# calculate average baseline on test sessions
FC_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(FC, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(FC_test_avg_val = FC_test_avg) %>% 
  mutate(FC_Bs_centred = test - FC_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(fc_diff = retest - test) %>%
  ungroup() 
#merge with clinical data
d_cl <- read_csv(
  "xxx",
  col_types = cols(Participant_ID = col_character())
) %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    treatment_group_code = as.factor(as.character(treatment_group_code)),
    SiteNumber = as.factor(as.character(SiteNumber)),
    treatment_group_code = factor(treatment_group_code, levels = c("2", "1"))
  )
d_eeg_cl <- dplyr::left_join(d_filt, d_cl, by = "Participant_ID")
#select variables for model
fixed_cols <- c(
  "test", "retest", "fc_diff", "treatment_group_code", "Condition", "Regions", "FC_Bs_centred", "Age_patient_v0"
)
data_al_fc <- d_eeg_cl[complete.cases(d_eeg_cl[, fixed_cols]), ]
# clear up
rm(d_cl, d_eeg, d_eeg_trt, d_eeg_cl, d_eeg_filt, d_filt, FC_test_avg, fixed_cols)

# Run statistical model
m_int <- lmerTest::lmer(
  fc_diff ~
    (treatment_group_code * FC_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_al_fc,
  control = lmerControl(optimizer = "bobyqa")
)
# model failed to converge, so removing nested effect
m_int <- lmerTest::lmer(
  fc_diff ~
    (treatment_group_code * FC_Bs_centred) + (Age_patient_v0) +
    (Condition) + (1 | Regions) + (1 | Participant_ID),
  data = data_al_fc,
  control = lmerControl(optimizer = "bobyqa")
)

summary(m_int)
# create table for reporting
fc_alpha_table <- lmer_FCresults_table(m_int)





# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(fc_diff ~ (Age_patient_v0) +
                       (Condition) + (1 + Regions | Participant_ID), 
                     data = data_al_fc, control = lmerControl(optimizer = "bobyqa")
)
# Now calculate partial residuals
data_al_fc$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_al_fc %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_al_fc$FC_Bs_centred, na.rm = TRUE), max(data_al_fc$FC_Bs_centred, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int, specs = ~ treatment_group_code, var = "FC_Bs_centred")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int, specs = ~ treatment_group_code, 
               at = list(FC_Bs_centred = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(FC_Bs_centred = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (FC_Bs_centred.trend * FC_Bs_centred))
# For individual points: average within participant
participant_averages <- data_al_fc %>%
  group_by(Participant_ID) %>%
  reframe(
    FC_Bs_centred_avg = mean(FC_Bs_centred, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = FC_Bs_centred_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = FC_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_alpha_FC_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_al_fc, 
             aes(x = FC_Bs_centred, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = FC_Bs_centred, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_alpha_FC_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)

# Check the model assumptions
plot_ass <- check_lmer_assumptions(m_int, data_al_fc, title = "Model assumptions for alpha connectivity")
ggsave("Alpha_FC_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
fc_alpha_VIFtable <- lmer_VIF_table(m_int)

# Clear up
rm(nb_vals, slopes, slopes_df, lines_df, plot, plot_ass)
rm(m_int)



# On Yeo-Johnson transformed data
# transform FC values using the Yeo-Johnson transform
# Estimate the optimal lambda for test data
pt <- powerTransform(data_al_fc$test, family = "yjPower")
# Apply the transformation to test data using that lambda
data_al_fc$test_yj <- yjPower(data_al_fc$test, pt$lambda)
# Apply the transformation to retest data using same lambda
data_al_fc$retest_yj <- yjPower(data_al_fc$retest, pt$lambda)
# calculate average baseline on test sessions
FC_test_avg_yj <- data_al_fc %>% 
  summarise(val = mean(test_yj, na.rm = TRUE))
# centre the baseline around the average
d_eeg_1 <- data_al_fc %>%
  mutate(FC_test_avg_val = FC_test_avg_yj) %>% 
  mutate(FC_Bs_centred_yj = test_yj - FC_test_avg_val$val)
# calculate change scores
data_al_fc_yj <-d_eeg_1 %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(fc_diff_yj = retest_yj - test_yj) %>%
  ungroup() 


par(mfrow = c(4, 2))
hist(data_al_fc$test, main = "Original Al FC test")
hist(data_al_fc_yj$test_yj, main = "Transformed Al FC test")
hist(data_al_fc$retest, main = "Original Al FC retest")
hist(data_al_fc_yj$retest_yj, main = "Transformed Al FC retest")
hist(data_al_fc$fc_diff, main = "Original Al FC change")
hist(data_al_fc_yj$fc_diff_yj, main = "Transformed Al FC change")
hist(data_al_fc$FC_Bs_centred, main = "Original Al FC centred baseline")
hist(data_al_fc_yj$FC_Bs_centred_yj, main = "Transformed Al FC centred baseline")


# clear up
rm(pt, FC_test_avg_yj, d_eeg_1)

# Run statistical model
m_int_yj <- lmerTest::lmer(
  fc_diff_yj ~
    (treatment_group_code * FC_Bs_centred_yj) + (Age_patient_v0) +
    (Condition) + (1 + Regions | Participant_ID),
  data = data_al_fc_yj,
  control = lmerControl(optimizer = "bobyqa")
)
# model failed to converge, so removing nested effect
m_int_yj <- lmerTest::lmer(
  fc_diff_yj ~
    (treatment_group_code * FC_Bs_centred_yj) + (Age_patient_v0) +
    (Condition) + (1 | Regions) + (1 | Participant_ID),
  data = data_al_fc_yj,
  control = lmerControl(optimizer = "bobyqa")
)

summary(m_int_yj)
# create table for reporting
fc_alpha_table_yj <- lmer_FCresults_table(m_int_yj)

# Check the sample size
ppts <- (rownames(coef(m_int_yj)$Participant_ID))
ppts_incl <- data_al_fc %>%
  filter(Participant_ID %in% ppts) %>%
  distinct(Participant_ID, treatment_group_code) %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl, ppts)


# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model without including the interaction term
mplot_no_int <- lmer(fc_diff_yj ~ (Age_patient_v0) + (Condition) + (1 | Regions) + (1 | Participant_ID), 
                     data = data_al_fc_yj, control = lmerControl(optimizer = "bobyqa")
)
# Now calculate partial residuals
data_al_fc_yj$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_al_fc_yj %>%
  mutate(treatment_group_code = case_match(as.character(treatment_group_code),
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_al_fc_yj$FC_Bs_centred_yj, na.rm = TRUE), max(data_al_fc_yj$FC_Bs_centred_yj, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int_yj, specs = ~ treatment_group_code, var = "FC_Bs_centred_yj")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Get estimated marginal mean at mean of covariates, per group
emm <- emmeans(m_int_yj, specs = ~ treatment_group_code, 
               at = list(FC_Bs_centred_yj = 0))  # centred = 0 means at the mean
emm_df <- as.data.frame(emm)
# Join intercept per group onto lines_df
lines_df <- slopes_df %>%
  tidyr::expand_grid(FC_Bs_centred_yj = nb_vals) %>%
  left_join(emm_df %>% select(treatment_group_code, emmean), 
            by = "treatment_group_code") %>%
  mutate(predicted = emmean + (FC_Bs_centred_yj.trend * FC_Bs_centred_yj))
# For individual points: average within participant
participant_averages <- data_al_fc_yj %>%
  group_by(Participant_ID) %>%
  reframe(
    FC_Bs_centred_yj_avg = mean(FC_Bs_centred_yj, na.rm = TRUE),
    partial_residuals_avg = mean(partial_residuals, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0)) %>%
  distinct()
# Plot with participant averages
ggplot() +
  geom_point(data = participant_averages, 
             aes(x = FC_Bs_centred_yj_avg, y = partial_residuals_avg, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = FC_Bs_centred_yj, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10_alpha_FCyj_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 

# Plot with all individual points
ggplot() +
  geom_point(data = data_al_fc_yj, 
             aes(x = FC_Bs_centred_yj, y = partial_residuals, 
                 color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, 
            aes(x = FC_Bs_centred_yj, y = predicted, 
                color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(x = "EEG at baseline",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V10sm_alpha_FCyj_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)
rm(emm, emm_df, participant_averages)





# Visualise the effects
# 1) Interaction effect: interaction model with partial residuals
# Fit model on complete data
mplot_no_int <- lmer(fc_diff_yj ~ (Age_patient_v0) + (Condition) + (1 | Regions) + (1 | Participant_ID), 
                     data = data_al_fc_yj, control = lmerControl(optimizer = "bobyqa")
)
# Now calculate partial residuals
data_al_fc_yj$partial_residuals <- residuals(mplot_no_int) + fitted(mplot_no_int)
data_al_fc_yj %>%
  mutate(treatment_group_code = case_match(treatment_group_code,
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Range of mean centred test values for plotting
nb_vals <- seq(min(data_al_fc_yj$FC_Bs_centred_yj, na.rm = TRUE), max(data_al_fc_yj$FC_Bs_centred_yj, na.rm = TRUE), length.out = 100)
# Get slope (trend) of Norm_Bs by group 
slopes <- emtrends(m_int_yj, specs = ~ treatment_group_code, var = "FC_Bs_centred_yj")
# Convert to dataframe
slopes_df <- as.data.frame(slopes)
# Expand over mean-centred test values
lines_df <- slopes_df %>%
  tidyr::expand_grid(FC_Bs_centred_yj = nb_vals) %>%
  mutate(predicted = FC_Bs_centred_yj.trend * FC_Bs_centred_yj)
lines_df %>%
  mutate(treatment_group_code = case_match(treatment_group_code,
                                           "1" ~ "Arbaclofen",
                                           "2" ~ "Placebo"))
# Plot the datasets together
ggplot() +
  geom_point(data = data_al_fc_yj, aes(x = FC_Bs_centred_yj, y = partial_residuals, color = treatment_group_code), 
             alpha = 0.6, size = 1.5, show.legend = FALSE) +
  geom_line(data = lines_df, aes(x = FC_Bs_centred_yj, y = predicted, color = treatment_group_code, linetype = treatment_group_code),
            linewidth = 0.5, alpha = 0.6) +
  labs(title = "Alpha frequency functional connectivity (Yeo-Johnson transformed)",
       x = "Baseline EEG",
       y = "Change in EEG (adjusted)",
       color = "Treatment Group",
       linetype = "Treatment Group"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold")) +
  scale_color_manual(
    name = "Treatment Group",
    values = c("1" = "cornflowerblue", "2" = "darkorange2"),
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  ) +
  scale_linetype_manual(
    name = "Treatment Group",  
    values = c("1" = "dashed", "2" = "solid"),  
    labels = c("1" = "Arbaclofen", "2" = "Placebo")
  )

ggsave("V8_Al_FCyj_by_Group.jpg", plot = last_plot(), width = 7, height = 5, dpi = 300, bg = "white") 
rm(mplot_no_int, slopes, slopes_df, lines_df)


# Check the model assumptions
plot_ass <- check_lmer_assumptions(m_int_yj, data_al_fc_yj, title = "Model assumptions for alpha connectivity (YJt)")
ggsave("Alpha_FCyj_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 
# VIF table to check collinearity
fc_alpha_VIFtable_yj <- lmer_VIF_table(m_int_yj)


# Clear up
rm(nb_vals, slopes, slopes_df, lines_df, plot, plot_ass)
rm(data_al_fc_yj, m_int_yj)



rm(check_lmer_assumptions)

#### Save tables into Word document ####
doc <- read_docx() %>%
  body_add_par("Table 1.1: Theta band connectivity (Yeo-Johnson transformed) - model output", style = "heading 1") %>%
  body_add_flextable(fc_theta_table_yj) %>%
  body_add_par("") %>%
  body_add_par("Table 1.2: Theta band connectivity (Yeo-Johnson transformed) - VIF table", style = "heading 2") %>%
  body_add_flextable(fc_theta_VIFtable_yj) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 2.1: Alpha band connectivity (Yeo-Johnson transformed) - model output", style = "heading 1") %>%
  body_add_flextable(fc_alpha_table_yj) %>%
  body_add_par("") %>%
  body_add_par("Table 2.2: Alpha band connectivity (Yeo-Johnson transformed) - VIF table", style = "heading 2") %>%
  body_add_flextable(fc_alpha_VIFtable_yj) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 3.1: Theta band connectivity - model output", style = "heading 1") %>%
  body_add_flextable(fc_theta_table) %>%
  body_add_par("") %>%
  body_add_par("Table 3.2: Theta band connectivity - VIF table", style = "heading 2") %>%
  body_add_flextable(fc_theta_VIFtable) %>%
  body_add_break() %>%  # Page break
  body_add_par("Table 4.1: Alpha band connectivity - model output", style = "heading 1") %>%
  body_add_flextable(fc_alpha_table) %>%
  body_add_par("") %>%
  body_add_par("Table 4.2: Alpha band connectivity - VIF table", style = "heading 2") %>%
  body_add_flextable(fc_alpha_VIFtable) 

print(doc, target = "FC_tables_ScientificNotation.docx")

# Clear up
rm(fc_theta_table, fc_theta_VIFtable)
rm(fc_alpha_table, fc_alpha_VIFtable)  
rm(doc)



#### 4) Clinical prediction ####
# # Clinical prediction data
clinical<-read_excel("xxx/rct1_clinical_20240619_triple_blind_ForEmily.xlsx")

clinical_sublong <- clinical[, c("IDcode", "Age_patient_v0", "gender_patient_v0", "ITT", "PERPROTOCOL", "treatment_group_code", "SiteNumber",
                                 "FSIQ", "ADOS2_CSS_Total", "vi3_v7_com_ss", "vi3_v1_soc_ss", "vi3_v1_com_ss", "vi3_v7_soc_ss", "SRSParentv1_T_Total", "SRSParentv7_T_Total", "AIM_v1_Total", "AIM_v7_Total", 
                                 "vi3_v1_ipr_gsv","vi3_v1_pla_gsv","vi3_v1_cop_gsv", "vi3_v7_ipr_gsv","vi3_v7_pla_gsv","vi3_v7_cop_gsv",
                                 "ABCC_Leth_Soc_Withdr_v1", "ABCC_Leth_Soc_Withdr_v7")]
clinical_sublong$IDcode<-as.character(clinical_sublong$IDcode)
clinical_sublong$vi3_v7_com_ss<-as.numeric(clinical_sublong$vi3_v7_com_ss)
clinical_sublong$vi3_v1_com_ss<-as.numeric(clinical_sublong$vi3_v1_com_ss)
clinical_sublong$vi3_v7_soc_ss<-as.numeric(clinical_sublong$vi3_v7_soc_ss)
clinical_sublong$vi3_v1_soc_ss<-as.numeric(clinical_sublong$vi3_v1_soc_ss)
clinical_sublong$SRSParentv1_T_Total<-as.numeric(clinical_sublong$SRSParentv1_T_Total)
clinical_sublong$SRSParentv7_T_Total<-as.numeric(clinical_sublong$SRSParentv7_T_Total)
clinical_sublong$AIM_v1_Total<-as.numeric(clinical_sublong$AIM_v1_Total)
clinical_sublong$AIM_v7_Total<-as.numeric(clinical_sublong$AIM_v7_Total)
clinical_sublong$ABCC_Leth_Soc_Withdr_v1<-as.numeric(clinical_sublong$ABCC_Leth_Soc_Withdr_v1)
clinical_sublong$ABCC_Leth_Soc_Withdr_v7<-as.numeric(clinical_sublong$ABCC_Leth_Soc_Withdr_v7)

clinical_sublong$vi3_v1_ipr_gsv<-as.numeric(clinical_sublong$vi3_v1_ipr_gsv)
clinical_sublong$vi3_v1_pla_gsv<-as.numeric(clinical_sublong$vi3_v1_pla_gsv)
clinical_sublong$vi3_v1_cop_gsv<-as.numeric(clinical_sublong$vi3_v1_cop_gsv)
clinical_sublong$vi3_v7_ipr_gsv<-as.numeric(clinical_sublong$vi3_v7_ipr_gsv)
clinical_sublong$vi3_v7_pla_gsv<-as.numeric(clinical_sublong$vi3_v7_pla_gsv)
clinical_sublong$vi3_v7_cop_gsv<-as.numeric(clinical_sublong$vi3_v7_cop_gsv)

# check values in histograms
clinical_sublong %>%
  select_if(is.numeric) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")

# Change 999 to NA 
clinical_sublong <- clinical_sublong %>% mutate(SRSParentv1_T_Total = na_if(SRSParentv1_T_Total, 999))
clinical_sublong <- clinical_sublong %>% mutate(SRSParentv7_T_Total = na_if(SRSParentv7_T_Total, 999))
clinical_sublong <- clinical_sublong %>% mutate(vi3_v7_com_ss = na_if(vi3_v7_com_ss, 999))
clinical_sublong <- clinical_sublong %>% mutate(vi3_v7_soc_ss = na_if(vi3_v7_soc_ss, 999))
clinical_sublong <- clinical_sublong %>% mutate(vi3_v7_cop_gsv = na_if(vi3_v7_cop_gsv, 999))
clinical_sublong <- clinical_sublong %>% mutate(vi3_v7_ipr_gsv = na_if(vi3_v7_ipr_gsv, 999))
clinical_sublong <- clinical_sublong %>% mutate(vi3_v7_pla_gsv = na_if(vi3_v7_pla_gsv, 999))
clinical_sublong <- clinical_sublong %>% mutate(AIM_v1_Total = na_if(AIM_v1_Total, 999))
clinical_sublong <- clinical_sublong %>% mutate(AIM_v7_Total = na_if(AIM_v7_Total, 999))
clinical_sublong <- clinical_sublong %>% mutate(ABCC_Leth_Soc_Withdr_v1 = na_if(ABCC_Leth_Soc_Withdr_v1, 999))
clinical_sublong <- clinical_sublong %>% mutate(ABCC_Leth_Soc_Withdr_v7 = na_if(ABCC_Leth_Soc_Withdr_v7, 999))


# check values in histograms to see if NA worked
clinical_sublong %>%
  select_if(is.numeric) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")

# Calculate difference scores
clinical_sublong$vi3_diff_com_ss<-clinical_sublong$vi3_v7_com_ss-clinical_sublong$vi3_v1_com_ss
clinical_sublong$vi3_diff_soc_ss<-clinical_sublong$vi3_v7_soc_ss-clinical_sublong$vi3_v1_soc_ss
clinical_sublong$AIM_diff_Total<-clinical_sublong$AIM_v7_Total-clinical_sublong$AIM_v1_Total
clinical_sublong$SRSParent_diff_T_Total<-clinical_sublong$SRSParentv7_T_Total-clinical_sublong$SRSParentv1_T_Total
clinical_sublong$vi3_diff_ipr_gsv<-clinical_sublong$vi3_v7_ipr_gsv-clinical_sublong$vi3_v1_ipr_gsv
clinical_sublong$vi3_diff_pla_gsv<-clinical_sublong$vi3_v7_pla_gsv-clinical_sublong$vi3_v1_pla_gsv
clinical_sublong$vi3_diff_cop_gsv<-clinical_sublong$vi3_v7_cop_gsv-clinical_sublong$vi3_v1_cop_gsv
clinical_sublong$ABCC_Leth_Soc_Withdr_diff<-clinical_sublong$ABCC_Leth_Soc_Withdr_v7-clinical_sublong$ABCC_Leth_Soc_Withdr_v1

# copy values from IDcode to new column named Participant_ID
clinical_sublong$Participant_ID <- clinical_sublong$IDcode


#### 4a) Alpha log power - EEG baseline as predictor ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Alpha") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(ap_diff = retest - test) %>%
  ungroup() 
# get ppts included in change EEG analyses
load("IDs_incl_LFreq_abspow.RData")
d_eeg_filt_ids <- d_filt %>%
  filter(Participant_ID %in% ppts)
#merge with clinical data
d_eeg_cl <- dplyr::left_join(d_eeg_filt_ids, clinical_sublong, by = "Participant_ID")
# clear up
rm(d_eeg, d_eeg_trt, d_eeg_filt, LPow_test_avg, d_eeg_filt_ids, d_filt)

#### 4a.1) AIM ####
# set AIM_v7 values to nan if there is APow_Bs_centred it nan
d_eeg_cl_AIM <- d_eeg_cl 
d_eeg_cl_AIM$AIM_v7_Total[is.na(d_eeg_cl_AIM$APow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "AIM_v7_Total", "AIM_v1_Total", "AIM_diff_Total", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_AIM2 <- d_eeg_cl_AIM[complete.cases(d_eeg_cl_AIM[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_AIM)

participant_averages <- d_eeg_cl_AIM2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    AIM_diff_Total=first(AIM_diff_Total),
    AIM_v1_Total=first(AIM_v1_Total),
    AIM_v7_Total=first(AIM_v7_Total)) %>%
  distinct()

# Run interaction model
m_int <- lm(AIM_diff_Total ~ AIM_v1_Total + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Al_lpow_AIM_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for alpha log power predicting AIM change")
ggsave("ClinCh_Alpha_lpow_AIM_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

# clean up
rm(m_int, d_eeg_cl_AIM2, participant_averages, plot_ass)

#### 4a.2) VABS-3 Soc ss ####
# set clinical values of interest to nan if there is APow_Bs_centred it nan
d_eeg_cl_VABS <- d_eeg_cl 
d_eeg_cl_VABS$vi3_v7_soc_ss[is.na(d_eeg_cl_VABS$APow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "vi3_v7_soc_ss", "vi3_v1_soc_ss", "vi3_diff_soc_ss", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_VABS2 <- d_eeg_cl_VABS[complete.cases(d_eeg_cl_VABS[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_VABS)
participant_averages <- d_eeg_cl_VABS2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    vi3_diff_soc_ss=first(vi3_diff_soc_ss),
    vi3_v1_soc_ss=first(vi3_v1_soc_ss),
    vi3_v7_soc_ss=first(vi3_v7_soc_ss)) %>%
  distinct()

# Run interaction model
m_int <- lm(vi3_diff_soc_ss ~ vi3_v1_soc_ss + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# Check the sample size
ppts_incl <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl)
# create table for reporting
Al_lpow_VABSsocss_table <- lm_results_table(m_int)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for alpha log power predicting VABS-3 Soc change")
ggsave("ClinCh_Alpha_lpow_VABS3Soc_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

rm(m_int, d_eeg_cl_VABS2, participant_averages, plot_ass)

#### 4a.3) SRS-2 Total, parent report ####
# set clinical values of interest to nan if there is APow_Bs_centred it nan
d_eeg_cl_SRS <- d_eeg_cl 
d_eeg_cl_SRS$SRSParentv7_T_Total[is.na(d_eeg_cl_SRS$APow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "SRSParentv7_T_Total", "SRSParentv1_T_Total", "SRSParent_diff_T_Total", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_SRS2 <- d_eeg_cl_SRS[complete.cases(d_eeg_cl_SRS[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_SRS)

participant_averages <- d_eeg_cl_SRS2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    SRSParent_diff_T_Total=first(SRSParent_diff_T_Total),
    SRSParentv1_T_Total=first(SRSParentv1_T_Total),
    SRSParentv7_T_Total=first(SRSParentv7_T_Total)) %>%
  distinct()

# Run interaction model
m_int <- lm(SRSParent_diff_T_Total ~ SRSParentv1_T_Total + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# Check the sample size
ppts_incl <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl)
# create table for reporting
Al_lpow_SRS2tot_table <- lm_results_table(m_int)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for alpha log power predicting SRS-2 Total change")
ggsave("ClinCh_Alpha_lpow_SRS2Total_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

rm(m_int, d_eeg_cl_SRS2, participant_averages, plot_ass)



#### 4b) Alpha log power - EEG change as predictor ####
#### 4b.1) AIM ####
# set AIM_v7 values to nan if there is ap_diff it nan
d_eeg_cl_AIM <- d_eeg_cl 
d_eeg_cl_AIM$AIM_v7_Total[is.na(d_eeg_cl_AIM$ap_diff)] <- NA
#select variables for model
fixed_cols <- c(
  "AIM_v7_Total", "AIM_v1_Total", "AIM_diff_Total", "treatment_group_code", "ap_diff", "Age_patient_v0"
)
d_eeg_cl_AIM2 <- d_eeg_cl_AIM[complete.cases(d_eeg_cl_AIM[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_AIM)

participant_averages <- d_eeg_cl_AIM2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Ap_change_avg = mean(ap_diff, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    AIM_diff_Total=first(AIM_diff_Total),
    AIM_v1_Total=first(AIM_v1_Total),
    AIM_v7_Total=first(AIM_v7_Total)) %>%
  distinct()

# Run interaction model
m_int <- lm(AIM_diff_Total ~ AIM_v1_Total + Ap_change_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Change_Al_lpow_AIM_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# clean up
rm(m_int, d_eeg_cl_AIM2, participant_averages)

#### 4b.2) VABS-3 Soc ss ####
# set clinical values of interest to nan if there is APow_Bs_centred it nan
d_eeg_cl_VABS <- d_eeg_cl 
d_eeg_cl_VABS$vi3_v7_soc_ss[is.na(d_eeg_cl_VABS$ap_diff)] <- NA
#select variables for model
fixed_cols <- c(
  "vi3_v7_soc_ss", "vi3_v1_soc_ss", "vi3_diff_soc_ss", "treatment_group_code", "ap_diff", "Age_patient_v0"
)
d_eeg_cl_VABS2 <- d_eeg_cl_VABS[complete.cases(d_eeg_cl_VABS[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_VABS)
participant_averages <- d_eeg_cl_VABS2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Ap_change_avg = mean(ap_diff, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    vi3_diff_soc_ss=first(vi3_diff_soc_ss),
    vi3_v1_soc_ss=first(vi3_v1_soc_ss),
    vi3_v7_soc_ss=first(vi3_v7_soc_ss)) %>%
  distinct()

# Run interaction model
m_int <- lm(vi3_diff_soc_ss ~ vi3_v1_soc_ss + Ap_change_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# Check the sample size
ppts_incl <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl)
# create table for reporting
Change_Al_lpow_VABSsocss_table <- lm_results_table(m_int)

rm(m_int, d_eeg_cl_VABS2, participant_averages)

#### 4a.3) SRS-2 Total, parent report ####
# set clinical values of interest to nan if there is APow_Bs_centred it nan
d_eeg_cl_SRS <- d_eeg_cl 
d_eeg_cl_SRS$SRSParentv7_T_Total[is.na(d_eeg_cl_SRS$ap_diff)] <- NA
#select variables for model
fixed_cols <- c(
  "SRSParentv7_T_Total", "SRSParentv1_T_Total", "SRSParent_diff_T_Total", "treatment_group_code", "ap_diff", "Age_patient_v0"
)
d_eeg_cl_SRS2 <- d_eeg_cl_SRS[complete.cases(d_eeg_cl_SRS[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_SRS)

participant_averages <- d_eeg_cl_SRS2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Ap_change_avg = mean(ap_diff, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    SRSParent_diff_T_Total=first(SRSParent_diff_T_Total),
    SRSParentv1_T_Total=first(SRSParentv1_T_Total),
    SRSParentv7_T_Total=first(SRSParentv7_T_Total)) %>%
  distinct()

# Run interaction model
m_int <- lm(SRSParent_diff_T_Total ~ SRSParentv1_T_Total + Ap_change_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# Check the sample size
ppts_incl <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl)
# create table for reporting
Change_Al_lpow_SRS2tot_table <- lm_results_table(m_int)

rm(m_int, d_eeg_cl_SRS2, participant_averages)












#### 4c) Beta log power - EEG baseline as predictor ####
# Get baseline values for beta power + average within ppts
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Beta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(bp_diff = retest - test) %>%
  ungroup() 
# get ppts included in change EEG analyses
load("IDs_incl_LFreq_abspow.RData")
d_eeg_filt_ids <- d_filt %>%
  filter(Participant_ID %in% ppts)
#merge with clinical data
d_eeg_cl <- dplyr::left_join(d_eeg_filt_ids, clinical_sublong, by = "Participant_ID")
# clear up
rm(d_eeg, d_eeg_trt, d_eeg_filt, LPow_test_avg, d_eeg_filt_ids, d_filt)

#### 4c.1) AIM ####
# set AIM_v7 values to nan if there is LPow_Bs_centred it nan
d_eeg_cl_AIM <- d_eeg_cl 
d_eeg_cl_AIM$AIM_v7_Total[is.na(d_eeg_cl_AIM$LPow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "AIM_v7_Total", "AIM_v1_Total", "AIM_diff_Total", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_AIM2 <- d_eeg_cl_AIM[complete.cases(d_eeg_cl_AIM[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_AIM)

participant_averages <- d_eeg_cl_AIM2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    AIM_diff_Total=first(AIM_diff_Total),
    AIM_v1_Total=first(AIM_v1_Total),
    AIM_v7_Total=first(AIM_v7_Total)) %>%
  distinct()

# Run interaction model
m_int <- lm(AIM_diff_Total ~ AIM_v1_Total + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)
# create table for reporting
B_lpow_AIM_table <- lm_results_table(m_int)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for beta log power predicting AIM change")
ggsave("ClinCh_B_lpow_AIM_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 




#### 4c.1.b) AIM follow-up in split groups ####
# Arbaclofen group
d_eeg_cl_AIM2_arb <- participant_averages %>%
  filter(treatment_group_code == 1)
# Run interaction model
m_int_fu <- lm(AIM_diff_Total ~ AIM_v1_Total + Lpow_Bs_centred_avg + scale(Age_patient_v0), 
               na.action = "na.omit", 
               data = d_eeg_cl_AIM2_arb)
summary(m_int_fu)
# create table for reporting
B_lpow_AIM_arb_table <- lm_results_table(m_int_fu)

# Placebo group
d_eeg_cl_AIM2_pla <- participant_averages %>%
  filter(treatment_group_code == 2)
# Run interaction model
m_int_fu2 <- lm(AIM_diff_Total ~ AIM_v1_Total + Lpow_Bs_centred_avg + scale(Age_patient_v0), 
                na.action = "na.omit", 
                data = d_eeg_cl_AIM2_pla)
summary(m_int_fu2)
# create table for reporting
B_lpow_AIM_pla_table <- lm_results_table(m_int_fu2)

# clean up
rm(m_int_fu, m_int_fu2, d_eeg_cl_AIM2, participant_averages)
rm(d_eeg_cl_AIM2_arb, d_eeg_cl_AIM2_pla)

#### 4c.2) VABS-3 Soc ss ####
d_eeg_cl_VABS <- d_eeg_cl 
#select variables for model
fixed_cols <- c(
  "vi3_v7_soc_ss", "vi3_v1_soc_ss", "vi3_diff_soc_ss", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_VABS2 <- d_eeg_cl_VABS[complete.cases(d_eeg_cl_VABS[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_VABS)

participant_averages <- d_eeg_cl_VABS2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    vi3_diff_soc_ss=first(vi3_diff_soc_ss),
    vi3_v1_soc_ss=first(vi3_v1_soc_ss),
    vi3_v7_soc_ss=first(vi3_v7_soc_ss)) %>%
  distinct()

# Run interaction model
m_int <- lm(vi3_diff_soc_ss ~ vi3_v1_soc_ss + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# Check the sample size
ppts_incl <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl)
# create table for reporting
B_lpow_VABSsocss_table <- lm_results_table(m_int)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for beta log power predicting VABS-3 Soc change")
ggsave("ClinCh_B_lpow_VABS3Soc_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

rm(m_int, d_eeg_cl_VABS2, participant_averages, plot_ass)

#### 4c.3 SRS-2 Total, parent report ####
d_eeg_cl_SRS <- d_eeg_cl 
#select variables for model
fixed_cols <- c(
  "SRSParentv7_T_Total", "SRSParentv1_T_Total", "SRSParent_diff_T_Total", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_SRS2 <- d_eeg_cl_SRS[complete.cases(d_eeg_cl_SRS[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_SRS)

participant_averages <- d_eeg_cl_SRS2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    SRSParent_diff_T_Total=first(SRSParent_diff_T_Total),
    SRSParentv1_T_Total=first(SRSParentv1_T_Total),
    SRSParentv7_T_Total=first(SRSParentv7_T_Total)) %>%
  distinct()

# Run interaction model
m_int <- lm(SRSParent_diff_T_Total ~ SRSParentv1_T_Total + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# Check the sample size
ppts_incl <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl)
# create table for reporting
B_lpow_SRS2tot_table <- lm_results_table(m_int)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for beta log power predicting SRS-2 Total change")
ggsave("ClinCh_B_lpow_SRS2Total_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

rm(m_int, d_eeg_cl_SRS2, participant_averages, plot_ass)




#### 4d) Beta log power - EEG change as predictor ####
#### 4d.1) AIM ####
d_eeg_cl_AIM <- d_eeg_cl 
d_eeg_cl_AIM$AIM_v7_Total[is.na(d_eeg_cl_AIM$bp_diff)] <- NA
#select variables for model
fixed_cols <- c(
  "AIM_v7_Total", "AIM_v1_Total", "AIM_diff_Total", "treatment_group_code", "bp_diff", "Age_patient_v0"
)
d_eeg_cl_AIM2 <- d_eeg_cl_AIM[complete.cases(d_eeg_cl_AIM[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_AIM)

participant_averages <- d_eeg_cl_AIM2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Bp_change_avg = mean(bp_diff, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    AIM_diff_Total=first(AIM_diff_Total),
    AIM_v1_Total=first(AIM_v1_Total),
    AIM_v7_Total=first(AIM_v7_Total)) %>%
  distinct()

# Run interaction model
m_int <- lm(AIM_diff_Total ~ AIM_v1_Total + Bp_change_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)
# create table for reporting
Change_B_lpow_AIM_table <- lm_results_table(m_int)

# clean up
rm(m_int, d_eeg_cl_AIM2, participant_averages, plot_ass)


#### 4d.2) VABS-3 Soc ss ####
d_eeg_cl_VABS <- d_eeg_cl 
#select variables for model
fixed_cols <- c(
  "vi3_v7_soc_ss", "vi3_v1_soc_ss", "vi3_diff_soc_ss", "treatment_group_code", "bp_diff", "Age_patient_v0"
)
d_eeg_cl_VABS2 <- d_eeg_cl_VABS[complete.cases(d_eeg_cl_VABS[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_VABS)

participant_averages <- d_eeg_cl_VABS2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Bp_change_avg = mean(bp_diff, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    vi3_diff_soc_ss=first(vi3_diff_soc_ss),
    vi3_v1_soc_ss=first(vi3_v1_soc_ss),
    vi3_v7_soc_ss=first(vi3_v7_soc_ss)) %>%
  distinct()

# Run interaction model
m_int <- lm(vi3_diff_soc_ss ~ vi3_v1_soc_ss + Bp_change_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# Check the sample size
ppts_incl <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl)
# create table for reporting
Change_B_lpow_VABSsocss_table <- lm_results_table(m_int)

rm(m_int, d_eeg_cl_VABS2, participant_averages)

#### 4d.3 SRS-2 Total, parent report ####
d_eeg_cl_SRS <- d_eeg_cl 
#select variables for model
fixed_cols <- c(
  "SRSParentv7_T_Total", "SRSParentv1_T_Total", "SRSParent_diff_T_Total", "treatment_group_code", "bp_diff", "Age_patient_v0"
)
d_eeg_cl_SRS2 <- d_eeg_cl_SRS[complete.cases(d_eeg_cl_SRS[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_SRS)

participant_averages <- d_eeg_cl_SRS2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Bp_change_avg = mean(bp_diff, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    SRSParent_diff_T_Total=first(SRSParent_diff_T_Total),
    SRSParentv1_T_Total=first(SRSParentv1_T_Total),
    SRSParentv7_T_Total=first(SRSParentv7_T_Total)) %>%
  distinct()

# Run interaction model
m_int <- lm(SRSParent_diff_T_Total ~ SRSParentv1_T_Total + Bp_change_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# Check the sample size
ppts_incl <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl
rm(ppts_incl)
# create table for reporting
Change_B_lpow_SRS2tot_table <- lm_results_table(m_int)

rm(m_int, d_eeg_cl_SRS2, participant_averages)

rm(d_eeg_cl, clinical, clinical_sublong)

#### Save tables into Word document ####
doc <- read_docx() %>%
  # Alpha log power - EEG baseline as predictor
  body_add_par("Table 1.1: Baseline alpha log power predicting change on AIM - model output", style = "heading 1") %>%
  body_add_flextable(Al_lpow_AIM_table) %>%
  body_add_par("") %>%
  body_add_par("Table 1.2: Baseline alpha log power predicting change on VABS Soc ss - model output", style = "heading 1") %>%
  body_add_flextable(Al_lpow_VABSsocss_table) %>%
  body_add_par("") %>%
  body_add_par("Table 1.3: Baseline alpha log power predicting change on SRS-2 Total - model output", style = "heading 1") %>%
  body_add_flextable(Al_lpow_SRS2tot_table) %>%
  body_add_break() %>%  
  # Alpha log power - EEG change as predictor
  body_add_par("Table 2.1: Change in alpha log power predicting change on AIM - model output", style = "heading 1") %>%
  body_add_flextable(Change_Al_lpow_AIM_table) %>%
  body_add_par("") %>%
  body_add_par("Table 2.2: Change in alpha log power predicting change on VABS Soc ss - model output", style = "heading 1") %>%
  body_add_flextable(Change_Al_lpow_VABSsocss_table) %>%
  body_add_par("") %>%
  body_add_par("Table 2.3: Change in alpha log power predicting change on SRS-2 Total - model output", style = "heading 1") %>%
  body_add_flextable(Change_Al_lpow_SRS2tot_table) %>%
  body_add_break() %>%  
  # Beta log power - EEG baseline as predictor
  body_add_par("Table 3.1: Baseline beta log power predicting change on AIM - model output", style = "heading 1") %>%
  body_add_flextable(B_lpow_AIM_table) %>%
  body_add_par("") %>%
  body_add_par("Table 3.1.a: Follow up in arbaclofen group: Baseline beta log power predicting change on AIM - model output", style = "heading 2") %>%
  body_add_flextable(B_lpow_AIM_arb_table) %>%
  body_add_par("") %>%
  body_add_par("Table 3.1.b: Follow up in arbaclofen group: Baseline beta log power predicting change on AIM - model output", style = "heading 2") %>%
  body_add_flextable(B_lpow_AIM_pla_table) %>%
  body_add_par("") %>%
  body_add_par("Table 3.2: Baseline beta log power predicting change on VABS Soc ss - model output", style = "heading 1") %>%
  body_add_flextable(B_lpow_VABSsocss_table) %>%
  body_add_par("") %>%
  body_add_par("Table 3.3: Baseline beta log power predicting change on SRS-2 Total - model output", style = "heading 1") %>%
  body_add_flextable(B_lpow_SRS2tot_table) %>%
  body_add_break() %>%
  # Beta log power - EEG baseline as predictor
  body_add_par("Table 4.1: Change in beta log power predicting change on AIM - model output", style = "heading 1") %>%
  body_add_flextable(Change_B_lpow_AIM_table) %>%
  body_add_par("") %>%
  body_add_par("Table 4.2: Change in beta log power predicting change on VABS Soc ss - model output", style = "heading 1") %>%
  body_add_flextable(Change_B_lpow_VABSsocss_table) %>%
  body_add_par("") %>%
  body_add_par("Table 4.3: Change in beta log power predicting change on SRS-2 Total - model output", style = "heading 1") %>%
  body_add_flextable(Change_B_lpow_SRS2tot_table)

print(doc, target = "ClinPred_tables.docx")

rm(Al_lpow_AIM_table, Al_lpow_VABSsocss_table, Al_lpow_SRS2tot_table)
rm(BLg_lpow_AIM_table, B_lpow_VABSsocss_table, B_lpow_SRS2tot_table)
rm(doc)





#### 5) Characterisation of the sample ####
# load data
clinical<-read_excel("xxx")

# get ppts included in change EEG analyses
load("IDs_incl_LFreq_abspow.RData")
# select needed variables
clinical <- clinical[, c("IDcode", "Age_patient_v0", "gender_patient_v0", "ITT", "PERPROTOCOL", "treatment_group_code", "SiteNumber",
                         "FSIQ", "NVIQ", "VIQ", "ADOS2_CSS_Total", "ADOS2_CSS_SA", "ADOS2_CSS_RRB")]
clinical$IDcode<-as.character(clinical$IDcode)
clinical$Participant_ID <- clinical$IDcode
clinical$Age_patient_v0<-as.numeric(clinical$Age_patient_v0)
clinical$FSIQ<-as.numeric(clinical$FSIQ)
clinical$NVIQ<-as.numeric(clinical$NVIQ)
clinical$VIQ<-as.numeric(clinical$VIQ)
clinical$ADOS2_CSS_SA<-as.numeric(clinical$ADOS2_CSS_SA)
clinical$ADOS2_CSS_RRB<-as.numeric(clinical$ADOS2_CSS_RRB)
clinical$ADOS2_CSS_Total<-as.numeric(clinical$ADOS2_CSS_Total)
# change 999 values in NA
clinical <- clinical %>% mutate(Age_patient_v0 = na_if(Age_patient_v0, 999))
clinical <- clinical %>% mutate(FSIQ = na_if(FSIQ, 999))
clinical <- clinical %>% mutate(NVIQ = na_if(NVIQ, 999))
clinical <- clinical %>% mutate(VIQ = na_if(VIQ, 999))
clinical <- clinical %>% mutate(ADOS2_CSS_SA = na_if(ADOS2_CSS_SA, 999))
clinical <- clinical %>% mutate(ADOS2_CSS_RRB = na_if(ADOS2_CSS_RRB, 999))
clinical <- clinical %>% mutate(ADOS2_CSS_Total = na_if(ADOS2_CSS_Total, 999))
# check values in histograms
clinical %>%
  select_if(is.numeric) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")
# included ppts
clinical <- clinical %>%
  filter(Participant_ID %in% ppts)
# recode group and gender from values to characters
clinical$treatment_group_code <- case_match(as.character(clinical$treatment_group_code),
                                         "1" ~ "Arbaclofen",
                                         "2" ~ "Placebo")
clinical$gender_patient_v0 <- case_match(as.character(clinical$gender_patient_v0),
                                            "1" ~ "Male",
                                            "2" ~ "Female")

# Get charactrisation values to report
# Group size
txgroup_incl <- clinical %>%
  filter(Participant_ID %in% ppts) %>%
  distinct(Participant_ID, treatment_group_code) %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
txgroup_incl
# Gender
gender_incl <- clinical %>%
  filter(Participant_ID %in% ppts) %>%
  distinct(Participant_ID, treatment_group_code, gender_patient_v0) %>%
  group_by(treatment_group_code, gender_patient_v0) %>%
  summarise(N=n())
gender_incl
# Age (years)
age_incl <- clinical %>%
  filter(Participant_ID %in% ppts) %>%
  group_by(treatment_group_code) %>%
  summarize(Age_patient_v0.mean = mean(Age_patient_v0, na.rm = TRUE),
            Age_patient_v0.sd = sd(Age_patient_v0, na.rm = TRUE),
            Age_patient_v0.min = min(Age_patient_v0, na.rm = TRUE),
            Age_patient_v0.max = max(Age_patient_v0, na.rm = TRUE),
            Age_patient_v0.N_NA = sum(is.na(Age_patient_v0)))
age_incl
# Full scale IQ
fsiq_incl <- clinical %>%
  filter(Participant_ID %in% ppts) %>%
  group_by(treatment_group_code) %>%
  summarize(FSIQ.mean = mean(FSIQ, na.rm = TRUE),
            FSIQ.sd = sd(FSIQ, na.rm = TRUE),
            FSIQ.min = min(FSIQ, na.rm = TRUE),
            FSIQ.max = max(FSIQ, na.rm = TRUE),
            FSIQ.N_NA = sum(is.na(FSIQ)))
fsiq_incl
# Verbal IQ
viq_incl <- clinical %>%
  filter(Participant_ID %in% ppts) %>%
  group_by(treatment_group_code) %>%
  summarize(VIQ.mean = mean(VIQ, na.rm = TRUE),
            VIQ.sd = sd(VIQ, na.rm = TRUE),
            VIQ.min = min(VIQ, na.rm = TRUE),
            VIQ.max = max(VIQ, na.rm = TRUE),
            VIQ.N_NA = sum(is.na(VIQ)))
viq_incl
# Non-verbal IQ
nviq_incl <- clinical %>%
  filter(Participant_ID %in% ppts) %>%
  group_by(treatment_group_code) %>%
  summarize(NVIQ.mean = mean(NVIQ, na.rm = TRUE),
            NVIQ.sd = sd(NVIQ, na.rm = TRUE),
            NVIQ.min = min(NVIQ, na.rm = TRUE),
            NVIQ.max = max(NVIQ, na.rm = TRUE),
            NVIQ.N_NA = sum(is.na(NVIQ)))
nviq_incl


# Updated ADOS file
# load data
clinical_ADOS<-read_excel("xxx.xlsx")
clinical_ADOS$Participant_ID<-as.character(clinical_ADOS$Participant_ID)
clinical_ADOS$ADOS2_CSS_SA<-as.numeric(clinical_ADOS$ADOS2_CSS_SA)
clinical_ADOS$ADOS2_CSS_RRB<-as.numeric(clinical_ADOS$ADOS2_CSS_RRB)
clinical_ADOS$ADOS2_CSS_Total<-as.numeric(clinical_ADOS$ADOS2_CSS_Total)
clinical_ADOS <- clinical_ADOS %>% mutate(ADOS2_CSS_SA = na_if(ADOS2_CSS_SA, 999))
clinical_ADOS <- clinical_ADOS %>% mutate(ADOS2_CSS_RRB = na_if(ADOS2_CSS_RRB, 999))
clinical_ADOS <- clinical_ADOS %>% mutate(ADOS2_CSS_Total = na_if(ADOS2_CSS_Total, 999))
clinical_joined <- dplyr::left_join(clinical, clinical_ADOS, by = "Participant_ID")

# ADOS - SA
ADOSSA_incl <- clinical_joined %>%
  filter(Participant_ID %in% ppts) %>%
  group_by(treatment_group_code) %>%
  summarize(ADOS2_CSS_SA.mean = mean(ADOS2_CSS_SA.y, na.rm = TRUE),
            ADOS2_CSS_SA.sd = sd(ADOS2_CSS_SA.y, na.rm = TRUE),
            ADOS2_CSS_SA.min = min(ADOS2_CSS_SA.y, na.rm = TRUE),
            ADOS2_CSS_SA.max = max(ADOS2_CSS_SA.y, na.rm = TRUE),
            ADOS2_CSS_SA.N_NA = sum(is.na(ADOS2_CSS_SA.y)))
ADOSSA_incl 
# ADOS - RRBs
ADOSRRBs_incl <- clinical_joined %>%
  filter(Participant_ID %in% ppts) %>%
  group_by(treatment_group_code) %>%
  summarize(ADOS2_CSS_RRB.mean = mean(ADOS2_CSS_RRB.y, na.rm = TRUE),
            ADOS2_CSS_RRB.sd = sd(ADOS2_CSS_RRB.y, na.rm = TRUE),
            ADOS2_CSS_RRB.min = min(ADOS2_CSS_RRB.y, na.rm = TRUE),
            ADOS2_CSS_RRB.max = max(ADOS2_CSS_RRB.y, na.rm = TRUE),
            ADOS2_CSS_RRB.N_NA = sum(is.na(ADOS2_CSS_RRB.y)))
ADOSRRBs_incl 
# ADOS - Total
ADOST_incl <- clinical_joined %>%
  filter(Participant_ID %in% ppts) %>%
  group_by(treatment_group_code) %>%
  summarize(ADOS2_CSS_Total.mean = mean(ADOS2_CSS_Total.y, na.rm = TRUE),
            ADOS2_CSS_Total.sd = sd(ADOS2_CSS_Total.y, na.rm = TRUE),
            ADOS2_CSS_Total.min = min(ADOS2_CSS_Total.y, na.rm = TRUE),
            ADOS2_CSS_Total.max = max(ADOS2_CSS_Total.y, na.rm = TRUE),
            ADOS2_CSS_Total.N_NA = sum(is.na(ADOS2_CSS_Total.y)))
ADOST_incl 

# SiteNumber
site_incl <- clinical %>%
  filter(Participant_ID %in% ppts) %>%
  distinct(Participant_ID, treatment_group_code, SiteNumber) %>%
  group_by(treatment_group_code, SiteNumber) %>%
  summarise(N=n())
site_incl

# Plots for inspection of distributions
# Create plots
# 0) Number
p0 <- ggplot(clinical, aes(x=as.factor(treatment_group_code), fill=as.factor(treatment_group_code) )) +  
  geom_bar(alpha = 0.5 ) +
  scale_fill_manual(values = c("cyan", "red"), 
                    labels = c("Arbaclofen", "Placebo")) +
  labs(title = "Treatment group",
       x = "Treatment goup",
       y = "Count",
       fill = "Treatment group") +
  theme_minimal()
# 1) Gender
p1 <- ggplot(clinical, aes(x=as.factor(treatment_group_code), fill=as.factor(gender_patient_v0) )) +  
  geom_bar(alpha = 0.9 ) +
  scale_fill_manual(values = c("cornflowerblue", "chartreuse3"), 
                    labels = c("Female", "Male")) +
  labs(title = "Gender",
       x = "Treatment group",
       y = "Count",
       fill = "Gender") +
  theme_minimal()
# 2) Age
p2 <- ggplot(clinical , aes(x = Age_patient_v0, fill = factor(treatment_group_code))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  scale_fill_manual(values = c("cyan", "red"), 
                    labels = c("Arbaclofen", "Placebo")) +
  labs(title = "Age",
       x = "Age (in years)",
       y = "Count",
       fill = "Treatment group") +
  theme_minimal()
# 3) FSIQ
p3 <- ggplot(clinical , aes(x = FSIQ, fill = factor(treatment_group_code))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  scale_fill_manual(values = c("cyan", "red"), 
                    labels = c("Arbaclofen", "Placebo")) +
  labs(title = "Full scale IQ",
       x = "FSIQ",
       y = "Count",
       fill = "Treatment group") +
  theme_minimal()
# 4) VIQ
p4 <- ggplot(clinical , aes(x = VIQ, fill = factor(treatment_group_code))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  scale_fill_manual(values = c("cyan", "red"), 
                    labels = c("Arbaclofen", "Placebo")) +
  labs(title = "Verbal IQ",
       x = "VIQ",
       y = "Count",
       fill = "Treatment group") +
  theme_minimal()
# 5) NVIQ
p5 <- ggplot(clinical , aes(x = NVIQ, fill = factor(treatment_group_code))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  scale_fill_manual(values = c("cyan", "red"), 
                    labels = c("Arbaclofen", "Placebo")) +
  labs(title = "Non-verbal IQ",
       x = "VIQ",
       y = "Count",
       fill = "Treatment group") +
  theme_minimal()
# 6) ADOS SA
p6 <- ggplot(clinical , aes(x = ADOS2_CSS_SA, fill = factor(treatment_group_code))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  scale_fill_manual(values = c("cyan", "red"), 
                    labels = c("Arbaclofen", "Placebo")) +
  labs(title = "ADOS-2 Social Affect",
       x = "SA CSS",
       y = "Count",
       fill = "Treatment group") +
  theme_minimal()
# 7) ADOS RRBs
p7 <- ggplot(clinical , aes(x = ADOS2_CSS_RRB, fill = factor(treatment_group_code))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  scale_fill_manual(values = c("cyan", "red"), 
                    labels = c("Arbaclofen", "Placebo")) +
  labs(title = "ADOS-2 Restricted & Repetitive Behaviours",
       x = "RRB CSS",
       y = "Count",
       fill = "Treatment group") +
  theme_minimal()
# 8) ADOS Total
p8 <- ggplot(clinical , aes(x = ADOS2_CSS_Total, fill = factor(treatment_group_code))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  scale_fill_manual(values = c("cyan", "red"), 
                    labels = c("Arbaclofen", "Placebo")) +
  labs(title = "ADOS-2 Total",
       x = "Total CSS",
       y = "Count",
       fill = "Treatment group") +
  theme_minimal()

# Combine plots
library(gridExtra)
library(grid)
plot_sampleinfo <- grid.arrange(p0, p1, p2, p3, p4, p5, p6, p7, p8, ncol = 3,
                                top = textGrob("EEG sample characteristics in AIMS-CT1 study", 
                                               gp = gpar(fontsize = 16, fontface = "bold")))
ggsave("EEGSampleCharacteristics.jpg", plot_sampleinfo, width = 12, height = 8, dpi = 300, bg = "white") 





#### 6) Exploratory analyses upon request from AIMS-2-TRIALS A-reps ####
# Check for associations with anxiety measures
# # Clinical prediction data
clinical<-read_excel("xxx")

clinical_sublong <- clinical[, c("IDcode", "Age_patient_v0", "gender_patient_v0", "ITT", "PERPROTOCOL", "treatment_group_code", "SiteNumber",
                                 "CBCL_AnxiousDepressed_v1_T", "CBCL_AnxiousDepressed_v7_T", 
                                 "CBCL_Internalization_v1_T","CBCL_Internalization_v7_T",
                                 "ABCC_Leth_Soc_Withdr_v1", "ABCC_Leth_Soc_Withdr_v7")]
clinical_sublong$IDcode<-as.character(clinical_sublong$IDcode)
clinical_sublong$CBCL_AnxiousDepressed_v1_T<-as.numeric(clinical_sublong$CBCL_AnxiousDepressed_v1_T)
clinical_sublong$CBCL_AnxiousDepressed_v7_T<-as.numeric(clinical_sublong$CBCL_AnxiousDepressed_v7_T)
clinical_sublong$CBCL_Internalization_v1_T<-as.numeric(clinical_sublong$CBCL_Internalization_v1_T)
clinical_sublong$CBCL_Internalization_v7_T<-as.numeric(clinical_sublong$CBCL_Internalization_v7_T)
clinical_sublong$ABCC_Leth_Soc_Withdr_v1<-as.numeric(clinical_sublong$ABCC_Leth_Soc_Withdr_v1)
clinical_sublong$ABCC_Leth_Soc_Withdr_v7<-as.numeric(clinical_sublong$ABCC_Leth_Soc_Withdr_v7)

# check values in histograms
clinical_sublong %>%
  select_if(is.numeric) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")

# Change 999 to NA 
clinical_sublong <- clinical_sublong %>% mutate(CBCL_AnxiousDepressed_v1_T = na_if(CBCL_AnxiousDepressed_v1_T, 999))
clinical_sublong <- clinical_sublong %>% mutate(CBCL_AnxiousDepressed_v7_T = na_if(CBCL_AnxiousDepressed_v7_T, 999))
clinical_sublong <- clinical_sublong %>% mutate(CBCL_Internalization_v1_T = na_if(CBCL_Internalization_v1_T, 999))
clinical_sublong <- clinical_sublong %>% mutate(CBCL_Internalization_v7_T = na_if(CBCL_Internalization_v7_T, 999))
clinical_sublong <- clinical_sublong %>% mutate(ABCC_Leth_Soc_Withdr_v1 = na_if(ABCC_Leth_Soc_Withdr_v1, 999))
clinical_sublong <- clinical_sublong %>% mutate(ABCC_Leth_Soc_Withdr_v7 = na_if(ABCC_Leth_Soc_Withdr_v7, 999))

# check values in histograms to see if NA worked
clinical_sublong %>%
  select_if(is.numeric) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")

# Calculate difference scores
clinical_sublong$CBCL_diff_anxdep<-clinical_sublong$CBCL_AnxiousDepressed_v7_T-clinical_sublong$CBCL_AnxiousDepressed_v1_T
clinical_sublong$CBCL_diff_int<-clinical_sublong$CBCL_Internalization_v7_T-clinical_sublong$CBCL_Internalization_v1_T
clinical_sublong$ABCC_Leth_Soc_Withdr_diff<-clinical_sublong$ABCC_Leth_Soc_Withdr_v7-clinical_sublong$ABCC_Leth_Soc_Withdr_v1

# copy values from IDcode to new column named Participant_ID
clinical_sublong$Participant_ID <- clinical_sublong$IDcode


#### 6a) Alpha log power - EEG baseline as predictor ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Alpha") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(ap_diff = retest - test) %>%
  ungroup() 
# get ppts included in change EEG analyses
load("IDs_incl_LFreq_abspow.RData")
d_eeg_filt_ids <- d_filt %>%
  filter(Participant_ID %in% ppts)
#merge with clinical data
d_eeg_cl <- dplyr::left_join(d_eeg_filt_ids, clinical_sublong, by = "Participant_ID")
# clear up
rm(d_eeg, d_eeg_trt, d_eeg_filt, LPow_test_avg, d_eeg_filt_ids, d_filt)

#### 6a.1) CBCL Anxiety / Depression ####
# set v7 values to nan if there is APow_Bs_centred it nan
d_eeg_cl_Cad <- d_eeg_cl 
d_eeg_cl_Cad$CBCL_AnxiousDepressed_v7_T[is.na(d_eeg_cl_Cad$LPow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "CBCL_AnxiousDepressed_v7_T", "CBCL_AnxiousDepressed_v1_T", "CBCL_diff_anxdep", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_Cad2 <- d_eeg_cl_Cad[complete.cases(d_eeg_cl_Cad[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_Cad)

participant_averages <- d_eeg_cl_Cad2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    CBCL_diff_anxdep=first(CBCL_diff_anxdep),
    CBCL_AnxiousDepressed_v1_T=first(CBCL_AnxiousDepressed_v1_T),
    CBCL_AnxiousDepressed_v7_T=first(CBCL_AnxiousDepressed_v7_T)) %>%
  distinct()

# Run interaction model
m_int <- lm(CBCL_diff_anxdep ~ CBCL_AnxiousDepressed_v1_T + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Al_lpow_Cad_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for alpha log power predicting CBCL Anx/Dep change")
ggsave("ClinCh_Alpha_lpow_CBCLAnxDep_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

# clean up
rm(m_int, d_eeg_cl_Cad2, participant_averages, plot_ass)

#### 6a.2) CBCL Internalizing ####
# set v7 values to nan if there is APow_Bs_centred it nan
d_eeg_cl_Cint <- d_eeg_cl 
d_eeg_cl_Cint$CBCL_AnxiousDepressed_v7_T[is.na(d_eeg_cl_Cint$LPow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "CBCL_Internalization_v7_T", "CBCL_Internalization_v1_T", "CBCL_diff_int", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_Cint2 <- d_eeg_cl_Cint[complete.cases(d_eeg_cl_Cint[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_Cint)

participant_averages <- d_eeg_cl_Cint2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    CBCL_diff_int=first(CBCL_diff_int),
    CBCL_Internalization_v1_T=first(CBCL_Internalization_v1_T),
    CBCL_Internalization_v7_T=first(CBCL_Internalization_v7_T)) %>%
  distinct()

# Run interaction model
m_int <- lm(CBCL_diff_int ~ CBCL_Internalization_v1_T + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Al_lpow_Cint_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for alpha log power predicting CBCL Int change")
ggsave("ClinCh_Alpha_lpow_CBCLInt_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

# clean up
rm(m_int, d_eeg_cl_Cint2, participant_averages, plot_ass)


#### 6b) Beta log power - EEG baseline as predictor ####
# Get baseline values for beta power + average within ppts
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Beta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(bp_diff = retest - test) %>%
  ungroup() 
# get ppts included in change EEG analyses
load("IDs_incl_LFreq_abspow.RData")
d_eeg_filt_ids <- d_filt %>%
  filter(Participant_ID %in% ppts)
#merge with clinical data
d_eeg_cl <- dplyr::left_join(d_eeg_filt_ids, clinical_sublong, by = "Participant_ID")
# clear up
rm(d_eeg, d_eeg_trt, d_eeg_filt, LPow_test_avg, d_eeg_filt_ids, d_filt)

#### 6b.1) CBCL Anxiety / Depression ####
# set v7 values to nan if there is APow_Bs_centred it nan
d_eeg_cl_Cad <- d_eeg_cl 
d_eeg_cl_Cad$CBCL_AnxiousDepressed_v7_T[is.na(d_eeg_cl_Cad$LPow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "CBCL_AnxiousDepressed_v7_T", "CBCL_AnxiousDepressed_v1_T", "CBCL_diff_anxdep", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_Cad2 <- d_eeg_cl_Cad[complete.cases(d_eeg_cl_Cad[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_Cad)

participant_averages <- d_eeg_cl_Cad2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    CBCL_diff_anxdep=first(CBCL_diff_anxdep),
    CBCL_AnxiousDepressed_v1_T=first(CBCL_AnxiousDepressed_v1_T),
    CBCL_AnxiousDepressed_v7_T=first(CBCL_AnxiousDepressed_v7_T)) %>%
  distinct()

# Run interaction model
m_int <- lm(CBCL_diff_anxdep ~ CBCL_AnxiousDepressed_v1_T + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Be_lpow_Cad_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for beta log power predicting CBCL Anx/Dep change")
ggsave("ClinCh_Beta_lpow_CBCLAnxDep_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

# clean up
rm(m_int, d_eeg_cl_Cad2, participant_averages, plot_ass)

#### 6b.2) CBCL Internalizing ####
# set v7 values to nan if there is APow_Bs_centred it nan
d_eeg_cl_Cint <- d_eeg_cl 
d_eeg_cl_Cint$CBCL_AnxiousDepressed_v7_T[is.na(d_eeg_cl_Cint$LPow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "CBCL_Internalization_v7_T", "CBCL_Internalization_v1_T", "CBCL_diff_int", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_Cint2 <- d_eeg_cl_Cint[complete.cases(d_eeg_cl_Cint[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_Cint)

participant_averages <- d_eeg_cl_Cint2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    CBCL_diff_int=first(CBCL_diff_int),
    CBCL_Internalization_v1_T=first(CBCL_Internalization_v1_T),
    CBCL_Internalization_v7_T=first(CBCL_Internalization_v7_T)) %>%
  distinct()

# Run interaction model
m_int <- lm(CBCL_diff_int ~ CBCL_Internalization_v1_T + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Be_lpow_Cint_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for beta log power predicting CBCL Int change")
ggsave("ClinCh_Beta_lpow_CBCLInt_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

# clean up
rm(m_int, d_eeg_cl_Cint2, participant_averages, plot_ass)
rm(d_eeg_cl, clinical, clinical_sublong)

#### 6.2) Check combining CBCL data from different age groups #### 
# Anxiety measures and EEG
clinical<-read_excel("xxx")

clinical_sublong <- clinical[, c("IDcode", "Age_patient_v0", "gender_patient_v0", "ITT", "PERPROTOCOL", "treatment_group_code", "SiteNumber",
                                 "CBCL_AnxiousDepressed_v1_T", "CBCL_AnxiousDepressed_v7_T", 
                                 "CBCL_Internalization_v1_T","CBCL_Internalization_v7_T",
                                 "ABCC_Leth_Soc_Withdr_v1", "ABCC_Leth_Soc_Withdr_v7")]
clinical_sublong$IDcode<-as.character(clinical_sublong$IDcode)
clinical_sublong$CBCL_AnxiousDepressed_v1_T<-as.numeric(clinical_sublong$CBCL_AnxiousDepressed_v1_T)
clinical_sublong$CBCL_AnxiousDepressed_v7_T<-as.numeric(clinical_sublong$CBCL_AnxiousDepressed_v7_T)
clinical_sublong$CBCL_Internalization_v1_T<-as.numeric(clinical_sublong$CBCL_Internalization_v1_T)
clinical_sublong$CBCL_Internalization_v7_T<-as.numeric(clinical_sublong$CBCL_Internalization_v7_T)
clinical_sublong$ABCC_Leth_Soc_Withdr_v1<-as.numeric(clinical_sublong$ABCC_Leth_Soc_Withdr_v1)
clinical_sublong$ABCC_Leth_Soc_Withdr_v7<-as.numeric(clinical_sublong$ABCC_Leth_Soc_Withdr_v7)

# check values in histograms
clinical_sublong %>%
  select_if(is.numeric) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")

# Change 999 to NA 
clinical_sublong <- clinical_sublong %>% mutate(CBCL_AnxiousDepressed_v1_T = na_if(CBCL_AnxiousDepressed_v1_T, 999))
clinical_sublong <- clinical_sublong %>% mutate(CBCL_AnxiousDepressed_v7_T = na_if(CBCL_AnxiousDepressed_v7_T, 999))
clinical_sublong <- clinical_sublong %>% mutate(CBCL_Internalization_v1_T = na_if(CBCL_Internalization_v1_T, 999))
clinical_sublong <- clinical_sublong %>% mutate(CBCL_Internalization_v7_T = na_if(CBCL_Internalization_v7_T, 999))
clinical_sublong <- clinical_sublong %>% mutate(ABCC_Leth_Soc_Withdr_v1 = na_if(ABCC_Leth_Soc_Withdr_v1, 999))
clinical_sublong <- clinical_sublong %>% mutate(ABCC_Leth_Soc_Withdr_v7 = na_if(ABCC_Leth_Soc_Withdr_v7, 999))

# check values in histograms to see if NA worked
clinical_sublong %>%
  select_if(is.numeric) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name, scales = "free")

# Calculate difference scores
clinical_sublong$CBCL_diff_anxdep<-clinical_sublong$CBCL_AnxiousDepressed_v7_T-clinical_sublong$CBCL_AnxiousDepressed_v1_T
clinical_sublong$CBCL_diff_int<-clinical_sublong$CBCL_Internalization_v7_T-clinical_sublong$CBCL_Internalization_v1_T
clinical_sublong$ABCC_Leth_Soc_Withdr_diff<-clinical_sublong$ABCC_Leth_Soc_Withdr_v7-clinical_sublong$ABCC_Leth_Soc_Withdr_v1

# copy values from IDcode to new column named Participant_ID
clinical_sublong$Participant_ID <- clinical_sublong$IDcode


#### 6.2a) Alpha log power - EEG baseline as predictor ####
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Alpha") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(ap_diff = retest - test) %>%
  ungroup() 
# get ppts included in change EEG analyses
load("IDs_incl_LFreq_abspow.RData")
d_eeg_filt_ids <- d_filt %>%
  filter(Participant_ID %in% ppts)
#merge with clinical data
d_eeg_cl <- dplyr::left_join(d_eeg_filt_ids, clinical_sublong, by = "Participant_ID")
# clear up
rm(d_eeg, d_eeg_trt, d_eeg_filt, LPow_test_avg, d_eeg_filt_ids, d_filt)

#### 6.2a.1) CBCL Anxiety / Depression ####
# set v7 values to nan if there is APow_Bs_centred it nan
d_eeg_cl_Cad <- d_eeg_cl 
d_eeg_cl_Cad$CBCL_AnxiousDepressed_v7_T[is.na(d_eeg_cl_Cad$LPow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "CBCL_AnxiousDepressed_v7_T", "CBCL_AnxiousDepressed_v1_T", "CBCL_diff_anxdep", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_Cad2 <- d_eeg_cl_Cad[complete.cases(d_eeg_cl_Cad[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_Cad)

participant_averages <- d_eeg_cl_Cad2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    CBCL_diff_anxdep=first(CBCL_diff_anxdep),
    CBCL_AnxiousDepressed_v1_T=first(CBCL_AnxiousDepressed_v1_T),
    CBCL_AnxiousDepressed_v7_T=first(CBCL_AnxiousDepressed_v7_T)) %>%
  distinct()

# Run interaction model
m_int <- lm(CBCL_diff_anxdep ~ CBCL_AnxiousDepressed_v1_T + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Al2_lpow_Cad_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for alpha log power predicting CBCL Anx/Dep change")
ggsave("ClinCh_Alpha2_lpow_CBCLAnxDep_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

# clean up
rm(m_int, d_eeg_cl_Cad2, participant_averages, plot_ass)

#### 6.2a.2) CBCL Internalizing ####
# set v7 values to nan if there is APow_Bs_centred it nan
d_eeg_cl_Cint <- d_eeg_cl 
d_eeg_cl_Cint$CBCL_AnxiousDepressed_v7_T[is.na(d_eeg_cl_Cint$LPow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "CBCL_Internalization_v7_T", "CBCL_Internalization_v1_T", "CBCL_diff_int", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_Cint2 <- d_eeg_cl_Cint[complete.cases(d_eeg_cl_Cint[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_Cint)

participant_averages <- d_eeg_cl_Cint2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    CBCL_diff_int=first(CBCL_diff_int),
    CBCL_Internalization_v1_T=first(CBCL_Internalization_v1_T),
    CBCL_Internalization_v7_T=first(CBCL_Internalization_v7_T)) %>%
  distinct()

# Run interaction model
m_int <- lm(CBCL_diff_int ~ CBCL_Internalization_v1_T + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Al2_lpow_Cint_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for alpha log power predicting CBCL Int change")
ggsave("ClinCh_Alpha2_lpow_CBCLInt_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

# clean up
rm(m_int, d_eeg_cl_Cint2, participant_averages, plot_ass)


#### 6.2b) Beta log power - EEG baseline as predictor ####
# Get baseline values for beta power + average within ppts
d_eeg <- read.csv(
  "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
  colClasses = c("character", rep(NA, ncol(
    read.csv(
      "RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.csv",
      nrows = 1
    )
  ) - 1))
) %>%
  rename(Participant_ID = ClinicalID) %>%
  filter(Frequencies == "Beta") %>%
  mutate(
    Participant_ID = as.factor(Participant_ID),
    Log_Power = ifelse(is.nan(Log_Power), NA, Log_Power),
    Condition = as.factor(Condition),
    Regions = as.factor(Regions),
    Sessions = as.factor(Sessions)
  ) %>%
  filter(N_trials > 0)
# generate test/retest columns 
d_eeg_trt <- d_eeg %>%
  select(Participant_ID, Sessions, Condition, Regions, Frequencies, Log_Power) %>%
  dcast(Participant_ID + Condition + Frequencies + Regions ~ Sessions, value.var = "Log_Power")
# calculate average baseline on test sessions
LPow_test_avg <- d_eeg %>%
  filter(Sessions == 'test') %>% 
  summarise(val = mean(Log_Power, na.rm = TRUE))
# centre the baseline around the average
d_eeg_filt <- d_eeg_trt %>%
  mutate(LPow_test_avg_val = LPow_test_avg) %>% 
  mutate(LPow_Bs_centred = test - LPow_test_avg_val$val)
# calculate change scores
d_filt <-d_eeg_filt %>%
  group_by(Participant_ID, Condition, Regions) %>%
  mutate(bp_diff = retest - test) %>%
  ungroup() 
# get ppts included in change EEG analyses
load("IDs_incl_LFreq_abspow.RData")
d_eeg_filt_ids <- d_filt %>%
  filter(Participant_ID %in% ppts)
#merge with clinical data
d_eeg_cl <- dplyr::left_join(d_eeg_filt_ids, clinical_sublong, by = "Participant_ID")
# clear up
rm(d_eeg, d_eeg_trt, d_eeg_filt, LPow_test_avg, d_eeg_filt_ids, d_filt)

#### 6.2b.1) CBCL Anxiety / Depression ####
# set v7 values to nan if there is APow_Bs_centred it nan
d_eeg_cl_Cad <- d_eeg_cl 
d_eeg_cl_Cad$CBCL_AnxiousDepressed_v7_T[is.na(d_eeg_cl_Cad$LPow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "CBCL_AnxiousDepressed_v7_T", "CBCL_AnxiousDepressed_v1_T", "CBCL_diff_anxdep", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_Cad2 <- d_eeg_cl_Cad[complete.cases(d_eeg_cl_Cad[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_Cad)

participant_averages <- d_eeg_cl_Cad2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    CBCL_diff_anxdep=first(CBCL_diff_anxdep),
    CBCL_AnxiousDepressed_v1_T=first(CBCL_AnxiousDepressed_v1_T),
    CBCL_AnxiousDepressed_v7_T=first(CBCL_AnxiousDepressed_v7_T)) %>%
  distinct()

# Run interaction model
m_int <- lm(CBCL_diff_anxdep ~ CBCL_AnxiousDepressed_v1_T + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Be2_lpow_Cad_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for beta log power predicting CBCL Anx/Dep change")
ggsave("ClinCh_Beta2_lpow_CBCLAnxDep_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

# clean up
rm(m_int, d_eeg_cl_Cad2, participant_averages, plot_ass)

#### 6.2b.2) CBCL Internalizing ####
# set v7 values to nan if there is APow_Bs_centred it nan
d_eeg_cl_Cint <- d_eeg_cl 
d_eeg_cl_Cint$CBCL_AnxiousDepressed_v7_T[is.na(d_eeg_cl_Cint$LPow_Bs_centred)] <- NA
#select variables for model
fixed_cols <- c(
  "CBCL_Internalization_v7_T", "CBCL_Internalization_v1_T", "CBCL_diff_int", "treatment_group_code", "LPow_Bs_centred", "Age_patient_v0"
)
d_eeg_cl_Cint2 <- d_eeg_cl_Cint[complete.cases(d_eeg_cl_Cint[, fixed_cols]), ]
rm(fixed_cols, d_eeg_cl_Cint)

participant_averages <- d_eeg_cl_Cint2 %>%
  group_by(Participant_ID) %>%
  reframe(
    Lpow_Bs_centred_avg = mean(LPow_Bs_centred, na.rm = TRUE),
    treatment_group_code = first(treatment_group_code),
    Age_patient_v0 = first(Age_patient_v0),
    CBCL_diff_int=first(CBCL_diff_int),
    CBCL_Internalization_v1_T=first(CBCL_Internalization_v1_T),
    CBCL_Internalization_v7_T=first(CBCL_Internalization_v7_T)) %>%
  distinct()

# Run interaction model
m_int <- lm(CBCL_diff_int ~ CBCL_Internalization_v1_T + Lpow_Bs_centred_avg * treatment_group_code + scale(Age_patient_v0), 
            na.action = "na.omit", 
            data = participant_averages)
summary(m_int)
# create table for reporting
Be2_lpow_Cint_table <- lm_results_table(m_int)

# Check the sample size
ppts_incl_curr <- participant_averages %>%
  group_by(treatment_group_code) %>%
  summarise(N=n())
ppts_incl_curr
rm(ppts_incl_curr)

# Check models for assumptions for lm
plot_ass <- check_lm_clinchange_assumptions(m_int, participant_averages, title = "Model assumptions for beta log power predicting CBCL Int change")
ggsave("ClinCh_Beta2_lpow_CBCLInt_Check_model.jpg", plot_ass, width = 8, height = 6, dpi = 300, bg = "white") 

# clean up
rm(m_int, d_eeg_cl_Cint2, participant_averages, plot_ass)

#### Save tables into Word document ####
doc <- read_docx() %>%
  # Alpha log power - EEG baseline as predictor
  body_add_par("Table 1.1: Baseline alpha log power predicting change on CBCL Anxious/ Depressed - model output", style = "heading 1") %>%
  body_add_flextable(Al2_lpow_Cad_table) %>%
  body_add_par("") %>%
  body_add_par("Table 1.2: Baseline alpha log power predicting change on CBCL Internalization - model output", style = "heading 1") %>%
  body_add_flextable(Al2_lpow_Cint_table) %>%
  body_add_break() %>%  
  # Beta log power - EEG baseline as predictor
  body_add_par("Table 2.1: Baseline beta log power predicting change on CBCL Anxious/ Depressed - model output", style = "heading 1") %>%
  body_add_flextable(Be2_lpow_Cad_table) %>%
  body_add_par("") %>%
  body_add_par("Table 2.2: Baseline beta log power predicting change on CBCL Internalization - model output", style = "heading 1") %>%
  body_add_flextable(Be2_lpow_Cint_table) %>%

print(doc, target = "ClinAnxietyPred_tables.docx")

rm(Al_lpow_AIM_table, Al_lpow_VABSsocss_table, Al_lpow_SRS2tot_table)
rm(BLg_lpow_AIM_table, B_lpow_VABSsocss_table, B_lpow_SRS2tot_table)
rm(doc)