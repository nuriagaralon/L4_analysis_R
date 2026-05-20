# Growth analysis for L4 Nagi data

#--------------------------
# Author: Núria Garriga Alonso
# GitHub: @nuriagaralon
# Repository: https://github.com/nuriagaralon/L4_analysis_R
#--------------------------

# Run everything IN ORDER except what is marked as [O], which means optional
# [CUSTOM] means that the line can be changed to suit conditions or controls
# in the experiment

# Load necessary libraries
library(ggfortify)
library(tidyverse)
library(readxl)
library(plotly)
library(pracma)
library(DescTools)
library(rstatix)
library(ggpubr)
library(ggsci)
#library(nlme) #[O], explained later

# Create results directory
dir.create("results/growth", showWarnings = FALSE, recursive = TRUE)

# Get growth data file
# Takes all files in the folder data/ which contain "growth_filtered_raw" in the name
path_gro <- list.files("data", pattern = "growth_filtered_raw", full.names = TRUE)

if(length(path_gro) < 1){
  stop("There is no suitable growth file.")
}

# Load excel file
# Selects only area, volume, length (and not eccentricity, orientation, etc.)
table_gro <- map(path_gro, function(path){
  read_excel(path, col_names = TRUE) |> select(1:11)
})

table_gro <- bind_rows(table_gro)

# Prepare data to analyze and plot
# [CUSTOM] Do we want to normalize? What is the biological replicate?
normalize <- FALSE # TRUE if we want to normalize data, FALSE if not
channels <- TRUE # TRUE if replicates are exp_id-chip-channel, FALSE if they are exp_id

# Set y labels for plot:
if(normalize){
  lab_amod <- "Area (A. U.)"
  lab_lmod <- "Length (A. U.)"
  lab_vmod <- "Volume (A. U.)"
} else {
  lab_amod <- "Area (mm\u00b2)" # mm²
  lab_lmod <- "Length (mm)"
  lab_vmod <- "Volume (mm\u00b3)" # mm³
}
lab_a <- "Area (\u03bcm\u00b2)"  # μm²
lab_l <- "Length (\u03bcm)"  # μm
lab_v <- "Volume (\u03bcm\u00b3)"  # μm³


# Takes the value at t0 of each condition for normalizing the data
# (separated by experiment ID, but channels together, which might need change)
# Normalization consists of dividing every value by its value at time 0
t0 <- table_gro |>
  filter(step < 7) |>
  select(-chip, -channel, -chamber, -step, -time, -worm_id) |>
  group_by(exp_id, condition) |>
  summarise_all(mean)

data_gro <- full_join(table_gro, t0, suffix = c("", ".t0avg"), by = join_by(exp_id, condition))

# Normalization, select needed columns, add replicate ID.
# If more than one file, replicate ID is the experiment ID+condition
# But if channels is TRUE, it is exp_id-chip-channel
if(length(path_gro) > 1){
  if(channels){
    rep_cols <- c("exp_id", "chip", "channel")
    replicate <- "Modelled using experiment-chip-channel as replicate."
  } else if (!channels) {
    rep_cols <- c("exp_id", "condition")
    replicate <- "Modelled using experiment as replicate."
  }
# If only one file, replicate ID is the chip_channel combination
} else if (length(path_gro) == 1) {
  rep_cols <- c("chip", "channel")
  replicate <- "Modelled using chip-channel as replicate."
}

# If we want normalization, the values are in variable_mod.
# If not, length_mod is in mm, area_mod in mm^2, volume_mod in (mm^3)
# If you want units in um, use length, area, volume
data_gro <- data_gro |>
  mutate(
    rep_id = pmap_chr(across(all_of(rep_cols)), ~ paste(..., sep = "_")),
    day = ceiling(time / (60 * 60 * 24)),
    hour = ceiling(time / (60 * 60)),
    h_nr = time / (60 * 60),
    length_mod = if (normalize) length / length.t0avg else length / 1000,
    area_mod = if (normalize) area / area.t0avg else area / 1e6,
    volume_mod = if (normalize) volume / volume.t0avg else volume / 1e9
  ) |>
  select(condition, step, time, day, hour, h_nr, rep_id,
      length, length_mod,
      area, area_mod,
      volume, volume_mod
  )

# [CUSTOM] Filter out FUdR data. Can be removed or changed for another condition
data_gro <- data_gro |> filter(!str_detect(condition, "FUdR"))

# [CUSTOM] Filter out from wiz file, where they starved for a weekend
data_gro <- data_gro |>
  filter(!(str_detect(rep_id, "wiz6YHzJUXFL") & hour >= 395))

if(length(path_gro) == 1 && str_detect(path_gro, "wiz6YHzJUXFL")) {
  data_gro <- data_gro |> filter(hour < 395)
}


# Sets condition and replicate ID as factors, needed for statistics
data_gro$condition <- factor(data_gro$condition)
data_gro$rep_id <- factor(data_gro$rep_id)


# Functions to fit logistic growth and plot
# We are using summarized data per each hour, as all data points would be
# Very confusing and very full

growth_sigmoid <- function(dataset, column){
  # 1. Group and summarize data
  summ_column <- dataset |>
      group_by(condition, hour) |>
      summarise(
          mean = mean({{column}}),
          sd = sd({{column}}),
          n = n(),
          .groups = "drop"
      )
  
  # 2. Fit model per condition
  fit_column <- summ_column |>
      group_by(condition) |>
      nest() |>
      mutate(
          # Set reasonable start values
          A_start = map_dbl(data, ~ max(.x$mean, na.rm = TRUE)),
          B_start = 0.01,
          C_start = map2_dbl(data, A_start, ~ .x$hour[which.min(abs(.x$mean - .y / 2))]),
          # Fit the model
          model = pmap(list(data, A_start, B_start, C_start), function(df, A, B, C){
              tryCatch(
                  nls(mean ~ A / (1 + exp(-B * (hour - C))),
                      data = df,
                      start = list(A = A, B = B, C = C)),
                  error = function(e) NULL
              )
          }))
  # 2.5. If there is a model with NULL, recalculate C_start (mean of C_starts that worked) and try again
  value_for_c <- fit_column |> filter(!map_lgl(model, is.null)) |> pull(C_start) |> mean(na.rm = TRUE)

  fit_column <- fit_column |>
    mutate(
      C_start = if_else(map_lgl(model, is.null), value_for_c, C_start),
      model = pmap(list(model, data, A_start, B_start, C_start), function(m, df, A, B, C){
              if (is.null(m)) {
              tryCatch(
                  nls(mean ~ A / (1 + exp(-B * (hour - C))),
                      data = df,
                      start = list(A = A, B = B, C = C)),
                  error = function(e) NULL
              )} else {m}
          }),
      # Calculate values for y for each x according to the model, for plotting
      fitted = map2(model, data, ~ if (!is.null(.x)) predict(.x, newdata = .y) else rep(NA, nrow(.y))),

      # We compute R2 to see if they fit similarly. It is not, however, a way to compare fits.
      r2 = map2(data, fitted, ~ 1 - sum((.x$mean - .y)^2)/sum((.x$mean-mean(.x$mean))^2))
    )
  fit_column
}

plot_sigmoid <- function(fit_column){
  # 3. Unnest and plot
  plot_column <- fit_column |>
    select(condition, data, fitted) |>
    unnest(c(data, fitted))
  
  column_ggplot <- ggplot(plot_column, aes(x = hour, y = mean, color = condition)) + 
    geom_point(size = 1) + #[CUSTOM] Size is 1 to have small points
    geom_line(aes(y = fitted), linewidth = 1) + #[CUSTOM] Linewidth is 1 to have a slightly thicker line
    xlab("Time (hour)") + # [CUSTOM] Change to change the x axis label
    theme_minimal() +
    labs(color = "Condition") +
    scale_color_igv() #[CUSTOM] Color scale can be changed.
  column_ggplot
}

# Use the functions to plot area, length, volume, and save plotly
# Stops you if any model did not converge.
# [CUSTOM] Change area_mod, length_mod, volume_mod to area, length, volume
# for plotting raw data in um. Also change ylabs below (when using plot_sigmoid function).

area_data <- growth_sigmoid(data_gro, area_mod)
if(any(map_lgl(area_data$model, is.null))){
  stop("One of the area models is NULL, please rerun with a different C_start (try value_for_c = 40)")
}

length_data <- growth_sigmoid(data_gro, length_mod)
if(any(map_lgl(length_data$model, is.null))){
  stop("One of the length models is NULL, please rerun with a different C_start (try value_for_c = 40)")
}

volume_data <- growth_sigmoid(data_gro, volume_mod)
if(any(map_lgl(volume_data$model, is.null))){
  stop("One of the volume models is NULL, please rerun with a different C_start (try value_for_c = 40)")
}

## [O] If the model did not converge, change placeholder (for example, to volume_data). 
## Can also change the value_for_c. It might also not converge.
#value_for_c <- 40
#placeholder <- placeholder |>
#    mutate(
#      C_start = if_else(map_lgl(model, is.null), value_for_c, C_start),
#      model = pmap(list(model, data, A_start, B_start, C_start), function(m, df, A, B, C){
#              if (is.null(m)) {
#              tryCatch(
#                  nls(mean ~ A / (1 + exp(-B * (hour - C))),
#                      data = df,
#                      start = list(A = A, B = B, C = C)),
#                  error = function(e) NULL
#              )} else {m}
#          }),
#      # Calculate values for y for each x according to the model, for plotting
#      fitted = map2(model, data, ~ if (!is.null(.x)) predict(.x, newdata = .y) else rep(NA, nrow(.y))),
#
#      # We compute R2 to see if they fit similarly. It is not, however, a way to compare fits.
#      r2 = map2(data, fitted, ~ 1 - sum((.x$mean - .y)^2)/sum((.x$mean-mean(.x$mean))^2))
#    )

# Plot and save plotly
# [CUSTOM] Here you can change ylabs. If plotting area, length, volume,
# change lab_amod, lab_lmod, lab_vmod to lab_a, lab_l, lab_v
area_ggplot <- plot_sigmoid(area_data) + ylab(lab_amod)
area_plotly <- ggplotly(area_ggplot)

htmlwidgets::saveWidget(as_widget(area_plotly), "results/growth/growth_area.html")

length_ggplot <- plot_sigmoid(length_data) + ylab(lab_lmod)
length_plotly <- ggplotly(length_ggplot)

htmlwidgets::saveWidget(as_widget(length_plotly), "results/growth/growth_length.html")

volume_ggplot <- plot_sigmoid(volume_data) + ylab(lab_vmod)
volume_plotly <- ggplotly(volume_ggplot)

htmlwidgets::saveWidget(as_widget(volume_plotly), "results/growth/growth_volume.html")

# Compare conditions statistically
# [CUSTOM] significant p-value to check normality, check for post-hoc tests
sig_pval <- 0.05

# [CUSTOM] Set control variable: Detects "Water"
condition_levels <- levels(data_gro$condition)
control <- condition_levels[str_detect(condition_levels, fixed("Water"))]

# Allows wider lines and enough rows when saving statistics results to text file
options(width = 1000)
options(max.print = 2000)

# Get data summary by *time*, not *hour* for more accuracy
# Use scaled time (h_nr) so the numbers are smaller and modeling finds them
# And separate by rep_id (for ANOVA or Kruskal-Wallis means)
growth_summary <- function(dataset, column){
  dataset |>
    group_by(condition, rep_id, time, h_nr) |>
    summarise(
      mean = mean({{column}}),
      sd = sd({{column}}),
      n = n(),
      .groups = "drop"
    )
}

# Comparison 1: AUC
# Calculate AUC values with trapezoidal function
growth_AUC <- function(growth_summ_data){
  growth_summ_data |>
    group_by(condition, rep_id) |>
    arrange(time, .by_group = TRUE) |>
    summarise(AUC = trapz(time, mean)) |>
    ungroup()
}

# Check normality: Are the residuals normal?
# Save normality plot
save_normality <- function(model, growth_var, param){
  qq <- ggqqplot(residuals(model))
  sl <- autoplot(model, which = 3)[[1]]
  growth_plot <- ggarrange(sl, qq)
  ggsave(filename = paste0("results/growth/", growth_var, "_normality_", param, ".png"), plot = growth_plot,
       width = 20, height = 15, dpi = 1000, units = "cm")
}

# Area, length, volume
# [CUSTOM] Change area_mod, length_mod, volume_mod to area, length, volume
# for analysing raw data in um. Also change ylabs below (when using plot_rep_sigmoid function)
area_summ <- growth_summary(data_gro, area_mod)
area_AUC <- growth_AUC(area_summ)

length_summ <- growth_summary(data_gro, length_mod)
length_AUC <- growth_AUC(length_summ)

volume_summ <- growth_summary(data_gro, volume_mod)
volume_AUC <- growth_AUC(volume_summ)

# Save statistical tests to file
growth_list <- list(
  area   = area_AUC,
  length = length_AUC,
  volume = volume_AUC
)

# For area, length, volume
for(growth_var in names(growth_list)){
  df <- growth_list[[growth_var]]

  # Normality plot: save to file
  model  <- lm(AUC ~ condition, data = df)
  save_normality(model, growth_var, "AUC")

  # Sink (save) results
  sink(paste0("results/growth/growth_AUC_", growth_var, ".txt"))
  # Print what kind of replicate we are using
  cat(replicate)
  cat("\n\nTest AUC:\n")
  cat("Check normality:\n")
  # Check normality of residuals with shapiro test. If not normal,
  # can check QQ plot. If it looks okay, ANOVA is quite robust.
  shap <- shapiro.test(residuals(model))

  if(shap$p.value < sig_pval){
    cat("Residuals distribution not normal. Please check assumptions at ")
    cat(growth_var)
    cat("_normality_AUC.png to use ANOVA results.\n\n")
    cat("Otherwise, here is a Kruskal-Wallis test:\n")
    # Kruskal-Wallis test for non-normal residuals
    AUC_kwt <- kruskal.test(AUC ~ condition, data = df)
    print(AUC_kwt)

    # If Kruskal-Wallis is significant, do post-hoc testing
    if (AUC_kwt$p.value < sig_pval){
      # Dunn's test, equivalent to TukeyHSD
      cat("\nResults from Dunn's test:\n")
      AUC_dnn <- dunn_test(df, AUC ~ condition, p.adjust.method = "holm") #Change method to "BH" for less stringency
      print(as.data.frame(AUC_dnn))
      # Might want to add only a comparison with the control (for better statistical power, equivalent to Dunnett's)
      # This is done by taking the p values of AUC_dnn (AUC_dnn$p), and the comparisons (AUC_dnn$comparisons)
      # And correcting only the p values (with p.adjust) where the comparison has the control.
    } else {
      cat("\nNo significant results from Kruskal-Wallis. No post-hoc test performed.\n")
    }
  } else {cat("Residuals distribution is normal.")}

  # One-way ANOVA
  # Do ANOVA test anyway, in case shapiro is significant, but plots look good.
  # As well as for non significant shapiro test
  cat("\nOne-way ANOVA summary:\n\n")
  AUC_anova <- aov(AUC ~ condition, data = df)
  print(summary(AUC_anova))
  
  # If ANOVA is significant, do post-hoc Dunnett's and Tukey's
  if (summary(AUC_anova)[[1]]$`Pr(>F)`[1] < sig_pval){
    cat("\nResults from Dunnett's test:\n")
    AUC_dnt <- DunnettTest(AUC ~ condition, data = df, control = control)
    print(AUC_dnt)
    cat("\nResults from Tukey's test:\n")
    AUC_thsd <- TukeyHSD(AUC_anova)
    AUC_thsd <- as.data.frame(AUC_thsd$condition) |>
      mutate(p.adj.signif = case_when(`p adj` < 0.001 ~ "***",
                                      `p adj` < 0.01 ~ "**",
                                      `p adj` < 0.05 ~ "*",
                                      `p adj` < 0.1 ~ ".",
                                      `p adj` >= 0.1 ~ "ns"))
    print(AUC_thsd)
  } else {
    cat("\nNo significant results from ANOVA. No post-hoc test performed.\n")
  }

  # Finish saving results
  sink()
}

# Comparison 2: Logistic growth parameters:
# Compare models by comparing model parameters
# A = maximum value = plateau
# B = growth rate = slope
# C = half growth = inflection point

# Fit NLS like it was done for plotting but for each replicate
growth_params <- function(growth_summ_data){
  fit_param <- growth_summ_data |>
    group_by(condition, rep_id) |>
    nest() |>
    mutate(
        # Set reasonable start values
        A_start = map_dbl(data, ~ max(.x$mean, na.rm = TRUE)),
        B_start = 0.01,
        C_start = map2_dbl(data, A_start, ~ .x$h_nr[which.min(abs(.x$mean - .y / 2))]),
        # Fit the model
        model = pmap(list(data, A_start, B_start, C_start), function(df, A, B, C){
            tryCatch(
                nls(mean ~ A / (1 + exp(-B * (h_nr - C))),
                    data = df,
                    start = list(A = A, B = B, C = C)),
                error = function(e) NULL
            )
        }))
  # If there is a model with NULL, recalculate C_start (mean of C_starts that worked) and try again
  value_for_c <- fit_param |>
  filter(!map_lgl(model, is.null)) |>
  pull(C_start) |>
  mean(na.rm = TRUE)
  
  fit_param <- fit_param |>
  mutate(
    C_start = if_else(map_lgl(model, is.null), value_for_c, C_start),
    model = pmap(list(model, data, A_start, B_start, C_start), function(m, df, A, B, C){
            if (is.null(m)) {
            tryCatch(
                nls(mean ~ A / (1 + exp(-B * (h_nr - C))),
                    data = df,
                    start = list(A = A, B = B, C = C)),
                error = function(e) NULL
            )} else {m}
        }),
    # Calculate fitted values
    fitted = map2(model, data, ~ if (!is.null(.x)) predict(.x, newdata = .y) else rep(NA, nrow(.y))),
    
    # We compute R2 to see if they fit similarly. It is not, however, a way to compare fits.
    r2 = map2(data, fitted, ~ 1 - sum((.x$mean - .y)^2)/sum((.x$mean-mean(.x$mean))^2))
  )
  # Set condition as factor for statistics
  fit_param$condition <- factor(fit_param$condition)
  fit_param
}

# Extract parameters (A, B, C) from fits
get_params <- function(growth_param_data){
  growth_param_data |>
    mutate(params = map(model, ~ as_tibble(as.list(coef(.x))))) |>
    select(condition, rep_id, params) |>
    unnest(params) |>
    ungroup()
}

# Use the functions to model area, length, volume, per replicate
# Stops you if any model did not converge.
area_params_fit <- growth_params(area_summ)
if(any(map_lgl(area_params_fit$model, is.null))){
  stop("One of the area models is NULL, please rerun with a different C_start (try value_for_c = 40)")
}

length_params_fit <- growth_params(length_summ)
if(any(map_lgl(length_params_fit$model, is.null))){
  stop("One of the length models is NULL, please rerun with a different C_start (try value_for_c = 40)")
}

volume_params_fit <- growth_params(volume_summ)
if(any(map_lgl(volume_params_fit$model, is.null))){
  stop("One of the volume models is NULL, please rerun with a different C_start (try value_for_c = 40)")
}

# Plot the replicates separately to see:
# NOT plotting data points (there are too many)
# Instead, plot only the modelled logistic growth curves
# We can see if there is a lot of difference between replicates
plot_rep_sigmoid <- function(params_fit){
  plot_column <- params_fit |>
    select(condition, data, fitted) |>
    unnest(c(data, fitted))
  
  column_ggplot <- ggplot(plot_column, aes(x = h_nr, y = fitted, color = rep_id)) + 
    geom_line() +
    xlab("Time (hour)") + # [CUSTOM] Change to change the x axis label
    labs(color = "Replicate") +
    theme_minimal()
  column_ggplot
}

# Plot and save plotly
# [CUSTOM] Here you can change ylabs. If plotting area, length, volume,
# change lab_amod, lab_lmod, lab_vmod to lab_a, lab_l, lab_v
area_rep_ggplot <- plot_rep_sigmoid(area_params_fit) + ylab(lab_amod)
area_rep_plotly <- ggplotly(area_rep_ggplot)

htmlwidgets::saveWidget(as_widget(area_rep_plotly), "results/growth/growth_area_rep.html")

length_rep_ggplot <- plot_rep_sigmoid(length_params_fit) + ylab(lab_lmod)
length_rep_plotly <- ggplotly(length_rep_ggplot)

htmlwidgets::saveWidget(as_widget(length_rep_plotly), "results/growth/growth_length_rep.html")

volume_rep_ggplot <- plot_rep_sigmoid(volume_params_fit) + ylab(lab_vmod)
volume_rep_plotly <- ggplotly(volume_rep_ggplot)

htmlwidgets::saveWidget(as_widget(volume_rep_plotly), "results/growth/growth_volume_rep.html")


# Extract parameters from model fit
area_params <- get_params(area_params_fit)
length_params <- get_params(length_params_fit)
volume_params <- get_params(volume_params_fit)

# Save statistical tests to file
growth_params_list <- list(
  area   = area_params,
  length = length_params,
  volume = volume_params
)

# For area, length, volume
for(growth_var in names(growth_params_list)){
  df <- growth_params_list[[growth_var]]
    
  # Sink (save) results, all parameters in one file
  sink(paste0("results/growth/growth_params_", growth_var, ".txt"))
  # Print what kind of replicate we are using
  cat(replicate)

  # For each parameter
  for(param in c("A", "B", "C")){

    # Normality plot: save to file
    model  <- lm(reformulate("condition", response = param), data = df)
    save_normality(model, growth_var, param)
    
    cat("\n\nTest ")
    cat(param)
    cat(":\nCheck normality:\n")
    # Check normality of residuals with shapiro test. If not normal,
    # can check QQ plot. If it looks okay, ANOVA is quite robust.
    shap <- shapiro.test(residuals(model))

    if(shap$p.value < sig_pval){
      cat("Residuals distribution not normal. Please check assumptions at ")
      cat(growth_var)
      cat("_normality_")
      cat(param)
      cat(".png to use ANOVA results.\n\n")
      cat("Otherwise, here is a Kruskal-Wallis test:\n")
      # Kruskal-Wallis test for non-normal residuals
      param_kwt <- kruskal.test(reformulate("condition", response = param), data = df)
      print(param_kwt)

      # If Kruskal-Wallis is significant, do post-hoc testing
      if (param_kwt$p.value < sig_pval){
        # Dunn's test, equivalent to TukeyHSD 
        cat("\nResults from Dunn's test:\n")
        param_dnn <- dunn_test(df, reformulate("condition", response = param), p.adjust.method = "holm") #Change method to "BH" for less stringency
        print(as.data.frame(param_dnn))
        # Might want to add only a comparison with the control (for better statistical power)
        # This is done by taking the p values of dnn (dnn$p), and the comparisons (dnn$comparisons)
        # And correcting only the p values (with p.adjust) where the comparison has the control.
      } else {
        cat("\nNo significant results from Kruskal-Wallis. No post-hoc test performed.\n")
      }
    } else {cat("Residuals distribution is normal.")}

    # One-way ANOVA
    # Do ANOVA test anyway, in case shapiro is significant, but plots look good.
    # As well as for non significant shapiro test
    cat("\nOne-way ANOVA summary:\n\n")
    param_anova <- aov(reformulate("condition", response = param), data = df)
    print(summary(param_anova))

    # If ANOVA is significant, do post-hoc Dunnett's and Tukey's
    if (summary(param_anova)[[1]]$`Pr(>F)`[1] < sig_pval){
      cat("\nResults from Dunnett's test:\n")
      param_dnt <- DunnettTest(reformulate("condition", response = param), data = df, control = control)
      print(param_dnt)
      cat("\nResults from Tukey's test:\n")
      param_thsd <- TukeyHSD(param_anova)
      param_thsd <- as.data.frame(param_thsd$condition) |>
        mutate(p.adj.signif = case_when(`p adj` < 0.001 ~ "***",
                                        `p adj` < 0.01 ~ "**",
                                        `p adj` < 0.05 ~ "*",
                                        `p adj` < 0.1 ~ ".",
                                        `p adj` >= 0.1 ~ "ns"))
      print(param_thsd)
    } else {
      cat("\nNo significant results from ANOVA. No post-hoc test performed.\n")
    }
  }
  
  # Finish saving results, after all parameters are done
  sink()
}

# Comparison 3: Compare a certain time point
timepoint <- 80 #Hour

# Summarise at timepoint (keep replicates separate)
data_timepoint <- data_gro |>
  filter(hour == timepoint) |>
  group_by(condition, rep_id) |>
  summarise(
    length = mean(length),
    area = mean(area),
    volume = mean(volume),
    length_mod = mean(length_mod),
    area_mod = mean(area_mod),
    volume_mod = mean(volume_mod)
  ) |> ungroup()

# Set condition as factor
data_timepoint$condition <- factor(data_timepoint$condition)

# For area, length, volume
# [CUSTOM] Change area_mod, length_mod, volume_mod to area, length, volume
# for plotting raw data in um.
for(growth_var in c("area_mod", "length_mod", "volume_mod")){
    
  # Sink (save) results
  sink(paste0("results/growth/growth_", timepoint, "_", growth_var, ".txt"))
  # Print what kind of replicate we are using
  cat(replicate)

  # Normality plot: save to file
  model  <- lm(reformulate("condition", response = growth_var), data = data_timepoint)
  save_normality(model, growth_var, timepoint)
  
  cat("\n\nTest at ")
  cat(timepoint)
  cat(" hours:\nCheck normality:\n")
  # Check normality of residuals with shapiro test. If not normal,
  # can check QQ plot. If it looks okay, ANOVA is quite robust.
  shap <- shapiro.test(residuals(model))

  if(shap$p.value < sig_pval){
    cat("Residuals distribution not normal. Please check assumptions at ")
    cat(growth_var)
    cat("_normality_")
    cat(timepoint)
    cat(".png to use ANOVA results.\n\n")
    cat("Otherwise, here is a Kruskal-Wallis test:\n")
    # Kruskal-Wallis test for non-normal residuals
    tp_kwt <- kruskal.test(reformulate("condition", response = growth_var), data = data_timepoint)
    print(tp_kwt)
    # If Kruskal-Wallis is significant, do post-hoc testing
    if (tp_kwt$p.value < sig_pval){
      # Dunn's test, equivalent to TukeyHSD
      cat("\nResults from Dunn's test:\n")
      tp_dnn <- dunn_test(data_timepoint, reformulate("condition", response = growth_var), p.adjust.method = "holm") #Change method to "BH" for less stringency
      print(as.data.frame(tp_dnn))
      assign(paste0("tp_dnn_", growth_var), as.data.frame(tp_dnn))
      # Might want to add only a comparison with the control (for better statistical power)
      # This is done by taking the p values of dnn (dnn$p), and the comparisons (dnn$comparisons)
      # And correcting only the p values (with p.adjust) where the comparison has the control.
    } else {
      cat("\nNo significant results from Kruskal-Wallis. No post-hoc test performed.\n")
    }
  } else {cat("Residuals distribution is normal.")}

  # One-way ANOVA
  # Do ANOVA test anyway, in case shapiro is significant, but plots look good.
  # As well as for non significant shapiro test
  cat("\nOne-way ANOVA summary:\n\n")
  tp_anova <- aov(reformulate("condition", response = growth_var), data = data_timepoint)
  print(summary(tp_anova))

  # If ANOVA is significant, do post-hoc Dunnett's and Tukey's
  if (summary(tp_anova)[[1]]$`Pr(>F)`[1] < sig_pval){
    cat("\nResults from Dunnett's test:\n")
    tp_dnt <- DunnettTest(reformulate("condition", response = growth_var), data = data_timepoint, control = control)
    assign(paste0("tp_dnt_", growth_var), as.data.frame(tp_dnt[[1]]))
    print(tp_dnt)
    cat("\nResults from Tukey's test:\n")
    tp_thsd <- TukeyHSD(tp_anova)
    tp_thsd <- as.data.frame(tp_thsd$condition) |>
      mutate(p.adj.signif = case_when(`p adj` < 0.001 ~ "***",
                                      `p adj` < 0.01 ~ "**",
                                      `p adj` < 0.05 ~ "*",
                                      `p adj` < 0.1 ~ ".",
                                      `p adj` >= 0.1 ~ "ns"))
    assign(paste0("tp_thsd_", growth_var), tp_thsd)
    print(tp_thsd)
  } else {
    cat("\nNo significant results from ANOVA. No post-hoc test performed.\n")
  }
  
  # Finish saving results
  sink()
}

# Plot at timepoint
# Organise data
summ_timepoint <- data_timepoint |>
  group_by(condition) |>
  summarise(
    mean_length = mean(length),
    mean_area = mean(area),
    mean_volume = mean(volume),
    mean_length_mod = mean(length_mod),
    mean_area_mod = mean(area_mod),
    mean_volume_mod = mean(volume_mod),
    n = n(),
    sem_length = sd(length)/sqrt(n),
    sem_area = sd(area)/sqrt(n),
    sem_volume = sd(volume)/sqrt(n),
    sem_length_mod = sd(length_mod)/sqrt(n),
    sem_area_mod = sd(area_mod)/sqrt(n),
    sem_volume_mod = sd(volume_mod)/sqrt(n)
    )

# Plotting function
plot_timepoint <- function(plot_var, error_var){
  ggplot(summ_timepoint, aes(x = condition, y = {{plot_var}})) +
    geom_point(size = 2) + # [CUSTOM] Change to change point size
    geom_errorbar(aes(ymin = {{plot_var}}-{{error_var}}, ymax = {{plot_var}}+{{error_var}})) +
    scale_x_discrete(guide = guide_axis(angle = 90)) + # [CUSTOM] Change to change angle of text of x axis
    xlab("") # [CUSTOM] Change to change the x axis label
}

# Plotting significant values
# Function to get comparison names
# We already have the variable condition_levels from when we extracted the control
extract_comp_from_cond <- function(comparison, condition_levels){
  matches <- condition_levels[map_lgl(condition_levels, ~ str_detect(comparison, fixed(.x)))]
  if(length(matches) == 2){
    return(as.list(matches))
  } else {
    warning(paste("Could not uniquely extract 2 groups from:", comparison))
    return(list(NA, NA))
  }
}

# Function to plot
plot_signif <- function(tp_posthoc_var, plot_var, error_var){
  # 1. Get the sig_data
  # If it comes from Dunnett's test, we need to add sig stars
  if("pval" %in% colnames(tp_posthoc_var)){
    sig_data <- tp_posthoc_var |>
      filter(pval < sig_pval) |>
      rownames_to_column(var = "condition") |>
      mutate(p.adj.signif = case_when(pval < 0.001 ~ "***",
                                      pval < 0.01 ~ "**",
                                      pval < 0.05 ~ "*",
                                      pval >= 0.05 ~ "ns"))
  # For Dunn just filter
  } else if ("group1" %in% colnames(tp_posthoc_var)) {
    sig_data <- tp_posthoc_var |>
      filter(!(p.adj.signif %in% c("ns", ".")))
  # For Tukey, filter and get row names in condition
  } else {
    sig_data <- tp_posthoc_var |>
      filter(!(p.adj.signif %in% c("ns", "."))) |>
      rownames_to_column(var = "condition")
  }

  # 2. Return NULL if no significant comparisons
  if (nrow(sig_data) == 0) {
    return(NULL)  # nothing to add to plot
  }
  
  # 3. Get comparison names 
  # If it comes from Dunn's test
  if("group1" %in% colnames(tp_posthoc_var)){
    comparisons <- sig_data |>
      mutate(pair = pmap(list(group1, group2), c)) |> pull(pair)
  # For Dunnett and Tukey
  } else {
    comparisons <- map(sig_data$condition, ~ extract_comp_from_cond(.x, condition_levels))
    comparisons <- map(comparisons, unlist)
  }

  # 4. Calculate position of significance stars in plot
  max_y <- max(summ_timepoint[[plot_var]] + summ_timepoint[[error_var]])
  min_y <- min(summ_timepoint[[plot_var]] - summ_timepoint[[error_var]])
  gap <- 0.1 * (max_y - min_y)
  y_positions <- seq(from = max_y + gap, by = gap, length.out = length(comparisons))

  # 5. Build geom_signif for ggplot
  plot_tp <- geom_signif(comparison = comparisons,
                         annotations = sig_data$p.adj.signif,
                         y_position = y_positions,
                         tip_length = 0.01)
  plot_tp
}

# Plot and save plots
# [CUSTOM] tp_dnt_area_mod plots Dunnet values and can be changed to 
# tp_thsd_area_mod for Tukey significant values
# tp_dnn_area_mod for Dunn significant values

# [CUSTOM] Change area_mod, length_mod, volume_mod to area, length, volume
# if statistical analysis was done with raw data in um.
# Also change ylabs. Change lab_amod, lab_lmod, lab_vmod to lab_a, lab_l, lab_v

tp_area <- plot_timepoint(mean_area_mod, sem_area_mod) + ylab(lab_amod)
if (exists("tp_dnt_area_mod")){
  tp_area <- tp_area + plot_signif(tp_dnt_area_mod, "mean_area_mod", "sem_area_mod")
}

ggsave(filename = paste0("results/growth/area_", timepoint, ".png"), plot = tp_area,
       width = 17, height = 17, dpi = 1000, units = "cm")

tp_length <- plot_timepoint(mean_length_mod, sem_length_mod) + ylab(lab_lmod)
if (exists("tp_dnt_length_mod")){
  tp_length <- tp_length + plot_signif(tp_dnt_length_mod, "mean_length_mod", "sem_length_mod")
}
ggsave(filename = paste0("results/growth/length_", timepoint, ".png"), plot = tp_length,
       width = 17, height = 17, dpi = 1000, units = "cm")

tp_volume <- plot_timepoint(mean_volume_mod, sem_volume_mod) + ylab(lab_vmod)
if (exists("tp_dnt_volume_mod")){
  tp_volume <- tp_volume + plot_signif(tp_dnt_volume_mod, "mean_volume_mod", "sem_volume_mod")
}
ggsave(filename = paste0("results/growth/volume_", timepoint, ".png"), plot = tp_volume,
       width = 17, height = 17, dpi = 1000, units = "cm")

# Arrange plots for publication	   
growth_argd <- ggarrange(length_ggplot, area_ggplot, volume_ggplot, tp_length, tp_area, tp_volume,
                         labels = c("A", "B", "C", "D", "E", "F"),
                         ncol = 3, nrow = 2, common.legend = TRUE)
ggsave(filename = "results/growth/growth_argd.pdf", plot = growth_argd,
       width = 10, height = 8)


# [UNFINISHED] NLME OR GNLS MODEL
# An nlme model or gnls model is needed for official LRT testing, but it is a lot of parameters
# and a lot of work, and it is difficult for it to converge. Therefore, it is difficult to set up
# generalized modeling for any and all kinds of experiments.
# Might be worth looking into it if the previous testing is not enough

#nlme
# data_gro <- data_gro |>
#   mutate(subject = map2_chr(condition, rep_id, ~ paste(.x, .y, sep = ".")))
# 
# start_values <- fit_area |>
#   select(condition, model) |>
#   ungroup() |>
#   mutate(start_vals = map(model, coef)) |>
#   mutate(diff = map(start_vals, ~ .x - baseline)) |>
#   mutate(named_vals = map2(start_vals, condition, ~ set_names(.x, paste(names(.x), .y, sep = ".")))) |>
#   pull(named_vals) |>
#   unlist()
# 
## Unfinished nlme
# model <- nlme(area_mod ~ A / (1 + exp(-B * (hour - C))),
#               fixed = A + B + C ~ condition,
#               random = A + B + C ~ 1 | subject,
#               data = data_gro,
#               start = start_values)

#gnls
#comp_area <- fit_area |>
#  select(condition, model) |>
#  ungroup() |>
#  mutate(start_vals = map(model, coef))
#
## Get start values
#control <- "Water"
#
#baseline <- comp_area |>
#  filter(str_detect(condition, control)) |>
#  pull(start_vals) |>
#  pluck(1)
#
#start_values <- comp_area |>
#  filter(!str_detect(condition, control)) |>
#  mutate(diff = map(start_vals, ~ .x - baseline)) |>
#  mutate(named_vals = map2(start_vals, condition, ~ set_names(.x, paste(names(.x), .y, sep = ".")))) |>
#  pull(named_vals) |>
#  unlist()
#
#start_values <- c(baseline, start_values)
#
## Get null model
#null_data <- data_gro |>
#    group_by(hour) |>
#    summarise(
#        mean = mean(area_mod),
#        sd = sd(area_mod),
#        n = n(),
#        .groups = "drop"
#    )
## Unfinished gnls
#area_null <- gnls(area_mod ~ A / (1 + exp(-B * (hour - C))),
#                    data = data_gro,
#                    start = baseline)