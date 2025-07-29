# Motility analysis for L4 Nagi data

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
library(afex)
library(DescTools)
library(emmeans)
library(ggsci)

# Create results directory
dir.create("results/motility", showWarnings = FALSE, recursive = TRUE)

# Get motility data file
# Takes all files in the folder data/ which contain "motility_filtered_raw" in the name
path_mot <- list.files("data", pattern = "motility_filtered_raw", full.names = TRUE)

if(length(path_mot) < 1){
  stop("There is no suitable motility file.")
}

# Load excel file
# Selects only area, volume, length (and not eccentricity, orientation, etc.)
table_mot <- map(path_mot, function(path){
  read_excel(path, col_names = TRUE) |> mutate(exp_id = str_extract(basename(path), "^[^_]+"))
})

table_mot <- bind_rows(table_mot)

# [CUSTOM] Filter out FUdR data. Can be removed or changed for another condition
table_mot <- table_mot |> filter(!str_detect(condition, "FUdR"))

# [CUSTOM] Do we want to use worms as replicates? What about channels?
worm_rep <- FALSE
channels <- TRUE

# [CUSTOM] Filter worms with less than 10 measurements (useful when using worm as replicate)
# Can change filter_number

# The problem with using worms as replicates is that, if they do not have all the values,
# there are problems afterwards with the ANOVA: it will remove the replicates that 
# don't have data for all timepoints.

if(worm_rep){
  filter_number <- 10

  valid_worms <- table_mot |>
      group_by(exp_id, chip, channel, chamber, worm_id) |>
      summarise(worm_id_count = n(), .groups = "drop") |>
      filter(worm_id_count >= filter_number)

  data_mot <- table_mot |> inner_join(valid_worms, by = c("exp_id", "chip", "channel", "chamber", "worm_id"))
  data_mot <- data_mot |> select(-worm_id_count)
} else {
  data_mot <- table_mot
}

# Set what is the replicate
if(length(path_mot) > 1){
  if(worm_rep){
    rep_cols <- c("exp_id", "chip", "channel", "chamber", "worm_id")
    replicate <- "Modelled using worm as replicate"
  } else if (channels) {
    rep_cols <- c("exp_id", "chip", "channel")
    replicate <- "Modelled using experiment-chip-channel as replicate."
  } else if (!channels) {
    rep_cols <- c("exp_id", "condition")
    replicate <- "Modelled using experiment as replicate."
  }
# If only one file, replicate ID is the chip_channel combination, or the worm
} else if (length(path_mot) == 1) {
  if(worm_rep){
    rep_cols <- c("chip", "channel", "chamber", "worm_id")
    replicate <- "Modelled using worm as replicate"
  } else {
  rep_cols <- c("chip", "channel")
  replicate <- "Modelled using chip-channel as replicate."
  }
}

# Sort data table
data_mot <- data_mot |>
  mutate(
    rep_id = pmap_chr(across(all_of(rep_cols)), ~ paste(..., sep = "_")),
    day = ceiling(time / (60 * 60 * 24)),
    hour = ceiling(time / (60 * 60)),
    h_nr = time / (60 * 60)
  ) |>
  select(exp_id, rep_id, condition, step, time, day, hour, h_nr,
    head_amplitude, tail_amplitude, displacement_speed, bodybends_frequency
  )

# Plot each variable head_amplitude, tail_amplitude, displacement_speed, bodybends_frequency
# Aggregated per day

# Function to summarise data per each variable
plot_motility <- function(data, variable){
  plotdata <- data |>
    group_by(day, condition) |>
    summarise(
      mean = mean({{variable}}),
      std = sd({{variable}}),
      n = n()
    ) |>
    mutate(
      sem = std / sqrt(n)
    )
  ggplot(plotdata, aes(x = day, y = mean, color = condition)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem)) +
    scale_color_igv()  + #[CUSTOM] Color scale can be changed.
    xlab("Time (days)") + # [CUSTOM] Change to change the x axis label
    labs(color = "Condition")
}

# Plot and save plotly
ha_plot <- plot_motility(data_mot, head_amplitude) + ylab("Head amplitude (mm)")
ha_plotly <- ggplotly(ha_plot)

htmlwidgets::saveWidget(as_widget(ha_plotly), "results/motility/motility_head_amplitude.html")

ta_plot <- plot_motility(data_mot, tail_amplitude) + ylab("Tail amplitude (mm)")
ta_plotly <- ggplotly(ta_plot)

htmlwidgets::saveWidget(as_widget(ta_plotly), "results/motility/motility_tail_amplitude.html")

ds_plot <- plot_motility(data_mot, displacement_speed) + ylab("Displacement speed (mm/s)")
ds_plotly <- ggplotly(ds_plot)

htmlwidgets::saveWidget(as_widget(ds_plotly), "results/motility/motility_displacement_speed.html")

bf_plot <- plot_motility(data_mot, bodybends_frequency) + ylab("Bodybends frequency (Hz)")
bf_plotly <- ggplotly(bf_plot)

htmlwidgets::saveWidget(as_widget(bf_plotly), "results/motility/motility_bodybends_frequency.html")

# Get data for statistical analysis

# Allows wider lines and enough rows when saving statistics results to text file
options(width = 1000)
options(max.print = 2000)

# [CUSTOM] Set significant p-value
sig_pval <- 0.05

# Summarize each replicate (for the four variables) grouped by condition and day
data_testmot <- data_mot |>
  group_by(condition, rep_id, day) |>
  summarise(
    ha = mean(head_amplitude),
    ta = mean(tail_amplitude),
    ds = mean(displacement_speed),
    bf = mean(bodybends_frequency)
  )

# This function fills empty values with 0, and can also filter for experiment if needed
from_data_testmot_fill <- function(data_testmot, id_pattern = NULL){
  if(!is.null(id_pattern)){
    smalldata <- data_testmot |> filter(str_detect(rep_id, id_pattern))
  } else {
    smalldata <- data_testmot
  }
  condition_map <- smalldata |> distinct(rep_id, condition)

  # Store original rows
  original_rows <- smalldata |> ungroup() |> select(rep_id, day)

  # Fill missing rows with 0
  smalldata <- smalldata |> ungroup() |> select(-condition) |>
    tidyr::complete(rep_id, day, fill = list(ha = 0, ta = 0, ds = 0, bf = 0)) |>
    mutate(was_missing = !paste(rep_id, day) %in% paste(original_rows$rep_id, original_rows$day))

  # Add condition to filled rows
  smalldata <- smalldata |>
    left_join(condition_map, by = "rep_id")
  
  return(smalldata)
}

# The mixed-measures two way anova deals with unbalanced data by erasing
# the row of data. So filling it avoids the data being deleted. So then, for example
# if one data is longer in time than the others, maybe we should either let it get cut
# Or fill the shorter one, or analyse them separately and then see if the conclusions are
# the same.

# So these are the options:
# 1. Use all data, filled with 0 when missing
data_all <- from_data_testmot_fill(data_testmot)

# 2. Use all data as is (it will remove incomplete replicates)
data_all_nofill <- data_testmot

# 3. Use data separated per experiment, with filled missing values
exp_ids <- unique(table_mot$exp_id)

data_by_exp <- lapply(exp_ids, function(id) {
  from_data_testmot_fill(data_testmot, id)
})
names(data_by_exp) <- exp_ids

# 4. Use data separated per experiment as is (it will remove incomplete replicates)
# This option is not recommended, as if we end up with only one replicate
# the ANOVA will not work.

data_by_exp_nofill <- lapply(exp_ids, function(id) {
    data_testmot |> filter(str_detect(rep_id, id))
})
names(data_by_exp_nofill) <- paste0(exp_ids, "_nofill")


# Save statistical tests to file

# Function to check normality
save_normality <- function(aov_ez_result, age, variable){
  # Save to PNG file
  png(paste0("results/motility/", age, "_", variable, "_normality.png"), width = 1200, height = 600)

  # Set up 1 row, 2 columns layout
  par(mfrow = c(1, 2))

  # Residuals vs Fitted
  fitted_vals <- fitted(aov_ez_result$lm)
  residuals_vals <- resid(aov_ez_result$lm)
  plot(fitted_vals, residuals_vals,
       xlab = "Fitted Values",
       ylab = "Residuals",
       main = "Residuals vs Fitted Values")
  abline(h = 0, col = "red", lty = 2)

  # Q-Q plot
  res <- residuals(aov_ez_result$lm)
  qqnorm(res, main = "Q-Q Plot of Residuals")
  qqline(res, col = "red", lty = 2)

  # Turn off plotting device
  dev.off()
}

# Function to run all normality tests, normality plots, anova, and posthoc tests.
run_anova_and_posthoc <- function(df_age, name, variable, worm_age, sig_pval){
  # Sink (save) results, say which
  sink(paste0("results/motility/motility_", name, "_", variable, "_", worm_age, ".txt"))

  cat(replicate)
  # First let's do young worms
  cat("\n\nTest ")
  cat(worm_age)
  cat(" worms:\n")
  # Mixed two-way ANOVA with DAY as within factor, CONDITION as between factor
  # Within means subjects are repeated (measuring same worms on day 1 and day 10)
  # Between means subjects are different (measuring different worms on treatment and control)
  aov_age <- aov_ez(id = "rep_id",
                     dv = variable,
                     data = df_age,
                     within = "day",
                     between = "condition")

  # Save normality plot
  save_normality(aov_age, worm_age, variable)

  # Check normality of residuals with shapiro test. If not normal,
  # can check QQ plot. If it looks okay, ANOVA is quite robust.
  # Also check equality of variances with Levene test. Can check
  # residuals vs fitted plot.
  cat("Check normality and equality of variances:\n")

  shap <- shapiro.test(residuals(aov_age$lm))
  lev <- LeveneTest(reformulate("condition", response = variable), data = df_age)

  if(shap$p.value < sig_pval || is.na(lev$`Pr(>F)`[1]) || lev$`Pr(>F)`[1] < sig_pval){
    cat("Residuals distribution not normal or unequal variances. Please check assumptions at ")
    cat(worm_age)
    cat("_")
    cat(variable)
    cat("_normality.png to use ANOVA results.\n\n")
    # Work in progress
    #cat("Otherwise, here is an ART ANOVA test:\n")
    # Use library(ARTool)
  } else {cat("Residuals distribution is normal, variances are equal.")}

  # Print Two-way ANOVA summary
  # Sphericity correction information:
  # Generally, you should use the Greenhouse-Geisser correction (more strict),
  # specially when epsilon < 0.75 (GG eps). If epsilon > 0.75,
  # some statisticians recommend the Huynh-Feldt correction (Girden 1992).
  # Here we will use GG correction.
  cat("\nMixed measures Two-way ANOVA summary:\n\n")
  print(summary(aov_age))

  # Post-hoc with emmeans
  if(any(aov_age$anova_table$`Pr(>F)` < sig_pval, na.rm = TRUE)){

    if(aov_age$anova_table["condition", "Pr(>F)"] < sig_pval){
      cat("\nPost hoc: condition main effect\n")
      print(emmeans(aov_age, pairwise ~ condition))
    }

    if (aov_age$anova_table["day", "Pr(>F)"] < sig_pval) {
      cat("\nPost hoc: day main effect\n")
      print(emmeans(aov_age, pairwise ~ day))
    }

    if (aov_age$anova_table["condition:day", "Pr(>F)"] < sig_pval) {
      # Compare conditions at each day
      # At each day, how do the conditions differ?
      cat("\nPost hoc: condition differences at each day\n")
      emm <- emmeans(aov_age, ~ condition | day)
      pwres <- pairs(emm, adjust = "holm")
      sig_pwres <- summary(pwres) |>
        as.data.frame() |>
        filter(p.value < sig_pval)
      print(sig_pwres)

      # Compare days within each condition
      # Does this condition improve or decline over time?
      cat("\nPost hoc: day differences at each condition\n")
      emm2 <- emmeans(aov_age, ~ day | condition)
      pwres2 <- pairs(emm2, adjust = "holm")
      sig_pwres2 <- summary(pwres2) |>
        as.data.frame() |>
        filter(p.value < sig_pval)
      print(sig_pwres2)
    }
  } else {
    cat("\nNo significant results from ANOVA. No post-hoc test performed.\n")
  }

  # Finish saving results
  sink()
}

# Significant results, first set variables to test
# (head amplitude, tail amplitude, displacement speed, bodybends frequency)
variables <- c("ha", "ta", "ds", "bf")

# Choose which option to use
# 1 = all data with fill, 2 = all data no fill,
# 3 = per-experiment with fill, 4 = per-experiment no fill
option <- 1

if (option == 1) {
  dataset_names <- "all"
  dataset_list <- list(all = data_all)

} else if (option == 2) {
  dataset_names <- "all_no_fill"
  dataset_list <- list(all_no_fill = data_all_nofill)

} else if (option == 3) {
  dataset_names <- names(data_by_exp)
  dataset_list <- data_by_exp

} else if (option == 4) {
  dataset_names <- names(data_by_exp_nofill)
  dataset_list <- data_by_exp_nofill

}

# Loop through datasets, variables, and do young and old results
for(name in dataset_names){
  df <- dataset_list[[name]]

  # Separate the data in young and old
  df_young <- df |> filter(day <= 10)
  df_old <- df |> filter(day > 10)

  for(variable in variables){
    run_anova_and_posthoc(df_young, name, variable, worm_age = "young", sig_pval)
    run_anova_and_posthoc(df_old, name, variable, worm_age = "old", sig_pval)
    
  }
}