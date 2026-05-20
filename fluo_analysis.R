# Fluorescence analysis for L4 Nagi data

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
library(ggpubr)
library(ggsci)

# Create results directory
dir.create("results/fluo", showWarnings = FALSE, recursive = TRUE)

# Get fluorescence data file
# Takes all files in the folder data/ which contain "fluo_raw" in the name
path_fluo <- list.files("data", pattern = "fluo_raw", full.names = TRUE)

if(length(path_fluo) < 1){
  stop("There is no suitable fluorescence file.")
}

# Load excel file
# Selects only area, volume, length (and not eccentricity, orientation, etc.)
table_fluo <- map(path_fluo, function(path){
  table <- read_excel(path, col_names = TRUE)
  table |> mutate(
    exp_id = str_extract(basename(path), "^[^_]+"),
    chamber_id = paste(chip, channel, chamber, sep = "-"))
})

table_fluo <- bind_rows(table_fluo)

# Prepare data to analyze and plot
# [CUSTOM] Write which chambers to use for analysis
# Remember to choose chambers with only one worm, at least 3 (better 5 or more)
# per condition, if possible the same amount per condition

if(file.exists("data/fluo_chambers_use.xlsx")){
  cham_to_use <- read_excel("data/fluo_chambers_use.xlsx")
  cham_to_use <- cham_to_use |>
    map(~ na.omit(.x)) |> # Removes NA (empty values)
    map(as.character)
} else {
  stop("There is no file that describes the chambers to use.")
}

data_all_fluo <- imap_dfr(cham_to_use, function(chambers, experiment) {
  table_fluo |> filter(exp_id == experiment, chamber_id %in% chambers)
})

# Set condition as factor
data_all_fluo$condition <- as.factor(data_all_fluo$condition)

# [CUSTOM] Set control variable: Detects "Water"
condition_levels <- levels(data_all_fluo$condition)
control <- condition_levels[str_detect(condition_levels, fixed("Water"))]

# Add day, hour, change intensity minus background to "fluo_intensity"
# In the current data, the variable name has a " at the end, I don't know if it
# is like this from the beginning (check)
data_all_fluo <- data_all_fluo |>
  mutate(
    day = ceiling(time / (60 * 60 * 24)),
    hour = ceiling(time / (60 * 60)),
    fluo_intensity = `tot_Intensity_minus_background"` 
  ) |>
  select(condition, step, time, day, hour,
      exp_id, chamber_id, fluo_intensity
  )

# Make a plot to check at which timepoint (which hour) the initial fluo intensity
# levels to real measurements
check_plot <- ggplot(data_all_fluo, aes(x = hour, y = fluo_intensity, color = condition)) +
  geom_point() + geom_line() + 
  labs(x = "Hour (h)", y = "Fluorescence intensity", color = "Condition")

check_plotly <- ggplotly(check_plot)
check_plotly

htmlwidgets::saveWidget(as_widget(check_plotly), "results/fluo/fluo_initial_check.html")

# [CUSTOM] Set start_hour to the hour we want to begin using the data
# Check fluo_initial_check.html or check_plotly to decide 
start_hour <- 483

data_fluo <- data_all_fluo |> filter(hour >= start_hour)

# Calculate fold change values
# Each experiment should use its own control (but the control condition
# should be the same, and it was set above at the variable control)
data_fluo <- data_fluo |>
  group_by(exp_id) |>
  group_modify(~ {
    # First we get a data frame with the data with only one experiment
    exp_data <- .

    # Then, we filter the control values
    # And average in the case two measurements are from the exact same timepoint
    # Which shouldn't happen, as we're dealing with timepoints in seconds
    data_control <- exp_data |>
      filter(condition == control) |>
      group_by(time) |>
      summarise(
        fluo_intensity_m = mean(fluo_intensity),
        .groups = "drop"
      )

    # Then we calculate the fold change, separate by each experiment_id
    exp_data <- exp_data |>
      rowwise() |>
      mutate(
        closest_control_index = which.min(abs(data_control$time - time)),
        closest_control_value = data_control$fluo_intensity_m[closest_control_index],
        fold_change = fluo_intensity / closest_control_value
      ) |>
      ungroup()

      #Then we pass the resulting dataframe out of the function
      return(exp_data)
  }) |>
  # And remove the exp_id groups
  ungroup()

# [CUSTOM] Change keep_control to FALSE to filter out
# the control values (which will all be FC = 1)
keep_control <- TRUE

if(!keep_control){
data_fluo <- data_fluo |> filter(condition != control)
}

# Summarise data by hour to plot (all chambers of the same
# condition are aggregated)
summ_fluo <- data_fluo |>
    group_by(condition, hour) |>
    summarise(
        fold_change_m = mean(fold_change),
        fold_change_sd = sd(fold_change),
        n_chambers = n_distinct(chamber_id)
    )

# Plot Fold Change values per hour
fluo_hour_plot <- ggplot(summ_fluo, aes(x = hour, y = fold_change_m, color = condition)) +
    geom_point() + geom_line() +
    geom_errorbar(aes(ymin = fold_change_m-fold_change_sd, ymax = fold_change_m+fold_change_sd)) + 
    labs(x = "Hour (h)", y = "Fluorescence Fold Change", color = "Condition")

fluo_hour_plotly <- ggplotly(fluo_hour_plot)

htmlwidgets::saveWidget(as_widget(fluo_hour_plotly), "results/fluo/fluo_FC_hour.html")

# [CUSTOM] Difference between conditions at a specific timepoint
timepoint <- 485 #Hour
sig_pval <- 0.05

# Allows wider lines and enough rows when saving statistics results to text file
options(width = 1000)
options(max.print = 2000)

# Summarise at timepoint (keep replicates separate)
# There should only be one point per hour, condition, exp_id, chamber_id
# but we will summarise just in case
data_fluotp <- data_fluo |>
  filter(hour == timepoint) |>
  group_by(condition, exp_id, chamber_id) |>
  summarise(
    fold_change_m = mean(fold_change),
    fold_change_sd = sd(fold_change),
    n = n()
  ) |> ungroup()

# Check normality: Are the residuals normal?
# Save normality plot
model  <- lm(fold_change_m ~ condition, data = data_fluotp)

qq <- ggqqplot(residuals(model))
sl <- autoplot(model, which = 3)[[1]]
fluotp_norm_plot <- ggarrange(sl, qq)

ggsave(filename = paste0("results/fluo/fluo_", timepoint, "_normality.png"),
       plot = fluotp_norm_plot, width = 20, height = 15,
       dpi = 1000, units = "cm")

# Check normality of residuals with shapiro test. If not normal,
# can check QQ plot. If it looks okay, ANOVA is quite robust.
shap <- shapiro.test(residuals(model))

# Sink (save) results
sink(paste0("results/fluo/fluo_", timepoint, ".txt"))
cat("\n\nTest at ")
cat(timepoint)
cat(" hours:\nCheck normality:\n")
# If residuals are not normally distributed
if(shap$p.value < sig_pval){
  cat("Residuals distribution not normal. Please check assumptions at fluo_")
  cat(timepoint)
  cat("_normality.png to use ANOVA results.\n\n")
  cat("Otherwise, here is a Kruskal-Wallis test:\n")
  # Kruskal-Wallis test for non-normal residuals
  ftp_kwt <- kruskal.test(fold_change_m ~ condition, data = data_fluotp)
  print(ftp_kwt)

  # If Kruskal-Wallis is significant, do post-hoc testing
  if (ftp_kwt$p.value < sig_pval){
    # Dunn's test, equivalent to TukeyHSD
    cat("\nResults from Dunn's test:\n")
    ftp_dnn <- dunn_test(data_fluotp, fold_change_m ~ condition, p.adjust.method = "holm") #Change method to "BH" for less stringency
    ftp_dnn <- as.data.frame(ftp_dnn)
    print(ftp_dnn)
    # Might want to add only a comparison with the control (for better statistical power, equivalent to Dunnett's)
  } else {
    cat("\nNo significant results from Kruskal-Wallis. No post-hoc test performed.\n")
  }
} else {cat("Residuals distribution is normal.")}

# One-way ANOVA
# Do ANOVA test anyway, in case shapiro is significant, but plots look good.
# As well as for non significant shapiro test
cat("\nOne-way ANOVA summary:\n\n")
ftp_anova <- aov(fold_change_m ~ condition, data = data_fluotp)
print(summary(ftp_anova))
  
# If ANOVA is significant, do post-hoc Dunnett's and Tukey's
if (summary(ftp_anova)[[1]]$`Pr(>F)`[1] < sig_pval){
  cat("\nResults from Dunnett's test:\n")
  ftp_dnt <- DunnettTest(fold_change_m ~ condition, data = data_fluotp, control = control)
  print(ftp_dnt)
  ftp_dnt <- as.data.frame(ftp_dnt[[1]])
  cat("\nResults from Tukey's test:\n")
  ftp_thsd <- TukeyHSD(ftp_anova)
  ftp_thsd <- as.data.frame(ftp_thsd$condition) |>
    mutate(p.adj.signif = case_when(`p adj` < 0.001 ~ "***",
                                    `p adj` < 0.01 ~ "**",
                                    `p adj` < 0.05 ~ "*",
                                    `p adj` < 0.1 ~ ".",
                                    `p adj` >= 0.1 ~ "ns"))
  print(ftp_thsd)
} else {
  cat("\nNo significant results from ANOVA. No post-hoc test performed.\n")
}

# Finish saving results
sink()

# Plot at timepoint
# Organise data
data_fluotp_plot <- data_fluotp |>
    ungroup() |>
    group_by(condition) |>
    summarise(
        fc_mean = mean(fold_change_m),
        fc_sd = sd(fold_change_m),
        n = n()
    ) |> mutate(
        fc_sem = fc_sd / sqrt(n)
    )

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
plot_signif <- function(posthoc_res, plot_var, error_var, point_var){
  # 1. Get the sig_data
  # If it comes from Dunnett's test, we need to add sig stars
  if("pval" %in% colnames(posthoc_res)){
    sig_data <- posthoc_res |>
      filter(pval < sig_pval) |>
      rownames_to_column(var = "condition") |>
      mutate(p.adj.signif = case_when(pval < 0.001 ~ "***",
                                      pval < 0.01 ~ "**",
                                      pval < 0.05 ~ "*",
                                      pval >= 0.05 ~ "ns"))
  # For Dunn just filter
  } else if ("group1" %in% colnames(posthoc_res)) {
    sig_data <- posthoc_res |>
      filter(!(p.adj.signif %in% c("ns", ".")))
  # For Tukey, filter and get row names in condition
  } else {
    sig_data <- posthoc_res |>
      filter(!(p.adj.signif %in% c("ns", "."))) |>
      rownames_to_column(var = "condition")
  }

  # 2. Return NULL if no significant comparisons
  if (nrow(sig_data) == 0) {
    return(NULL)  # nothing to add to plot
  }
  
  # 3. Get comparison names 
  # If it comes from Dunn's test
  if("group1" %in% colnames(posthoc_res)){
    comparisons <- sig_data |>
      mutate(pair = pmap(list(group1, group2), c)) |> pull(pair)
  # For Dunnett and Tukey
  } else {
    comparisons <- map(sig_data$condition, ~ extract_comp_from_cond(.x, condition_levels))
    comparisons <- map(comparisons, unlist)
  }
  # 4. Calculate position of significance stars in plot
  max_y <- max(data_fluotp_plot[[plot_var]] + data_fluotp_plot[[error_var]], data_fluotp[[point_var]], na.rm = TRUE)
  min_y <- min(data_fluotp_plot[[plot_var]] - data_fluotp_plot[[error_var]], data_fluotp[[point_var]], na.rm = TRUE)
  gap <- 0.1 * (max_y - min_y)
  y_positions <- seq(from = max_y + gap, by = gap, length.out = length(comparisons))

  # 5. Build geom_signif for ggplot
  plot_tp <- geom_signif(comparison = comparisons,
                         annotations = sig_data$p.adj.signif,
                         y_position = y_positions,
                         tip_length = 0.01,
                         color = "black")
  plot_tp
}

# Plot values and significance
# General plot
plot_fluotp <- ggplot(data_fluotp_plot, aes(x = condition, y = fc_mean, color = condition)) +
    geom_boxplot(show.legend = FALSE) +
    geom_jitter(aes(x = condition, y = fold_change_m), data = data_fluotp, show.legend = FALSE) +
    geom_errorbar(aes(ymin = fc_mean-fc_sem, ymax = fc_mean+fc_sem),  show.legend = FALSE) +
    #scale_color_igv() + #[CUSTOM] Color scale can be changed.
    scale_x_discrete(guide = guide_axis(angle = 90)) + # [CUSTOM] Change to change angle of labels of x axis
    labs(x = "", y = "Fluorescence (Fold Change vs control)") # [CUSTOM] Change to change the x, y axis labels

# Use ftp_dnt, ftp_dnn, ftp_thsd
if (exists("ftp_dnn")){
  plot_fluotp <- plot_fluotp + plot_signif(ftp_dnn, "fc_mean", "fc_sem", "fold_change_m")
}

ggsave(filename = paste0("results/fluo/fluo_", timepoint, ".png"), plot = plot_fluotp,
       width = 17, height = 15, dpi = 1000, units = "cm")

# Save plot for publication
ggsave(filename = paste0("results/fluo/fluo_", timepoint, ".pdf"), plot = plot_fluotp,
       width = 5, height = 4.5)