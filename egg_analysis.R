# Egg analysis for L4 Nagi data
# Run everything IN ORDER except what is marked as [O], which means optional
# [CUSTOM] means that the line can be changed to suit conditions or controls
# in the experiment

# Load necessary libraries
library(tidyverse)
library(readxl)
library(DescTools)
library(ggpubr)
library(ggfortify)
library(plotly)
library(ggsci)

# General raw data is confusing so we are using aggregated data
# [CUSTOM] Do we want to use channel as a replicate?
# If channel is true, it will take "eggs_channel" file/s
# If condition is true, it will take "eggs_condition" files
# If we only have one file to analyse, we need it to be "eggs_channel"
channel <- TRUE

if(channel){
  path_egg <- list.files("data", pattern = "eggs_channel", full.names = TRUE)
  replicate <- "Using experiment-chip-channel as replicate."
} else {
  path_egg <- list.files("data", pattern = "eggs_condition", full.names = TRUE)
  if(length(path_egg) == 1){
    stop("There is only one file. Please use the exp-id_eggs_channel.xlsx file.")
  } 
  replicate <- "Using experiment as replicate."
}
if(length(path_egg) < 1){
  stop("There is no suitable egg file. Please check if your channel flag is set up correctly.")
}



# FIRST: EGG EMERGENCE
# Load tables, fix headers
start_tables <- map(path_egg, function(path){
  table <- read_excel(path, col_names = FALSE, sheet = 1)

  # Gets condition, summary type (mean, std, n) and replicate ID
  # (either experiment ID or experiment ID + chip-channel) for each column
  nanum <- sum(is.na(table[[1]]))
  
  headers_x <- table[1, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
  headers_y <- table[nanum - 1, ] |> unlist(use.names = FALSE) |> discard(is.na)
  headers_z <- str_extract(basename(path), "^[^_]+")

  if(nanum == 4){
    headers_z2 <- table[nanum - 2, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
    headers_z <- paste(headers_z, headers_z2, sep = "-")
  }
  headers <- paste(headers_x, headers_z, headers_y, sep = "_")

  # Strips the table of "header" rows, then adds correct ones: Type for
  # Emergence step and Emergence hour column, then the previously extracted headers
  table <- table |> slice((nanum + 1):nrow(table))
  names(table) <- c("Type", headers)

  # Then shifts the table from horizontal to vertical
  table <- table |>
    pivot_longer(
      cols = -Type,
      names_to = c("condition", "Summary"),
      names_pattern = "(.*)_(mean|std|N_chambers)",
      values_to = "value"
    ) |>
    pivot_wider(
      names_from = c("Type", "Summary"),
      values_from = value,
    ) |>
    mutate(across(-condition, as.numeric))
  table
})

# Gets data from all files in one table
egg_table <- bind_rows(start_tables) |>
  separate(
    col = "condition",
    into = c("condition", "rep_id"),
    sep = "_"
  )

# [CUSTOM] Filter out FUdR data. Can be removed or changed for another condition
egg_table <- egg_table |> filter(!str_detect(condition, "FUdR"))

# Set condition as factor
egg_table$condition <- as.factor(egg_table$condition)

# Compare egg_emergence statistically
# [CUSTOM] significant p-value to check normality, check for post-hoc tests
sig_pval <- 0.05

# [CUSTOM] Set control variable: Detects "Water"
condition_levels <- levels(egg_table$condition)
control <- condition_levels[str_detect(condition_levels, fixed("Water"))]

# Normality plot function
save_normality <- function(model, egg_var){
  qq <- ggqqplot(residuals(model))
  sl <- autoplot(model, which = 3)[[1]]
  egg_plot <- ggarrange(sl, qq)
  ggsave(filename = paste0("results/egg/", egg_var, "_normality.png"), plot = egg_plot,
       width = 20, height = 15, dpi = 1000, units = "cm")
}


# Check normality: Are the residuals normal?
# Save normality plot
model  <- lm(`Emergence hour_mean` ~ condition, data = egg_table)
save_normality(model, "egg_emergence")

# Check normality of residuals with shapiro test. If not normal,
# can check QQ plot. If it looks okay, ANOVA is quite robust.
shap <- shapiro.test(residuals(model))

# Sink (save) results
sink("results/egg/egg_emergence.txt")
# Print what kind of replicate we are using
cat(replicate)
cat("\n\nTest egg emergence hour:\n\n")
cat("Check normality:\n")
# If residuals are not normally distributed
if(shap$p.value < sig_pval){
  cat("Residuals distribution not normal. Please check assumptions at ")
  cat("egg_emergence_normality.png to use ANOVA results.\n\n")
  cat("Otherwise, here is a Kruskal-Wallis test:\n")
  # Kruskal-Wallis test for non-normal residuals
  egge_kwt <- kruskal.test(`Emergence hour_mean` ~ condition, data = egg_table)
  print(egge_kwt)

  # If Kruskal-Wallis is significant, do post-hoc testing
  if (egge_kwt$p.value < sig_pval){
    # Dunn's test, equivalent to TukeyHSD
    cat("\nResults from Dunn's test:\n")
    egge_dnn <- dunn_test(egg_table, `Emergence hour_mean` ~ condition, p.adjust.method = "holm") #Change method to "BH" for less stringency
    egge_dnn <- as.data.frame(egge_dnn)
    print(egge_dnn)
    # Might want to add only a comparison with the control (for better statistical power, equivalent to Dunnett's)
  } else {
    cat("\nNo significant results from Kruskal-Wallis. No post-hoc test performed.\n")
  }
} else {cat("Residuals distribution is normal.")}

# One-way ANOVA
# Do ANOVA test anyway, in case shapiro is significant, but plots look good.
# As well as for non significant shapiro test
cat("\nOne-way ANOVA summary:\n\n")
egge_anova <- aov(`Emergence hour_mean` ~ condition, data = egg_table)
print(summary(egge_anova))
  
# If ANOVA is significant, do post-hoc Dunnett's and Tukey's
if (summary(egge_anova)[[1]]$`Pr(>F)`[1] < sig_pval){
  cat("\nResults from Dunnett's test:\n")
  egge_dnt <- DunnettTest(`Emergence hour_mean` ~ condition, data = egg_table, control = control)
  print(egge_dnt)
  egge_dnt <- as.data.frame(egge_dnt[[1]])
  cat("\nResults from Tukey's test:\n")
  egge_thsd <- TukeyHSD(egge_anova)
  egge_thsd <- as.data.frame(egge_thsd$condition) |>
    mutate(p.adj.signif = case_when(`p adj` < 0.001 ~ "***",
                                    `p adj` < 0.01 ~ "**",
                                    `p adj` < 0.05 ~ "*",
                                    `p adj` < 0.1 ~ ".",
                                    `p adj` >= 0.1 ~ "ns"))
  print(egge_thsd)
} else {
  cat("\nNo significant results from ANOVA. No post-hoc test performed.\n")
}

# Finish saving results
sink()

# Plot egg emergence hour
# Mean, std, n of replicates
egg_eh <- egg_table |>
  group_by(condition) |>
  summarise(
    mean = mean(`Emergence hour_mean`),
    std = sd(`Emergence hour_mean`),
    n = n()
  ) |>
  mutate(
    sem = std / sqrt(n)
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
  max_y <- max(egg_eh[[plot_var]] + egg_eh[[error_var]], egg_table[[point_var]], na.rm = TRUE)
  min_y <- min(egg_eh[[plot_var]] - egg_eh[[error_var]], egg_table[[point_var]], na.rm = TRUE)
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
egg_eh_plot <- ggplot(egg_eh, aes(x = condition, y = mean, color = condition)) +
    geom_boxplot(show.legend = FALSE) +
    geom_jitter(aes(x = condition, y = `Emergence hour_mean`), data = egg_table, show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem),  show.legend = FALSE) +
    #scale_color_igv() + #[CUSTOM] Color scale can be changed.
    scale_x_discrete(guide = guide_axis(angle = 90)) + # [CUSTOM] Change to change angle of labels of x axis
    labs(x = "", y = "Egg emergence (hours)") # [CUSTOM] Change to change the x, y axis labels

if (exists("egge_dnt")){
  egg_eh_plot <- egg_eh_plot + plot_signif(egge_dnt, "mean", "sem", "Emergence hour_mean")
}

ggsave(filename = "results/egg/egg_emergence_plot.png", plot = egg_eh_plot,
       width = 17, height = 15, dpi = 1000, units = "cm")


# SECOND: EGG COUNT
# Load tables, fix headers
count_tables <- map(path_egg, function(path){
  table <- read_excel(path, col_names = FALSE, sheet = 2) |> select(-1)

  # Gets condition, summary type (mean, std, n) and replicate ID
  # (either experiment ID or experiment ID + chip-channel) for each column
  nanum <- sum(is.na(table[[1]]))
  
  headers_x <- table[1, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
  headers_y <- table[nanum, ] |> unlist(use.names = FALSE) |> discard(is.na)
  headers_z <- str_extract(basename(path), "^[^_]+")

  if(nanum == 3){
    headers_z2 <- table[nanum - 1, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
    headers_z <- paste(headers_z, headers_z2, sep = "-")
  }
  headers <- c(headers_y[1:2], paste(headers_x, headers_z, headers_y[-c(1:2)], sep = "_"))

  # Strips the table of "header" rows, then adds correct ones: Type for
  # Emergence step and Emergence hour column, then the previously extracted headers
  table <- table |> slice((nanum + 2):nrow(table))
  names(table) <- headers

  # Then shifts the table from horizontal to vertical
  # With columns Step, Hour, condition, mean, std, N
  table <- table |>
    pivot_longer(
      cols = -c("Step", "Hour"),
      names_to = c("condition", "Summary"),
      names_pattern = "(.*)_(mean_.*|std_.*|N_.*)",
      values_to = "value"
    ) |>
    pivot_wider(
      names_from = Summary,
      values_from = value,
    ) |>
    mutate(across(-condition, as.numeric))
  table
})

# Gets data from all files in one table
egg_count_data <- bind_rows(count_tables) |>
  separate(
    col = "condition",
    into = c("condition", "rep_id"),
    sep = "_"
  ) |>
  rename_with(~ "mean_eggs_per_worm", .cols = matches("mean")) |>
  rename_with(~ "std_eggs_per_worm", .cols = matches("std")) |>
  rename_with(~ "n_chambers", .cols = matches("N_chambers"))

# [CUSTOM] Filter out FUdR data. Can be removed or changed for another condition
egg_count_data <- egg_count_data |> filter(!str_detect(condition, "FUdR"))

# Set condition as factor
egg_count_data$condition <- as.factor(egg_count_data$condition)

# Get Total egg count
egg_total <- egg_count_data |>
  group_by(rep_id) |>
  summarise(eggs_total = sum(mean_eggs_per_worm))

# Add normalized egg count
egg_count_table <- full_join(egg_count_data, egg_total, by = join_by(rep_id))

egg_count_table <- egg_count_table |>
  mutate(
    egg_norm = mean_eggs_per_worm / eggs_total,
    day_hour = ceiling(Hour / 24) * 24,
    day = ceiling(Hour / 24)
  )

# Plot of egg_counts
# Summarized data per day
egg_counts_day <- egg_count_table |>
  group_by(day, condition) |>
  summarise(
    mean = mean(mean_eggs_per_worm),
    std = sd(mean_eggs_per_worm),
    mean_norm = mean(egg_norm),
    std_norm = sd(egg_norm),
    n = n()
  ) |>
  mutate(
    sem = std / sqrt(n),
    sem_norm = std_norm / sqrt(n)
  )

# Plot and save plotly
egg_count_plot <- ggplot(egg_counts_day, aes(x = day, y = mean_norm, color = condition)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_norm-sem_norm, ymax = mean_norm+sem_norm)) +
  scale_color_igv()  + #[CUSTOM] Color scale can be changed.
  xlab("Time (days)") + # [CUSTOM] Change to change the x axis label
  ylab("Normalized egg count (A. U.)") + # [CUSTOM] Change to change the y axis label
  labs(color = "Condition")

egg_count_plotly <- ggplotly(egg_count_plot)

htmlwidgets::saveWidget(as_widget(egg_count_plotly), "results/egg/egg_count.html")


# Statistical analysis
# Get data: maximum values of egg laying
# As normalization is per rep_id, this contains maximums
# of both normalized and raw data
data_max <- egg_count_table |>
  group_by(rep_id) |>
  filter(mean_eggs_per_worm == max(mean_eggs_per_worm))

# Testing Maximum eggs laid, Maximum normalized eggs laid, and time at maximum
egg_list <- list(
  mean_eggs_per_worm = "Maximum egg laying value",
  egg_norm = "Maximum normalized egg laying value",
  Hour = "Time at maximum egg laying"
) 

# For mean_eggs_per_worm, egg_norm (normalized eggs), Hour
for(egg_var in names(egg_list)){

  # Normality plot: save to file
  model  <- lm(reformulate("condition", response = egg_var), data = data_max)
  save_normality(model, egg_var)

  # Sink (save) results
  sink(paste0("results/egg/egg_count_", egg_var, ".txt"))
  # Print what kind of replicate we are using
  cat(replicate)
  cat("\n\nTest ")
  cat(egg_list[[egg_var]])
  cat(":\nCheck normality:\n")
  # Check normality of residuals with shapiro test. If not normal,
  # can check QQ plot. If it looks okay, ANOVA is quite robust.
  shap <- shapiro.test(residuals(model))

  if(shap$p.value < sig_pval){
    cat("Residuals distribution not normal. Please check assumptions at ")
    cat(egg_var)
    cat("_normality.png to use ANOVA results.\n\n")
    cat("Otherwise, here is a Kruskal-Wallis test:\n")
    # Kruskal-Wallis test for non-normal residuals
    eggc_kwt <- kruskal.test(reformulate("condition", response = egg_var), data = data_max)
    print(eggc_kwt)

    # If Kruskal-Wallis is significant, do post-hoc testing
    if (eggc_kwt$p.value < sig_pval){
      # Dunn's test, equivalent to TukeyHSD
      cat("\nResults from Dunn's test:\n")
      eggc_dnn <- dunn_test(data_max, reformulate("condition", response = egg_var), p.adjust.method = "holm") #Change method to "BH" for less stringency
      print(as.data.frame(eggc_dnn))
      assign(paste0("eggc_dnn_", egg_var), as.data.frame(eggc_dnn))
      # Might want to add only a comparison with the control (for better statistical power, equivalent to Dunnett's)
      # This is done by taking the p values of eggc_dnn (eggc_dnn$p), and the comparisons (eggc_dnn$comparisons)
      # And correcting only the p values (with p.adjust) where the comparison has the control.
    } else {
      cat("\nNo significant results from Kruskal-Wallis. No post-hoc test performed.\n")
    }
  } else {cat("Residuals distribution is normal.")}

  # One-way ANOVA
  # Do ANOVA test anyway, in case shapiro is significant, but plots look good.
  # As well as for non significant shapiro test
  cat("\nOne-way ANOVA summary:\n\n")
  eggc_anova <- aov(reformulate("condition", response = egg_var), data = data_max)
  print(summary(eggc_anova))
  
  # If ANOVA is significant, do post-hoc Dunnett's and Tukey's
  if (summary(eggc_anova)[[1]]$`Pr(>F)`[1] < sig_pval){
    cat("\nResults from Dunnett's test:\n")
    eggc_dnt <- DunnettTest(reformulate("condition", response = egg_var), data = data_max, control = control)
    assign(paste0("eggc_dnt_", egg_var), as.data.frame(eggc_dnt[[1]]))
    print(eggc_dnt)
    cat("\nResults from Tukey's test:\n")
    eggc_thsd <- TukeyHSD(eggc_anova)
    eggc_thsd <- as.data.frame(eggc_thsd$condition) |>
      mutate(p.adj.signif = case_when(`p adj` < 0.001 ~ "***",
                                      `p adj` < 0.01 ~ "**",
                                      `p adj` < 0.05 ~ "*",
                                      `p adj` < 0.1 ~ ".",
                                      `p adj` >= 0.1 ~ "ns"))
    assign(paste0("eggc_thsd_", egg_var), eggc_thsd)
    print(eggc_thsd)
  } else {
    cat("\nNo significant results from ANOVA. No post-hoc test performed.\n")
  }

  # Finish saving results
  sink()
}
