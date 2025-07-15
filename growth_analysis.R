# Growth analysis for L4 Nagi data
# Run everything IN ORDER except what is marked as [O], which means optional
# [CUSTOM] means that the line can be changed to suit conditions or controls
# in the experiment

# Load necessary libraries
library(tidyverse)
library(readxl)
library(plotly)
library(pracma)
library(rstatix)
library(DescTools)
library(ggpubr)
library(ggsci)
#library(nlme) #[O], explained later

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

# Sort data to plot
# Takes the value at t0 of each condition for normalizing the data
# (separated by experiment ID, but channels together, which might need change)
t0 <- table_gro |>
  filter(step < 7) |>
  select(-chip, -channel, -chamber, -step, -time, -worm_id) |>
  group_by(exp_id, condition) |>
  summarise_all(mean)

data_gro <- full_join(table_gro, t0, suffix = c("", ".t0avg"), by = join_by(exp_id, condition))

# Normalization, select needed columns, add replicate ID.
# If more than one file, replicate ID is the experiment ID+condition
if(length(path_gro) > 1){
  data_gro <- data_gro |>
  mutate(rep_id = map2_chr(exp_id, condition, ~ paste(.x, .y, sep="_"))) |>
  mutate(day = ceiling(time / (60 * 60 * 24))) |>
  mutate(hour = ceiling(time / (60 * 60))) |>
  mutate(h_nr = time / (60 * 60)) |>
  mutate(length.t0norm = length / length.t0avg,
         area.t0norm = area / area.t0avg,
         volume.t0norm = volume / volume.t0avg
  ) |>
  select(condition, step, time, day, hour, h_nr, rep_id,
    length, length.t0norm,
    area, area.t0norm,
    volume, volume.t0norm
    )
  replicate <- "Modelled using experiment as replicate."
# If only one file, replicate ID is the chip_channel combination
}else if (length(path_gro) == 1) {
  data_gro <- data_gro |>
  mutate(rep_id = map2_chr(chip, channel, ~ paste(.x, .y, sep="_"))) |>
  mutate(day = ceiling(time / (60 * 60 * 24))) |>
  mutate(hour = ceiling(time / (60 * 60))) |>
  mutate(h_nr = time / (60 * 60)) |>
  mutate(length.t0norm = length / length.t0avg,
         area.t0norm = area / area.t0avg,
         volume.t0norm = volume / volume.t0avg
  ) |>
  select(condition, step, time, day, hour, h_nr, rep_id,
    length, length.t0norm,
    area, area.t0norm,
    volume, volume.t0norm
    )
  replicate <- "Modelled using chip-channel as replicate."
}

# [CUSTOM] Filter out FUdR data. Can be removed or changed for another condition
data_gro <- data_gro |> filter(!str_detect(condition, "FUdR"))

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
          A_start = map_dbl(data, ~ max(.x$mean, na.rm = TRUE)),
          B_start = 0.01,
          C_start = map2_dbl(data, A_start, ~ .x$hour[which.min(abs(.x$mean - .y / 2))]),
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
    geom_point(size = 1) +
    geom_line(aes(y = fitted), linewidth = 1) +
    xlab("Time (hour)") +
    theme_minimal() +
    scale_color_igv()
  column_ggplot
}

# Use the functions to plot area, length, volume, and save plotly
area_data <- growth_sigmoid(data_gro, area.t0norm)
area_ggplot <- plot_sigmoid(area_data) + ylab("Area (A. U.)")
area_plotly <- ggplotly(area_ggplot)

htmlwidgets::saveWidget(as_widget(area_plotly), "results/area.html")

length_data <- growth_sigmoid(data_gro, length.t0norm)
length_ggplot <- plot_sigmoid(length_data) + ylab("Length (A. U.)")
length_plotly <- ggplotly(length_ggplot)

htmlwidgets::saveWidget(as_widget(length_plotly), "results/length.html")

volume_data <- growth_sigmoid(data_gro, volume.t0norm)
volume_ggplot <- plot_sigmoid(volume_data) + ylab("Volume (A. U.)")
volume_plotly <- ggplotly(volume_ggplot)

htmlwidgets::saveWidget(as_widget(volume_plotly), "results/volume.html")

# Compare conditions statistically
sig_pval <- 0.05
control_name <- str_subset(unique(data_gro$condition), "Water")

# Get data summary by *time*, not *hour* for more accuracy
# Use scaled time so the numbers are smaller and modeling finds them
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
## Calculate trapezoidal function
growth_AUC <- function(growth_summ_data){
  growth_summ_data |>
    group_by(condition, rep_id) |>
    arrange(time, .by_group = TRUE) |>
    summarise(AUC = trapz(time, mean)) |>
    ungroup()
}

## Check normality and see
normality_AUC <- function(growth_AUC_data){
  mu <- mean(growth_AUC_data$AUC)
  sigma <- sd(growth_AUC_data$AUC)

  hi <- ggplot(growth_AUC_data, aes(x = AUC)) +
    geom_histogram(bins = 15, fill = "skyblue", color = "black", aes(y = after_stat(density))) +
    stat_function(fun = dnorm, args = list(mean = mu, sd = sigma), color = "red", linewidth = 1) +
    theme_minimal()

  qq <- ggplot(growth_AUC_data, aes(sample = AUC)) +
       stat_qq() +
       stat_qq_line() +
       theme_minimal()

  ggarrange(hi, qq)
}

# Area
area_summ <- growth_summary(data_gro, area.t0norm)
area_AUC <- growth_AUC(area_summ)
area_normality_AUC <- normality_AUC(area_AUC)

ggsave(filename = "results/area_normality_AUC.png", plot = area_normality_AUC,
       width = 17, height = 15, dpi = 1000, units = "cm")

sink("results/growth_area.txt")
cat(replicate)
cat("\n\nTest AUC:\n")
cat("Check normality:\n")
a <- shapiro.test(area_AUC$AUC)

if(a$p.value < sig_pval){
  cat("Data distribution not normal. Please check area_normality_AUC.png to use ANOVA results.\n")
}

AUC_welchaov <- area_AUC |> welch_anova_test(AUC ~ condition)

print(AUC_welchaov, n = Inf)

  AUC_gh <- area_AUC |> games_howell_test(AUC ~ condition)

print(AUC_gh, n = Inf) #maybe csv

if (AUC_welchaov$p < sig_pval){
    AUC_gh <- area_AUC |> games_howell_test(AUC ~ condition)

} else {
  cat("\nNo significant results from Welch's ANOVA. No post-hoc test performed.\n")
}

sink()

# Length
length_summ <- growth_summary(data_gro, length.t0norm)

# Volume
volume_summ <- growth_summary(data_gro, volume.t0norm)


# 4. ANOVA and Dunnett test
AUC_anova <- aov(AUC ~ condition, data = AUC_area)
summary(AUC_anova)

condition_levels <- levels(data_gro$condition)
control <- condition_levels[str_detect(condition_levels, "Water")]

AUC_dnt <- DunnettTest(AUC ~ condition, data = AUC_area, control = control)
AUC_dnt

# Compare models by comparing model parameters (treating channels as biological replicates)
# 1. Fit NLS for each replicate
fit_charea <- summ_gro |>
    group_by(condition, rep_id) |>
    nest() |>
    mutate(
        A_start = map_dbl(data, ~ max(.x$mean, na.rm = TRUE)),
        B_start = 0.01,
        C_start = map2_dbl(data, A_start, ~ .x$h_nr[which.min(abs(.x$mean - .y / 2))]),
        model = pmap(list(data, A_start, B_start, C_start), function(df, A, B, C){
            tryCatch(
                nls(mean ~ A / (1 + exp(-B * (h_nr - C))),
                    data = df,
                    start = list(A = A, B = B, C = C)),
                error = function(e) NULL
            )
        }))

# REDO NULL MODELS, IF IT DOES NOT WORK, INSTEAD OF C_START = 40,
# CHECK ALL VALUES AND PUT AN AVERAGE
value_for_c <- fit_charea |>
  filter(!map_lgl(model, is.null)) |>
  pull(C_start) |>
  mean(na.rm = TRUE)

fit_charea <- fit_charea |>
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
        })
  )

params_area <- fit_charea |>
  mutate(params = map(model, ~ as_tibble(as.list(coef(.x))))) |>
  select(condition, rep_id, params) |>
  unnest(params)

# A = maximum value = plateau
A_area_anova <- aov(A ~ condition, data = params_area)
summary(A_area_anova)

A_area_dnt <- DunnettTest(A ~ condition, data = params_area, control = control)
A_area_dnt

# B = growth rate = slope
B_area_anova <- aov(B ~ condition, data = params_area)
summary(B_area_anova)

B_area_dnt <- DunnettTest(B ~ condition, data = params_area, control = control)
B_area_dnt

# C = half growth = inflection point
C_area_anova <- aov(C ~ condition, data = params_area)
summary(C_area_anova)

C_area_dnt <- DunnettTest(C ~ condition, data = params_area, control = control)
C_area_dnt

# An nlme model or gnls model is needed for official LRT testing, but it is a lot of parameters
# and a lot of work, and it is difficult for it to converge.
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
# model <- nlme(area.t0norm ~ A / (1 + exp(-B * (hour - C))),
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
#        mean = mean(area.t0norm),
#        sd = sd(area.t0norm),
#        n = n(),
#        .groups = "drop"
#    )
#
#area_null <- gnls(area.t0norm ~ A / (1 + exp(-B * (hour - C))),
#                    data = data_gro,
#                    start = baseline)