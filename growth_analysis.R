library(tidyverse)
library(readxl)
library(plotly)
library(pracma)
library(DescTools)
#library(nlme)

# Get growth data file
setwd("C:\\Users\\ngarriga\\Documents\\SydLab-One\\L4_analysis_R")

path_gro <- list.files("data", pattern = "growth_filtered_raw", full.names = TRUE)

# Load excel
table_gro <- read_excel(path_gro, col_names = TRUE) |> select(1:11)


# Sort data to plot
t0 <- table_gro |>
  filter(step < 7) |>
  select(-chip, -channel, -chamber, -step, -time, -worm_id) |>
  group_by(exp_id, condition) |>
  summarise_all(mean)

data_gro <- full_join(table_gro, t0, suffix = c("", ".t0avg"), by = join_by(exp_id, condition))

# From here it is treating chip-channel as a replicate. Needs to be checked.
data_gro <- data_gro |>
  mutate(chip_channel = map2_chr(chip, channel, ~ paste(.x, .y, sep="_"))) |>
  mutate(day = ceiling(time / (60 * 60 * 24))) |>
  mutate(hour = ceiling(time / (60 * 60))) |>
  mutate(h_nr = time / (60 * 60)) |>
  mutate(length.t0norm = length / length.t0avg,
         area.t0norm = area / area.t0avg,
         volume.t0norm = volume / volume.t0avg
  ) |>
  select(exp_id, condition, step, time, day, hour, h_nr, chip_channel,
    length, length.t0norm,
    area, area.t0norm,
    volume, volume.t0norm
    )

data_gro$condition <- factor(data_gro$condition)
#data_gro$condition <- substr(data_gro$condition, 11, 100)
#data_gro$condition <- str_replace(data_gro$condition, "\\?", "\U03BC")

# Function to fit logistic growth and plot

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
          }),
          fitted = map2(model, data, ~ if (!is.null(.x)) predict(.x, newdata = .y) else rep(NA, nrow(.y))),
          r2 = map2(data, fitted, ~ 1 - sum((.x$mean - .y)^2)/sum((.x$mean-mean(.x$mean))^2))
      )
  
  # We compute R2 to see if they fit similarly. It is not, however, a way to compare fits.
  
  # 3. Unnest and plot
  plot_column <- fit_column |>
    select(condition, data, fitted) |>
    unnest(c(data, fitted))
  
  column_ggplot <- ggplot(plot_column, aes(x = hour, y = mean, color = condition)) + 
    geom_point(size = 1) +
    geom_line(aes(y = fitted))
    #theme_minimal() +
    #geom_errorbar(aes(x=hour,ymin=mean-sd, ymax=mean+sd, color = condition))
  column_ggplot
}

# Use the function to plot area, length, volume
area_ggplot <- growth_sigmoid(data_gro, area.t0norm)
area_plotly <- ggplotly(area_ggplot)

length_ggplot <- growth_sigmoid(data_gro, length.t0norm)
length_plotly <- ggplotly(length_ggplot)

volume_ggplot <- growth_sigmoid(data_gro, volume.t0norm)
volume_plotly <- ggplotly(volume_ggplot)

htmlwidgets::saveWidget(as_widget(area_plotly), "results/area.html")

# Compare models by comparing AUC (treating channels as biological replicates)
# 1. Get data summary by *time*, not *hour* for more accuracy
AUC_summ <- data_gro |>
    group_by(condition, chip_channel, time, h_nr) |>
    summarise(
        mean = mean(area.t0norm),
        sd = sd(area.t0norm),
        n = n(),
        .groups = "drop"
    )

# 2. Calculate trapezoidal function
AUC_area <- AUC_summ |>
  group_by(condition, chip_channel) |>
  arrange(time, .by_group = TRUE) |>
  summarise(
    AUC = trapz(time, mean)
  )

# 3. Check normality and see
ggplot(AUC_area, aes(x = AUC)) +
  geom_histogram(bins = 15, fill = "skyblue", color = "black") +
  theme_minimal()

ggplot(AUC_area, aes(sample = AUC)) +
     stat_qq() +
     stat_qq_line() +
     theme_minimal()

shapiro.test(AUC_area$AUC)

# 4. ANOVA and Dunnett test
AUC_anova <- aov(AUC ~ condition, data = AUC_area)
summary(AUC_anova)

condition_levels <- levels(data_gro$condition)
control <- condition_levels[str_detect(condition_levels, "Water")]

AUC_dnt <- DunnettTest(AUC ~ condition, data = AUC_area, control = control)
AUC_dnt

# Compare models by comparing model parameters (treating channels as biological replicates)
# 1. Fit NLS for each replicate
fit_charea <- AUC_summ |>
    group_by(condition, chip_channel) |>
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
value_for_c <- 40

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
  select(condition, chip_channel, params) |>
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
#   mutate(subject = map2_chr(condition, chip_channel, ~ paste(.x, .y, sep = ".")))
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