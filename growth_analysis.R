library(tidyverse)
library(readxl)

# Get growth data file
setwd("C:\\Users\\ngarriga\\Documents\\SydLab-One\\L4_analysis_R")

path_gro <- list.files("data", pattern = "growth_filtered_raw", full.names = TRUE)

# Load excel
table_gro <- read_excel(path_gro, col_names = TRUE) |> select(2:11)


# Sort data to plot
t0 <- table_gro |>
  filter(step < 7) |>
  select(-chip, -channel, -chamber, -step, -time, -worm_id) |>
  group_by(condition) |>
  summarise_all(mean)

data_gro <- full_join(table_gro, t0, suffix = c("", ".t0avg"), by = join_by(condition))
data_gro <- data_gro |>
  mutate(day = ceiling(time / (60 * 60 * 24))) |>
  mutate(hour = ceiling(time / (60 * 60))) |>
  mutate(length.t0norm = length / length.t0avg,
         area.t0norm = area / area.t0avg,
         volume.t0norm = volume / volume.t0avg
  ) |>
  select(condition, step, time, day, hour, worm_id,
    length, length.t0norm,
    area, area.t0norm,
    volume, volume.t0norm
    )

#data_gro$condition <- substr(data_gro$condition, 11, 100)
#data_gro$condition <- str_replace(data_gro$condition, "\\?", "\U03BC")

# Fit sigmoid and plot
# 1. Group and summarize data
summ_area <- data_gro |>
    group_by(condition, hour) |>
    summarise(
        mean = mean(area.t0norm),
        sd = sd(area.t0norm),
        n = n(),
        .groups = "drop"
    )

# 2. Fit model per condition

fit_area <- summ_area |>
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
    )

# We could compute R2 to see if they fit similarly. It is not, however, goodness of fit. 

# 3. Unnest and plot
plot_area <- fit_area |>
  select(condition, data, fitted) |>
  unnest(c(data, fitted))

ggplot(plot_area, aes(x = hour, y = mean, color = condition)) + 
  geom_point(size = 1) +
  geom_line(aes(y = fitted))
  #theme_minimal() +
  #geom_errorbar(aes(x=hour,ymin=mean-sd, ymax=mean+sd, color = condition))

# 4. Compare models
comp_area <- fit_area |>
  select(condition, model) |>
  mutate(start_vals = map(model, coef))

# Get start values
control <- "Water"

a <- comp_area |>
  filter(str_detect(condition, control)) |>
  pull(start_vals) |>
  pluck(1)

b <- comp_area |>
  filter(!str_detect(condition, control)) |>
  mutate(diff = map(start_vals, ~ .x - a)) |>
  mutate(named_vals = map2(start_vals, condition, ~ set_names(.x, paste(names(.x), .y, sep = ".")))) |>
  select(named_vals)
# condition is grouping, maybe remove.

set_names(a, paste(names(a), comp_area$condition[1], sep = "."))

data_surv_set <- data_surv |> filter(str_detect(condition, conds))

start_vals <- A

# Get null model
area_null <- gnls()

# Get variable models