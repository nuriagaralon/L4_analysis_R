library(tidyverse)

# Get motility data file
path_mot <- list.files("old_data", pattern = "motility_analysis_raw_scaled_per_day_only-alive", full.names = TRUE)

# Load table
table_mot <- read.csv(file = path_mot, header = TRUE) |> select(-alive)

# Sort data to plot
t0 <- table_mot |>
  filter(step < 7) |>
  select(-chip, -channel, -chamber, -step, -time, -worm_id) |>
  group_by(condition) |>
  summarise_all(mean)

data_mot <- full_join(table_mot, t0, suffix = c("", ".t0avg"), by = join_by(condition))
data_mot <- data_mot |>
  mutate(days = ceiling(time / (60 * 60 * 24))) |>
  mutate(ampl_head.t0norm = ampl_head / ampl_head.t0avg,
         ampl_mid.t0norm = ampl_mid / ampl_mid.t0avg,
         ampl_tail.t0norm = ampl_tail / ampl_tail.t0avg,
         vel.t0norm = vel / vel.t0avg,
         freq_m.t0norm = freq_m / freq_m.t0avg,
         curvature.t0norm = curvature / curvature.t0avg,
         area.t0norm = area / area.t0avg
  ) |>
  select(condition, step, time, days, worm_id,
    ampl_head, ampl_head.t0norm,
    ampl_mid, ampl_mid.t0norm,
    ampl_tail, ampl_tail.t0norm,
    vel, vel.t0norm,
    freq_m, freq_m.t0norm,
    curvature, curvature.t0norm,
    area, area.t0norm
    )

data_mot$condition <- substr(data_mot$condition, 11, 100)
data_mot$condition <- str_replace(data_mot$condition, "\\?", "\U03BC")

# Autoplots: are the residuals normal?
test_lm <- lm(formula = area.t0norm ~ days * condition, data=data_mot)
test_lm2 <- lm(formula = ampl_head.t0norm ~ days * condition, data=data_mot)
test_lm3 <- lm(formula = ampl_mid.t0norm ~ days * condition, data=data_mot)
test_lm4 <- lm(formula = ampl_tail.t0norm ~ days * condition, data=data_mot)
test_lm5 <- lm(formula = vel.t0norm ~ days * condition, data=data_mot)
test_lm6 <- lm(formula = freq_m.t0norm ~ days * condition, data=data_mot)
test_lm7 <- lm(formula = curvature.t0norm ~ days * condition, data=data_mot)

# Fit sigmoid
plot_mot <- data_mot |> filter(str_detect(condition, "FUdR")) |> group_by(step) |> summarise(mean = mean(area.t0norm), stdev = sd(area.t0norm), n = n())

fit_mot <- nls(mean ~ A / (1 + exp(-B * (step - C))), data = plot_mot, start = list(A = 1, B = 1, C = 2))
summary(fit_mot)

plot_mot$y_pred <- predict(fit_mot)

# Plot original data and fitted sigmoidal curve
ggplot(plot_mot, aes(x = step, y = mean)) +
  geom_point() +
  geom_line(aes(x = step, y = y_pred))








fits <- data_mot |>
  group_by(condition, step) |>
  summarise(
    mean = mean(area.t0norm),
    stdev = sd(area.t0norm),
    n = n()
  ) |>
  group_by(condition) |>
  nest() |>
  mutate(
    fit = map(data, ~ {
      # Fit the sigmoidal model to each condition's data
      nls(mean ~ A / (1 + exp(-B * (step - C))),
          data = .x, 
          start = list(A = max(.x$mean), B = 0.1, C = median(.x$step)))
    }),
    predictions = map2(fit, data, ~ {
      data_extended <- data.frame(step = seq(min(.x$step), max(.x$step), length.out = 100))
      data_extended$y_pred <- predict(.y, newdata = data_extended)
      data_extended
    })
  ) |>
  unnest(predictions)

fits_good <- fits |> 
  mutate(model = map(data, function(df) nls(mean ~ A / (1 + exp(-B * (step - C))), 
  data = df, start = list(A = max(df$mean), B = 0.1, C = median(df$step)))))

nls(mean ~ A / (1 + exp(-B * (step - C))), data = .x, 
          start = )

# Plot the observed data and the fitted curves
ggplot(fits, aes(x = step, y = mean, color = condition)) +
  geom_point() +
  geom_line(aes(x = step, y = y_pred), size = 1) +
  labs(title = "Sigmoidal Fits for Each Condition",
       x = "Step", y = "Mean Area") +
  theme_minimal() +
  theme(legend.title = element_blank())















# Area plots, need fixing

plot_100 <- plot_mot |> filter(grepl("100|Water", condition))


plot_mot <- data_mot |>
  select(condition, days, area.t0norm) |>
  group_by(condition, days) |>
  summarise(mean = mean(area.t0norm), stdev = sd(area.t0norm), n = n()) |> mutate(SEM = stdev/sqrt(n))


ggplot(plot_mot, aes(days, mean, colour = condition)) +
  geom_point() + geom_line() +
  geom_ribbon(aes(x = days, ymin=mean-stdev, ymax=mean+stdev, color = condition), alpha=.01) +
  scale_x_continuous(breaks = 0:max(plot_mot$days))

ggplot(plot_mot, aes(days, mean, colour = condition)) + geom_point() + geom_line() + geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev)) + scale_x_continuous(breaks = 0:max(plot_mot$days))


# Ribbon plot
ribbon_data <- plot_mot |> group_by(days) |> summarise(ymax = max(mean+stdev), ymin = min(mean-stdev))

ggplot(plot_mot, aes(x = days)) +
  geom_point(aes(y = mean, colour = condition)) + geom_line(aes(y = mean, colour = condition)) +
  geom_ribbon(data = ribbon_data, aes(ymin = ymin, ymax = ymax), alpha = .3) +
  scale_x_continuous(breaks = 0:max(plot_mot$days))
  # but should i plot ribbon with stdev or sderror? 2*stdv/sqrt(n)

aov_result <- aov(area.t0norm ~ condition * days + Error(worm_id/days), data = data_mot)





dnt_result <- DunnettTest(area.t0norm ~ condition, data = data_mot, control = "Water 10 %")



aov_result <- aov(area.t0norm ~ condition * days, data = data_mot)

dunnett_result <- glht(aov_result, linfct = mcp(condition = "Dunnett", days = "Dunnett"))

# maybe some need to be set as factor and the control releveled.
data_mot$condition <- relevel(data_mot$condition, ref = "Water 10 %")
data_mot$days = as.factor(data_mot$days)
data_mot$condition = as.factor(data_mot$condition)


# Amplitude plots: head

plot_mot <- data_mot |> select(condition, days, ampl_head.t0norm) |> group_by(condition, days) |> summarise(mean = mean(ampl_head.t0norm), stdev = sd(ampl_head.t0norm), n = n()) 


ggplot(plot_mot, aes(x = days, y = condition, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of Head Amplitude by Time and Condition", x = "Time", y = "Condition") +
  theme_minimal()

# Velocity plot

# Frequence plot

# Curvature plot