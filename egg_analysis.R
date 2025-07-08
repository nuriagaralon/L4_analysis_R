library(tidyverse)
library(readxl)
library(DescTools)
library(ggsignif)

#
setwd("C:\\Users\\ngarriga\\Documents\\SydLab-One\\L4_analysis_R")

# General raw data is confusing so we are using the data
# per channel (when using channels as replicates),
# per experiment (when using experiments as replicates)
path_egg <- list.files("data", pattern = "eggs_channel", full.names = TRUE)

# Load tables, fix headers

table_egg_start <- read_excel(path_egg, sheet = )

all_tables <- map(path_surv, function(path){
  table_surv <- read_excel(path, col_names = FALSE) |> select(-1)

  headers_x <- table_surv[1, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
  headers_y <- rep(c("time", "event"), times = length(headers_x) / 2)
  headers <- paste(headers_x, headers_y, sep = "_")

  table_surv <- table_surv |> slice(4:nrow(table_surv))
  names(table_surv) <- headers

  table_surv
})

# Start and end are not showing up in the general raw data
# So we need to combine all conditions

# Get all egg files in one table, with condition column added
path_egg <- list.files("old_data/egg_stats", full.names = TRUE)

# Read all files into a list of data frames
table_egg_list <- lapply(path_egg, function(file){
  df <- read_excel(file) |>
    mutate(condition = str_replace_all(file, c("^.*_" = "", "\\.xls." = "")))
  return(df)
})

# Combine all data frames into one table, set condition as factor
table_egg <- do.call(rbind, table_egg_list) |>
  mutate(condition = as.factor(condition))

# New eggs:
path_egg <- list.files("data", pattern = "eggs_condition", full.names = TRUE)
path_egg_raw <- list.files("data", pattern = "eggs_raw", full.names = TRUE)

table_egg_emerg <- read_excel(path_egg, col_names = FALSE, sheet = 1)
table_egg_count <- read_excel(path_egg, col_names = FALSE, sheet = 2)
table_egg_raw <- read_excel(path_egg_raw, col_names = TRUE)

# Fix egg_emerg
headers_x <- table_egg_emerg[1,] |> select(-1) |> unlist(use.names = FALSE) |> zoo::na.locf()
headers_y <- table_egg_emerg[2,] |> select(-1) |> unlist(use.names = FALSE)
headers <- paste(headers_x, headers_y, sep = "_")

table_egg_emerg <- table_egg_emerg |> slice(4:nrow(table_egg_emerg)) |> column_to_rownames("...1")
names(table_egg_emerg) <- headers

# slice(table_egg_emerg, 2) to get the hours, and slice 1 to get the steps.

# Multiple One-way ANOVAS
# But do we want to check interactions between start/end/eggs_per_worm?
# If so we might want two or three-way ANOVAS?

# + post-hoc Dunnett's Test (Water 10% as control)

sig_pval <- 0.05

sink("results/anova_egg_start.txt")
cat("One-way ANOVA summary for 'start' (first step with egg presence):\n\n")

fit_start <- aov(start ~ condition, data = table_egg)
summary(fit_start)

if (summary(fit_start)[[1]]$`Pr(>F)`[1] < sig_pval){
  dnt_start <- DunnettTest(start ~ condition, data = table_egg, control = "Water 10 %")
  dnt_start_res <- as.data.frame(dnt_start[[1]])
  dnt_start_res <- dnt_start_res |>
    filter(pval < sig_pval) |>
    rownames_to_column(var = "condition") |>
    mutate(sigval = case_when(pval < 0.001 ~ "***",
                              pval < 0.01 ~ "**",
                              pval < 0.05 ~ "*",
                              pval >= 0.05 ~ "NS"))
  print(dnt_start)
  cat("\n")
} else {
  cat("\nNo significant results from ANOVA. No Dunnett Test performed.\n")
}
sink()

sink("results/anova_egg_end.txt")
cat("One-way ANOVA summary for 'end' (last step with egg presence):\n\n")

fit_end <- aov(end ~ condition, data = table_egg)
summary(fit_end)

if (summary(fit_end)[[1]]$`Pr(>F)`[1] < sig_pval){
  dnt_end <- DunnettTest(end ~ condition, data = table_egg, control = "Water 10 %")
  dnt_end_res <- as.data.frame(dnt_end[[1]])
  dnt_end_res <- dnt_end_res |>
    filter(pval < sig_pval) |>
    rownames_to_column(var = "condition") |>
    mutate(sigval = case_when(pval < 0.001 ~ "***",
                              pval < 0.01 ~ "**",
                              pval < 0.05 ~ "*",
                              pval >= 0.05 ~ "NS"))
  print(dnt_end)
  cat("\n")
} else {
  cat("\nNo significant results from ANOVA. No Dunnett Test performed.\n")
}
sink()

sink("results/anova_egg_eggs.txt")
cat("One-way ANOVA summary for 'eggs_per_worm' (eggs/worms in chamber):\n\n")

fit_eggs <- aov(eggs_per_worm ~ condition, data = table_egg)
summary(fit_eggs)

if (summary(fit_eggs)[[1]]$`Pr(>F)`[1] < sig_pval){
  dnt_eggs <- DunnettTest(eggs_per_worm ~ condition, data = table_egg, control = "Water 10 %")
  dnt_eggs_res <- as.data.frame(dnt_eggs[[1]])
  dnt_eggs_res <- dnt_eggs_res |>
    filter(pval < sig_pval) |>
    rownames_to_column(var = "condition") |>
    mutate(sigval = case_when(pval < 0.001 ~ "***",
                              pval < 0.01 ~ "**",
                              pval < 0.05 ~ "*",
                              pval >= 0.05 ~ "NS"))
  print(dnt_eggs)
  cat("\n")
} else {
  cat("\nNo significant results from ANOVA. No Dunnett Test performed.\n")
}
sink()

# Plots
plot_data <- table_egg |>
  select(condition, start, end, eggs_per_worm) |>
  group_by(condition) |>
  summarise_all(list(mean = "mean",  sd = "sd", n = "length", sem = ~sd(.)/sqrt(length(.)))) |>
  select(-end_n, -eggs_per_worm_n) |>
  rename(n = start_n) |>
  select(condition, start_mean, start_sd, start_sem, end_mean, end_sd, end_sem,
         eggs_per_worm_mean, eggs_per_worm_sd, eggs_per_worm_sem, n)

write.csv(plot_data, "results/eggs_plot_data.csv")

# Start
start_plot <- ggplot(plot_data, aes(condition, start_mean)) +
  geom_col(fill = "deepskyblue") +
  geom_errorbar(aes(ymin = start_mean - start_sd,
                    ymax = start_mean + start_sd)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "", y = "First step with eggs")

if (exists("dnt_start_res")) {
  comparisons <- str_split(dnt_start_res$condition, pattern = regex("-(?=[^-]+$)"))
  y_positions <- seq(from = max(plot_data$start_mean+plot_data$start_sd) + 1, by = 50, length.out = length(comparisons))
  start_plot <- start_plot + geom_signif(comparison = comparisons,
                                         annotations = dnt_start_res$sigval,
                                         y_position = y_positions,
                                         tip_length = 0)
}

ggsave(filename = "results/eggs_start_plot.png", plot = start_plot,
       width = 17, height = 15, dpi = 1000, units = "cm")

# Alternate start plot

start_plot <- ggplot(plot_data, aes(x = condition, y = start_mean)) +
  geom_jitter(aes(x = condition, y = start, color = condition), data = table_egg, show.legend = FALSE) +
  geom_crossbar(aes(ymin = start_mean-start_sem, ymax = start_mean+start_sem, color= condition),  show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = "", y = "First step with eggs")

if (exists("dnt_start_res")) {
  comparisons <- str_split(dnt_start_res$condition, pattern = regex("-(?=[^-]+$)"))
  y_positions <- seq(from = max(table_egg$start) + 1, by = 50, length.out = length(comparisons))
  start_plot <- start_plot + geom_signif(comparison = comparisons,
                                         annotations = dnt_start_res$sigval,
                                         y_position = y_positions,
                                         tip_length = 0.01)
}

# End

end_plot <- ggplot(plot_data, aes(condition, end_mean)) +
  geom_col(fill = "deepskyblue") +
  geom_errorbar(aes(ymin = end_mean - end_sd,
                    ymax = end_mean + end_sd)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "", y = "Last step with eggs")

if (exists("dnt_end_res")) {
  comparisons <- str_split(dnt_end_res$condition, pattern = regex("-(?=[^-]+$)"))
  y_positions <- seq(from = max(plot_data$end_mean+plot_data$end_sd) + 1, by = 50, length.out = length(comparisons))
  end_plot <- end_plot + geom_signif(comparison = comparisons,
                                     annotations = dnt_end_res$sigval,
                                     y_position = y_positions,
                                     tip_length = 0)
}

ggsave(filename = "results/eggs_end_plot.png", plot = end_plot,
       width = 17, height = 15, dpi = 1000, units = "cm")

# Alternate end plot

end_plot <- ggplot(plot_data, aes(x = condition, y = end_mean)) +
  geom_jitter(aes(x = condition, y = end, color = condition), data = table_egg, show.legend = FALSE) +
  geom_crossbar(aes(ymin = end_mean-end_sem, ymax = end_mean+end_sem, color= condition),  show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = "", y = "Last step with eggs")

if (exists("dnt_end_res")) {
  comparisons <- str_split(dnt_end_res$condition, pattern = regex("-(?=[^-]+$)"))
  y_positions <- seq(from = max(table_egg$end) + 1, by = 50, length.out = length(comparisons))
  end_plot <- end_plot + geom_signif(comparison = comparisons,
                                         annotations = dnt_end_res$sigval,
                                         y_position = y_positions,
                                         tip_length = 0.01)
}

# Eggs per worm

egg_plot <- ggplot(plot_data, aes(condition, eggs_per_worm_mean)) +
  geom_col(fill = "deepskyblue") +
  geom_errorbar(aes(ymin = eggs_per_worm_mean - eggs_per_worm_sd,
                    ymax = eggs_per_worm_mean + eggs_per_worm_sd)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "", y = "Eggs per worm")

if (exists("dnt_eggs_res")) {
  comparisons <- str_split(dnt_eggs_res$condition, pattern = regex("-(?=[^-]+$)"))
  y_positions <- seq(from = max(plot_data$eggs_per_worm_mean+plot_data$eggs_per_worm_sd) + 1, by = 50, length.out = length(comparisons))
  egg_plot <- egg_plot + geom_signif(comparison = comparisons,
                                     annotations = dnt_eggs_res$sigval,
                                     y_position = y_positions,
                                     tip_length = 0)
}

ggsave(filename = "results/eggs_eggs_per_worm_plot.png", plot = egg_plot,
       width = 17, height = 15, dpi = 1000, units = "cm")

# Alternate eggs per worm plot

egg_plot <- ggplot(plot_data, aes(x = condition, y = eggs_per_worm_mean)) +
  geom_jitter(aes(x = condition, y = eggs_per_worm, color = condition), data = table_egg, show.legend = FALSE) +
  geom_crossbar(aes(ymin = eggs_per_worm_mean-eggs_per_worm_sem, ymax = eggs_per_worm_mean+eggs_per_worm_sem, color= condition),  show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = "", y = "Eggs per worm")

if (exists("dnt_eggs_res")) {
  comparisons <- str_split(dnt_eggs_res$condition, pattern = regex("-(?=[^-]+$)"))
  y_positions <- seq(from = max(table_egg$eggs_per_worm) + 1, by = 50, length.out = length(comparisons))
  egg_plot <- egg_plot + geom_signif(comparison = comparisons,
                                         annotations = dnt_eggs_res$sigval,
                                         y_position = y_positions,
                                         tip_length = 0.01)
}
