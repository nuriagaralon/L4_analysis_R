library(tidyverse)
library(readxl)
library(DescTools)
library(ggsignif)

# General raw data is confusing so we are using the data
# per channel (when using channels as replicates),
# per experiment (when using experiments as replicates)
path_egg <- list.files("data", pattern = "eggs_condition", full.names = TRUE)

# FIRST: EGG EMERGENCE
# Load tables, fix headers
start_tables <- map(path_egg, function(path){
  table <- read_excel(path, col_names = FALSE, sheet = 1)

  nanum <- sum(is.na(table[[1]]))
  
  headers_x <- table[1, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
  headers_y <- table[nanum - 1, ] |> unlist(use.names = FALSE) |> discard(is.na)

  if(nanum == 4){
    headers_z <- table[nanum - 2, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
  }else{
    headers_z <- str_extract(basename(path), "^[^_]+")
  }
  headers <- paste(headers_x, headers_z, headers_y, sep = "_")

  table <- table |> slice((nanum + 1):nrow(table))
  names(table) <- c("Type", headers)

  table <- table |>
    pivot_longer(
      cols = -Type,
      names_to = c("Condition", "Summary"),
      names_pattern = "(.*)_(mean|std|N_chambers)",
      values_to = "value"
    ) |>
    pivot_wider(
      names_from = c("Type", "Summary"),
      values_from = value,
    ) |>
    mutate(across(-Condition, as.numeric)) |>
    filter(!str_detect(Condition, "FUdR"))
  table
})

egg_table <- bind_rows(start_tables) |>
  separate(
    col = "Condition",
    into = c("Condition", "rep_id"),
    sep = "_"
  )

egg_table$Condition <- as.factor(egg_table$Condition)

# Check for normality
ggplot(egg_table, aes(sample = `Emergence hour_mean`)) +
     stat_qq() +
     stat_qq_line() +
     theme_minimal()

shapiro.test(egg_table$`Emergence hour_mean`)

# One-way ANOVA for egg emergence
sig_pval <- 0.05
control_name <- str_subset(unique(egg_table$Condition), "Water")

sink("results/anova_egg_emergence.txt")
cat("One-way ANOVA summary for egg emergence hour:\n\n")

fit_start <- aov(`Emergence hour_mean` ~ Condition, data = egg_table)
summary(fit_start)

if (summary(fit_start)[[1]]$`Pr(>F)`[1] < sig_pval){
  dnt_start <- DunnettTest(`Emergence hour_mean` ~ Condition, data = egg_table, control = control_name)
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

# Plot
# Mean, std, n of replicates
egg_eh <- egg_table |>
  group_by(Condition) |>
  summarise(
    mean = mean(`Emergence hour_mean`),
    std = sd(`Emergence hour_mean`),
    n = n()
  ) |>
  mutate(
    sem = std / sqrt(n)
  )

egg_eh_plot <- ggplot(egg_eh, aes(Condition, mean)) +
  geom_boxplot(fill = "deepskyblue") +
  geom_errorbar(aes(ymin = mean - sem,
                    ymax = mean + sem)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "", y = "Egg emergence (hours)")

if (exists("dnt_start_res")) {
  comparisons <- str_split(dnt_start_res$Condition, pattern = regex("-(?=[^-]+$)"))
  y_positions <- seq(from = max(egg_eh$mean+egg_eh$sem) + 1, by = 50, length.out = length(comparisons))
  egg_eh_plot <- egg_eh_plot + geom_signif(comparison = comparisons,
                                         annotations = dnt_start_res$sigval,
                                         y_position = y_positions,
                                         tip_length = 0)
}

ggsave(filename = "results/egg_emergence_plot.png", plot = egg_eh_plot,
       width = 17, height = 15, dpi = 1000, units = "cm")

# Alternate plot

egg_eh_plot <- ggplot(egg_eh, aes(x = Condition, y = mean)) +
  geom_point(aes(x = Condition, y = `Emergence hour_mean`, color = Condition), data = egg_table, show.legend = FALSE) +
  geom_crossbar(aes(ymin = mean-sem, ymax = mean+sem, color= Condition),  show.legend = FALSE) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = "", y = "Egg emergence (hours)")

if (exists("dnt_start_res")) {
  comparisons <- str_split(dnt_start_res$condition, pattern = regex("-(?=[^-]+$)"))
  y_positions <- seq(from = max(table_egg$start) + 1, by = 50, length.out = length(comparisons))
  egg_eh_plot <- egg_eh_plot + geom_signif(comparison = comparisons,
                                         annotations = dnt_start_res$sigval,
                                         y_position = y_positions,
                                         tip_length = 0.01)
}

# SECOND: EGG COUNT
# Load tables, fix headers
count_tables <- map(path_egg, function(path){
  table <- read_excel(path, col_names = FALSE, sheet = 2) |> select(-1)

  nanum <- sum(is.na(table[[1]]))
  
  headers_x <- table[1, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
  headers_y <- table[nanum, ] |> unlist(use.names = FALSE) |> discard(is.na)

  if(nanum == 3){
    headers_z <- table[nanum - 1, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
  }else{
    headers_z <- str_extract(basename(path), "^[^_]+")
  }
  headers <- c(headers_y[1:2], paste(headers_x, headers_z, headers_y[-c(1:2)], sep = "_"))

  table <- table |> slice((nanum + 2):nrow(table))
  names(table) <- headers

  table <- table |>
    pivot_longer(
      cols = -c("Step", "Hour"),
      names_to = c("Condition", "Summary"),
      names_pattern = "(.*)_(mean_E_per_worm_cond|std_E_per_worm_cond|N_chambers_cond)",
      values_to = "value"
    ) |>
    pivot_wider(
      names_from = Summary,
      values_from = value,
    ) |>
    mutate(across(-Condition, as.numeric)) |>
    filter(!str_detect(Condition, "FUdR"))
  table
})

egg_count_table <- bind_rows(count_tables) |>
  separate(
    col = "Condition",
    into = c("Condition", "rep_id"),
    sep = "_"
  )

egg_count_table$Condition <- as.factor(egg_count_table$Condition)

# Plot
# Mean, std, n of replicates
egg_counts <- egg_count_table |>
  group_by(Hour, Condition, rep_id) |>
  summarise(
    mean = mean(mean_E_per_worm_cond),
    std = sd(mean_E_per_worm_cond),
    n = n()
  ) |>
  mutate(
    sem = std / sqrt(n)
  )

# However, egg counts look quite different between replicates, maybe should just anova test
# the maximum value (max(mean)) or the nematode age (time at max(mean))
# like this paper: DOI: 10.1186/1472-6785-9-14
