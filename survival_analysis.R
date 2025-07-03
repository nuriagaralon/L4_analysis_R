library(tidyverse)
library(readxl)
library(survival)
library(survminer)

# Get surv data file
setwd("C:\\Users\\ngarriga\\Documents\\SydLab-One\\L4_analysis_R")

path_surv <- list.files("data", pattern = "survival", full.names = TRUE)

# Load old table, add headers
# headers <- read.csv(path_surv, nrows = 2, header = FALSE)
# headers <- sapply(headers, paste, collapse = "_")
# table_surv <- read.csv(file = path_surv, skip = 2, header = FALSE)
# names(table_surv) <- substr(headers, 11, 100)

# Load new table, fix headers

all_tables <- map(path_surv, function(path){
  table_surv <- read_excel(path, col_names = FALSE) |> select(-1)

  headers_x <- table_surv[1, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
  headers_y <- rep(c("time", "event"), times = length(headers_x) / 2)
  headers <- paste(headers_x, headers_y, sep = "_")

  table_surv <- table_surv |> slice(4:nrow(table_surv))
  names(table_surv) <- headers

  table_surv
})

table_surv <- bind_rows(all_tables)

# Sort data to plot
data_surv <- table_surv |>
  pivot_longer(
    cols = everything(),
    names_to = c("condition", "type"),
    names_pattern = "(.*)_(time|event)",
    values_to = "value"
  ) |>
  pivot_wider(
    names_from = type,
    values_from = value,
    values_fn = list
  ) |>
  unnest(c(time, event)) |>
  drop_na()

data_surv$event <- as.numeric(data_surv$event)
data_surv$time <- as.numeric(data_surv$time)

# Plot survival
conds <- paste(c("Water", "S-medium"), collapse = "|")

data_surv_set <- data_surv |> filter(str_detect(condition, conds))

km_fit <- survfit(Surv(time, event) ~ condition,
                  data = data_surv, conf.type = "log")

km_fit_set <- survfit(Surv(time, event) ~ condition,
                  data = data_surv_set, conf.type = "log")

# This calculates p-values and plots ONLY what is filtered through conds
surv_plot_set <- ggsurvplot(
  km_fit_set,
  data = data_surv_set,
  conf.int = TRUE,
  pval = TRUE,
  xlab = "Time (hour)",
  legend.title = "",
  legend.labs = levels(factor(data_surv_set$condition))
)

# This calculates p-values with everything, and plots conds
surv_plot_all_filtered <- ggsurvplot(
  km_fit_set,
  data = data_surv,
  conf.int = TRUE,
  pval = TRUE,
  xlab = "Time (hour)",
  legend.title = "",
  legend.labs = levels(factor(data_surv_set$condition))
)

# p-values with everything, and plots everything.
surv_plot_all <- ggsurvplot(
  km_fit,
  data = data_surv,
  conf.int = FALSE,
  conf.int.style = "step",
  pval = TRUE,
  xlab = "Time (hour)",
  legend.title = "",
  legend.labs = levels(factor(data_surv$condition))
)

# Save plot
ggsave_workaround <- function(g){
  survminer:::.build_ggsurvplot(x = g,
                                surv.plot.height = NULL,
                                risk.table.height = NULL,
                                ncensor.plot.height = NULL)
}

g_to_save <- ggsave_workaround(surv_plot)

ggsave(filename = "results/survival_all.png", plot = g_to_save,
       width = 17, height = 15, dpi = 1000, units = "cm")