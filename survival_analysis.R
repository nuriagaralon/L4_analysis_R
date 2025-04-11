library(tidyverse)
library(survival)
library(survminer)

# Get surv data file
setwd("C:\\Users\\ngarriga\\Documents\\SydLab-One\\L4_analysis_R")

path_surv <- list.files("data", pattern = "pred_deaths_worm", full.names = TRUE)

# Load table with proper headers
headers <- read.csv(path_surv, nrows = 2, header = FALSE)
headers <- sapply(headers, paste, collapse = "_")
table_surv <- read.csv(file = path_surv, skip = 2, header = FALSE)
names(table_surv) <- substr(headers, 11, 100)

# Sort data to plot
data_surv <- table_surv |>
  pivot_longer(
    cols = everything(),
    names_to = c("condition", "type"),
    names_pattern = "(.*)_(deaths|event)",
    values_to = "value"
  ) |>
  pivot_wider(
    names_from = type,
    values_from = value
  ) |>
  unnest(c(deaths, event)) |>
  drop_na()

# Plot survival
km_fit <- survfit(Surv(deaths, event) ~ condition, data = data_surv)

surv_plot <- ggsurvplot(
    km_fit,
    data = data_surv,
    pval = TRUE,
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