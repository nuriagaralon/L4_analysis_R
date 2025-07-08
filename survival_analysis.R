library(tidyverse)
library(readxl)
library(survival)
library(survminer)
library(plotly)

# Get surv data file
setwd("C:\\Users\\ngarriga\\Documents\\SydLab-One\\L4_analysis_R")

path_surv <- list.files("data", pattern = "survival", full.names = TRUE)

# Pivot function
surv_pivot <- function(surv_data){
  surv_table <- surv_data |>
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
  surv_table$event <- as.numeric(surv_table$event)
  surv_table$time <- as.numeric(surv_table$time)
  surv_table
}

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

# Check if the replicates (files in path_surv) are pool-able
# by a survival log-rank test (if p > 0.05, they are similar enough)

if(length(path_surv) > 1){

# Always first what comes first. Since control is
# N2 | OP50 100 % | Water 10 %, we match N2 first
  control <- paste(c("N2", "Water"), collapse = ".*")

  cont_tab <- map2(all_tables, str_extract(basename(path_surv), "^[^_]+_"),
                   function(table, id){
                     tab <- table |> select(matches(control))
                     tab <- tab |> setNames(paste0(id, names(tab)))
                     tab
                   })

  cont_tab <- bind_rows(cont_tab)
  data_con <- surv_pivot(cont_tab)
  cont_lrt <- survdiff(Surv(time, event) ~ condition, data = data_con)

  if(cont_lrt$pvalue < 0.05){
    control_plot <- ggsurvplot(
      survfit(Surv(time, event) ~ condition, data = data_con),
      data = data_con,
      conf.int = TRUE,
      pval = TRUE,
      xlab = "Time (hour)",
      legend.title = ""
    )
    print(control_plot)
    stop(paste("Replicates are too different, LRT pvalue =", cont_lrt$pvalue))
  }
} else if (length(path_surv < 1)) {
   stop("There is no suitable survival file.")
}

# bind rows
table_surv <- bind_rows(all_tables)

# Plot survival
data_surv <- surv_pivot(table_surv)

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

# p-values with everything, and plots everything.
surv_plot_all <- ggsurvplot(
  km_fit,
  data = data_surv,
  conf.int = FALSE,
  conf.int.style = "step",
  pval = TRUE,
  xlab = "Time (hour)",
  legend.title = "",
  legend.labs = levels(factor(data_surv$condition)),
  palette = "viridis"
)

surv_plotly <- ggplotly(surv_plot_all[[1]])

htmlwidgets::saveWidget(as_widget(surv_plotly), "results/survival.html")


# Save plot
ggsave_workaround <- function(g){
  survminer:::.build_ggsurvplot(x = g,
                                surv.plot.height = NULL,
                                risk.table.height = NULL,
                                ncensor.plot.height = NULL)
}

g_to_save <- ggsave_workaround(surv_plot_set)

ggsave(filename = "results/survival_set.png", plot = g_to_save,
       width = 17, height = 15, dpi = 1000, units = "cm")