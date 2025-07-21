# Survival analysis for L4 Nagi data
# Run everything IN ORDER except what is marked as [O], which means optional
# [CUSTOM] means that the line can be changed to suit conditions or controls
# in the experiment

# Load necessary libraries
library(tidyverse)
library(readxl)
library(plotly)
library(survival)
library(survminer)

# Get surv data file
# Takes all files in the folder data/ which contain "survival" in the name
path_surv <- list.files("data", pattern = "survival", full.names = TRUE)

if(length(path_surv) < 1){
  stop("There is no suitable survival file.")
}

# Pivot function
# Transforms data from [condition1_time, condition1_event, condition2_time, ...]
# to [condition, time, event]
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

# Load excel files, fix headers
# Loads each survival file in a table, contained as a list in all_tables
# Cleans data so it has headers and it is in the format
# [condition1_time, condition1_event, condition2_time, condition2_event, ...]

all_tables <- map(path_surv, function(path){
  table_surv <- read_excel(path, col_names = FALSE) |> select(-1)

  headers_x <- table_surv[1, ] |> unlist(use.names = FALSE) |> zoo::na.locf()
  headers_y <- rep(c("time", "event"), times = length(headers_x) / 2)
  headers <- paste(headers_x, headers_y, sep = "_")

  table_surv <- table_surv |> slice(4:nrow(table_surv))
  names(table_surv) <- headers

  table_surv
})

# Check if the biological replicates (files in path_surv) are pool-able
# by a survival log-rank test (if p > 0.05, they are similar enough)

# Only check if there is more than one file
if(length(path_surv) > 1){

  # Set the control (or condition to be compared)
  # Using paste with collapse ".*" we match a regular expression.
  # We need to match first what comes first. Since control is
  # N2 | OP50 100 % | Water 10 %, we match N2, OP50 100, and Water, in that order
  control <- paste(c("N2", "OP50 100", "Water"), collapse = ".*") # [CUSTOM]

  # Takes the experiment ID from the file name (expID_survival...)
  # from the beginning up until and including the first _
  cont_tab <- map2(all_tables, str_extract(basename(path_surv), "^[^_]+_"),
                   function(table, id){
                     tab <- table |> select(matches(control))
                     tab <- tab |> setNames(paste0(id, names(tab)))
                     tab
                   })

  # Combines all tables from cont_tab into one, uses surv_pivot function
  # explained above, and conducts Log-Rank test
  cont_tab <- bind_rows(cont_tab)
  data_con <- surv_pivot(cont_tab)
  cont_lrt <- survdiff(Surv(time, event) ~ condition, data = data_con)

  # Stop if log-rank test is significant
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

# [O] If log rank test was significant and it did not print the control_plot
# It can be printed by uncommenting and running the following line
# control_plot

# Combines all tables from all_tables into one, uses surv_pivot function
table_surv <- bind_rows(all_tables)
data_surv <- surv_pivot(table_surv)

# Plot survival: Kaplan-Meier and Log-rank test
# Test and plot ALL conditions
km_fit <- survfit(Surv(time, event) ~ condition,
                  data = data_surv, conf.type = "log")

#km_lrt <- survdiff(Surv(time, event) ~ condition, data = data_surv)

surv_plot_all <- ggsurvplot(
  km_fit,
  data = data_surv,
  conf.int = FALSE, # [CUSTOM] Set to TRUE for confidence intervals
  conf.int.style = "step",
  pval = TRUE, # [CUSTOM] Change to FALSE to not display p-value
  xlab = "Time (hour)", # [CUSTOM] Change to change the x axis label
  legend.title = "", # [CUSTOM] Change to change the legend title
  legend.labs = levels(factor(data_surv$condition)),
  #palette = "igv" # [CUSTOM] Change colors
)

surv_plotly <- ggplotly(surv_plot_all[[1]])

# Fix labels messed up by the palette
legend_labels <- levels(factor(data_surv$condition))

# Rename traces (this assumes one trace per group, typical in KM plots)
for (i in seq_along(legend_labels)) {
  surv_plotly$x$data[[i]]$name <- legend_labels[i]
  surv_plotly$x$data[[i]]$legendgroup <- legend_labels[i]
  surv_plotly$x$data[[i]]$hovertemplate <- sub("^.*?<extra>", paste0(legend_labels[i], "<extra>"), surv_plotly$x$data[[i]]$hovertemplate)
}

# [O] To see the plots, execute the following lines
surv_plot_all # For the regular ggplot
surv_plotly # For the plotly plot

# Save plotly plot as html
htmlwidgets::saveWidget(as_widget(surv_plotly), "results/survival/survival.html")


# [O] This calculates p-values and plots ONLY what is filtered through conds
# conds must include all conditions we want in the plot, order does not matter
# If we have "PLA 100" and "PLA 200" and we add PLA to conds, it takes both
conds <- paste(c("Water", "S-medium", "PLA"), collapse = "|")

data_surv_set <- data_surv |> filter(str_detect(condition, conds))

km_fit_set <- survfit(Surv(time, event) ~ condition,
                      data = data_surv_set, conf.type = "log")

#km_set_lrt <- survdiff(Surv(time, event) ~ condition, data = data_surv_set)

surv_plot_set <- ggsurvplot(
  km_fit_set,
  data = data_surv_set,
  conf.int = TRUE, # [CUSTOM] Set to FALSE for no confidence intervals
  pval = TRUE, # [CUSTOM] Change to FALSE to not display p-value
  xlab = "Time (hour)", # [CUSTOM] Change to change the x axis label
  legend.title = "", # [CUSTOM] Change to change the legend title
  legend = "right", # [CUSTOM] Legend position: if top it gets very long
  legend.labs = levels(factor(data_surv_set$condition))
)

# [O] To see the plot
surv_plot_set

# [O] To get a plotly
# Note: Plotly does not work well with confidence intervals,
# so it will throw warnings and will not plot them
surv_set_plotly <- ggplotly(surv_plot_set[[1]])

# Save plotly plot as html
htmlwidgets::saveWidget(as_widget(surv_set_plotly), "results/survival/survival_set.html")

# [O] Save plot of the set of conds
ggsave_workaround <- function(g){
  survminer:::.build_ggsurvplot(x = g,
                                surv.plot.height = NULL,
                                risk.table.height = NULL,
                                ncensor.plot.height = NULL)
}

g_to_save <- ggsave_workaround(surv_plot_set)

ggsave(filename = "results/survival/survival_set.png", plot = g_to_save,
       width = 25, height = 13, dpi = 1000, units = "cm")