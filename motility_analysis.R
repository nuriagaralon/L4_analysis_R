# Growth analysis for L4 Nagi data
# Run everything IN ORDER except what is marked as [O], which means optional
# [CUSTOM] means that the line can be changed to suit conditions or controls
# in the experiment

# Load necessary libraries
library(ggfortify)
library(tidyverse)
library(readxl)
library(plotly)
library(afex)


library(pracma)
library(DescTools)
library(rstatix)
library(ggpubr)
library(ggsci)
#library(nlme) #[O], explained later

# Get motility data file
# Takes all files in the folder data/ which contain "motility_filtered_raw" in the name
path_mot <- list.files("data", pattern = "motility_filtered_raw", full.names = TRUE)

if(length(path_mot) < 1){
  stop("There is no suitable motility file.")
}

# Load excel file
# Selects only area, volume, length (and not eccentricity, orientation, etc.)
table_mot <- map(path_mot, function(path){
  read_excel(path, col_names = TRUE) |> mutate(exp_id = str_extract(basename(path), "^[^_]+"))
})

table_mot <- bind_rows(table_mot)

# [CUSTOM] Filter out FUdR data. Can be removed or changed for another condition
table_mot <- table_mot |> filter(!str_detect(condition, "FUdR"))

## [CUSTOM] Filter worms with less than 10 measurements
#filter_number <- 10
#
#valid_worms <- table_mot |>
#    group_by(exp_id, chip, channel, chamber, worm_id) |>
#    summarise(worm_id_count = n(), .groups = "drop") |>
#    filter(worm_id_count >= filter_number)
#
#data_mot <- table_mot |> inner_join(valid_worms, by = c("exp_id", "chip", "channel", "chamber", "worm_id"))
#data_mot <- data_mot |> select(-worm_id_count)

data_mot <- table_mot

# Do we want to use worms as replicates? What about channels?
worm_rep <- FALSE
channels <- TRUE

if(length(path_mot) > 1){
  if(worm_rep){
    rep_cols <- c("exp_id", "chip", "channel", "chamber", "worm_id")
    replicate <- "Modelled using worm as replicate"
  } else if (channels) {
    rep_cols <- c("exp_id", "chip", "channel")
    replicate <- "Modelled using experiment-chip-channel as replicate."
  } else if (!channels) {
    rep_cols <- c("exp_id", "condition")
    replicate <- "Modelled using experiment as replicate."
  }
# If only one file, replicate ID is the chip_channel combination
} else if (length(path_mot) == 1) {
  if(worm_rep){
    rep_cols <- c("chip", "channel", "chamber", "worm_id")
    replicate <- "Modelled using worm as replicate"
  } else {
  rep_cols <- c("chip", "channel")
  replicate <- "Modelled using chip-channel as replicate."
  }
}

# Sort data table
data_mot <- data_mot |>
  mutate(
    rep_id = pmap_chr(across(all_of(rep_cols)), ~ paste(..., sep = "_")),
    day = ceiling(time / (60 * 60 * 24)),
    hour = ceiling(time / (60 * 60)),
    h_nr = time / (60 * 60)
  ) |>
  select(exp_id, rep_id, condition, step, time, day, hour, h_nr,
    head_amplitude, tail_amplitude, displacement_speed, bodybends_frequency
  )

#
plot_motility <- function(data, variable){
  plotdata <- data |>
    group_by(day, condition) |>
    summarise(
      mean = mean({{variable}}),
      std = sd({{variable}}),
      n = n()
    ) |>
    mutate(
      sem = std / sqrt(n)
    )
  ggplot(plotdata, aes(x = day, y = mean, color = condition)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem)) +
    scale_color_igv()  + #[CUSTOM] Color scale can be changed.
    xlab("Time (days)") + # [CUSTOM] Change to change the x axis label
    labs(color = "Condition")
}


ha_plot <- plot_motility(data_mot, head_amplitude) + ylab("Head amplitude (mm)")
ha_plotly <- ggplotly(ha_plot)

htmlwidgets::saveWidget(as_widget(ha_plotly), "results/motility/motility_head_amplitude.html")

ta_plot <- plot_motility(data_mot, tail_amplitude) + ylab("Tail amplitude (mm)")
ta_plotly <- ggplotly(ta_plot)

htmlwidgets::saveWidget(as_widget(ta_plotly), "results/motility/motility_tail_amplitude.html")

ds_plot <- plot_motility(data_mot, displacement_speed) + ylab("Displacement speed (mm/s)")
ds_plotly <- ggplotly(ds_plot)

htmlwidgets::saveWidget(as_widget(ds_plotly), "results/motility/motility_displacement_speed.html")

bf_plot <- plot_motility(data_mot, bodybends_frequency) + ylab("Bodybends frequency (Hz)")
bf_plotly <- ggplotly(bf_plot)

htmlwidgets::saveWidget(as_widget(bf_plotly), "results/motility/motility_bodybends_frequency.html")

# Get data for statistical analysis
# This is to analyze the full dataset
data_testmot <- data_mot |>
  group_by(condition, rep_id, day) |>
  summarise(
    ha = mean(head_amplitude),
    ta = mean(tail_amplitude),
    ds = mean(displacement_speed),
    bf = mean(bodybends_frequency)
  )

# This function fills empty values with 0, and can also filter for experiment
from_data_testmot_fill <- function(data_testmot, id_pattern = NULL){
  if(!is.null(id_pattern)){
    smalldata <- data_testmot |> filter(str_detect(rep_id, id_pattern))
  } else {
    smalldata <- data_testmot
  }
  condition_map <- smalldata |> distinct(rep_id, condition)
  smalldata <- smalldata |> ungroup() |> select(-condition) |>
    tidyr::complete(rep_id, day, fill = list(ha = 0, ta = 0, ds = 0, bf = 0))

  smalldata <- smalldata |>
    left_join(condition_map, by = "rep_id")
}

# So if we want to use all data as is, use data_testmot
# If we want to use data with missing values set to 0, separated by experiment
data_wiz <- from_data_testmot_fill(data_testmot, "wiz")
data_CGa <- from_data_testmot_fill(data_testmot, "CGa")
data_fgo <- from_data_testmot_fill(data_testmot, "fgo")
# If we want to use all data with missing values set to 0
data_all <- from_data_testmot_fill(data_testmot)

# the mixed-measures two way anova deals with unbalanced data by erasing
# the row of data. So filling it avoids the data being deleted. So then, for example
# if one data is longer in time than the others, maybe we should either let it get cut
# Or fill the shorter one, or analyse them separately and then see if the conclusions are
# the same.

# Save statistical tests to file
dataset_list <- list(
  wiz = data_wiz,
  CGa = data_CGa,
  fgo = data_fgo,
  all = data_all
)

for(replicate in names(dataset_list)){
  df <- dataset_list[[replicate]]

  for(variable in c("ha", "ta", "ds", "bf"){
  # Sink (save) results, say which 
  sink(paste0("results/motility/motility_", replicate, "_", variable, ".txt"))
  




  }
}




aov_result <- aov_ez(id = "rep_id",
                      dv = "ha",
                      data = placeh,
                      within = "day",
                      between = "condition")

placeh_young <- placeh |> filter(day <=10)
placeh_old <- placeh |> filter(day >10)

aov_young <- aov_ez(id = "rep_id",
                     dv = "ha",
                     data = placeh_young,
                     within = "day",
                     between = "condition")

aov_old <- aov_ez(id = "rep_id",
                     dv = "ha",
                     data = placeh_old,
                     within = "day",
                     between = "condition")

> emm <- emmeans(aov_young, ~ condition | day)
> emm2 <- pairs(emm, adjust = "holm") 
> str(emm2)
'emmGrid' object with variables:
    contrast = 1, 2, 3, 4, 5, 6, 7, 8, 9, ..., 120
    day = multivariate response levels: 1, 2, 3, 4, 5, 6, 7, 8, 9, ..., 10
> emm2 <- as.data.frame(emm2)
> View(emm2)
> emm2 <- emm2 |> filter(p.value < 0.01)