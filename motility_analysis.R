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

# [CUSTOM] Do we want to use worms as replicates? What about channels?
worm_rep <- FALSE
channels <- TRUE

# [CUSTOM] Filter worms with less than 10 measurements (useful when using worm as replicate)
# Can change filter_number

# The problem with using worms as replicates is that, if they do not have all the values,
# there are problems afterwards with the ANOVA: it will remove the replicates that 
# don't have data for all timepoints.

if(worm_rep){
  filter_number <- 10

  valid_worms <- table_mot |>
      group_by(exp_id, chip, channel, chamber, worm_id) |>
      summarise(worm_id_count = n(), .groups = "drop") |>
      filter(worm_id_count >= filter_number)

  data_mot <- table_mot |> inner_join(valid_worms, by = c("exp_id", "chip", "channel", "chamber", "worm_id"))
  data_mot <- data_mot |> select(-worm_id_count)
} else {
  data_mot <- table_mot
}

# Set what is the replicate
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
# If only one file, replicate ID is the chip_channel combination, or the worm
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

# Plot each variable head_amplitude, tail_amplitude, displacement_speed, bodybends_frequency
# Aggregated per day

# Function to summarise data per each variable
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

# Plot and save plotly
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
# Summarize each replicate (for the four variables) grouped by condition and day
data_testmot <- data_mot |>
  group_by(condition, rep_id, day) |>
  summarise(
    ha = mean(head_amplitude),
    ta = mean(tail_amplitude),
    ds = mean(displacement_speed),
    bf = mean(bodybends_frequency)
  )

# This function fills empty values with 0, and can also filter for experiment if needed
from_data_testmot_fill <- function(data_testmot, id_pattern = NULL){
  if(!is.null(id_pattern)){
    smalldata <- data_testmot |> filter(str_detect(rep_id, id_pattern))
  } else {
    smalldata <- data_testmot
  }
  condition_map <- smalldata |> distinct(rep_id, condition)

  # Store original rows
  original_rows <- smalldata |> ungroup() |> select(rep_id, day)

  # Fill missing rows with 0
  smalldata <- smalldata |> ungroup() |> select(-condition) |>
    tidyr::complete(rep_id, day, fill = list(ha = 0, ta = 0, ds = 0, bf = 0)) |>
    mutate(was_missing = !paste(rep_id, day) %in% paste(original_rows$rep_id, original_rows$day))

  # Add condition to filled rows
  smalldata <- smalldata |>
    left_join(condition_map, by = "rep_id")
  
  return(smalldata)
}

# So 
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

# So these are the options:
# 1. Use all data, filled with 0 when missing
data_all <- from_data_testmot_fill(data_testmot)

# 2. Use all data as is (it will remove incomplete replicates)
data_all_nofill <- data_testmot

# 3. Use data separated per experiment, with filled missing values
exp_ids <- unique(table_mot$exp_id)

data_by_exp <- lapply(exp_ids, function(id) {
  from_data_testmot_fill(data_testmot, id)
})
names(data_by_exp) <- exp_ids

# 4. Use data separated per experiment as is (it will remove incomplete replicates)
# This option is not recommended, as if we end up with only one replicate
# the ANOVA will not work.

data_by_exp_nofill <- lapply(exp_ids, function(id) {
    data_testmot |> filter(str_detect(rep_id, id))
})
names(data_by_exp_nofill) <- paste0(exp_ids, "_nofill")


# Save statistical tests to file
variables <- c("ha", "ta", "ds", "bf")

# Choose which option to use
# 1 = all data with fill, 2 = all data no fill,
# 3 = per-experiment with fill, 4 = per-experiment no fill
option <- 2

if (option == 1) {
  dataset_names <- "all"
  dataset_list <- list(all = data_all)

} else if (option == 2) {
  dataset_names <- "all_no_fill"
  dataset_list <- list(all_no_fill = data_all_nofill)

} else if (option == 3) {
  dataset_names <- names(data_by_exp)
  dataset_list <- data_by_exp

} else if (option == 4) {
  dataset_names <- names(data_by_exp_nofill)
  dataset_list <- data_by_exp_nofill

}

for(name in dataset_names){
  df <- dataset_list[[name]]

  for(variable in variables){
  # Sink (save) results, say which 
  sink(paste0("results/motility/motility_", name, "_", variable, ".txt"))
  




  }
}




aov_result <- aov_ez(id = "rep_id",
                      dv = "ha",
                      data = df,
                      within = "day",
                      between = "condition")

data_wiz_young <- data_wiz |> filter(day <=10)
data_wiz_old <- data_wiz |> filter(day >10)

aov_young <- aov_ez(id = "rep_id",
                     dv = "ha",
                     data = data_wiz_young,
                     within = "day",
                     between = "condition")

aov_old <- aov_ez(id = "rep_id",
                     dv = "ha",
                     data = data_wiz_old,
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