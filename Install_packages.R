# Install required packages

#--------------------------
# Author: Núria Garriga Alonso
# GitHub: @nuriagaralon
# Repository: https://github.com/nuriagaralon/L4_analysis_R
#--------------------------

# List of required packages for all scripts
packages <- c(
  "tidyverse",
  "readxl",
  "survival",
  "survminer",
  "plotly",
  "ggfortify",
  "pracma",
  "DescTools",
  "rstatix",
  "ggpubr",
  "ggsci",
  "htmlwidgets",
  "zoo",
  "afex",
  "emmeans"
)

# Function: install only the ones that are not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Apply the function to each package
invisible(lapply(packages, install_if_missing))

# Remove variable
rm(packages)