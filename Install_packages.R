# Reproducible package environment managed by renv

#--------------------------
# Author: Núria Garriga Alonso
# GitHub: @nuriagaralon
# Repository: https://github.com/nuriagaralon/L4_analysis_R
#--------------------------

# If the project does not automatically prompt for restore:
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::restore()