# L4_analysis_R: An R project for SydLab-One phenotypic data

## Overview
R project for *Caenorhabditis elegans* phenotypic analysis of SydLab One microfluidics system data. It accepts `.xlsx` files and generates analyses for growth, survival, fertility, motility and fluorescence data.

The following is a tutorial on the usage of this R project.
- If you need to install it, start from **Installation**.
- If you are analyzing data, start from **Growth analysis**, **Survival analysis**, etc.
- Remember to add the data as detailed in **Data files**.

## Installation
### 1. Install R and RStudio
Install R from [CRAN](https://cran.rstudio.com/). Choose the installer for your operating system  (Windows, Linux, or macOS). The analyses were tested with version 4.5.2, but the latest release will likely also work.

Install RStudio from [posit](https://posit.co/download/rstudio-desktop) for Windows, Linux, or macOS.  

### 2. Clone the repository
In RStudio:

1. Click **File -> New project**
2. Select **Version Control -> Git**
3. Paste the repository URL:

```
https://github.com/nuriagaralon/L4_analysis_R
```
4. Click **Create Project**



## Data files
All analysis scripts look for the relevant data files in the `L4_analysis_R/data` folder. Copy the following files there:

- Growth data: `[experimentID]_growth_filtered_raw.xlsx`
- Survival data: `[experimentID]_survival_condition.xlsx`
- Fertility data: : `[experimentID]_eggs_condition.xlsx` or `[experimentID]_eggs_channel.xlsx`
- Motility data: `[experimentID]_motility_filtered_raw.xlsx`
- Fluorescence data: `[experimentID]_fluo_raw_[light].xlsx`

## Growth analysis

## Survival analysis

