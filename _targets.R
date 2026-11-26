# Created by use_targets().

# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)

# Set target options:
tar_option_set(
  packages = c("tibble", "broom.mixed", "purrr", "patchwork", "performance", "tidyverse", "fixr", "glmmTMB", "ggplot2") # Packages that your targets need for their tasks.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Power analysis for goby~cormorant + year + (1|year) + (1|site):
list(
  tar_target(goby_data, get_goby_data("Data/CPUE_BaitedOnly.csv")), #get data
  tar_target(cormorant_data, get_cormorant_data("Data/CormorantHWRaw.csv")),#get data
  tar_target(comb_data, combine_data(goby_data, cormorant_data)),#put in one dataset
  tar_target(sites, get_sites(comb_data)),#extract sites as levels
  tar_target(nest_data, fix_nest_data(comb_data)),#separate cormorant data set 2002-2022 for easier simulation
  tar_target(cormorant_residual_sd, model_cormorant_trend(nest_data)), #modelling cormorant data to simulate future years
  tar_target(designs, create_four_designs(4000, cormorant_residual_sd, sites, nest_data)), #simulate future cormorant years with a estimated average of 4000 and normally distributed deviations
  tar_target(power_analysis, run_power_analysis(designs, nsim=1000)),#Run power analysis can change nsim (should be 1000)
  tar_target(power_plots, plot_power_analysis(power_analysis[[2]]))#power analysis plots 
)
