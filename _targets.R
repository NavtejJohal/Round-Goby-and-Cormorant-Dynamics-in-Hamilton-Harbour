# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tibble", "broom.mixed", "purrr", "patchwork", "performance", "tidyverse", "fixr", "glmmTMB", "ggplot2") # Packages that your targets need for their tasks.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(goby_data, get_goby_data("Data/CPUE_BaitedOnly.csv")),
  tar_target(cormorant_data, get_cormorant_data("Data/CormorantHWRaw.csv")),
  tar_target(comb_data, combine_data(goby_data, cormorant_data)),
  tar_target(sites, get_sites(comb_data)),
  tar_target(nest_data, fix_nest_data(comb_data)),
  tar_target(cormorant_residual_sd, model_cormorant_trend(nest_data)),
  tar_target(designs, create_four_designs(4000, cormorant_residual_sd, sites, nest_data)),
  tar_target(power_analysis, run_power_analysis(designs)),
  tar_target(power_plots, plot_power_analysis(power_analysis[[2]]))
)
