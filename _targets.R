# Setup
library(targets)
library(tarchetypes)
source("R/functions.R")
source("R/figures.R")
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(tidyverse.quiet = TRUE)
# Run targets
tar_option_set(packages = c(
  "tidyverse", "readxl", "rlang", "png", "grid", "viridis", "viridisLite",
  "scales", "ggh4x", "patchwork", "rstan", "brms", "DHARMa", "tidybayes",
  "ggdist", "insight", "knitr", "quarto"
))
list(
  tar_target(data_file, "data/data.csv", format = "file"),
  tar_target(raw_data, read.csv(data_file, check.names = FALSE)),
  tar_target(dens_data, filter_densities_data(raw_data)),
  tar_target(mod_data, filter_model_data(raw_data)),
  tar_target(model, fit_metreg(mod_data)),
  tar_target(fig_2, make_fig_2(mod_data)),
  tar_target(sim_data, simulate_densities(1e3, data = dens_data)),
  tar_target(dens_plot_data, make_densities_plot_data(sim_data)),
  tar_target(fig_3, make_fig_3(model)),
  tar_target(fig_4, make_fig_4(model, mod_data)),
  tar_target(fig_5, make_fig_5(dens_plot_data)),
  tar_target(fig_s1, make_fig_s1(model)),
  tar_target(fig_s2, make_fig_s2(model)),
  tar_target(null_model, fit_null_metreg(mod_data)),
  tar_target(mod_comp, brms::loo_compare(null_model, model))
)
