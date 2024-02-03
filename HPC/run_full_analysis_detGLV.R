# Load packages
# devtools::install_github("KCEvers/bifurcationEWS", force = TRUE,
# upgrade = c("default", "ask", "always", "never")[4])
# # # # #


# Setup parallelization
cl <- parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
# cl <- parallel::makeCluster(1)
# doParallel::registerDoParallel(cl)

library(bifurcationEWS)
library(foreach)
library(dplyr)
library(ggplot2)

# Define template parameters
pars_general_template = setup_pars(model_name = "detGLV",
                                   pars_add = list(nr_timesteps_full = 1000,
                                                   nr_timesteps_trans = 1000,
                                                   min_length_regime_full = 10,
                                                   min_length_regime_trans = 10,
                                                   nr_required_models = 25,
                                                   data_idxs = 1:(18*2),#1:(18),
                                                   pre_steps = 100,
                                                   baseline_steps = 100,#c(25,100),#c(25,50,100),
                                                   default_transition_steps = 100,
                                                   default_baseline_steps = 100,
                                                   transition_steps = 100,#c(25,100), #c(25,50,100),
                                                   post_steps = 100,
                                                   nr_smooth_full = 0,
                                                   nr_smooth_trans = 0,
                                                   factor_k = .1,
                                                   variable_name = "X1",
                                                   keep_nr_timesteps = "all",
                                                   thresh_node = .001,
                                                   thresh_coord_spread_full = .025,
                                                   thresh_coord_spread_trans = .025,
                                                   thresh_peak_idx_spread=2,
                                                   thresh_full_band_full = .9,
                                                   thresh_full_band_trans = .7,
                                                   select_regime_switches = c("PD_1to1", "PD_1to2", "PD_2to4","PD_4to8","PD_8to16",
                                                                              "PD_Mixed-Periodic_to_Chaotic1",
                                                                              "SUBD_Chaotic_to_Mixed-Periodic1",
                                                                              "PH_16to8", "PH_8to4", "PH_4to2","PH_2to1",
                                                                              "PH_Chaotic_to_Mixed-Periodic1",
                                                                              "SUBD_Mixed-Periodic_to_Chaotic1",
                                                                              "Interior-Crisis-Separation", "Interior-Crisis-Merging",
                                                                              "Boundary-Crisis"
                                                   ),
                                                   sigma_crit_step = .01,

                                                   thresh_max_sigma_crit = 150
                                                   # sigmas_crit = seq(.25, 10, by = .25),

                                   ))

regime_switch_names = c(
  # "Saddle-Node" # near check
  # "Hopf" # near check
  # "PD_to_Chaos" # doing
  # "PH_from_Chaos" # doing
  "Interior-Crisis-Separation"# doing 10 1 .1
  # "Interior-Crisis-Merging" # TO DO
  # "Boundary-Crisis" # check 10 1 .1
  #
  # "complete_forwards"
  # "complete_backwards"
)
# regime_switch_name = "Interior-Crisis-Separation"
# regime_switch_name = "Boundary-Crisis"
# regime_switch_name = "Interior-Crisis-Merging"
# regime_switch_name =  "PD_to_Chaos"
# regime_switch_name =  "PH_from_Chaos"
# regime_switch_name =  "Hopf"
# regime_switch_name =  "Saddle-Node"

# Find regime switch specific adjustments of template
pars_template_adjust = setup_bifpars()


for (regime_switch_name in regime_switch_names){

  # Adjust template to match regime switch parameters
  pars_template = modify_list(pars_general_template,
                                pars_template_adjust[[regime_switch_name]])

  # File parameters
  pars_template$type_output = 'data'
  pars_template$analysis_type = pars_template$s_string
  # pars_template$bifpar_list = do.call(get_bifurcation_range, pars_template$bifpar_pars)
  # start_s = pars_template$bifpar_list[[1]]$s
  # second_s = pars_template$bifpar_list[[2]]$s
  # length_s = length(pars_template$bifpar_list)
  # pars_template$analysis_type = sprintf(
  #   "s%.05f-s%.05f-by%.05f",
  #   start_s, pars_template$bifpar_list[[length_s]]$s, diff(c(start_s, second_s)))

  print(c(pars_template$bifpar_pars$bifpar_start, pars_template$bifpar_pars$bifpar_end))
  print(pars_template$nr_timesteps_full)
  print(pars_template$downsample_fs)
  print(pars_template$sigma_obs_noise)

  # Run scripts
  # source('generate_full_GLV.R')
  # source('generate_transitions_GLV.R')
  # source('compute_EWS_GLV.R')
  source('eval_performance_EWS_GLV.R')

}
parallel::stopCluster(cl)
