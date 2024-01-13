# Load packages
# devtools::install_github("KCEvers/bifurcationEWS", force = TRUE,
# upgrade = c("default", "ask", "always", "never")[4])
# # # # #


# Setup parallelization
# cl <- parallel::makeCluster(parallel::detectCores())
# doParallel::registerDoParallel(cl)
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
<<<<<<< HEAD
                                                   nr_required_models = 25,
                                                   data_idxs = 1:(18*2),#1:(18),
=======
                                                   deSolve_method = c("lsoda", "rk4", "euler")[2],
                                                   timestep = .01,
                                                   nr_required_models = 25,
                                                   data_idxs = 1:(18*2),#1:(18),
                                                   nr_noise_iters = 10,#10,
                                                   downsample_fs = c(10, 1, .1), #10,#c(10, 5, 1),
>>>>>>> 61ac03fb56642a658567af132c0bdcd19b6a0497
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
<<<<<<< HEAD
                                                   thresh_full_band = .85,#.9,
=======
                                                   thresh_full_band = .9,
>>>>>>> 61ac03fb56642a658567af132c0bdcd19b6a0497
                                                   select_regime_switches = c("PD_1to1", "PD_1to2", "PD_2to4","PD_4to8","PD_8to16",
                                                                              "PD_Mixed-Periodic_to_Chaotic1",
                                                                              "SUBD_Chaotic_to_Mixed-Periodic1",
                                                                              "PH_16to8", "PH_8to4", "PH_4to2","PH_2to1",
                                                                              "PH_Chaotic_to_Mixed-Periodic1",
                                                                              "SUBD_Mixed-Periodic_to_Chaotic1",
                                                                              "Interior-Crisis-Separation", "Interior-Crisis-Merging",
                                                                              "Boundary-Crisis"
<<<<<<< HEAD
                                                   )
                                                   # sigmas_crit = seq(.25, 10, by = .25),

                                   ))
=======
                                                   ),
                                                   fs = 10,
                                                   # downsample_pars = list(type = c("average", "one_sample")[1],
                                                   #                        win_size = 10,
                                                   #                        which_X = c(10, "first", "middle", "last", "random")[1],
                                                   #                        seed = 123),
                                                   sigma_obs_noise = c(0.0001, .02, .04),
                                                   noise_constant = 0,
                                                   sigmas_crit = seq(.25, 10, by = .25),
                                                   sigma_crit_step = .25,
                                                   nr_consecutive_warnings = 1
                                   ))
# Nested lists that already exist in the original have to be updated separately

# # Choose stop function
# # Avoid nodes
# pars_template$stopifregime = function(out, X_names=pars_template$X_names){
# all(apply(apply(scale(out[,X_names], center = TRUE, scale = FALSE), 2, range), 2, diff) < .01)
# }
# Avoid nodes in the beginning
stopifregime_no_nodes_in_beg = function(out, X_names=pars_template$X_names){
  all(apply(apply(scale(out[,X_names], center = TRUE, scale = FALSE), 2, range), 2, diff) < .01) & (unique(out[,"bifpar_idx"]) < 50)
}
#
#
# Only nodes
stopifregime_only_nodes = function(out, X_names=pars_template$X_names){
  all(apply(apply(scale(out[,X_names], center = TRUE, scale = FALSE), 2, range), 2, diff) > .1) & (unique(out[,"bifpar_idx"])) > 2 & (unique(out[,"bifpar_idx"]) < 50)
}
# pars_template$bifpar_pars = list(bifpar_start = 1.4, bifpar_end = .9, pre_steps = 0, baseline_steps = 0, transition_steps = 501, post_steps = 0) # Steps of .0001
# pars_template$nr_timesteps =500
#
# Complete range
# pars_template$nr_timesteps =500
# pars_template$bifpar_pars = list(bifpar_start = .7, bifpar_end = 1.4, pre_steps = 0, baseline_steps = 0, transition_steps = 7001, post_steps = 0) # Step of .0001
# pars_template$bifpar_pars = list(bifpar_start = 1.4, bifpar_end = .75, pre_steps = 0, baseline_steps = 0, transition_steps = 651, post_steps = 0) # Step of .001
# pars_template$bifpar_pars = list(bifpar_start = 1.275, bifpar_end = .75, pre_steps = 0, baseline_steps = 0, transition_steps = 526, post_steps = 0) # Step of .001
# pars_template$timestep=.1
# pars_template$bifpar_pars = list(bifpar_start = .9, bifpar_end = 1.4, pre_steps = 0, baseline_steps = 0, transition_steps = 401, post_steps = 0) # Step of .001
# pars_template$factor_k = 1

>>>>>>> 61ac03fb56642a658567af132c0bdcd19b6a0497

regime_switch_names = c(
  # "Saddle-Node" # near check
  # "Hopf" # near check
  "PD_to_Chaos" # doing
  # "PH_from_Chaos" # doing
  # "Interior-Crisis-Separation"# doing 10 1 .1
  # "Interior-Crisis-Merging" # TO DO
  # "Boundary-Crisis" # check 10 1 .1
  #
  # "complete_forwards"
  # "complete_backwards"
)

<<<<<<< HEAD
# Find regime switch specific adjustments of template
pars_template_adjust = setup_bifpars()

for (regime_switch_name in regime_switch_names){

  # Adjust template to match regime switch parameters
  pars_template = utils::modifyList(pars_general_template,
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
=======
for (regime_switch_name in regime_switch_names){
  ##### SADDLE NODE
  if (regime_switch_name == "Saddle-Node"){
    pars_template = utils::modifyList(pars_general_template,
                                      list(data_idxs = 1,
                                           factor_k = 1,
                                           nr_timesteps_full = 100,
                                           nr_noise_iters = 250,
                                           bifpar_pars = list(bifpar_start = 1.15, bifpar_end = .9, pre_steps = 0, baseline_steps = 0, transition_steps = 2501, post_steps = 0),
                                           select_regime_switches = c("PD_1to1"),
                                           # max_iter = 100,
                                           stopifregime = stopifregime_only_nodes
                                      ))
    ##### HOPF
  } else if (regime_switch_name == "Hopf"){
    pars_template = utils::modifyList(pars_general_template,
                                      list(data_idxs = 1, factor_k = 1,
                                           nr_timesteps_full = 200,
                                           nr_noise_iters = 250,
                                           bifpar_pars = list(bifpar_start = .6, bifpar_end = .85, pre_steps = 0, baseline_steps = 0, transition_steps = 251, post_steps = 0),
                                           select_regime_switches = c("PD_1to2")
                                      ))
    ##### PERIOD-DOUBLING TO CHAOS
  } else if (regime_switch_name == "PD_to_Chaos"){
    pars_template = utils::modifyList(pars_general_template,
                                      list(
                                        # nr_timesteps_full = 2000,
                                        # thresh_coord_spread_trans = .01,
                                        min_length_regime_full = 8,
                                        bifpar_pars = list(bifpar_start = .85, bifpar_end = .985, pre_steps = 0, baseline_steps = 0,
                                                           transition_steps = 1351,
                                                           post_steps = 0),
                                        select_regime_switches = c(
                                          "PD_2to4","PD_4to8",
                                          "PD_8to16",
                                          "PD_Mixed-Periodic_to_Chaotic1",
                                          "SUBD_Chaotic_to_Mixed-Periodic1"
                                        )
                                      ))
    ##### PERIOD-HALVING FROM CHAOS
  } else if (regime_switch_name == "PH_from_Chaos"){
    pars_template = utils::modifyList(pars_general_template,
                                      list(
                                        thresh_coord_spread_trans = .01,
                                        factor_k = .01,
                                        min_length_regime_trans = 5,
                                        bifpar_pars = list(bifpar_start = .985, bifpar_end = .8,
                                                           pre_steps = 0, baseline_steps = 0, transition_steps = 1851,
                                                           post_steps = 0),
                                        select_regime_switches = c("PH_16to8", "PH_8to4", "PH_4to2","PH_2to1",
                                                                   "PH_Chaotic_to_Mixed-Periodic1",
                                                                   "SUBD_Mixed-Periodic_to_Chaotic1"),
                                        stopifregime = stopifregime_no_nodes_in_beg
                                      ))
    ##### EXPANSION CHAOS
  } else if (regime_switch_name == "Interior-Crisis-Merging"){

    pars_template = utils::modifyList(pars_general_template,
                                      list(
                                        nr_timesteps_full = 2000,
                                        # factor_k = 1,
                                        nr_smooth_full = 5,
                                        nr_smooth_trans = 0,
                                        min_length_regime_full = 5,
                                        min_length_regime_trans = 1,
                                        # thresh_coord_spread_full=.01,
                                        thresh_full_band = .85,
                                        # bifpar_pars = list(bifpar_start = 1.01, bifpar_end = 1.03, pre_steps = 0, baseline_steps = 0, transition_steps = 201, post_steps = 0),
                                        bifpar_pars = list(bifpar_start = 1.01, bifpar_end = 1.02, pre_steps = 0, baseline_steps = 0, transition_steps = 201, post_steps = 0),
                                        select_regime_switches = c("Interior-Crisis-Merging"),
                                        stopifregime = stopifregime_no_nodes_in_beg
                                      ))
    ##### REDUCTION CHAOS
  } else if (regime_switch_name == "Interior-Crisis-Separation"){

    pars_template = utils::modifyList(pars_general_template,
                                      list(
                                        nr_timesteps_full = 2000,
                                        factor_k = 1, nr_smooth_full = 10, nr_smooth_trans = 10,
                                        thresh_coord_spread_full=.01,
                                        thresh_full_band = .85,
                                        bifpar_pars = list(bifpar_start = 1.03, bifpar_end = 1.01, pre_steps = 0, baseline_steps = 0, transition_steps = 201, post_steps = 0),
                                        select_regime_switches = c("Interior-Crisis-Separation"),
                                        stopifregime = stopifregime_no_nodes_in_beg
                                      ))

    ##### BOUNDARY CRISIS
  } else if (regime_switch_name == "Boundary-Crisis"){
    pars_template = utils::modifyList(pars_general_template,
                                      list(factor_k = 1,
                                           deSolve_method = "euler",
                                           timestep = .1,
                                           nr_timesteps_full = 10000,
                                           keep_nr_timesteps = 1000,
                                           bifpar_pars = list(bifpar_start = 1.1, bifpar_end = 1.3, pre_steps = 0,baseline_steps = 0, transition_steps = 201, post_steps = 0),
                                           select_regime_switches = c("Boundary-Crisis"),
                                           stopifregime = stopifregime_no_nodes_in_beg,
                                           min_length_regime_full = 20,
                                           min_length_regime_trans = 20,
                                           variable_name = "X4",
                                           thresh_coord_spread_full = .1

                                      ))
  } else if (regime_switch_name == "complete_forwards"){
    pars_template = utils::modifyList(pars_general_template,
                                      list(factor_k = .1,
                                           data_idxs = 1,
                                           # deSolve_method = "euler",
                                           # timestep = .1,
                                           nr_timesteps_full = 2000,
                                           # keep_nr_timesteps = 1000,
                                          bifpar_pars = list(bifpar_start = .6, bifpar_end = 1.3, pre_steps = 0,baseline_steps = 0, transition_steps = 1401, post_steps = 0),
                                           # select_regime_switches = c("Boundary-Crisis"),
                                           min_length_regime_full = 20,
                                           min_length_regime_trans = 20,
                                           # variable_name = "X4",
                                           thresh_coord_spread_full = .1

                                      ))
  } else if (regime_switch_name == "complete_backwards"){

    pars_template = utils::modifyList(pars_general_template,
                                      list(factor_k = .1,
                                           data_idxs=1:18,
                                           # deSolve_method = "euler",
                                           # timestep = .1,
                                           nr_timesteps_full = 2000,
                                           # keep_nr_timesteps = 1000,
                                           bifpar_pars = list(bifpar_start = 1.25, bifpar_end = .6, pre_steps = 0,baseline_steps = 0, transition_steps = 1301, post_steps = 0),
                                           # select_regime_switches = c("Boundary-Crisis"),
                                           # stopifregime = stopifregime_no_nodes_in_beg,
                                           min_length_regime_full = 10,
                                           # min_length_regime_trans = 10,
                                           # variable_name = "X4",
                                           thresh_coord_spread_full = .1
                                      ))
  }

  print(c(pars_template$bifpar_pars$bifpar_start, pars_template$bifpar_pars$bifpar_end))

  # File parameters
  pars_template$type_output = 'data'
  pars_template$bifpar_list = do.call(get_bifurcation_range, pars_template$bifpar_pars)
  start_s = pars_template$bifpar_list[[1]]$s
  second_s = pars_template$bifpar_list[[2]]$s
  length_s = length(pars_template$bifpar_list)
  pars_template$analysis_type = sprintf(
    "s%.05f-s%.05f-by%.05f",
    start_s, pars_template$bifpar_list[[length_s]]$s, diff(c(start_s, second_s)))

>>>>>>> 61ac03fb56642a658567af132c0bdcd19b6a0497
  print(pars_template$nr_timesteps_full)
  print(pars_template$downsample_fs)
  print(pars_template$sigma_obs_noise)

<<<<<<< HEAD
  # Run scripts
  # source('generate_full_GLV.R')
  #
=======

  # source('generate_full_GLV.R')
>>>>>>> 61ac03fb56642a658567af132c0bdcd19b6a0497
  # source('generate_transitions_GLV.R')
  # source('compute_EWS_GLV.R')
  source('eval_performance_EWS_GLV.R')

}
parallel::stopCluster(cl)
