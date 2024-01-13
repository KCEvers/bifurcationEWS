## code to prepare `GLV_Hopf_trans` dataset goes here

## Code to prepare GLV datasets for package by cutting timeseries length

# # Read full data set
# filepath_stem = "C:/Users/kever/Documents/PhD_demos/scripts/firstpkg/detGLV"
# s_IC_forward = "s1.01000-s1.02000-by0.00005" # "s1.01000-s1.03000-by0.00010"
# s_IC_backward = "s1.03000-s1.01000-by-0.00010"
# s_BC = "s1.10000-s1.30000-by0.00100"
# s_SN = "s1.15000-s0.90000-by-0.00010"
# s_Hopf = "s0.60000-s0.85000-by0.00100"
# s_PD = "s0.85000-s0.98500-by0.00010"
# s_PH = "s0.98500-s0.80000-by-0.00010"
# s_forward = "s0.60000-s1.30000-by0.00050"
# s_backward = "s1.25000-s0.60000-by-0.00050"
# desired_regime_switch = "PD_1to2"
# pre_steps = 100
# baseline_steps = c(100)
# transition_steps = 100
# data_idx = 1
# trans_model = sprintf("transition_%dtransSteps", transition_steps)
# null_model = "null"
# factor_k = 1
# min_length_regime = 5
#
# GLV_Hopf_trans <- readRDS(file.path(filepath_stem, "data", s_Hopf, "PD_1to2", sprintf("nr%d_rk4_T1000_ts0.010_fs10.00_Xsigma0.00000_%s.RDS", data_idx, trans_model)))
# GLV_Hopf_null <- readRDS(file.path(filepath_stem, "data", s_Hopf, "PD_1to2", sprintf("nr%d_rk4_T1000_ts0.010_fs10.00_Xsigma0.00000_%s.RDS", data_idx, null_model)))
# format(object.size(GLV_Hopf_trans), "Mb")
# format(object.size(GLV_Hopf_null), "Mb")
#
# # Find regimes
# regime_list_trans = find_regimes(GLV_Hopf_trans, factor_k = factor_k)
# regime_list_null = find_regimes(GLV_Hopf_null, factor_k = factor_k)
#
# # Define desired regimes
# regime_switch_list = make_filter_regime_switches(min_length_regime, default_baseline_steps, transition_steps)[desired_regime_switch]
#
# # Match transition and null model
# regime_bounds_successful = match_trans_null_model(
#     regime_bounds_trans = regime_list_trans$regime_bounds %>% mutate(trans_or_null = "transition", data_idx = !!data_idx),
#     regime_bounds_null = regime_list_null$regime_bounds %>% mutate(trans_or_null = "null", data_idx = !!data_idx),
#     min_length_regime = min_length_regime, regime_switch_list = regime_switch_list,
#     pre_steps = pre_steps, baseline_steps = baseline_steps, transition_steps = transition_steps)
# regime_bounds_successful %>% as.data.frame()
#
# # Cut length of data set, resulting in 100 baseline steps and 100 transition steps
# GLV_Hopf_trans$df = GLV_Hopf_trans$df  %>%
#   filter(bifpar_idx %in% seq(unique(regime_bounds_successful$baseline_start_idx), unique(regime_bounds_successful$transition_end_idx))) %>%
#   select(-c(.data$time, .data$win_nr))
# GLV_Hopf_null$df = GLV_Hopf_null$df  %>%
#   filter(bifpar_idx %in% seq(unique(regime_bounds_successful$baseline_start_idx), unique(regime_bounds_successful$transition_end_idx))) %>%
#   select(-c(.data$time, .data$win_nr))
# format(object.size(GLV_Hopf_trans), "Mb")
# format(object.size(GLV_Hopf_null), "Mb")

# GLV_Hopf_trans = readRDS("C:/Users/kever/Documents/PhD_demos/scripts/firstpkg/detGLV/data/s0.60000-s0.85000-by0.00100/PD_1to2/nr1_rk4_T1000_ts0.010_fs10.00_Xsigma0.00000_transition_100transSteps_cut.RDS")
usethis::use_data(GLV_Hopf_trans, overwrite = TRUE)


