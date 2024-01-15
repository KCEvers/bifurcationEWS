## Code to prepare GLV datasets for package by cutting timeseries length
library(dplyr)
library(foreach)

# Define parameters
filepath_stem = "C:/Users/kever/Documents/PhD_demos/scripts/firstpkg/detGLV"
pre_steps = 100
baseline_steps = 100
transition_steps = 100

# Helpers
trans_model = sprintf("transition_%dtransSteps", transition_steps)
null_model = "null"
deSolve_method = "rk4"
timestep = 0.01
nr_timesteps = 1000

# Find regime switch specific adjustments of template
pars_template_adjust = setup_bifpars() %>%
  # Remove bifpar_list
  purrr::map(., function(x){x[names(x) != "bifpar_list"]}) %>%
  # Set parameters with substring "_trans" as standard parameters by removing substring
  purrr::map(., function(x){modify_list(x, x[grepl("_trans", names(x))] %>% stats::setNames(stringr::str_replace(names(.), "_trans", "")))}) %>%
  # Unwrap template: Repeat when there are multiple selected regime switches for s
  purrr::map(., function(x){rep(list(x), length(x$select_regime_switches)) %>%
      purrr::imap(., function(y, idx){modify_list(y,
                                                        list(regime_switch_name = y$select_regime_switches[idx],
                                                             data_idx = 1, # Add data_idx
                                                             GLV_string = y$GLV_string[idx],
                                                             select_regime_switches = y$select_regime_switches[idx])
                                                        )
        } )
    }) %>% unlist(recursive = F, use.names = F) %>%
  stats::setNames(unlist(purrr::map(., "regime_switch_name")))

# Update data idx for some regime switches
pars_template_adjust[["SUBD_Mixed-Periodic_to_Chaotic1"]]$data_idx = 19
pars_template_adjust[["Interior-Crisis-Merging"]]$data_idx = 3
pars_template_adjust[["Interior-Crisis-Separation"]]$data_idx = 3
pars_template_adjust[["Boundary-Crisis"]]$data_idx = 3


# for_par = pars_template_adjust[["Boundary-Crisis"]]
# for_par = pars_template_adjust[["Interior-Crisis-Merging"]]
for_par = pars_template_adjust[["Interior-Crisis-Separation"]]
for_par = pars_template_adjust[["Interior-Crisis-Separation"]]

for_par=pars_template_adjust[[5]]
s_string=for_par[["s_string"]]
select_regime_switches=for_par[["select_regime_switches"]]
min_length_regime = for_par[["min_length_regime"]]
for_par[["thresh_coord_spread"]]
data_idx=for_par$data_idx

# Loop through desired regime switches
filepaths_data_cut = foreach(for_par = pars_template_adjust) %do% {

  with(for_par, {
      # Read transition and null model
    filepath_trans = file.path(filepath_stem, "data", s_string, select_regime_switches,
              sprintf("nr%d_%s_T%d_ts%.03f_fs10.00_Xsigma0.00000_%s.RDS", data_idx, deSolve_method, nr_timesteps, timestep, trans_model))
    filepath_null = file.path(filepath_stem, "data", s_string, select_regime_switches,
                              sprintf("nr%d_%s_T%d_ts%.03f_fs10.00_Xsigma0.00000_%s.RDS", data_idx, deSolve_method, nr_timesteps, timestep, null_model))
    filepath_trans_cut = stringr::str_replace(filepath_trans, ".RDS", "_cut.RDS")
    filepath_null_cut = stringr::str_replace(filepath_null, ".RDS", "_cut.RDS")

    if (!file.exists(filepath_trans_cut) | !file.exists(filepath_null_cut)){
      print(filepath_trans)
      print(filepath_null)
      GLV_trans <- readRDS(filepath_trans)
      GLV_null <- readRDS(filepath_null)
      format(object.size(GLV_trans), "Mb")
      format(object.size(GLV_null), "Mb")

      # Select names in for_par that match arguments required in find_regimes
      for_par_find_regimes = for_par[intersect(names(for_par), names(formals(find_regimes)))]
      # Find all regimes in transition and null model
      regime_list_trans = do.call(find_regimes,
              modify_list(for_par_find_regimes, list(GLV=GLV_trans)))
      regime_list_null = do.call(find_regimes,
                                  modify_list(for_par_find_regimes, list(GLV=GLV_null)))

      # Define desired regime filter
      regime_switch_list = make_filter_regime_switches(min_length_regime,
                                                       baseline_steps, transition_steps)[select_regime_switches]

      # Match transition and null model
      regime_bounds_successful = match_trans_null_model(
        regime_bounds_trans = regime_list_trans$regime_bounds %>% mutate(trans_or_null = "transition", data_idx = data_idx),
        regime_bounds_null = regime_list_null$regime_bounds %>% mutate(trans_or_null = "null", data_idx = data_idx),
        min_length_regime = min_length_regime, regime_switch_list = regime_switch_list,
        pre_steps = pre_steps, baseline_steps = baseline_steps, transition_steps = transition_steps)

      # Cut length of data set, resulting in 100 baseline steps and 100 transition steps
      seq_bifpar = seq(unique(regime_bounds_successful$baseline_start_idx),
                       max(regime_bounds_successful$transition_end_idx))
      GLV_trans$df = GLV_trans$df  %>%
        filter(bifpar_idx %in% seq_bifpar) %>%
        select(-c(.data$time, .data$win_nr))
      GLV_null$df = GLV_null$df  %>%
        filter(bifpar_idx %in% seq_bifpar) %>%
        select(-c(.data$time, .data$win_nr))

      format(object.size(GLV_trans), "Mb")
      format(object.size(GLV_null), "Mb")

      # Save
      saveRDS(GLV_trans, filepath_trans_cut)
      saveRDS(GLV_null, filepath_null_cut)

    }
    return(list("trans" = filepath_trans_cut, "null" = filepath_null_cut))
  })
}

# j="trans"

# Finally, load cut data sets as R objects
# pars_template_adjust %>%
  # purrr::imap(., function(x, idx){
for (idx in seq_along(pars_template_adjust)){
   for (j in c("trans", "null")){
      GLV_name = sprintf("%s_%s", pars_template_adjust[[idx]][["GLV_string"]], j)
      print(GLV_name)
      print(idx)

    # Create R object
    assign(GLV_name, readRDS(filepaths_data_cut[[idx]][[j]]))

    # # Add to raw data
    # usethis::use_data_raw(GLV_name)
    # # Add: ## Run prepGLV.R to prepare data set
    # # Run script
    #
    # # Create documentation
    # file.create(sprintf("R/%s.R", GLV_name))
}
  }


