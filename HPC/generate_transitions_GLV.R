print("Start simulating selected regime switches!")
rerun = F
pars_template$nr_timesteps=pars_template$nr_timesteps_full
pars_template$min_length_regime=pars_template$min_length_regime_full

# Select regime switches
filepath_all_regime_bounds = format_path(format_pars(utils::modifyList(
  pars_template,
  list(type_output = "regimes", filename = "all_regime_bounds")
)))
regime_bounds_df = readRDS(filepath_all_regime_bounds)


# Find selected regimes in full GLVs
regime_switch_list = make_filter_regime_switches(pars_template$min_length_regime, pars_template$default_baseline_steps, pars_template$default_transition_steps)
selected_regime_bounds = apply_filter_regime_switches(regime_bounds_df, regime_switch_list[purrr::map(regime_switch_list, "regime_switch") %in% pars_template$select_regime_switches])

print("SELECTED")
print(selected_regime_bounds %>% group_by(data_idx, regime_switch) %>% dplyr::summarise(n=n(), .groups = 'drop') %>% as.data.frame())

# Switch to transition models
pars_template$nr_timesteps=pars_template$nr_timesteps_trans
pars_template$min_length_regime=pars_template$min_length_regime_trans
filepath_successful_regime_bounds = format_path(format_pars(utils::modifyList(
  pars_template,
  list(type_output = "regimes", filename = "all_trans_regime_bounds")
)))
filepath_unfiltered_trans_regime_bounds = format_path(format_pars(utils::modifyList(
  pars_template,
  list(type_output = "regimes", filename = "unfiltered_trans_regime_bounds")
)))
regime_switch_list = make_filter_regime_switches(pars_template$min_length_regime, pars_template$default_baseline_steps, pars_template$default_transition_steps)


forloop_ = c(
  get_forloop(
    data_idx = pars_template$data_idxs,
    regime_switch = pars_template$select_regime_switches,
    transition_steps = pars_template$transition_steps,
    trans_or_null = c("transition", "null")
  )
  # Hopf bifurcation only needs to be simulated once (as it starts with a node)
) %>% purrr::map(function(x){
  if (x$regime_switch == "PD_1to2" & x$data_idx != 1) {
    return(NULL)
  } else if (x$trans_or_null == "null" & x$transition_steps != pars_template$default_transition_steps) {
    return(NULL)
    } else {
    return(x)
    }
  })
forloop = forloop_[!unlist(purrr::map(forloop_, is.null))]

#
# filepaths = foreach(
#   for_par = forloop,
#   .packages = c("bifurcationEWS", "dplyr", "ggplot2"),
#   .export = c("pars_template")
# ) %dopar% {
#   .GlobalEnv$pars_template <- pars_template # Don't ask me why this is necessary...
#
#   pars <- utils::modifyList(pars_template, for_par)
#   pars$subfolder1 = for_par$regime_switch
#
#   # Get details of regime switch
#   desired_switch = selected_regime_bounds %>% dplyr::filter(data_idx == pars$data_idx,
#                                                             regime_switch == pars$regime_switch)
#   if (nrow(desired_switch) > 0) {
#     if (pars$trans_or_null == "transition") {
#       bifpar_end = pars$bifpar_list[[desired_switch$regime2_start_idx]]$s
#       pars$filename = sprintf("%s_%dtransSteps",
#                               pars$trans_or_null,
#                               pars$transition_steps)
#       transition_steps = pars$transition_steps
#     } else if (pars$trans_or_null == "null") {
#       bifpar_end = NA
#       pars$filename = for_par$trans_or_null
#       transition_steps = pars$default_transition_steps
#     }
#
#     filepath_GLV = format_path(format_pars(pars))
#     filepath_regimes = format_path(format_pars(utils::modifyList(pars, list(type_output = "regimes"))))
#
#     print(filepath_GLV)
#     if (file.exists(filepath_GLV)){
#           err1 = "try-error" %in% class(try(readRDS(filepath_GLV)))
#     } else {
#       err1 = FALSE
#     }
#     if (file.exists(filepath_regimes)){
#       err2 = "try-error" %in% class(try(readRDS(filepath_regimes)))
#     } else {
#       err2 = FALSE
#     }
#     print(err1)
#     print(err2)
#     if (err1){
#       file.remove(filepath_GLV)
#     }
#     if (err2){
#       file.remove(filepath_regimes)
#     }
#     # return(file.exists(filepath_GLV))
#   return(c(err1,err2))
#   }
# }
#

start_t = Sys.time()
foreach(
  for_par = forloop,
  .combine = 'cfun',
  .packages = c("bifurcationEWS", "dplyr", "ggplot2"),
  .export = c("pars_template")
) %dopar% {
  .GlobalEnv$pars_template <- pars_template # Don't ask me why this is necessary...

  pars <- utils::modifyList(pars_template, for_par)
  pars$subfolder1 = for_par$regime_switch

  # Get details of regime switch
  desired_switch = selected_regime_bounds %>% dplyr::filter(data_idx == pars$data_idx,
                                                            regime_switch == pars$regime_switch)
  if (nrow(desired_switch) > 0) {
    if (pars$trans_or_null == "transition") {
      bifpar_end = pars$bifpar_list[[desired_switch$regime2_start_idx]]$s
      pars$filename = sprintf("%s_%dtransSteps",
                              pars$trans_or_null,
                              pars$transition_steps)
      transition_steps = pars$transition_steps
    } else if (pars$trans_or_null == "null") {
      bifpar_end = NA
      pars$filename = for_par$trans_or_null
      transition_steps = pars$default_transition_steps
    }
    s_seq = get_bifurcation_range(
      bifpar_start = pars$bifpar_list[[desired_switch$regime1_halfway_idx]]$s,
      bifpar_end = bifpar_end,
      pre_steps = pars$pre_steps,
      baseline_steps = pars$default_baseline_steps,
      transition_steps = transition_steps,
      post_steps = pars$post_steps
    )
    X0 = unlist(desired_switch[, pars$X_names])

    filepath_GLV = format_path(format_pars(pars))
    filepath_regimes = format_path(format_pars(utils::modifyList(pars, list(type_output = "regimes"))))

    print(filepath_GLV)
    print(as.data.frame(desired_switch))
    if (!file.exists(filepath_regimes) | rerun) {
      if (!file.exists(filepath_GLV)) {
        # Generate timeseries
        GLV = bifurcation_ts(
          model = GLV_model,
          model_pars = pars$model_pars,
          bifpar_list = s_seq,
          timestep = pars$timestep,
          X0 = X0,
          nr_timesteps = pars$nr_timesteps_trans,
          fs = pars$fs,
          deSolve_method = pars$deSolve_method,
          X_names = pars$X_names
        )
        # df = GLV$df
        # head(df)

        saveRDS(GLV, filepath_GLV)
      }
      GLV = readRDS(filepath_GLV)
      # GLV$X0s = GLV$X0s %>% magrittr::set_colnames(c("bifpar_idx", pars$X_names))
      # saveRDS(GLV, filepath_GLV)

      # Find regimes
      regime_list = find_regimes(
        GLV,
        thresh_node = pars$thresh_node,
        thresh_coord_spread = pars$thresh_coord_spread_trans,
        thresh_peak_idx_spread = pars$thresh_peak_idx_spread,
        thresh_full_band = pars$thresh_full_band,
        min_length_regime = pars$min_length_regime,
        nr_smooth = pars$nr_smooth_trans,
        factor_k = pars$factor_k,
        variable_name = pars$variable_name,
        keep_nr_timesteps = pars$keep_nr_timesteps
      )

      #   print(regime_list$periods)
        # print(regime_list$regimes)
        # print(regime_list$regime_bounds)

      saveRDS(regime_list, filepath_regimes)
    }

    # if (TRUE & pars$data_idx == 1 & file.exists(filepath_regimes)){
    if (F & file.exists(filepath_regimes)){
    regime_list = readRDS(filepath_regimes)

    # Plot regime results
    peaks_df = regime_list$peaks_df %>%
      dplyr::filter(bifpar_idx > 1)# %>%
    # group_by(bifpar_idx, variable) %>%
    # arrange(.data$time_idx, .by_group=TRUE) %>%
    # Only transient
    # dplyr::filter(dplyr::row_number() <= round(n()*.2) ) %>% ungroup()
    # Discard transient
    # dplyr::filter(dplyr::row_number() >= round(n()*.8) ) %>% ungroup()

    point_size = 10 / length(unique(peaks_df$bifpar_idx))

    pl_peaks = peaks_df %>% dplyr::filter(variable=="X1") %>% ggplot() + geom_point(aes(x = bifpar_idx, y = X), size = point_size) +
      labs(x = 'Bifurcation parameter', y = "X1", title = "Transient bifurcation diagram")

    filepath_image = format_path(format_pars(utils::modifyList(
      pars,
      list(
        type_output = "figs",
        subfolder1 = "regimes",
        subfolder2 = pars$subfolder1,
        filename = paste0("peaks-x1", "_", pars$filename),
        # filename = paste0("only-transient_peaks-x1", "_", pars$filename),
        # filename = paste0("discard-transient_peaks-x1", "_", pars$filename),
        file_ext = ".png"
      )
    )))
    save_plot(
      style_plot(pl_peaks),
      filepath_image,
      w = 30,
      h = 20,
      formats = ".png"
    )

    pl_regimes = peaks_df %>% ggplot() +
      geom_rect(
        data = regime_list$regimes,
        aes(
          xmin = start_bifpar_idx - .5,
          xmax = end_bifpar_idx + .5,
          ymin = -Inf,
          ymax = Inf,
          fill = regime
        ),
        linewidth = .5,
        alpha = .4
      ) +
      geom_vline(
        data = regime_list$regime_bounds,
        aes(xintercept = regime1_end_idx,
            col = regime_switch_type),
        linewidth = 1.3,
        linetype = 'dashed'
      ) +
      geom_vline(
        data = regime_list$regime_bounds,
        aes(xintercept = regime2_start_idx,
            col = regime_switch_type),
        linewidth = 1.3,
        linetype = 'dashed'
      ) +
      geom_point(
        aes(x = bifpar_idx, y = X),
        col = 'grey30',
        alpha = .8,
        size = point_size
      ) +

      viridis::scale_color_viridis(
        discrete = TRUE,
        option = 'magma',
        begin = 0,
        end = .8,
        name = "Regime switch (end of first and start of second regime)"
      ) +
      viridis::scale_fill_viridis(
        discrete = TRUE,
        option = 'turbo',
        begin = 0,
        end = 1,
        name = "Period",
        guide = guide_legend(override.aes = list(alpha = 1,
                                                 size = 3))
      ) +
      guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) + # x rows in legend
      ggh4x::facet_grid2(variable ~ .) +
      labs(x = 'Bifurcation parameter', y = "", title = "Transient bifurcation diagram")

    save_plot(
      style_plot(pl_regimes),
      stringr::str_replace(filepath_image, "peaks-x1", "peaks"),
      w = 40,
      h = 50,
      formats = ".png"
    )
    print(filepath_image)
}

       return(NULL)
  }
  return(NULL)

}

end_t = Sys.time()
print(end_t - start_t)

# Compile regimes of all simulations
regimes_list = foreach(for_par = forloop) %dopar% {
                       pars <- utils::modifyList(pars_template, for_par)
                       pars$subfolder1 = for_par$regime_switch

                       # Get details of regime switch
                       desired_switch = selected_regime_bounds %>% dplyr::filter(data_idx == pars$data_idx,
                                                                                 regime_switch == pars$regime_switch)
                       if (nrow(desired_switch) > 0){
                         if (pars$trans_or_null == "transition") {
                           pars$filename = sprintf("%s_%dtransSteps",
                                                   pars$trans_or_null,
                                                   pars$transition_steps)
                         } else if (pars$trans_or_null == "null") {
                           pars$filename = for_par$trans_or_null
                         }
                         filepath_GLV = format_path(format_pars(pars))
                         filepath_regimes = format_path(format_pars(utils::modifyList(pars, list(type_output = "regimes"))))

                         if (file.exists(filepath_regimes)) {
                           regime_list = readRDS(filepath_regimes)[c("regime_bounds")]
                           return(regime_list)
                         } else {
                           return(NULL)
                         }
}
                       }
regimes_list_idx = which(!unlist(purrr::map(regimes_list, is.null)))


regime_bounds_trans_df = plyr::llply(regimes_list_idx,
                                    function(x) {
                                        regimes_list[[x]]$regime_bounds %>% add_par_as_cols(forloop[[x]])
                                    }) %>%
  do.call(rbind, .) %>% as.data.frame() %>%
  tidyr::unnest(cols = colnames(.), keep_empty = TRUE) %>%
  utils::type.convert(as.is = TRUE)
saveRDS(regime_bounds_trans_df, filepath_unfiltered_trans_regime_bounds)
regime_bounds_trans_df = readRDS(filepath_unfiltered_trans_regime_bounds)

print("On to successful!")
# Check if simulation went correctly - sufficient length, null and trans
regime_bounds_successful_ = regime_bounds_trans_df %>%
  # filter(data_idx==1) %>%
  group_by(.data$data_idx, .data$regime_switch, .data$transition_steps) %>%
  # Match transition and null model
  dplyr::group_modify(~ match_trans_null_model(
    regime_bounds_trans = .x %>% dplyr::filter(trans_or_null == "transition"),
    regime_bounds_null =  .x %>% dplyr::filter(trans_or_null == "null"),
    min_length_regime = pars_template$min_length_regime,
    # Define desired regime filter
    regime_switch_list = make_filter_regime_switches(pars_template$min_length_regime,
                                                     pars_template$baseline_steps,
                                                     .y$transition_steps)[.y$regime_switch],
    pre_steps = pars_template$pre_steps, baseline_steps = pars_template$baseline_steps,
    transition_steps = .y$transition_steps) %>% select(-any_of(names(.y))), .keep = T)

if (nrow(regime_bounds_successful_) > 0){
  regime_bounds_successful = regime_bounds_successful_ %>%
  # Check for matching transition and null model per condition -> 2 timeseries per condition
  group_by(.data$data_idx, .data$regime_switch, .data$transition_steps, .data$baseline_steps) %>%
  # dplyr::summarise(n = n()) %>% as.data.frame()
  dplyr::filter(n() == 2) %>% ungroup() %>%
  # Only need X amount of models
  group_by(.data$trans_or_null, .data$regime_switch, .data$transition_steps, .data$baseline_steps) %>%
  dplyr::slice(1:pars_template$nr_required_models) %>% ungroup() %>%
  # Make sure there are enough transition steps
  arrange(.data$data_idx, .data$regime_switch, .data$transition_steps,
          .data$baseline_steps,
          .data$trans_or_null)

# regime_bounds_successful = dplyr::bind_rows(
#   # Null models should have correct first regime & regime switch as well as be of sufficient length
#   regime_bounds_trans_df  %>%
#     group_by(.data$regime_switch) %>%
#     group_modify(~ apply_filter_regime_switches(.x,
#                                         regime_switch_list[unlist(purrr::map(regime_switch_list, "regime_switch")) == .y$regime_switch],
#                                         trans_or_null = "null") %>%
#                    # Remove grouping variable from returned result as it will be appended anyway
#                    select(-.data$regime_switch), .keep = TRUE) %>%
#     # Add transition conditions
#     mutate(count = rep(length(pars_template$transition_steps), n())) %>%
#     tidyr::uncount(count) %>% dplyr::mutate(transition_steps = rep(pars_template$transition_steps, n() / length(pars_template$transition_steps))),
#   regime_bounds_trans_df  %>%
#     group_by(.data$regime_switch, .data$transition_steps) %>%
#     group_modify(~ apply_filter_regime_switches(.x, regime_switch_list[unlist(purrr::map(regime_switch_list, "regime_switch")) == .y$regime_switch], trans_or_null = "transition") %>%
#                    # Remove grouping variable from returned result as it will be appended anyway
#                    select(-c(.data$regime_switch, .data$transition_steps)), .keep = TRUE)
# ) %>% distinct() %>%
#   # Add baseline condition
#   mutate(count = rep(length(pars_template$baseline_steps), n())) %>%
#   tidyr::uncount(count) %>%
#   dplyr::mutate(baseline_steps = rep(pars_template$baseline_steps, n() / length(pars_template$baseline_steps))) %>%
#   # Check for matching transition and null model per condition -> 2 timeseries per condition
#   group_by(regime_switch, data_idx, transition_steps, baseline_steps) %>%
#   # dplyr::summarise(n = n()) %>% as.data.frame()
#   dplyr::filter(n() == 2) %>%
#   ## Null models do not have an index where regime 2 started, so copy from corresponding transition condition
#   # Copy start regime 2 from transition condition
#   dplyr::mutate(transition_end_idx = regime2_start_idx[trans_or_null == "transition"] -1, seed_nr = cur_group_id()) %>%
#   ungroup() %>%
      # trans_or_null == "transition",
#       pars_template$pre_steps + pars_template$default_baseline_steps,
#       NA
#     ),
#     par_change_end_idx = ifelse(
#       trans_or_null == "transition",
#       pars_template$pre_steps + pars_template$default_baseline_steps + transition_steps,
#       NA
#     )
#   ) %>%
#   # Make sure there is enough baseline time for each model
#   dplyr::filter_at(c("transition_start_idx", "transition_end_idx",
#                      "baseline_start_idx", "baseline_end_idx"), ~ .x > 0) %>%
#   # Only need X amount of models
#   group_by(regime_switch, trans_or_null, transition_steps, baseline_steps) %>%
#   dplyr::slice(1:pars_template$nr_required_models) %>% ungroup() %>%
#   # Make sure there are enough transition steps
#   # dplyr::filter(transition_end_idx <= (pars_template$pre_steps + pars_template$default_baseline_steps + transition_steps)) %>%
#   arrange(data_idx,
#           regime_switch,
#           transition_steps,
#           baseline_steps,
#           trans_or_null)


# Check
regime_bounds_successful %>%
  # Should be n = 2, one transition and one null model
  group_by(data_idx, regime_switch, transition_steps, baseline_steps) %>%
  summarise(n = n()) %>% as.data.frame()

regime_bounds_successful %>%
  group_by(regime_switch, transition_steps, baseline_steps) %>%
  summarise(n = n()) %>% as.data.frame()
#
saveRDS(regime_bounds_successful, filepath_successful_regime_bounds)
regime_bounds_successful = readRDS(filepath_successful_regime_bounds)

}

# Per condition, check if simulation was successfull
forloop_cond = unique(purrr::map(forloop, function(x) {
  x[setdiff(names(x), "data_idx")]
}))
forloop_inner = get_forloop(data_idx = pars_template$data_idxs)

foreach(
  for_par_cond = rev(forloop_cond),
  .combine = 'cfun',
  .packages = c("bifurcationEWS", "dplyr", "foreach", "ggplot2")
) %dopar% {
  pars_cond <- utils::modifyList(pars_template, for_par_cond)
  pars_cond$subfolder1 = pars_cond$regime_switch
  if (pars_cond$trans_or_null == "transition") {
    pars_cond$filename = sprintf("%s_%dtransSteps",
                                 pars_cond$trans_or_null,
                                 pars_cond$transition_steps)
  } else if (pars_cond$trans_or_null == "null") {
    pars_cond$filename = pars_cond$trans_or_null
  }

  regime_bounds = regime_bounds_successful %>% merge(as.data.frame(for_par_cond))
  # add_par_as_cols(data.frame(NA), for_par_cond) %>% merge(regime_bounds_successful)

  if (nrow(regime_bounds) > 0){
    # Get peaks from all simulations of this condition
    all_peaks_df = foreach(for_par_inner = forloop_inner, .combine = 'rbind') %do% {
      pars = utils::modifyList(pars_cond, for_par_inner)
      # Get details of regime switch
      desired_switch = selected_regime_bounds %>% dplyr::filter(data_idx == pars$data_idx,
                                                                regime_switch == pars$regime_switch)

      if (nrow(desired_switch) == 1) {
        filepath_GLV = format_path(format_pars(pars))
        filepath_regimes = format_path(format_pars(utils::modifyList(pars, list(type_output = "regimes"))))

        print(filepath_GLV)
        if (file.exists(filepath_regimes)) {
          return(
            readRDS(filepath_regimes)$peaks_df %>%
              dplyr::rename(timestep_idx = bifpar_idx) %>%
              dplyr::mutate(
                data_idx = pars$data_idx,
                bifpar_value = unlist(readRDS(filepath_regimes)$bifpar_list)[timestep_idx]
              )
          )
        }
      } else {
        return(NULL)
      }
    }

    range_bifpar_value = range(all_peaks_df$bifpar_value) + (range(all_peaks_df$bifpar_value) * c(-0.05, 0.05))
    range_X = range(all_peaks_df$X) + (range(all_peaks_df$X) * c(-0.05, 0.05))

    all_peaks_df_X1 = all_peaks_df %>%
      dplyr::filter(variable == "X1") %>%
      merge(
        regime_bounds %>% select(
          data_idx,
          par_change_start_idx,
          par_change_end_idx,
          transition_start_idx,
          transition_end_idx,
          regime_switch_type
        ),
        all = TRUE
      )

    for (range_to_plot in c("full", "partial")){
      if (range_to_plot == "partial"){
        x_range_to_plot = c(min(all_peaks_df_X1 %>% pull(transition_end_idx), na.rm = T) - 12,
                            max(all_peaks_df_X1 %>% pull(transition_end_idx), na.rm = T) + 6)
        # all_peaks_df_X1_ = all_peaks_df_X1 %>% dplyr::filter(timestep_idx > (transition_end_idx - 10) & timestep_idx < (transition_end_idx + 10))
      } else {
        x_range_to_plot = c(min(all_peaks_df_X1 %>% pull(timestep_idx)), max(all_peaks_df_X1 %>% pull(timestep_idx)))
        # all_peaks_df_X1_ = all_peaks_df_X1
      }
        # Per regime switch, plot all simulations together
      pl_peaks = all_peaks_df_X1 %>%
        ggplot() +
        # Bifurcation parameter value
        geom_line(aes(
          x  = timestep_idx,
          y = scale_range(bifpar_value, a = range_X[1], b = range_X[2])
        ),
        linewidth = .005,
        col = 'orange2') +
        geom_point(aes(x  = timestep_idx, y = X, col = minmax), size = .005) +
        viridis::scale_color_viridis(
          discrete = T,
          option = "magma",
          name = "Minimum or maximum of peak",
          begin = .1,
          end = .7
        ) +
        ggnewscale::new_scale_colour() +
        # Intended regime bounds
        geom_vline(
          aes(xintercept = par_change_start_idx,
              col = regime_switch_type),
          alpha = .5,
          linewidth = 1.1,
          linetype = 'dotted'
        ) +
        geom_vline(
          aes(xintercept = par_change_end_idx,
              col = regime_switch_type),
          alpha = .5,
          linewidth = 1.1,
          linetype = 'dotted'
        ) +
        viridis::scale_color_viridis(
          discrete = T,
          option = "plasma",
          name = "Parameter change",
          na.translate = F,
          begin = .5,
          end = .6
        ) +
        ggnewscale::new_scale_colour() +
        # Empirical regime bounds
        geom_vline(
          aes(xintercept = transition_start_idx, # regime1_end_idx,
              col = regime_switch_type),
          alpha = .5,
          linewidth = 1.1,
          linetype = 'dashed'
        ) +
        geom_vline(
          aes(xintercept = transition_end_idx, # regime2_start_idx,
              col = regime_switch_type),
          alpha = .5,
          linewidth = 1.1,
          linetype = 'dashed'
        ) +
        viridis::scale_color_viridis(
          discrete = T,
          option = "viridis",
          name = "Empirical regime bounds",
          na.translate = F,
          begin = .1,
          end = .2
        ) +
        scale_y_continuous(
          name = "X",
          limits = range_X,
          sec.axis = sec_axis(
            trans =  ~ rev_scale_range(
              .,
              a = range_X[1],
              b = range_X[2],
              min_x = range_bifpar_value[1],
              max_x = range_bifpar_value[2]
            ),
            name = "s"
          )
        ) +
        scale_x_continuous(name = "Index bifurcation parameter", limits = x_range_to_plot) +
        coord_cartesian(clip = 'off') +
        ggh4x::facet_wrap2(. ~ data_idx, ncol = 6, axes='all')

      filepath_image = format_path(format_pars(utils::modifyList(
        pars_cond,
        list(
          type_output = "figs",
          subfolder1 = "regimes",
          subfolder2 = pars$subfolder1,
          filename = sprintf("sim-comp-peaks-x1_%s_%s", range_to_plot, pars_cond$filename),
          file_ext = ".png"
        )
      )))
      save_plot(
        style_plot(pl_peaks),
        filepath_image,
        w = 50,
        h = 30,
        formats = ".png"
      )
  }
}
  return(NULL)
}

