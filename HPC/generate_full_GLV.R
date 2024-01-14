rerun = T
print("Start simulating full GLV!")
pars_template$nr_timesteps=pars_template$nr_timesteps_full
pars_template$min_length_regime=pars_template$min_length_regime_full

# Define foreach loop parameters
forloop = get_forloop(data_idx = pars_template$data_idxs)

start_t = Sys.time()
foreach(for_par = forloop,
        .combine = 'cfun',
        .packages = c(
          "bifurcationEWS", "dplyr", "ggplot2"
        ),
        .export = c("pars_template")
) %dopar% {
  .GlobalEnv$pars_template <- pars_template # Don't ask me why this is necessary...

  pars <- modify_list(pars_template, for_par)
  filepath_GLV = format_path(format_pars(pars))
  filepath_regimes = format_path(format_pars(modify_list(pars, list(type_output = "regimes"))))
  print(filepath_GLV)
  if (!file.exists(filepath_regimes) | rerun){
  if (!file.exists(filepath_GLV)){
    seed_nr = for_par$data_idx * 1000

    # Generate timeseries
    GLV = bifurcation_ts(model = GLV_model,
                         model_pars = pars$model_pars,
                         bifpar_pars = pars$bifpar_pars,
                         timestep = pars$timestep,
                         nr_timesteps = pars$nr_timesteps_full,
                         stopifregime = pars$stopifregime,
                         seed_nr = seed_nr,
                         deSolve_method = pars$deSolve_method,
                         fs = pars$fs,
                         X_names = pars$X_names,
                         max_iter = pars$max_iter)
    df = GLV$df
    head(df)

    saveRDS(GLV, filepath_GLV)
  }
  GLV = readRDS(filepath_GLV)
  # GLV$X0s = GLV$X0s %>% magrittr::set_colnames(c("bifpar_idx", pars$X_names))
  # saveRDS(GLV, filepath_GLV)

  # Find regimes
  regime_list = find_regimes(GLV,
                             thresh_node = pars$thresh_node,
                             thresh_coord_spread = pars$thresh_coord_spread_full,
                             thresh_peak_idx_spread=pars$thresh_peak_idx_spread,
                             thresh_full_band=pars$thresh_full_band,
                             min_length_regime = pars$min_length_regime,
                             nr_smooth = pars$nr_smooth_full,
                              factor_k = pars$factor_k,
                             variable_name = pars$variable_name,
                             keep_nr_timesteps = pars$keep_nr_timesteps
  )

  print(regime_list$periods)
  print(as.data.frame(regime_list$regimes))
  print(as.data.frame(regime_list$regime_bounds))

  saveRDS(regime_list, filepath_regimes)
  rm(regime_list)
  rm(GLV)
  }

  if (T){
    regime_list = readRDS(filepath_regimes)

    # Plot regime results
    peaks_df = regime_list$peaks_df %>%
      dplyr::filter(bifpar_idx > 1) # %>%
      # group_by(bifpar_idx, variable) %>%
      # arrange(.data$time_idx, .by_group=TRUE) %>%
      # Only transient
    # dplyr::filter(dplyr::row_number() <= round(n()*.2) ) %>% ungroup()
    # Discard transient
    # dplyr::filter(dplyr::row_number() >= round(n()*.8) ) %>% ungroup()

    point_size = 40 / length(unique(peaks_df$bifpar_idx))
    pl_peaks = peaks_df %>% dplyr::filter(variable=="X1") %>%
      ggplot() + geom_point(aes(x = bifpar_idx, y = X, col = minmax), size = point_size) +
      labs(x = 'Bifurcation parameter', y = "X1", title = "Transient bifurcation diagram")

    filepath_image = format_path(format_pars(modify_list(pars, list(type_output = "figs",
                                                                                   subfolder1 = "regimes",
                                                                          # filename = "only-transient_peaks-x1",
                                                                          # filename = "discard-transient_peaks-x1",
                                                                          filename = "peaks-x1",
                                                                          file_ext = ".png"))))
    save_plot(style_plot(pl_peaks), filepath_image, w = 40, h = 20, formats = ".png")

    # # Periodicity
    # pl_periods = peaks_df %>%
    #   ggplot() +
    #   geom_vline(data = regime_list$periods, aes(
    #     xintercept = bifpar_idx,
    #     col = period_bifpar), linewidth = .4, alpha=.8) +
    #   geom_point(aes(x = bifpar_idx, y = X),col='grey30',alpha = .8, size = .5) +
    #
    #   viridis::scale_color_viridis(
    #     discrete=TRUE,
    #     option = 'magma',
    #     begin = 0, end = .8, name = "Period",
    #     guide = guide_legend(override.aes = list(alpha = 1, linewidth = 3,                                                                                                            size = 5))) +
    #
    #   guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) + # x rows in legend
    #   ggh4x::facet_grid2(variable ~ .) +
    #   labs(x = 'Bifurcation parameter', y = "", title = "Transient bifurcation diagram")
    #
    # save_plot(style_plot(pl_periods), stringr::str_replace(filepath_image, "peaks-x1", "periods"), w = 40, h = 50, formats = ".png")


    # Regimes and regime boundaries
    pl_regimes = peaks_df %>%
      dplyr::filter(variable=="X1") %>%ggplot() +
      geom_rect(data=regime_list$regimes,aes(xmin=start_bifpar_idx,xmax=end_bifpar_idx,
                                             ymin=-Inf, ymax=Inf, fill = regime), linewidth = .5, alpha = .4)
    if (any(!is.na(regime_list$regime_bounds$regime1))){
      pl_regimes = pl_regimes +
        geom_vline(data = regime_list$regime_bounds, aes(
          xintercept = regime1_end_idx,
          col = regime_switch_type), linewidth = 1.3, linetype = 'dashed') +
        geom_vline(data = regime_list$regime_bounds, aes(
          xintercept = regime2_start_idx,
          col = regime_switch_type), linewidth = 1.3, linetype = 'dashed')  +
        viridis::scale_color_viridis(
            discrete=TRUE,
            option = 'magma',
            begin = 0, end = .8, name = "Regime switch (end of first and start of second regime)") +
        viridis::scale_fill_viridis(
          discrete=TRUE,
          option = 'turbo',
          begin = 0, end = 1, name = "Period",
          guide = guide_legend(override.aes = list(alpha = 1,
                                                   size = 3)))
    }
    pl_regimes = pl_regimes +
      geom_point(aes(x = bifpar_idx, y = X),col='grey30',alpha = .8, size = point_size*2) +
      guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) + # x rows in legend
      ggh4x::facet_grid2(variable ~ .) +
      labs(x = 'Bifurcation parameter', y = "", title = "Transient bifurcation diagram")
    save_plot(style_plot(pl_regimes), stringr::str_replace(filepath_image, "peaks-x1", "regimes-x1"), w = 40, h = 20, formats = ".png")

      pl_regimes = peaks_df %>% ggplot() +
      geom_rect(data=regime_list$regimes,aes(xmin=start_bifpar_idx,xmax=end_bifpar_idx,
                                             ymin=-Inf, ymax=Inf, fill = regime), linewidth = .5, alpha = .4)
      if (any(!is.na(regime_list$regime_bounds$regime1))){
        pl_regimes = pl_regimes +
      geom_vline(data = regime_list$regime_bounds, aes(
        xintercept = regime1_end_idx,
        col = regime_switch_type), linewidth = 1.3, linetype = 'dashed') +
      geom_vline(data = regime_list$regime_bounds, aes(
        xintercept = regime2_start_idx,
        col = regime_switch_type), linewidth = 1.3, linetype = 'dashed') +
      viridis::scale_color_viridis(
        discrete=TRUE,
        option = 'magma',
        begin = 0, end = .8, name = "Regime switch (end of first and start of second regime)") +
      viridis::scale_fill_viridis(
        discrete=TRUE,
        option = 'turbo',
        begin = 0, end = 1, name = "Period",
        guide = guide_legend(override.aes = list(alpha = 1,
                                                 size = 3)))
        }

      pl_regimes = pl_regimes +
        geom_point(aes(x = bifpar_idx, y = X),col='grey30',alpha = .8, size = point_size*2) +
        guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) + # x rows in legend
      ggh4x::facet_grid2(variable ~ .) +
      labs(x = 'Bifurcation parameter', y = "", title = "Transient bifurcation diagram")

    save_plot(style_plot(pl_regimes), stringr::str_replace(filepath_image, "peaks-x1", "regimes"), w = 40, h = 50, formats = ".png")
    print(filepath_image)
}
  return(NULL)
}

end_t = Sys.time()
print(end_t - start_t)


# Because of the period doubling sequence, we know which period we expect after period k, namely k*2. We check with max_k whether we had enough points in the timeseries to detect this period.
# regime_list$period_per_var %>% group_by(bifpar_idx, variable) %>% dplyr::summarise(m = paste0(unique(max_k), collapse = ","), .groups = 'drop') %>% head(n=30) %>% as.data.frame()
#
# regime_list$period_per_var %>% group_by(variable) %>% dplyr::summarise(m = max(k), .groups = 'drop') %>% head(n=30) %>% as.data.frame()


# Compile regimes of all simulations
regimes_list = foreach(for_par = forloop,
                       .packages = c(
                         "dplyr"
                       )
                       # .combine = 'rbind'
) %dopar% {

  pars <- modify_list(pars_template, for_par)
  filepath_regimes = format_path(format_pars(modify_list(pars, list(type_output = "regimes"))))
  if (file.exists(filepath_regimes)){
    regime_list = readRDS(filepath_regimes)
  }
  return(regime_list)

}

regimes_df = purrr::map(seq_along(regimes_list),
                        function(x){regimes_list[[x]]$regimes %>% add_par_as_cols(forloop[[x]])
                          }) %>%
  do.call(rbind, .) %>% as.data.frame()
filepath_all_regimes = format_path(format_pars(modify_list(pars_template, list(type_output = "regimes", filename = "all_regimes"))))
saveRDS(regimes_df, filepath_all_regimes)

regime_bounds_df = purrr::map(seq_along(regimes_list),
                              function(x){regimes_list[[x]]$regime_bounds %>% add_par_as_cols(forloop[[x]])}) %>%
  do.call(rbind, .) %>% as.data.frame()
filepath_all_regime_bounds = format_path(format_pars(modify_list(pars_template, list(type_output = "regimes", filename = "all_regime_bounds"))))
saveRDS(regime_bounds_df, filepath_all_regime_bounds)
regime_bounds_df = readRDS(filepath_all_regime_bounds)


# Plot summary of regimes
idxs_x = seq(1, length(pars_template$bifpar_list), length.out = 8) #c(1,25,50, 75,100)
pl_regime_comp = regimes_df %>% arrange(data_idx, start_bifpar_idx) %>%
  dplyr::filter(length_region >= pars_template$min_length_regime) %>%
  ggplot() +
  geom_segment(data = regime_bounds_df, aes(
    x = regime1_end_idx,
    xend = regime2_start_idx,
    y = data_idx-.25,  yend = data_idx+.25,
    col = regime_switch_type), linewidth = 1, alpha = 1) +
  # geom_segment(data = regime_bounds_df, aes(
    # x = regime2_start_idx,
    # xend = regime2_start_idx,
    # y = data_idx,  yend = data_idx,
    # col = regime_switch_type), linewidth = .25, linetype = 'dashed') +
  viridis::scale_color_viridis(discrete=TRUE,name = "Regime switch", option = "magma", end = .9) +
  ggnewscale::new_scale_colour() +
  geom_segment(aes(x = start_bifpar_idx, xend = end_bifpar_idx,
                   y = data_idx, yend = data_idx,
                   col = regime), linewidth = 3) +

  labs(y = "Dataset number", x = "Bifurcation parameter value", title = "Comparison of regimes across simulations") +
  viridis::scale_color_viridis(discrete=TRUE,name = "Regime", option = "turbo") +
  scale_x_continuous(breaks = idxs_x, labels = unlist(purrr::map(pars_template$bifpar_list, "s")[idxs_x]))

# print(style_plot(pl_regime_comp))

filepath_image = format_path(format_pars(modify_list(pars_template, list(type_output = "figs",
                                                                               subfolder1 = "regimes",
                                                                               filename = "compare_all_regimes",
                                                                               file_ext = ".png"))))
save_plot(style_plot(pl_regime_comp), filepath_image, w = 30, h = 30, formats = ".png")

