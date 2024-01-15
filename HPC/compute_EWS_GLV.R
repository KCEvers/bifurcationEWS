print("Start computing EWS!")
rerun = T
pars_template$nr_timesteps=pars_template$nr_timesteps_trans

# Define EWS functions
uni_metrics = c(
  "mean" = mean,
  "variance" = var,
  "COV" = get_COV,
  "autocorrelation" = get_autocorr,
  "skewness" = moments::skewness,
  "kurtosis" = moments::kurtosis,
  "spectral_exp" = get_spectral_exp,
  "spectral_ratio" = get_spectral_ratio
  # "Smax" = get_Smax
  # "Hurst" = get_Hurst_exp
)
multi_metrics = c(
  "meanAbsCrossCorr" = get_conn,
  "largestLambdaCovMatrix" = eigenvalue,
  "spatial_variance" = spatial_variance,
  "spatial_skewness" = spatial_skewness,
  "spatial_kurtosis" = spatial_kurtosis,
  "Smax" = get_Smax
  # "RQA" = runRQA
)


# "RQA" = list(
#   emDim = 1,
#   emLag = 1,
#   theiler = 1,
#   distNorm = "max",
#   targetValue = .05
# )
# )

filepath_successful_regime_bounds = format_path(format_pars(modify_list(
  pars_template,
  list(type_output = "regimes", filename = "all_trans_regime_bounds")
)))
regime_bounds_successful = readRDS(filepath_successful_regime_bounds)

regime_bounds_successful %>%
  group_by(regime_switch) %>%
  dplyr::summarise(n = length(unique(data_idx))) %>% as.data.frame()

forloop_ = c(
  get_forloop(
    data_idx = pars_template$data_idxs,
    regime_switch = pars_template$select_regime_switches,
    downsample_fs = pars_template$downsample_fs,
    transition_steps = pars_template$transition_steps,
    trans_or_null = c("transition", "null"),
    sigma_obs_noise = pars_template$sigma_obs_noise,
    noise_iter = 1:pars_template$nr_noise_iters
  )
) %>% purrr::map(function(x){
  if (x$trans_or_null == "null" & x$transition_steps != pars_template$default_transition_steps) {
    return(NULL)
  } else if (x$sigma_obs_noise == 0 & x$noise_iter != 1) {
    return(NULL)
  } else {
    return(x)
  }
})
forloop = forloop_[!unlist(purrr::map(forloop_, is.null))]

label_metrics = c(
  "meanAbsCrossCorr" = "Mean Abs. Cross-Correlation",
  "largestLambdaCovMatrix" = "Max. Eigenvalue of Cov. Matrix",
  "spatial_variance" = "Spatial Variance",
  "spatial_skewness" = "Spatial Skewness",
  "spatial_kurtosis" = "Spatial Kurtosis",
  "DET" = "Determinism",
  "MAX_dl" = "Max. Diagonal Line Length",
  "DIV_dl" = "Divergence Diagonal Lines",
  "N_dl" = "Number Diagonal Lines",
  "MEAN_dl" = "Mean Diagonal Line",
  "ENTrel_dl" = "Entropy Diagonal Lines (Relative)",
  "CoV_dl" = "Diagonal Coefficient of Variation",
  "LAM_vl" = "Laminarity",
  "MAX_vl" = "Max. Vertical Line Length",
  "DIV_vl" = "Divergence Vertical Lines",
  "N_vl" = "Number Vertical Lines",
  "TT_vl" = "Mean Vertical Line",
  "ENTrel_vl" = "Entropy Vertical Lines (Relative)",
  "CoV_vl" = "Vertical Coefficient of Variation",
  "LAM_wl" = "Laminarity (White Lines)",
  "MAX_wl" = "Max. White Line Length",
  "DIV_wl" = "Divergence White Lines",
  "N_wl" = "Number White Lines",
  "MEAN_wl" = "Mean White Line",
  "ENTrel_wl" = "Entropy White Lines (Relative)",
  "CoV_wl" = "White Coefficient of Variation",
  "RR" = "Recurrence Rate",
  "emRad" = "Radius",
  "DET_RR" = "Determinism / Recurrence Rate",
  "LAM_DET"  = "Laminarity / Determinism",
  "SING_rate" = "Singularity Rate",
  "SING_N" = "Number of Isolated Points",
  "REP_av" = "Repetitiveness (Average)"
)


print(sprintf("%d conditions", length(forloop)))


# filepaths = foreach(
#   for_par =forloop,
#   .combine = 'cfun',
#   .packages = c("bifurcationEWS", "dplyr", "ggplot2"),
#   .export = c("pars_template")
# ) %dopar% {
#   .GlobalEnv$pars_template <- pars_template # Don't ask me why this is necessary...
#
#   pars <- modify_list(pars_template, for_par)
#   pars$subfolder1 = for_par$regime_switch
#
#   if (pars$trans_or_null == "transition") {
#     filename_GLV = sprintf("%s_%dtransSteps",
#                            pars$trans_or_null,
#                            pars$transition_steps)
#
#     filename_EWS = sprintf(
#       "%s_%dtransSteps_sigmaObs%.4f_iter%04d",
#       pars$trans_or_null,
#       pars$transition_steps,
#       pars$sigma_obs_noise, pars$noise_iter
#     )
#   } else if (pars$trans_or_null == "null") {
#     filename_GLV = sprintf("%s", pars$trans_or_null)
#     filename_EWS = sprintf("%s_sigmaObs%.4f_iter%04d",
#                            pars$trans_or_null,
#                            pars$sigma_obs_noise, pars$noise_iter)
#   }
#
#   filepath_GLV = format_path(
#     format_pars(modify_list(pars,
#                                   list(filename = filename_GLV))))
#   filepath_regimes = format_path(format_pars(modify_list(
#     pars, list(type_output = "regimes", filename = filename_GLV)
#   )))
#
#   # Update downsampling parameters
#   pars$win_size = round(pars_template$fs / pars$downsample_fs)
#
#   filepath_EWS = format_path(format_pars(modify_list(
#     pars, list(type_output = "EWS", filename = filename_EWS, fs = pars$downsample_fs)
#   )))
#
#   if (file.exists(filepath_EWS)){
#     err1 = "try-error" %in% class(try(readRDS(filepath_EWS)))
#   } else {
#     err1 = FALSE
#   }
#   print(err1)
#   if (err1){
#     file.remove(filepath_EWS)
#   }
# }


start_t = Sys.time()
foreach(
  for_par = forloop[221:1000],
  .combine = 'cfun',
  .packages = c("bifurcationEWS", "dplyr", "ggplot2"),
  .export = c("pars_template")
) %dopar% {
  .GlobalEnv$pars_template <- pars_template # Don't ask me why this is necessary...

  pars <- modify_list(pars_template, for_par)
  pars$subfolder1 = for_par$regime_switch

  if (pars$trans_or_null == "transition") {
    filename_GLV = sprintf("%s_%dtransSteps",
                           pars$trans_or_null,
                           pars$transition_steps)

    filename_EWS = sprintf(
      "%s_%dtransSteps_sigmaObs%.4f_iter%04d",
      pars$trans_or_null,
      pars$transition_steps,
      pars$sigma_obs_noise, pars$noise_iter
    )
  } else if (pars$trans_or_null == "null") {
    filename_GLV = sprintf("%s", pars$trans_or_null)
    filename_EWS = sprintf("%s_sigmaObs%.4f_iter%04d",
                           pars$trans_or_null,
                           pars$sigma_obs_noise, pars$noise_iter)
  }

  filepath_GLV = format_path(
    format_pars(modify_list(pars,
                                  list(filename = filename_GLV))))
  filepath_regimes = format_path(format_pars(modify_list(
    pars, list(type_output = "regimes", filename = filename_GLV)
  )))

  # Update downsampling parameters
  pars$win_size = round(pars_template$fs / pars$downsample_fs)

  filepath_EWS = format_path(format_pars(modify_list(
    pars, list(type_output = "EWS", filename = filename_EWS, fs = pars$downsample_fs)
  )))


  # Check if simulation was successful
  regime_bounds = data.frame(for_par) %>% merge(regime_bounds_successful)

  if (nrow(regime_bounds) > 0) { # There may be multiple rows because of the baseline condition

    # Count baseline back
    baseline_start = min(regime_bounds$baseline_start_idx)
    transition_end = max(regime_bounds$transition_end_idx)
    bifpar_idx_ =  seq(baseline_start, transition_end)
    print(filepath_EWS)

    # Check if all desired EWS are in there
    if (file.exists(filepath_EWS) & !rerun){
      split_df_EWS_old = readRDS(filepath_EWS) %>%
        filter(bifpar_idx %in% bifpar_idx_)
      bifpar_idx_to_do = setdiff(bifpar_idx_, unique(split_df_EWS_old$bifpar_idx))

      # # Find which metrics aren't in the dataframe
      # uni_metrics_todo = uni_metrics[!unlist(lapply(names(uni_metrics), function(i){any(grepl(i, unique(split_df_EWS_old$metric), fixed = T)) }))]
      # multi_metrics_todo = multi_metrics[!unlist(lapply(names(multi_metrics), function(i){any(grepl(i, unique(split_df_EWS_old$metric), fixed = T)) }))]
      uni_metrics_todo = uni_metrics
      multi_metrics_todo = multi_metrics
    } else {
      bifpar_idx_to_do = bifpar_idx_
      uni_metrics_todo = uni_metrics
      multi_metrics_todo = multi_metrics
    }

    if (file.exists(filepath_GLV) &
        (((length(uni_metrics_todo) > 0) | length(multi_metrics_todo) > 0) | rerun)) {

      df = readRDS(filepath_GLV)$df %>%
        filter(bifpar_idx %in% bifpar_idx_)

      # Downsample with new sampling frequency
      if (pars$win_size != 1){
         df = downsample(df, pars$X_names, win_size = pars$win_size)
      }

      # Add observational noise
      noisy_df = add_obs_noise(
        df,
        pars$X_names,
        noise_mean = 0,
        noise_sigma = pars$sigma_obs_noise,
        noise_constant = pars$noise_constant,
        # Set seed constant across comparing transition and null models; adjust as the seed number cannot exceed the maximum integer
        seed_nr = round(as.numeric(regime_bounds$seed_nr) * as.numeric(pars$noise_iter) * as.numeric(pars$downsample_fs) * 100) # pars$noise_iter
      ) %>%
        filter(bifpar_idx %in% bifpar_idx_to_do)

      rm(df)

      # Compute EWS
      EWS_args = list(
        "spectral_exp" = list(
          fs = pars$downsample_fs,
          nr_timesteps = pars$nr_timesteps
        ),
        "spectral_ratio" = list(
          fs = pars$downsample_fs,
          nr_timesteps = pars$nr_timesteps,
          f_min_to_f_max = list(c(0.005, .5), c(.05, .5))
        ),
        "Smax" = list(
          fs = pars$downsample_fs,
          nr_timesteps = pars$nr_timesteps
        )
        # "Hurst" = list(
        #   fs = pars$downsample_fs,
        #   nr_timesteps = pars$nr_timesteps
        # )
        )
      start_t = Sys.time()

      split_df_EWS = run_bifEWS(noisy_df,# %>% dplyr::filter(bifpar_idx == 150),
                                # %>% dplyr::filter(bifpar_idx >= 150, bifpar_idx <= 159),
                                pars$X_names,
                                uni_metrics_todo,
                                multi_metrics_todo,
                                EWS_args = EWS_args,
                                save_intermediate = F
                                # do_parallel = TRUE
                                )
      print(Sys.time() - start_t)

      if (file.exists(filepath_EWS) & !rerun){
        split_df_EWS = dplyr::bind_rows(split_df_EWS_old, split_df_EWS) %>%
          filter(bifpar_idx %in% bifpar_idx_) %>%
          arrange(bifpar_idx, metric)
      }
      saveRDS(split_df_EWS, filepath_EWS)
      rm(noisy_df)
      rm(split_df_EWS)
}
      # Plot
      if (T & (pars$noise_iter == 1) & (pars$data_idx <= 5) & file.exists(filepath_EWS)){
        print("Plot EWS!")
        regime_list = readRDS(filepath_regimes)
        split_df_EWS = readRDS(filepath_EWS)
        df = readRDS(filepath_GLV)$df %>%
          dplyr::filter(bifpar_idx >= baseline_start, bifpar_idx <= transition_end)

        # Downsample with new sampling frequency
        if (pars$win_size != 1){
          df = downsample(df, pars$X_names, win_size = pars$win_size)
        }

        # Add observational noise
        noisy_df = add_obs_noise(
          df,
          pars$X_names,
          noise_mean = 0,
          noise_sigma = pars$sigma_obs_noise,
          noise_constant = pars$noise_constant, seed_nr = pars$noise_iter
        )
        rm(df)

         noisy_peaks_df = noisy_df %>% select(all_of(pars$X_names), bifpar_idx, time_idx) %>%
           tidyr::gather(variable, X, -setdiff(colnames(.), pars$X_names)) %>%
          merge(regime_list$peaks_df %>% select(variable, minmax, bifpar_idx,time_idx))
        #peaks_bifdiag(noisy_df, pars$X_names)

        regime_bounds_long = regime_bounds %>% select(
          par_change_start_idx,
          par_change_end_idx,
          transition_start_idx,
          transition_end_idx,
          regime1_end_idx,
          regime2_start_idx
          # regime_switch_type
        ) %>% distinct() %>% tidyr::gather(moment, idx)

        pl_regimes = regime_list$peaks_df %>%
          dplyr::filter(bifpar_idx > 1) %>%
          # merge(regime_bounds, all = TRUE) %>%
          ggplot() +
          # geom_rect(
          #   data = regime_list$regimes,
          #   aes(
          #     xmin = start_bifpar_idx - .5,
          #     xmax = end_bifpar_idx + .5,
          #     ymin = -Inf,
          #     ymax = Inf,
          #     fill = regime
          #   ),
          #   linewidth = .5,
          #   alpha = .4
        # ) +
        geom_vline(
          data = regime_bounds_long,
          aes(xintercept = idx,
              col = moment),
          linewidth = 1.3, alpha = .75,
          linetype = 'dashed'
        ) +
          viridis::scale_color_viridis(
            discrete = T,
            option = "viridis",
            name = "",
            na.translate = F,
            begin = .1,
            end = .7
          ) +
          geom_point(data = noisy_peaks_df %>% dplyr::filter(bifpar_idx > 1),
                     aes(x = bifpar_idx, y = X),
                     col = 'orange4',
                     alpha = .8,
                     size = 1
          ) +
          geom_point(
            aes(x = bifpar_idx, y = X),
            col = 'grey30',
            alpha = .8,
            size = 1
          ) +

          guides(color = guide_legend(ncol = 1)) + # x rows in legend
          ggh4x::facet_grid2(variable ~ .) +
          labs(x = 'Bifurcation parameter',
               y = "",
               title = "Transient bifurcation diagram")

        pl_EWS <- split_df_EWS %>%
          # merge(regime_bounds, all = TRUE) %>%
          ggplot() +
          geom_line(aes(x = bifpar_idx, y = value), linewidth = .45, alpha = .5, color ='grey50') +
          geom_point(aes(x = bifpar_idx, y = value), size = .6, color ='grey20') +
          geom_vline(
            data = regime_bounds_long,
            aes(xintercept = idx,
                col = moment),
            linewidth = 1.3, alpha = .75,
            linetype = 'dashed'
          ) +
          viridis::scale_color_viridis(
            discrete = T,
            option = "viridis",
            name = "",
            na.translate = F,
            begin = .1,
            end = .7
          ) +
          ggh4x::facet_wrap2(. ~ metric, ncol = 4, scales = "free_y") +
          labs(x = 'Bifurcation parameter',
               y = "",
               title = "Early Warning Signals")

        combo = cowplot::plot_grid(style_plot(pl_regimes) + theme(legend.position = 'none'),
                                   style_plot(pl_EWS) + theme(legend.position = 'bottom'),
                                   ncol = 2,
                                   rel_widths = c(1, 3))
        filepath_image = format_path(format_pars(modify_list(
          pars,
          list(
            type_output = "figs",
            subfolder1 = "EWS",
            subfolder2 = pars$subfolder1,
            filename =
              paste0("bifpar-EWS", filename_EWS),
            fs = pars$downsample_fs,
            file_ext = ".png"
          )
        )))
        save_plot(
          combo,
          filepath_image,
          w = 50,
          h = 40,
          formats = ".png"
        )
        print(filepath_image)
        rm(regime_list)
      }

  }
  return(NULL)

}

end_t = Sys.time()
print(end_t - start_t)


