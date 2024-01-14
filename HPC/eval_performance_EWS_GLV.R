print("Start evaluating performance of EWS!")
rerun = F

pars_template$nr_timesteps=pars_template$nr_timesteps_trans
filepath_successful_regime_bounds = format_path(format_pars(modify_list(
  pars_template,
  list(type_output = "regimes", filename = "all_trans_regime_bounds")
)))
regime_bounds_successful = readRDS(filepath_successful_regime_bounds)

filepath_all_warnings = format_path(format_pars(modify_list(pars_template, list(type_output = "warnings",
                                                                                      filename = "all_warnings"))))
filepath_EWS_warnings_ROC = format_path(format_pars(modify_list(pars_template, list(type_output = "warnings",
                                                                                          filename = "ROC"))))
filepath_EWS_warnings_AUC = format_path(format_pars(modify_list(pars_template, list(type_output = "warnings",
                                                                                          filename = "AUC"))))

forloop = get_forloop(data_idx = pars_template$data_idxs,
                      regime_switch = pars_template$select_regime_switches,
                      baseline_steps = pars_template$baseline_steps,
                      transition_steps = pars_template$transition_steps,
                      trans_or_null = c("transition", "null"),
                      downsample_fs = pars_template$downsample_fs,
                      sigma_obs_noise = pars_template$sigma_obs_noise,
                      noise_iter = 1:pars_template$nr_noise_iters
)

print(sprintf("%d conditions", length(forloop)))


start_t = Sys.time()

filepaths_warnings = foreach(for_par = forloop,
                             .packages = c(
                               "bifurcationEWS", "dplyr"
                             )
) %dopar% {

  pars <- modify_list(pars_template, for_par)
  pars$subfolder1 = for_par$regime_switch

  if (pars$trans_or_null == "transition"){
    filename_GLV = sprintf("%s_%dtransSteps", pars$trans_or_null, pars$transition_steps)
    filename_EWS = sprintf("%s_%dtransSteps_sigmaObs%.4f_iter%04d", pars$trans_or_null, pars$transition_steps, pars$sigma_obs_noise, pars$noise_iter)
  } else if (pars$trans_or_null == "null"){
    filename_GLV = sprintf("%s", pars$trans_or_null)
    filename_EWS = sprintf("%s_sigmaObs%.4f_iter%04d", pars$trans_or_null, pars$sigma_obs_noise, pars$noise_iter)
  }
  filename_warnings = sprintf("%s_%dbaseSteps_%dtransSteps_sigmaObs%.4f_iter%04d", pars$trans_or_null,
                              pars$baseline_steps, pars$transition_steps, pars$sigma_obs_noise, pars$noise_iter)

  filepath_GLV = format_path(format_pars(modify_list(pars, list(filename = filename_GLV))))
  filepath_EWS = format_path(format_pars(modify_list(pars, list(type_output = "EWS", filename=filename_EWS, fs = pars$downsample_fs))))
  filepath_warnings = format_path(format_pars(modify_list(pars, list(type_output = "warnings", filename=filename_warnings, fs = pars$downsample_fs))))

  # Check if simulation was successful
  regime_bounds = as.data.frame(for_par) %>% merge(regime_bounds_successful)

  if (nrow(regime_bounds) == 1){
    if (file.exists(filepath_GLV) & file.exists(filepath_EWS) & (!file.exists(filepath_warnings) | rerun)){
      print(filepath_warnings)

      # Count baseline back
      baseline_idx = round(seq(regime_bounds$baseline_start_idx, regime_bounds$baseline_end_idx))
      transition_idx = round(seq(regime_bounds$transition_start_idx, regime_bounds$transition_end_idx))

      # Read EWS
      split_df_EWS = readRDS(filepath_EWS)

      # Compute at which critical sigma value a warning occurs
      warning_df = get_warnings(split_df_EWS,
                                baseline_idx,
                                transition_idx,
                                sigma_crit_step = pars$sigma_crit_step,
                                # sigmas_crit =pars$sigmas_crit,
                                nr_consecutive_warnings = pars$nr_consecutive_warnings)$warning_df %>%
        dplyr::mutate(first_warning_bifpar_idx = first_warning_bifpar_idx - regime_bounds$transition_start_idx + 1)

      saveRDS(list(warning_df = warning_df %>% add_par_as_cols(for_par),
                   regime_bounds = regime_bounds %>% select(-c(all_of(pars$X_names), regime1, regime2))
      ), filepath_warnings)
    }

  }
  if (file.exists(filepath_warnings)){
    return(filepath_warnings)
  } else {
    return(NULL)
  }
}

end_t = Sys.time()
print(end_t - start_t)

EWS_warnings = foreach(filepath_warnings = filepaths_warnings,
                       for_par = forloop,
                       .combine='rbind') %dopar% {
                         if(!is.null(filepath_warnings)){
                           return(
                             readRDS(filepath_warnings)$warning_df
                           )
                         }
                       }
saveRDS(EWS_warnings, filepath_all_warnings)
print(head(EWS_warnings) %>% as.data.frame())

# Compute ROC
grouping_vars = setdiff(names(forloop[[1]]), c("data_idx","noise_iter", "trans_or_null"))
EWS_warnings_ROC = warnings_to_ROC(EWS_warnings, grouping_vars)
EWS_warnings_AUC = ROC_to_AUC(EWS_warnings_ROC, grouping_vars)

# Save
saveRDS(EWS_warnings_ROC, filepath_EWS_warnings_ROC)
saveRDS(EWS_warnings_AUC, filepath_EWS_warnings_AUC)


EWS_warnings = readRDS(filepath_all_warnings)
EWS_warnings_ROC = readRDS(filepath_EWS_warnings_ROC)
EWS_warnings_AUC = readRDS(filepath_EWS_warnings_AUC)

##### ROC
# if (FALSE){
plot_ROC <- function(pars_template, EWS_warnings_ROC_sub, grouping_keys){
  filename = grouping_keys %>% t() %>% as.data.frame()%>%
    tibble::rownames_to_column() %>%
    dplyr::mutate_all (~ trimws(.x))%>%
    apply(1, paste0, collapse = "_") %>% paste0(collapse="_")

  pl_ROC = EWS_warnings_ROC_sub %>%
    select(-c("sum_tp", "sum_fp", "sum_tn", "sum_fn", "acc")) %>%
    tidyr::gather(rate, sum_freq, -setdiff(colnames(.), c("fnr", "fpr", "tnr", "tpr"))) %>%
    dplyr::mutate(trans_or_null = ifelse(rate %in% c("fnr", "tpr"), "Transition", "Null" )) %>%
    mutate(
      rate = recode(
        rate,
        tpr = "% True Positive",
        fpr = "% False Positive",
        fnr = "% False Negative",
        tnr = "% True Negative"
      )
    ) %>%
    # dplyr::filter(metric == "mean_var1") %>%
    ggplot(aes(
      x = sigma_crit,
      y = sum_freq,
      colour = rate,
      fill = rate
    )) +
    viridis::scale_color_viridis(
      discrete = TRUE,
      option = "magma",
      begin = 0,
      end = 0.9,
      name = ""
    ) +
    viridis::scale_fill_viridis(
      discrete = TRUE,
      option = "magma",
      begin = 0,
      end = 0.9,
      alpha = .5,
      name = ""
    ) +
    geom_bar(position = "stack", stat = "identity") +
    ggh4x::facet_nested( metric ~ regime_switch + trans_or_null,
                         # scales = "free_y",
                         # switch = "y",
                         # axes='all',
                         labeller = label_parsed) +
    labs(y = "Rate", x = latex2exp::TeX("$\\sigma_{crit}$"), title = filename) +
    scale_x_continuous(expand = c(0, 0)
                       # limits = c(min(EWS_warnings_ROC_long$sigma_crit), max(EWS_warnings_ROC_long$sigma_crit))
    ) +
    scale_y_continuous(n.breaks = 2, limits = c(0,1))

  filepath_image = format_path(format_pars(modify_list(
    pars_template,
    list(
      type_output = "figs",
      subfolder1 = "ROC",
      # subfolder2 = pars$subfolder1,
      filename =filename,
      file_ext = ".png"
    )
  )))
  save_plot(style_plot(pl_ROC) + theme(legend.position = "bottom") +  theme(strip.text.y = element_text(angle = 0)),
            filepath_image,
            w = length(unique(EWS_warnings_ROC_sub$regime_switch)) * 10 + 15,
            h = 50,
            formats = ".png")
  print(filepath_image)
}

EWS_warnings_ROC_grouped=EWS_warnings_ROC%>%
  dplyr::group_by_at(c(setdiff(names(forloop[[1]]), c("data_idx","noise_iter", "regime_switch","trans_or_null")) ) )%>%# dplyr::group_keys()
  # dplyr::group_split() %>%
  dplyr::group_map(~ plot_ROC(pars_template,
                              .x, .y))
# }
##### AUC
source('visualise_helpers.R')

plot_AUC <- function(pars_template, EWS_warnings_AUC_sub, grouping_keys,
                     width = 14, height = 7.5
                     # width = 19, height = 8
){
  filename = grouping_keys %>% t() %>% as.data.frame()%>%
    tibble::rownames_to_column() %>%
    dplyr::mutate_all (~ trimws(.x))%>%
    apply(1, paste0, collapse = "_") %>% paste0(collapse="_")

  pl_summ_discr = EWS_warnings_AUC_sub %>% rowwise() %>%
    dplyr::mutate(AUC_class = get_AUC_class(AUC)) %>% ungroup() %>%
    mutate(regime_switch_label =
             recode_factor(regime_switch, !!!purrr::flatten(regimes_switch_labels), .default = NULL, .ordered=T)) %>%
    # mutate(regime_switch_class = dplyr::recode_factor(regime_switch, !!!regime_switch_to_class) %>% factor(levels = c("Fixed-Point", "Period-Doubling", "Period-Halving", "Chaotic")))  %>%
    mutate(metric_label =
             recode_factor(metric, !!!purrr::flatten(metric_labels), .default = NULL, .ordered=T)) %>%
    mutate(metric_class = dplyr::recode_factor(metric, !!!metric_to_class) %>% factor(levels = c("Generic", "Multivariate", "Spectral")))  %>%
    mutate(downsample_fs_label = latex2exp::TeX(sprintf("$f_{s} = %.2f$", downsample_fs), output = 'character')) %>%
    mutate(sigma_obs_noise_label = latex2exp::TeX(sprintf("$sigma_{obs} = $%.4f", sigma_obs_noise), output = 'character')) %>%
    ggplot() +
    geom_tile(aes(x = metric_label, y = 1,
                  fill=AUC_class
    )) +
    viridis::scale_fill_viridis(discrete = TRUE, name = "AUC class", drop = F,
                                labels = scales::parse_format(), option = "inferno"
    ) +
    scale_y_discrete(expand = c(0,0), position = 'right')+
    scale_x_discrete(expand = c(0,0),
                     labels = scales::parse_format())+
    ggh4x::facet_nested(sigma_obs_noise_label + downsample_fs_label ~ regime_switch_label + metric_class,
                        strip = ggh4x:: strip_nested(size = "variable"),
                        scale = "free",
                        space ='free_x',
                        switch = "y",
                        labeller = labeller(.cols = label_wrap_gen(width = 50), #label_parsed,
                                            .rows = label_parsed #label_wrap_gen(width = 12)
                        )
    )  +
    labs(x="", y = "", title = "")  +
    guides(
      fill = guide_legend(title.position = "top", ncol = 2, order=1,drop=F,
                          override.aes = list(alpha = 1
                          ))
    )

  filepath_image = format_path(format_pars(modify_list(
    pars_template,
    list(
      type_output = "figs",
      subfolder1 = "AUC",
      filename = sprintf("AUC-all-conditions_%s", filename),
      file_ext = ".pdf"# ".png"
    )
  )))

  save_plot(
    style_plot(pl_summ_discr) +
      theme(legend.position = "bottom") +
      theme(axis.text.x = element_text(angle = 45,
                                       size = 12,
                                       hjust=1
      )) +
      theme(strip.text = element_text(colour = 'gray10',
                                      margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
      )) +
      theme(strip.background = element_rect(color = col_facet_labels, fill = NA, linewidth=1.5),
            panel.border = element_rect(color = NA, fill = NA, linewidth = 1.5),
            strip.text.y.left = element_text(
              angle = 0,
              margin = margin(0.75, 0.5, 0.75, 0.5, "cm")
            ),
      ),
    filepath_image,
    w = width,
    h = height,
    formats = ".pdf")
  print(filepath_image)

}

EWS_warnings_AUC%>%
  # dplyr::group_by(regime_switch, baseline_steps, transition_steps)%>%# dplyr::group_keys()
  dplyr::group_by(regime_switch)%>%# dplyr::group_keys()
  dplyr::group_map(~ plot_AUC(pars_template,
                              .x, .y), .keep = T)

#
# # Summary AUC across conditions
# pl_summ_discr = EWS_warnings_AUC %>%
#   group_by(regime_switch, metric) %>%
#   dplyr::summarise(median_AUC = median(AUC, na.rm=TRUE), sd_AUC = sd(AUC, na.rm=TRUE),
#                    mad_AUC = mad(AUC, na.rm= TRUE),
#                    mean_AUC = mean(AUC, na.rm=TRUE), .groups='drop') %>%
#   rowwise() %>%
#   dplyr::mutate(median_AUC_class = get_AUC_class(median_AUC)) %>% ungroup() %>%
#   ggplot() +
#   geom_tile(aes(x = metric, y = regime_switch,
#                 fill=median_AUC_class
#   )) +
#   viridis::scale_fill_viridis(discrete = TRUE, name = "AUC class", labels = scales::parse_format(),
#                               option = "inferno", drop = F
#   ) +
#   # geom_point(aes(x = metric_full, y = regime_switch_name, size = mad_AUC_class),
#   #            color='bisque4', fill= 'white',# 'white',
#   #            # color='lemonchiffon',  fill = 'lemonchiffon',#'lemonchiffon',#'papayawhip' ,#'white',
#   #            shape = 21, alpha = .95, stroke = .45) +
#   # scale_size_manual(breaks = c("Low", "Medium", "High", "Very high"), values = c(.25, 1, 2, 3.25),
#   # labels = c(sprintf("Low (< %.2f)", low_sd),
#   # sprintf("Medium (< %.2f)", medium_sd),
#   # sprintf("High (< %.2f)", high_sd), sprintf("Very high (> %.2f)", high_sd)),
#   # drop=FALSE,
#   # name = latex2exp::TeX("MAD ${AUC}$"),
#   ## name = latex2exp::TeX("$\\sigma_{AUC}$") #"SD AUC"
# # ) +
# # colorspace::scale_fill_discrete_sequential(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
# #                                            drop=F,
# #                                            labels = labels,
# #                                            name = latex2exp::TeX("Median ${AUC}$"), #"SD AUC"
# #                                            nmax = 10, palette= "Greens 3",
# #                                            # palette='Heat 2',
# #                                            # palette='ag_GrnYl',
# #                                            # palette="Blue-Yellow",
# #                                            # palette="Tropic",
# #                                            # rev=TRUE,
# #                                            # order = rev(1:10),
# #                                            # palette = "YlOrBr", nmax = 15
# # )   + #
# scale_y_discrete(expand = c(0,0),
#                  labels = scales::parse_format())+
#   scale_x_discrete(expand = c(0,0),
#                    labels = scales::parse_format())+
#   # ggh4x::facet_nested(. ~ regime_switch ,
#   #                     # scale = "free",
#   #                     # space ='free_x',
#   #                     # switch = "y",
#   #                     labeller = label_parsed
#   # )  +
#   labs(x="", y = "")  +
#   guides(
#     # size = guide_legend(title.position = "top", ncol = 1, order=2, drop=F,
#                              # override.aes = list(alpha = 1
#         #                                         # color = 'gray30'
#                              # )),
#          fill = guide_legend(title.position = "top",
#                              ncol = 2,
#                              # order=1,drop=F,
#                              override.aes = list(alpha = 1))
#   )
#
# #
# filepath_image = format_path(format_pars(modify_list(
#   pars_template,
#   list(
#     type_output = "figs",
#     subfolder1 = "AUC",
#     filename = "AUC-median-across-conditions",
#     file_ext = ".png"
#   )
# )))
# save_plot(
#   style_plot(pl_summ_discr) +
#     theme(legend.position = "bottom") +
#     theme(axis.text.x = element_text(angle = 45,
#                                      # vjust = 1,
#                                      hjust=1
#     )) +
#     # theme(panel.grid.major = element_blank(),
#     #       panel.grid.minor = element_blank()) +
#     theme(strip.text = element_text(colour = 'gray10',
#                                     margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
#     )) +
#     theme(strip.background = element_rect(color = main_color, fill = NA, linewidth=1.5),
#           panel.border = element_rect(color = NA, fill = NA, linewidth = 1.5),
#           strip.text.y.left = element_text(
#             angle = 0,
#             margin = margin(0.75, 0.5, 0.75, 0.5, "cm")
#           ),
#     ),
#   stringr::str_replace(filepath_image, "AUC_per_condition", "AUC_across_conditions"),
#   w = 20,
#   h = 15,
#   formats = ".png")
#
