print("Start summarizing EWS across all bifurcations!")

library(bifurcationEWS)
library(foreach)
library(dplyr)
library(ggplot2)
source('visualise_helpers.R')

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
                                   ))


# Read all AUC and ROC ----------------------------------------------------
filepaths_AUC = list.files(
  path = "detGLV/AUC",
  pattern = "AUC*.RDS", full.names = TRUE, recursive = TRUE)
EWS_warnings_AUC = plyr::ldply(filepaths_AUC, function(i){readRDS(i)})

EWS_warnings_AUC %>% group_by(regime_switch, transition_steps, baseline_steps, sigma_obs_noise, downsample_fs) %>%
  dplyr::summarise(n = n()) %>% as.data.frame()
EWS_warnings_AUC %>% group_by(regime_switch, metric, transition_steps, baseline_steps, sigma_obs_noise, downsample_fs) %>%
  mutate(n = n(), id = cur_group_id() ) %>% as.data.frame() %>% filter(n!=1)

# ROC
filepaths_ROC = list.files(
  path = "detGLV/ROC",
  pattern = "ROC*.RDS", full.names = TRUE, recursive = TRUE)
EWS_warnings_ROC = plyr::ldply(filepaths_ROC, function(i){readRDS(i)})

# Inspect ROC
EWS_warnings_ROC %>% group_by(regime_switch, downsample_fs, sigma_obs_noise, metric) %>% dplyr::summarise(n = n(), .groups = 'drop') %>% arrange(n) %>% as.data.frame() %>% head()
EWS_warnings_ROC %>% pull(regime_switch) %>% unique
EWS_warnings_ROC %>% group_by(regime_switch) %>% dplyr::summarise(max(max_sigma_crit), min(max_sigma_crit), .groups = 'drop') %>% as.data.frame()
test = EWS_warnings_ROC %>% group_by(regime_switch, downsample_fs, sigma_obs_noise, metric) %>% dplyr::mutate(n = n()) %>% ungroup()
test %>% arrange(n, regime_switch, downsample_fs, sigma_obs_noise, metric) %>% as.data.frame() %>% head(20)

# Check how many did not "max out" sigma_crit -> fpr is still too high while tpr is also high
EWS_warnings_ROC %>%
  group_by(regime_switch, metric, transition_steps, baseline_steps, sigma_obs_noise, downsample_fs) %>%
  as.data.frame() %>%
  filter(fpr > .01)


EWS_warnings_ROC %>% group_by(regime_switch, metric, transition_steps, baseline_steps, sigma_obs_noise, downsample_fs) %>%
  mutate(n = n(), id = cur_group_id() ) %>% as.data.frame() %>% filter(n!=1) %>% head


# Direction of EWS for sigma_crit corresponding to best cut-off sigma_crit > --------

# Compute best threshold using Youden's J statistic (https://machinelearningmastery.com/threshold-moving-for-imbalanced-classification/)
criterion_df = EWS_warnings_ROC %>% dplyr::rename(sigma_crit = max_sigma_crit) %>%
group_by(regime_switch, metric, baseline_steps,
         transition_steps, downsample_fs, sigma_obs_noise) %>%
  group_modify(~ YoudensJ(.x)) %>%
  filter(!is.na(sigma_crit))
# %>%  filter(!grepl("LF0.05", metric, fixed = T))
criterion_df %>% as.data.frame %>% head()
criterion_df %>%  group_by(regime_switch, downsample_fs, sigma_obs_noise) %>% dplyr::summarise(n=n()) %>% as.data.frame

# Find corresponding warnings
filepaths_warnings = list.files(
  pattern = ".*nr.*RDS",
  path = "detGLV/warnings",
  full.names = TRUE, recursive = TRUE)

# Filter out null models
idxs = unlist(lapply(filepaths_warnings, function(i){grepl("null", i, fixed = T)}))
filepaths_warnings = filepaths_warnings[!idxs]

# filepaths_warnings = filepaths_warnings[1500:1530]
# Get warning per file
filepath_warnings_YoudensJ = format_path(format_pars(modify_list(pars_general_template, list(type_output = "warnings",
                                                                                     filename = "warnings_YoudensJ"))))

# Split into batches
nr_batches = 50
filepath_warnings_YoudensJ_batch = lapply(1:nr_batches, function(batch){stringr::str_replace(filepath_warnings_YoudensJ, ".RDS", paste0("_", batch, ".RDS"))})
batch_idxs = split(1:length(filepaths_warnings), cut(seq_along(1:length(filepaths_warnings)),nr_batches,labels = FALSE))

foreach(batch = 1:nr_batches) %do% {

  if (!file.exists(filepath_warnings_YoudensJ_batch[[batch]])){
    warnings_YoudensJ_batch = foreach(i = batch_idxs[[batch]],
                                                      filepath_warnings = filepaths_warnings[batch_idxs[[batch]]],
                                .packages = "dplyr",
                                                      .combine = 'rbind') %dopar% {
                                                        if(i %% 100 == 0){print(sprintf("%s i = %d %.4f %%", Sys.time(), i, i / length(filepaths_warnings)))}
                                                        warning_df = readRDS(filepath_warnings)$warning_df %>% select(-c("trans_or_null", "nr_patches", "nr_warnings", "baseline_steps", "transition_steps", "warning_signal"))
      return(merge(criterion_df, warning_df) %>% select(-c("sum_tp", "sum_fp", "sum_tn", "sum_fn")))
                                                      }
    saveRDS(warnings_YoudensJ_batch, filepath_warnings_YoudensJ_batch[[batch]])
    rm(warnings_YoudensJ_batch)
  }
  return(NULL)

}
# Compile batches
warnings_YoudensJ = foreach(batch = 1:nr_batches, .combine = 'rbind') %do% {
  return(readRDS(filepath_warnings_YoudensJ_batch[[batch]]))
}
saveRDS(warnings_YoudensJ, filepath_warnings_YoudensJ)

warnings_YoudensJ = readRDS(filepath_warnings_YoudensJ)

# Summarise
summ_dir_timing = function(x){
  summ_df =   x %>%
    dplyr::summarise(nr_positive = sum(score > 0, na.rm = T),
                     nr_negative = sum(score < 0, na.rm = T),
                     nr_zero = sum(score == 0, na.rm = T),
                     no_warning = sum(is.na(score)),
                     mean_score = mean(score, na.rm = T),
                     median_score = median(score, na.rm = T),
                     median_first_warning_bifpar_idx = median(first_warning_bifpar_idx, na.rm= T),
                     .groups = 'drop'
    )  %>%
    mutate(direction_summary = ifelse(nr_positive > 0 & nr_negative == 0 & nr_zero == 0, "Positive", ifelse(nr_positive == 0 & nr_negative > 0 & nr_zero == 0, "Negative", ifelse(nr_positive > 0 & nr_negative > 0 & nr_zero == 0, "Mixed", "?")))) %>%
    mutate(regime_switch_label =
             recode_factor(regime_switch, !!!purrr::flatten(regimes_switch_labels), .default = NULL, .ordered=T)) %>%
    mutate(regime_switch_class = dplyr::recode_factor(regime_switch, !!!regime_switch_to_class) %>% stringr::str_replace("-", "- ") %>% factor(levels = c("Fixed- Point", "Period- Doubling", "Period- Halving", "Chaotic")))  %>%
    mutate(metric_label =
             recode_factor(metric, !!!purrr::flatten(metric_labels), .default = NULL, .ordered=T)) %>%
    mutate(metric_class = dplyr::recode_factor(metric, !!!metric_to_class) %>% factor(levels = c("Generic", "Multivariate", "Spectral")))
  return(summ_df)
}

warnings_YoudensJ_per_cond = warnings_YoudensJ %>%
  group_by(metric, regime_switch, sigma_obs_noise, downsample_fs) %>%
  summ_dir_timing()

warnings_YoudensJ_summ = warnings_YoudensJ %>%
  group_by(metric, regime_switch) %>%
  summ_dir_timing()

warnings_YoudensJ_per_cond  %>%
  group_by(metric, regime_switch, sigma_obs_noise, downsample_fs) %>%
  dplyr::summarise(n = n()) %>% as.data.frame()

warnings_YoudensJ  %>%
  group_by(metric, regime_switch, sigma_obs_noise, downsample_fs) %>%
  dplyr::summarise(n = n()) %>% as.data.frame()

##
# filepaths_warnings = list.files(
#   pattern = "warnings*.RDS", full.names = TRUE, recursive = TRUE)
# filepaths_warnings_merged = stringr::str_replace(filepaths_warnings, ".RDS", sprintf("_merged_YoudensJ.RDS"))
# filepaths_warnings_cond_summ = stringr::str_replace(filepaths_warnings, ".RDS", sprintf("_cond_summ_YoudensJ.RDS"))
# filepaths_warnings_summ = stringr::str_replace(filepaths_warnings, ".RDS", sprintf("_summ_YoudensJ.RDS"))
#
# # Loop through all warning dataframes and find the row corresponding to the optimal sigma_crit; save result
# foreach(i = filepaths_warnings,
#         i_merged = filepaths_warnings_merged,
#         i_cond_summ = filepaths_warnings_cond_summ,
#         i_summ = filepaths_warnings_summ) %do% {
#   print(i)
#   warning_dfs = readRDS(i) %>% filter(trans_or_null == "transition",
#                                       sigma_crit <= max(criterion_df$sigma_crit, na.rm=T)) %>%
#     select(c(sigma_crit, metric, regime_switch,
#              first_warning_bifpar_idx, warning_signal,
#              score, data_idx, noise_iter, downsample_fs, sigma_obs_noise)) %>%
#     group_by(regime_switch, sigma_obs_noise, downsample_fs, metric) %>% group_split()
#
#    merged_df = foreach(warning_df = warning_dfs,.combine = 'rbind') %do% {
#
#      # nrow(warning_df %>% select(-c(sigma_crit,score,first_warning_bifpar_idx,warning_signal)) %>% unique())
#           merge(criterion_df, warning_df)
#    }
#
#    rm(warning_dfs)
#    saveRDS(merged_df, i_merged)
#
#    # Summarise direction and timing across condition and aggregating across conditions
#    merged_df_cond = merged_df %>%
#      group_by(metric, regime_switch, sigma_obs_noise, downsample_fs)
#    merged_df_aggr = merged_df %>%
#      group_by(metric, regime_switch)
#
#    summ_dir_timing = function(x){
#      summ_df =   x %>%
#        # # Remove outliers to not mess up the scaling of the plot
#        # filter((score > quantile(score, .01, na.rm = T, names = F)) & (score < quantile(score, .99, na.rm= T, names = F))) %>%
#        dplyr::summarise(nr_positive = sum(score > 0, na.rm = T),
#                         nr_negative = sum(score < 0, na.rm = T),
#                         nr_zero = sum(score == 0, na.rm = T),
#                         no_warning = sum(is.na(score)),
#                         mean_score = mean(score, na.rm = T),
#                         median_score = median(score, na.rm = T),
#                         median_first_warning_bifpar_idx = median(first_warning_bifpar_idx, na.rm= T),
#                         .groups = 'drop'
#        )  %>%
#        mutate(direction_summary = ifelse(nr_positive > 0 & nr_negative == 0 & nr_zero == 0, "Positive", ifelse(nr_positive == 0 & nr_negative > 0 & nr_zero == 0, "Negative", ifelse(nr_positive > 0 & nr_negative > 0 & nr_zero == 0, "Mixed", "?")))) %>%
#        mutate(regime_switch_label =
#                 recode_factor(regime_switch, !!!purrr::flatten(regimes_switch_labels), .default = NULL, .ordered=T)) %>%
#        mutate(regime_switch_class = dplyr::recode_factor(regime_switch, !!!regime_switch_to_class) %>% stringr::str_replace("-", "- ") %>% factor(levels = c("Fixed- Point", "Period- Doubling", "Period- Halving", "Chaotic")))  %>%
#        mutate(metric_label =
#                 recode_factor(metric, !!!purrr::flatten(metric_labels), .default = NULL, .ordered=T)) %>%
#        mutate(metric_class = dplyr::recode_factor(metric, !!!metric_to_class) %>% factor(levels = c("Generic", "Multivariate", "Spectral")))
#        return(summ_df)
#    }
#
#    summ_merged_df_cond = summ_dir_timing(merged_df_cond)
#    summ_merged_df = summ_dir_timing(merged_df_aggr)
#
#   rm(merged_df)
#   saveRDS(summ_merged_df_cond, i_cond_summ)
#   saveRDS(summ_merged_df, i_summ)
#   rm(summ_merged_df)
#   return(NULL)
# }
#
# # Read dataframes with selected warnings
# warnings_YoudensJ_per_cond = foreach(i =filepaths_warnings_cond_summ,
#                                              .combine = 'rbind') %do% {
#                                                return(readRDS(i))
#                                              }
# warnings_YoudensJ_summ = foreach(i =filepaths_warnings_summ,
#                                              .combine = 'rbind') %do% {
#                                                return(readRDS(i))
#                                              }
# format(object.size(warnings_YoudensJ_per_cond), "Mb")
# format(object.size(warnings_YoudensJ_summ), "Mb")
#
# warnings_YoudensJ_summ %>% as.data.frame() %>% head

# Direction EWS heatmap ---------------------------------------------------
# Direction aggregated across conditions
pl_dir_across_cond = warnings_YoudensJ_summ  %>%
  rowwise() %>%
  mutate(perc_pos = nr_positive / (nr_positive + nr_negative + nr_zero) * 100) %>% ungroup() %>%
  ggplot() +
  geom_tile(aes(x = metric_label,
                y = 1,
                fill=perc_pos
  )) +
  scico::scale_fill_scico(palette="roma",
                          name = "Direction of warnings",
                          limits = c(0, 100),
                          breaks = c(0, 25, 50, 75, 100),
                          labels = c("All negative", "", "Equally mixed", "", "All positive"),
                          direction = -1,
                          guide = guide_colorbar(
                            barheight = 1.15,
                            barwidth = 20,
                            title.position = "top",
                            title.vjust = 0.5,
                            title.hjust = 0.5
                          )) +

  scale_y_discrete(expand = c(0,0), position = 'right',
                   labels = scales::parse_format())+
  scale_x_discrete(expand = c(0,0),
                   labels = scales::parse_format())+
  ggh4x::facet_nested(regime_switch_class + regime_switch_label ~ metric_class,
                      strip = ggh4x:: strip_nested(size = "variable"),
                      scale = "free",
                      space ='free_x',
                      switch = "y",
                      labeller = labeller(
                        regime_switch_class = label_wrap_gen(width = 50),
                        regime_switch_label = label_wrap_gen(width = 50),
                        metric_class = label_parsed
                      )
  )  +
  labs(x="", y = "",  title = "")


filepath_image = format_path(format_pars(modify_list(
  pars_general_template,
  list(
    type_output = "figs",
    subfolder1 = "summary",
    filename = sprintf("direction-heatmap_EWS-across-cond_YoudensJ"),
    file_ext = ".pdf"
  )
)))

save_plot(
  style_plot(pl_dir_across_cond) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 22)
    ) +
    theme(axis.text.x = element_text(angle = 45,
                                     # vjust = 1,
                                     hjust=1, size = 13
    )) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.text = element_text(colour = 'gray10',
                                    # margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    )) +
    theme(plot.title = element_text(vjust=1)) +
    theme(strip.background = element_rect(color = col_facet_labels, fill = NA, linewidth=1.5),
          panel.border = element_rect(color = NA, fill = NA, linewidth = 1.5),
          strip.text.y.left = element_text(
            angle = 0,
          ),
    ),
  filepath_image,
  w = 15.5, h = 9.5,
  formats = ".pdf")
print(filepath_image)

# Direction per condition
pl_dir_per_cond = warnings_YoudensJ_per_cond  %>%
  mutate(downsample_fs_label = latex2exp::TeX(sprintf("$f_{s} = %.2f$", downsample_fs), output = 'character')) %>%
  mutate(sigma_obs_noise_label = latex2exp::TeX(sprintf("$sigma_{obs} = $%.4f", sigma_obs_noise), output = 'character')) %>%
  ggplot() +
  geom_tile(aes(x = metric_label,
                y = 1,
                fill=direction_summary
  )) +
  scale_fill_manual(name = "Direction EWS", labels=c("Increase", "Decrease", "Mixed"), breaks=c("Positive", "Negative", "Mixed"),
                    values =scales::brewer_pal(palette="Spectral")(11)[c(9, 2, 5)]
                    ) +
  scale_y_discrete(expand = c(0,0), position = 'right',
                   labels = scales::parse_format())+
  scale_x_discrete(expand = c(0,0),
                   labels = scales::parse_format())+
   ggh4x::facet_nested(regime_switch_class + regime_switch_label + sigma_obs_noise_label + downsample_fs_label ~ metric_class,
                      strip = ggh4x:: strip_nested(size = "variable"),
                      scale = "free",
                      space ='free_x',
                      switch = "y",
                      labeller = labeller(
                        regime_switch_class = label_wrap_gen(width = 50),
                        regime_switch_label = label_wrap_gen(width = 50),
                        sigma_obs_noise_label = label_parsed,
                        downsample_fs_label = label_parsed,
                        metric_class = label_parsed
                        )
  )  +
  labs(x="", y = "", title = "") +
  guides(
         fill = guide_legend(title.position = "top", ncol = 1, order=1,drop=F,
                             override.aes = list(alpha = 1))
  )


filepath_image = format_path(format_pars(modify_list(
  pars_general_template,
  list(
    type_output = "figs",
    subfolder1 = "summary",
    filename = sprintf("direction-heatmap_EWS-per-cond_YoudensJ"),
        file_ext = ".pdf"
  )
)))

save_plot(
  style_plot(pl_dir_per_cond) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 22)
    ) +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust=1, size = 13
    )) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(strip.text = element_text(colour = 'gray10',
    )) +
    theme(plot.title = element_text(vjust=1)) +
    theme(strip.background = element_rect(color = col_facet_labels, fill = NA, linewidth=1.5),
          panel.border = element_rect(color = NA, fill = NA, linewidth = 1.5),
          strip.text.y.left = element_text(
            angle = 0,
          ),
    ),
  filepath_image,
  w = 18, h = 50,
  formats = ".pdf")
print(filepath_image)


#  Timing EWS -------------------------------------------------------------
# Summary timing
pl_timing = warnings_YoudensJ_summ  %>%
  # Set transition at t = 0
  mutate(median_first_warning_bifpar_idx = (median_first_warning_bifpar_idx-pars_general_template$transition_steps)) %>%
  # Manipulate yintercept to fit all EWS (not meaningful placement, for illustration purposes)
  mutate(yintercept = as.numeric(metric_label)*.1) %>%
  mutate(yintercept = scale(yintercept, scale = F)) %>%
  ggplot() +
  geom_hline(aes(yintercept = median(unique(yintercept))), col = col_facet_labels, linewidth = .75, alpha = .75) +
  geom_jitter(aes(x = median_first_warning_bifpar_idx,
                  y = yintercept,
                  color = metric_label),
              height = 0,width = 0,
              size = 1.5, alpha = .95) +
  ggh4x::facet_nested(regime_switch_class + regime_switch_label ~ .,
                      switch = 'y',
                      strip = ggh4x:: strip_nested(size = "variable"),
                      remove_labels = T, axes = "all",
                      labeller = labeller(regime_switch_class = label_wrap_gen(width = 12),
                                          regime_switch_label = label_wrap_gen(width = 20))
                      # labeller = label_wrap_gen(width = 12)
  ) +
  viridis::scale_color_viridis(name="Metric", discrete = T,
                               labels = scales::parse_format(),
                               option = "turbo", begin = 0, end = 1) +
  scale_x_continuous(n.breaks = 5, limits = c(-100,0), expand = c(0, 0)) +
  labs(x = "Steps before transition",
       y= "", title = ""  ) +
  scale_y_continuous(position = 'right', expand = c(0.2,0.2)) +
  coord_cartesian(clip = 'off') #+

filepath_image = format_path(format_pars(modify_list(
  pars_general_template,
  list(
    type_output = "figs",
    subfolder1 = sprintf("summary"),
    filename =sprintf("median-timing_EWS-across-cond_YoudensJ"),
    file_ext = ".pdf"
  )
)))


save_plot(
  style_plot(pl_timing)  +
    theme(
      strip.background = element_rect(fill = NA, linewidth=1),
      panel.border = element_rect(color = NA, fill = NA, linewidth = 1.5),
      strip.text.y.left = element_text(
        angle = 0, size = 10,
        margin = margin(l = 0.1, b = 0.01, r = 0.1, t = 0.01, "cm")

      ),
      strip.text.x = element_text(
        margin = margin(0.01, 0.01, 0.01, 0.01, "cm")
      )
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    theme(strip.text = element_text(colour = 'gray10')) +
    theme(
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 12),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank())  +
    theme(legend.position = "right",
          legend.key.size = unit(.4, "cm"),
          legend.spacing.y = unit(.01, 'cm')) +
    guides(
      color = guide_legend(title.position = "top", ncol = 1,
                           byrow = T, reverse = T,
                           override.aes = list(alpha = 1, size = 2.1)
      )
    ) +
    theme(panel.spacing = unit(0.225, "cm")),

  filepath_image,
  w = 7,
  h = 10,
  formats = ".pdf")

# AUC ---------------------------------------------------------------------
pl_summ_discr = EWS_warnings_AUC %>%
    group_by(regime_switch, metric) %>%
    dplyr::summarise(median_AUC = median(AUC, na.rm=TRUE), sd_AUC = sd(AUC, na.rm=TRUE),
                     mad_AUC = mad(AUC, na.rm= TRUE),
                     mean_AUC = mean(AUC, na.rm=TRUE), .groups='drop') %>%
    rowwise() %>%
    mutate(AUC_class = get_AUC_class(median_AUC)) %>% ungroup() %>%
    mutate(regime_switch_label =
             recode_factor(regime_switch, !!!purrr::flatten(regimes_switch_labels), .default = NULL, .ordered=T)) %>%
    mutate(regime_switch_class = dplyr::recode_factor(regime_switch, !!!regime_switch_to_class) %>% stringr::str_replace("-", "- ") %>% factor(levels = c("Fixed- Point", "Period- Doubling", "Period- Halving", "Chaotic")))  %>%
    mutate(metric_label =
             recode_factor(metric, !!!purrr::flatten(metric_labels), .default = NULL, .ordered=T)) %>%
    mutate(metric_class = dplyr::recode_factor(metric, !!!metric_to_class) %>% factor(levels = c("Generic", "Multivariate", "Spectral")))  %>%
    ggplot() +
    geom_tile(aes(x = metric_label,
                  y = 1,
                  fill=AUC_class
    )) +
    viridis::scale_fill_viridis(discrete = TRUE, name = "AUC class", drop = F,
                                labels = scales::parse_format(), option = "inferno"
    ) +
  scale_y_discrete(expand = c(0,0), position = 'right',
                   labels = scales::parse_format())+
    scale_x_discrete(expand = c(0,0),
                     labels = scales::parse_format())+
    ggh4x::facet_nested(regime_switch_class + regime_switch_label ~ metric_class,
                        strip = ggh4x:: strip_nested(size = "variable"),
                        scale = "free",
                        space ='free_x',
                        switch = "y",
                        labeller = labeller(metric_class = label_parsed,
                                            regime_switch_class = label_wrap_gen(width =9),
                                            regime_switch_label = label_wrap_gen(width =100))
    )  +
    labs(x="", y = "")  +
    guides(size = guide_legend(title.position = "top", ncol = 1, order=2, drop=F,
                               override.aes = list(alpha = 1
                                                   # color = 'gray30'
                               )),
           fill = guide_legend(title.position = "top", ncol = 2, order=1,drop=F,
                               override.aes = list(alpha = 1
                               ))
    )

  filepath_image = format_path(format_pars(modify_list(
    pars_general_template,
    list(
      type_output = "figs",
      subfolder1 = "summary",
      filename ="AUC",
      file_ext = ".pdf"# ".png"
    )
  )))
  save_plot(
    style_plot(pl_summ_discr) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 19),
            legend.title = element_text(size = 22)) +
      theme(axis.text.x = element_text(angle = 45,
                                       # vjust = 1,
                                       hjust=1, size = 13
      )) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      theme(strip.text = element_text(colour = 'gray10',
                                      # margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
      )) +
      theme(strip.background = element_rect(color = col_facet_labels, fill = NA, linewidth=1.5),
            panel.border = element_rect(color = NA, fill = NA, linewidth = 1.5),
            strip.text.y.left = element_text(
              angle = 0,
            ),
      ),
    filepath_image,
    w = 15,h = 10,
    formats = ".pdf")
  print(filepath_image)

