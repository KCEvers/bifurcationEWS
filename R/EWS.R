#' Compute early warning signals (EWS)
#'
#' @param x Dataframe
#' @inheritParams run_bifEWS
#'
#' @return Dataframe with early warning signals
#' @importFrom dplyr .data
#' @export
#' @examples
run_EWS <- function(x, uni_metrics, multi_metrics, EWS_args = list()){

  if (length(uni_metrics) > 0){
  uni_EWS = plyr::ldply(names(uni_metrics), function(j){
    plyr::ldply(1:ncol(x), function(i){
      if (j %in% names(EWS_args)){
        EWS_arg = EWS_args[[j]]
      } else {
        EWS_arg = list()
      }

      uni_x = do.call(uni_metrics[[j]], modify_list(list(x = x[,i]), EWS_arg))

      # In case of multiple outputs
      if (length(uni_x) > 1){
        uni_df = data.frame(metric = paste0(j, "_", names(uni_x), "_", "var", i),
                            value = matrix(unlist(uni_x), ncol = 1))
      } else {
        # Single output
        uni_df = data.frame(metric = paste0(j, "_", "var", i),
                            value = matrix(unlist(uni_x), ncol = 1))
      }

      return(uni_df)
    })
    # apply(x, 2, uni_metric)
  })
  } else {
  uni_EWS = data.frame()
  }

  if (length(multi_metrics) > 0){
  multi_EWS = plyr::ldply(names(multi_metrics), function(j){
    # do.call(multi_metrics[[j]], modify_list(list(x = x), EWS_args[j]))

    if (j %in% names(EWS_args)){
      EWS_arg = EWS_args[[j]]
    } else {
      EWS_arg = list()
    }

    multi_x = do.call(multi_metrics[[j]], modify_list(list(x = x), EWS_arg))

    # In case of multiple outputs
    if (length(multi_x) > 1){
      multi_df = data.frame(metric = paste0(j, "_", names(multi_x)),
                          value = matrix(unlist(multi_x), ncol = 1))
    } else {
      # Single output
      multi_df = data.frame(metric = paste0(j),
                          value = matrix(unlist(multi_x), ncol = 1))
    }

    return(multi_df)
  })
  } else {
    multi_EWS = data.frame()
  }
  # %>%
  #   do.call(rbind, .) %>% magrittr::set_rownames(NULL) %>%
  #   as.data.frame() %>% cbind(metric = names(multi_metrics)) %>%
  #   dplyr::rename(value = .data$V1)

  return(rbind(uni_EWS, multi_EWS))
}

#' Run EWS for all bifurcation parameter values
#'
#' @param df Complete dataframe with timeseries for all bifurcation parameter values
#' @param X_names Names of variables in model
#' @param uni_metrics List of univariate EWS
#' @param multi_metrics List of multivariate EWS
#' @param EWS_args List of parameters to pass to EWS functions, empty list if no parameters need to be passed
#' @param save_intermediate Logical; should intermediate results be saved and combined at the end or should they be kept in memory?
#' @param do_parallel Logical; should the analysis be run in parallel?
#'
#' @return Dataframe with EWS
#' @importFrom dplyr .data arrange group_by group_split all_of select mutate
#' @importFrom foreach foreach `%do%` `%dopar%`
#' @export
#'
#' @examples
run_bifEWS <- function(df, X_names, uni_metrics = c("Smax" = get_Smax),
                       multi_metrics = c("spatialSkewness" = spatial_skewness),
                       EWS_args = list("Smax" = list(fs = 1, nr_timesteps = 100)),
                       save_intermediate = FALSE, do_parallel = FALSE){

  if (do_parallel){
    save_intermediate = TRUE
  }

  # Split dataframe per bifurcation parameter
  split_df = df %>% as.data.frame() %>% group_by(.data$bifpar_idx) %>% group_split()

  if (!save_intermediate){
    # Keep all results in memory; don't save intermediate results but combine in one go
    split_df_EWS = split_df %>%
      lapply(function(df_){run_EWS(df_ %>% arrange(.data$time_idx) %>%
                                     select(all_of(X_names)) %>% as.matrix(),
                                   uni_metrics, multi_metrics, EWS_args = EWS_args) %>%
          mutate(bifpar_idx = unique(df_$bifpar_idx))
        }) %>%
      do.call(rbind, .) %>% as.data.frame()

  } else if (save_intermediate){
    # Save each separately and combine after forloop to save memory
    tmps = lapply(seq_along(split_df), function(x){tempfile(fileext = ".RDS")})

    if (!do_parallel){
      foreach(df_ = split_df,
                            tmp = tmps,
                            .combine = function(...){NULL}
            ) %do% {
                              EWS_df = run_EWS(df_ %>% arrange(.data$time_idx) %>%
                                        select(all_of(X_names)) %>% as.matrix(),
                                      uni_metrics, multi_metrics, EWS_args = EWS_args) %>%
                                mutate(bifpar_idx = unique(df_$bifpar_idx))

                              saveRDS(EWS_df, tmp)
                              return(NULL)
            }

    } else if (do_parallel){
    foreach(df_ = split_df,
            tmp = tmps,
            .packages = c("bifurcationEWS", "dplyr"),
            .combine = function(...){NULL}
    ) %dopar% {
      EWS_df = run_EWS(df_ %>% arrange(.data$time_idx) %>%
                         select(all_of(X_names)) %>% as.matrix(),
                       uni_metrics, multi_metrics, EWS_args = EWS_args) %>%
        mutate(bifpar_idx = unique(df_$bifpar_idx))

      saveRDS(EWS_df, tmp)
      return(NULL)
    }
    }

    # Compile results
    split_df_EWS = foreach(tmp = tmps, .combine = 'rbind') %do% {
      return(readRDS(tmp))
    }
    # Delete intermediate files
    lapply(tmps, unlink)
  }
  return(split_df_EWS)
}


#' Get warnings per critical value of sigma
#'
#' @inheritParams get_warnings
#' @param y Vector with EWS
#' @param bifpar_idx Bifurcation parameter index
#' @param z_score Z-score
#'
#' @return Dataframe with warnings
#' @export
#'
#' @examples
get_warnings_per_sigma <- function(y, bifpar_idx, z_score, sigma_crit_step,
                                   thresh_max_sigma_crit = 500,
                                   nr_consecutive_warnings = 1){

  # Get sequence of critical cut-off values
  max_sigma_crit = sigma_crit_step * ceiling(max(abs(z_score), na.rm = T)/sigma_crit_step)
  if (is.na(max_sigma_crit) | is.infinite(max_sigma_crit)){
    max_sigma_crit = sigma_crit_step
  } else if (max_sigma_crit > thresh_max_sigma_crit){
    max_sigma_crit = thresh_max_sigma_crit
  }
  sigmas_crit =  seq(sigma_crit_step, max_sigma_crit, by=sigma_crit_step)

  lapply(sigmas_crit, function(sigma_crit){

    warnings_one_sigma = get_warnings_one_sigma(bifpar_idx, z_score, sigma_crit, nr_consecutive_warnings)

    # if (nrow(warnings_one_sigma) == 0){
    #   # No warning
    #   first_warning_bifpar_idx = NA
    #   score = NA
    # } else {
      # First warning
      # first_warning_bifpar_idx = warnings_one_sigma$bifpar_idx[1]
      # score = warnings_one_sigma$z_score[1]
    # }
    return(data.frame(sigma_crit = sigma_crit,
                      first_warning_bifpar_idx = warnings_one_sigma$bifpar_idx[1],
                      score = warnings_one_sigma$z_score[1],
                      nr_warnings = nrow(warnings_one_sigma),
                      nr_patches = warnings_one_sigma$nr_patches[1]))
  }) %>% do.call(rbind, .) %>% as.data.frame() %>% return()

}



#' Find warnings for cut-off value sigma-crit
#'
#' @inheritParams get_warnings
#' @inheritParams get_warnings_one_sigma
#'
#' @return Dataframe with raw values of EWS as well as the confidence bands corresponding to sigma_crit; ready for plotting
#' @export
#' @importFrom dplyr group_by group_modify ungroup mutate mutate_at first last .data
#'
#' @examples
get_warnings_raw <- function(split_df_EWS, baseline_idx, transition_idx, sigma_crit = 2, nr_consecutive_warnings = 1){
  # Compute baseline mean and standard deviation
  EWS_z_scores = get_zscore_EWS(split_df_EWS, baseline_idx)

  # EWS_z_scores = merge(split_df_EWS, EWS_df_CI) %>% arrange(.data$bifpar_idx) %>%
  #   group_by(.data$metric) %>%
  #   # Compute z-scores (normalized): Convert z-score so that it can be compared to any sigma (i.e. you don't have to rerun the procedure for every sigma) using the following formula:
  #   # current > mu_baseline + (sigma_baseline * sigma_crit)
  #   # (current - mu_baseline) > (sigma_baseline * sigma_crit)
  #   # ((current - mu_baseline) / sigma_baseline) > (sigma_crit)
  #   # z = (current - mu_baseline) / sigma_baseline
  #   mutate(z_score = ((.data$value - .data$mean_w0) / .data$sd_w0)) %>% ungroup() %>%
  #   apply(2, unlist) %>% as.data.frame %>%
  #   mutate_at(c("bifpar_idx", "value", "z_score", "mean_w0", "sd_w0", "quantile_000", "quantile_100"), ~ as.numeric(as.character(.x))) %>%
  #   group_by(.data$metric) %>%
  #   mutate(bifpar_idx = round(as.numeric(as.character(.data$bifpar_idx)))) %>%
  #   ungroup()

  # Get warnings for ONE critical sigma
  warnings_raw_value = EWS_z_scores %>%
    # filter(.data$bifpar_idx %in% transition_idx) %>%
    group_by(.data$metric) %>%
    group_modify(~ get_warnings_one_sigma(bifpar_idx = .x$bifpar_idx, z_score = .x$z_score,
                                          sigma_crit= sigma_crit,
                                          nr_consecutive_warnings = nr_consecutive_warnings)) %>%
    ungroup() %>%
    # Get indices of the rows corresponding to these warnings
    merge(EWS_z_scores %>% tibble::rownames_to_column(), all.x=TRUE) %>% ungroup()

  # Indicate the warning signals
  EWS_z_scores$warning_signal = FALSE
  EWS_z_scores[warnings_raw_value$rowname, "warning_signal"] = TRUE

  # Add confidence bands and whether the warnings were on time
  winEWS_df = EWS_z_scores %>%
    group_by(.data$metric) %>%
    mutate(ymin = .data$mean_w0 - sigma_crit * .data$sd_w0,
                  ymax = .data$mean_w0 + sigma_crit * .data$sd_w0,
                  warning_signal = (.data$value <= .data$ymin) | (.data$value >= .data$ymax),
                  value_no_warning = ifelse((.data$value > .data$ymin) & (.data$value < .data$ymax), .data$value, NA),
                  value_warning_on_time = ifelse(.data$warning_signal & .data$bifpar_idx %in% transition_idx, .data$value, NA),
                  value_warning_baseline = ifelse(.data$warning_signal & .data$bifpar_idx %in% baseline_idx, .data$value, NA),
                  value_warning_too_late = ifelse(.data$warning_signal & .data$bifpar_idx > max(transition_idx), .data$value, NA)
    ) %>%
    dplyr::mutate_at(c("ymin", "ymax", "warning_signal", "value_no_warning",
                       "value_warning_on_time", "value_warning_baseline", "value_warning_too_late"), ~suppressWarnings(as.numeric(as.character(.)))) %>%
                         mutate(baseline_start_idx = first(baseline_idx),
                                baseline_end_idx = last(baseline_idx),
                                transition_start_idx = first(transition_idx),
                                transition_end_idx = last(transition_idx))

  return(winEWS_df)
}


#' Detect warnings based on a critical cut-off value
#'
#' @inheritParams get_warnings_per_sigma
#' @param sigma_crit Critical cut-off value
#' @return Dataframe with warning indices and number of warning patches
#' @export
#'
#' @examples
get_warnings_one_sigma <- function(bifpar_idx, z_score, sigma_crit, nr_consecutive_warnings){
  # Find all warnings
  idx_all_warnings = which(abs(z_score) >= sigma_crit)

  # Find consecutive warnings
  list_conseq_seq = split(idx_all_warnings, cumsum(c(1, diff(idx_all_warnings) != 1)))
  conseq_seq = list_conseq_seq[unlist(purrr::map(list_conseq_seq, length)) >= nr_consecutive_warnings]
  idxs_warnings = unlist(unname(conseq_seq))
  nr_patches = ifelse(length(conseq_seq) == 0, numeric(0), length(conseq_seq))

  return(data.frame(bifpar_idx = bifpar_idx[idxs_warnings], z_score = z_score[idxs_warnings]) %>% dplyr::mutate(nr_patches = nr_patches))
  # return(list(idxs_warnings = idxs_warnings, nr_patches = nr_patches))
}

#' Get z-scores of EWS
#'
#' @inheritParams get_warnings
#'
#' @return Dataframe with the mean and standard deviation per EWS in the baseline period
#' @export
#' @importFrom dplyr arrange group_by filter summarise mutate mutate_at ungroup .data
#' @examples
get_zscore_EWS <- function(split_df_EWS, baseline_idx){
  # Compute baseline mean and standard deviation
  EWS_df_CI = split_df_EWS %>% arrange(.data$bifpar_idx) %>% group_by(.data$metric) %>%
    filter(.data$bifpar_idx %in% baseline_idx) %>%
    summarise(mean_w0 = mean(.data$value), sd_w0 = stats::sd(.data$value),
              quantile_000 = as.numeric(stats::quantile(.data$value, 0, na.rm=TRUE)),
              quantile_100 = as.numeric(stats::quantile(.data$value, 1, na.rm=TRUE)),
              .groups = 'drop')


  EWS_z_scores = merge(split_df_EWS, EWS_df_CI) %>% arrange(.data$bifpar_idx) %>%
    group_by(.data$metric) %>%
    # Compute z-scores (normalized): Convert z-score so that it can be compared to any sigma (i.e. you don't have to rerun the procedure for every sigma) using the following formula:
    # current > mu_baseline + (sigma_baseline * sigma_crit)
    # (current - mu_baseline) > (sigma_baseline * sigma_crit)
    # ((current - mu_baseline) / sigma_baseline) > (sigma_crit)
    # z = (current - mu_baseline) / sigma_baseline
    mutate(z_score = ((.data$value - .data$mean_w0) / .data$sd_w0)) %>% ungroup() %>%
    apply(2, unlist) %>% as.data.frame %>%
    mutate_at(c("bifpar_idx", "value", "z_score", "mean_w0", "sd_w0", "quantile_000", "quantile_100"), ~ as.numeric(as.character(.x))) %>%
    group_by(.data$metric) %>%
    mutate(bifpar_idx = round(as.numeric(as.character(.data$bifpar_idx)))) %>%
    ungroup()

 return(EWS_z_scores)
}


#' Get warnings per EWS metric
#'
#' @param split_df_EWS Dataframe of EWS split into lists, with one entry per value of the bifurcation parameter
#' @param baseline_idx Indices of baseline steps in bifurcation parameter
#' @param transition_idx Indices of transition steps in bifurcation parameter
#' @param sigma_crit_step Step size in sequence of critical values of sigma
#' @param nr_consecutive_warnings Number of consecutive warnings to look for
#' @param thresh_max_sigma_crit Maximum critical cut-off value to look for
#'
#' @return Dataframe with warnings per EWS metric
#' @importFrom dplyr arrange ungroup filter group_by mutate mutate_at summarise group_modify .data
#' @export
#'
#' @examples
get_warnings <- function(split_df_EWS, baseline_idx, transition_idx,
                         sigma_crit_step = .25,
                         thresh_max_sigma_crit = 500,
                         # sigmas_crit= seq(.25, 6, by = .25),
                         nr_consecutive_warnings = 1){
  # Compute baseline mean and standard deviation
  # EWS_df_CI = split_df_EWS %>% arrange(.data$bifpar_idx) %>% group_by(.data$metric) %>%
  #   filter(.data$bifpar_idx %in% baseline_idx) %>%
  #   summarise(mean_w0 = mean(.data$value), sd_w0 = stats::sd(.data$value),
  #                    quantile_000 = as.numeric(stats::quantile(.data$value, 0, na.rm=TRUE)),
  #                    quantile_100 = as.numeric(stats::quantile(.data$value, 1, na.rm=TRUE)),
  #                    .groups = 'drop')
  EWS_z_scores = get_zscore_EWS(split_df_EWS, baseline_idx)

  # # Compute z-scores - we're interested in absolute deviations, so either below or above mu + alpha_crit*sd
  # EWS_z_scores = merge(split_df_EWS, EWS_df_CI) %>% arrange(.data$bifpar_idx) %>%
  #   group_by(.data$metric) %>%
  #   # Compute z-scores (normalized): Convert z-score so that it can be compared to any sigma (i.e. you don't have to rerun the procedure for every sigma) using the following formula:
  #   # current > mu_baseline + (sigma_baseline * sigma_crit)
  #   # (current - mu_baseline) > (sigma_baseline * sigma_crit)
  #   # ((current - mu_baseline) / sigma_baseline) > (sigma_crit)
  #   # z = (current - mu_baseline) / sigma_baseline
  #   mutate(z_score = ((.data$value - .data$mean_w0) / .data$sd_w0)) %>% ungroup() %>%
  #   apply(2, unlist) %>% as.data.frame %>%
  #   mutate_at(c("bifpar_idx", "value", "z_score", "mean_w0", "sd_w0", "quantile_000", "quantile_100"), ~ as.numeric(as.character(.x))) %>%
  #   group_by(.data$metric) %>%
  #   mutate(bifpar_idx = round(as.numeric(as.character(.data$bifpar_idx)))) %>%
  #   ungroup()

  # Get warnings per critical sigma
  warning_df = EWS_z_scores %>%
    filter(.data$bifpar_idx %in% transition_idx) %>%
    group_by(.data$metric) %>%
    arrange(.data$bifpar_idx, .by_group = TRUE) %>%
    group_modify(~ get_warnings_per_sigma(y = .y,
                                          bifpar_idx = .x$bifpar_idx,
                                          z_score = .x$z_score,
                                          sigma_crit_step = sigma_crit_step,
                                          thresh_max_sigma_crit = thresh_max_sigma_crit,
                                          # sigmas_crit= sigmas_crit,
                                          nr_consecutive_warnings = nr_consecutive_warnings)) %>% ungroup() %>%
    rowwise() %>%
    dplyr::mutate(warning_signal = sum(.data$nr_warnings != 0)) %>% ungroup()
  # dplyr::mutate(warning_signal = sum(.data$nr_warnings != 0), no_warning_signal = sum(.data$nr_warnings == 0)) %>% ungroup()

  return(list(winEWS_df = EWS_z_scores,
              warning_df = warning_df))
}


#' Convert warnings to Receiver Operating Curve (ROC)
#'
#' @param EWS_warnings Dataframe with warnings
#' @param grouping_vars Names of grouping variables
#'
#' @return Dataframe with true positive rate, true negative rate, false positive rate, and false negative rate per critical value
#' @export
#' @importFrom dplyr mutate mutate_at select summarise rowwise ungroup group_by_at .data
#'
#' @examples
warnings_to_ROC <- function(EWS_warnings, grouping_vars){

  default_grouping_vars <- c("sigma_crit", "metric")
  grouping_vars = c(default_grouping_vars, grouping_vars)

  # Make sure the same sigma_crit is present for each model
  complete_sigma_crit = EWS_warnings %>% group_by_at(setdiff(grouping_vars, "sigma_crit")) %>%
  tidyr::complete(.data$noise_iter, .data$data_idx, .data$trans_or_null, .data$sigma_crit,
                  fill = list( first_warning_bifpar_idx=NA,
                               score = NA, nr_warnings = 0,
                               nr_patches = NA, warning_signal = 0)) %>% ungroup()

  # complete_sigma_crit %>% group_by(trans_or_null, metric, regime_switch) %>%
  # summarise(n = n(), .groups = 'drop') %>% as.data.frame() %>% head(n=100)
  # complete_sigma_crit %>% group_by_at(c(setdiff(grouping_vars, "sigma_crit"), "trans_or_null")) %>%
  # summarise(n = n(), .groups = 'drop') %>% as.data.frame() %>% head(n=100)


  # Add number of true positives, false negatives, true negatives, and false positives
  complete_sigma_crit_ = complete_sigma_crit %>%
   # rowwise() %>%
    # mutate(warning_signal = sum(.data$nr_warnings != 0), no_warning_signal = sum(.data$nr_warnings == 0)) %>%
    # ungroup() %>% rowwise() %>%
    mutate(nr_tp = ifelse(.data$trans_or_null == "transition", (.data$warning_signal == 1)*1,
                                 ifelse(.data$trans_or_null == "null", NA, "?")),
                  nr_fn = ifelse(.data$trans_or_null == "transition", (.data$warning_signal == 0)*1,
                                 ifelse(.data$trans_or_null == "null", NA, "?")),
                  nr_tn = ifelse(.data$trans_or_null == "transition", NA,
                                 ifelse(.data$trans_or_null == "null",(.data$warning_signal == 0)*1, "?")),
                  nr_fp = ifelse(.data$trans_or_null == "transition", NA,
                                 ifelse(.data$trans_or_null == "null", (.data$warning_signal == 1)*1, "?"))
    )

  # Compute false/true positive/negative rate
  EWS_warnings_ROC = complete_sigma_crit_ %>%
    mutate_at(c("sigma_crit", "nr_tp", "nr_fp", "nr_tn", "nr_fn"), ~as.numeric(as.character(.))) %>%
    group_by_at(grouping_vars) %>%
    summarise(acc = (sum(.data$nr_tp, na.rm = TRUE) + sum(.data$nr_tn, na.rm = TRUE)) / (sum(.data$nr_tp, na.rm = TRUE) + sum(.data$nr_tn, na.rm = TRUE) + sum(.data$nr_fp, na.rm = TRUE) + sum(.data$nr_fn, na.rm = TRUE)),
                     sum_tp = sum(.data$nr_tp, na.rm = TRUE),
                     sum_fp = sum(.data$nr_fp, na.rm = TRUE),
                     sum_tn = sum(.data$nr_tn, na.rm = TRUE),
                     sum_fn = sum(.data$nr_fn, na.rm = TRUE),
                     fnr = ifelse(sum(.data$nr_fn, na.rm = TRUE) == 0, 0, sum(.data$nr_fn, na.rm = TRUE) / (sum(.data$nr_fn, na.rm = TRUE) + sum(.data$nr_tp, na.rm = TRUE))), # miss rate
                     tnr = ifelse(sum(.data$nr_tn, na.rm = TRUE) == 0, 0, sum(.data$nr_tn, na.rm = TRUE) / (sum(.data$nr_tn, na.rm = TRUE) + sum(.data$nr_fp, na.rm = TRUE))), # specificity
                     fpr = ifelse(sum(.data$nr_fp, na.rm = TRUE) == 0, 0, sum(.data$nr_fp, na.rm = TRUE) / (sum(.data$nr_fp, na.rm = TRUE) + sum(.data$nr_tn, na.rm = TRUE))), # false alarm
                     tpr = ifelse(sum(.data$nr_tp, na.rm = TRUE) == 0, 0, sum(.data$nr_tp, na.rm = TRUE) / (sum(.data$nr_tp, na.rm = TRUE) + sum(.data$nr_fn, na.rm = TRUE))),  # sensitivity
                     .groups = 'drop')

  return(EWS_warnings_ROC)
}



#' Integrate area under receiver operator curve (ROC) to Area Under the Curve (AUC)
#'
#' @param EWS_warnings_ROC Dataframe with true positive rate, true negative rate, false positive rate, and false negative rate per critical value
#' @param grouping_vars Names of grouping variables
#'
#' @return AUC dataframe
#' @export
#' @importFrom dplyr summarise group_by_at .data
#'
#' @examples
ROC_to_AUC <- function(EWS_warnings_ROC, grouping_vars){

  default_grouping_vars <- c("metric")
  grouping_vars = c(default_grouping_vars, grouping_vars)

  # Add number of true positives, false negatives, true negatives, and false positives
  EWS_warnings_AUC = EWS_warnings_ROC %>%
    group_by_at(grouping_vars) %>%
    summarise(AUC = get_AUC(.data$fpr, .data$tpr) %>% round(6), .groups = 'drop')

  return(EWS_warnings_AUC)
}


#' Label AUC
#'
#' @param AUC Vector with Area Under the Curve
#' @param nbins Number of bins to divide AUC into
#'
#' @return Labelled AUC
#' @export
#'
#' @examples
get_AUC_class <- function(AUC, nbins = 10){
  breaks_AUC = seq(0, 1, length.out = nbins + 1)
  labels_AUC = plyr::laply(1:nbins,
                           function(i){
                             if (i==nbins){
                               comp_oper = "\\leq"
                             } else {
                               comp_oper = "<"
                             }
                             latex2exp::TeX(sprintf("$%.2f \\leq AUC %s %.2f$", breaks_AUC[i], comp_oper, breaks_AUC[i+1]), output = 'character')})

  return(factor(cut(AUC,
      breaks = breaks_AUC,
      labels = labels_AUC,
      include.lowest=TRUE, # Don't exclude AUC = 1
      right = FALSE # Brackets NOT closed on the right but closed on the left
      ), levels = labels_AUC))
}


#' Get Area Under the Curve (AUC)
#'
#' @param fpr False Positive Rate
#' @param tpr True Positive Rate
#'
#' @return AUC
#' @export
#'
#' @examples
get_AUC <- function(fpr, tpr){

  xy = cbind(fpr = fpr, tpr = tpr) %>%
    as.data.frame() %>%
    dplyr::distinct() %>% dplyr::arrange(fpr, tpr)
  x <- xy$fpr
  y <- xy$tpr

  # Only true positives
  if (all(x == 0) & all(y == 1)){
    return(1)
    # Only false negatives
  } else if (all(y == 0) & all(x == 1)){
    return(0)
  } else if (all(x == 0) & all(y == 0)){
    return(0)
  } else if (all(x == 1) & all(y == 1)){
    return(0.5)
  } else {
    return(pracma::trapz(x,y))
  }
}

#' Compute best threshold using Youden's J statistic
#'
#' @param x Dataframe with true/false positive/negative rates for each critical cut-off value
#'
#' @return Dataframe with only the row corresponding to the maximal Youden's J statistic
#' @export
#'
#' @examples
YoudensJ = function(x){
  x %>% dplyr::arrange(.data$sigma_crit) %>%
    dplyr::summarise(
      idx_sigma_crit = which.max(.data$tpr - .data$fpr),
      sigma_crit = .data$sigma_crit[.data$idx_sigma_crit],
      tpr = .data$tpr[.data$idx_sigma_crit],
      fpr = .data$fpr[.data$idx_sigma_crit],
      fnr = .data$fnr[.data$idx_sigma_crit],
      tnr = .data$tnr[.data$idx_sigma_crit],
      sum_tp = .data$sum_tp[.data$idx_sigma_crit],
      sum_tn = .data$sum_tn[.data$idx_sigma_crit],
      sum_fp = .data$sum_fp[.data$idx_sigma_crit],
      sum_fn = .data$sum_fn[.data$idx_sigma_crit]
    ) %>% return()
}

## Multivariate EWS Metrics
#' Eigenvalue
#'
#' @param x Matrix
#'
#' @return Eigenvalue
#' @export
#'
#' @examples
eigenvalue <- function(x) {
  eigen(stats::cov(x))$values[1]
}

#' Mean absolute cross-correlation
#'
#' @param x Matrix
#'
#' @return Mean absolute cross-correlation
#' @export
#'
#' @examples
get_conn <- function(x) {
  mean(abs(stats::cor(x)[upper.tri(stats::cor(x))]))
}


#' Spatial variance
#'
#' @param x Dataframe or matrix
#'
#' @return Spatial variance
#' @export
#'
#' @examples
spatial_variance <- function(x) {
  x <- as.matrix(x)
  return(1/(ncol(x)*nrow(x)) * sum((x - mean(x))**2))
}


#' Spatial skewness
#'
#' @param x Dataframe or matrix
#'
#' @return Spatial skewness
#' @export
#'
#' @examples
spatial_skewness <- function(x) {
  x <- as.matrix(x)
  sigma3 = spatial_variance(x)**(1.5)
  return(1/(ncol(x)*nrow(x)*sigma3) * sum((x - mean(x))**3) )
}


#' Spatial kurtosis
#'
#' @param x Dataframe or matrix
#'
#' @return Spatial kurtosis
#' @export
#'
#' @examples
spatial_kurtosis <- function(x) {
  x <- as.matrix(x)
  sigma4 = spatial_variance(x)**(2)
  return(1/(ncol(x)*nrow(x)*sigma4) * sum((x - mean(x))**4) )
}




## Univariate EWS Metrics

#' Coefficient of Variation (COV)
#'
#' @param x Timeseries
#'
#' @return COV
#' @export
#'
#' @examples
get_COV <- function(x){
  return(mean(x, na.rm = TRUE) / stats::sd(x, na.rm=TRUE))
}

#' Skewness
#'
#' @param x Vector
#'
#' @return Skewness
#' @export
#'
#' @examples
skewness <- function(x){moments::skewness(x)}

#' Kurtosis
#'
#' @param x Vector
#'
#' @return Kurtosis
#' @export
#'
#' @examples
kurtosis <- function(x){moments::kurtosis(x)}



#' Lag-1 Autocorrelation
#'
#' @param x Vector
#'
#' @return Lag-1 Autocorrelation
#' @export
#'
#' @examples
get_autocorr <- function(x) {
  stats::acf(x, plot = FALSE, lag = 1)$acf[[2]]
}




#' Compute Maximum Spectral Density (Bury, 2021)
#'
#' @param x Signal
#' @param fs Sampling frequency
#' @param nr_timesteps Number of timesteps
#'
#' @return Spectral EWS from Welch's PSD
#' @export
#'
#' @examples
get_Smax = function(x, fs, nr_timesteps){
  x = as.matrix(x)
  if (any(apply(x, 2, stats::sd) == 0)){
    Smax_df = matrix(0, ncol = ncol(x), nrow = 1) %>% magrittr::set_colnames(colnames(x)) %>% as.data.frame()

    # spectral_exp = 0
    # spec_ratio_LF_HF = 0
    # spec_ratio_LF_HF_perc = 0
  } else {
    lx <- fs * nr_timesteps
    pw <- gsignal::pwelch(x, window = lx, fs = fs,
                          # Remove mean so that the spectral peak isn't at zero
                          detrend = c("long-mean", "short-mean", "long-linear", "short-linear", "none")[1],
                          range = "half")
    Smax_df = as.data.frame(t(apply(pw$spec, 2, max)))
  }
  return(Smax_df)
}


#' Hurst exponent
#'
#'
#' @inheritParams get_Smax
#' @param scaleMin Minimum scale (in data points) to use for log-log regression
#' @param scaleMax Maximum scale (in data points) to use for log-log regression
#' @param polyOrderSegment The DFA order, the order of polynomial trend to remove from the bin
#'
#' @return Fitted linear slope of scale vs. detrended fluctuation (log-log)
#' @export
#'
#' @examples
get_Hurst_exp <- function(x, fs, nr_timesteps, scaleMin = 10, scaleMax = 100,
                          polyOrderSegment = 2){

  if (identical(all.equal(stats::sd(x), 0), TRUE)){
    Hurst_exp = 0
  } else {
  DFA = casnet::fd_dfa(
    stats::ts(x, frequency = fs, end = nr_timesteps),
    removeTrend = c("no", "poly", "adaptive", "bridge")[1], # Prettyman (2020)
    # polyOrder = 2, # Prettyman (2020)
    removeTrendSegment = c("no", "poly", "adaptive", "bridge")[2], # Prettyman (2020)
    polyOrderSegment = polyOrderSegment, # Prettyman (2020)
    scaleMin = scaleMin, # Prettyman (2020)
    scaleMax = scaleMax, # Prettyman (2020)
    fs = fs, doPlot = FALSE, silent = TRUE
  )
  Hurst_exp = unname(DFA$fullRange$sap)
  }
  return(Hurst_exp)
}



#' Compute Spectral Exponent (Prettyman, 2020; Wijnants, 2013)
#'
#' @inheritParams get_Smax
#' @param method Method of selecting frequency range in which to compute spectral exponent
#'
#' @return Slope of power spectral density
#' @export
#'
#' @examples
get_spectral_exp <- function(x, fs, nr_timesteps,
                             method = c("Prettyman2020", "Wijnants2013")[2]
                           ){

  if (identical(all.equal(stats::sd(x), 0), TRUE)){
    pse_value = 0
  } else {
    N = length(x)
    xdft = stats::fft(x)
    xdft = xdft[1:floor(N/2)+1] # Only keep positive frequencies
    psdx = (1/(fs*N)) * abs(xdft)**2
    psdx[2:(length(psdx)-1)] = 2*psdx[2:(length(psdx)-1)] # Don't multiply zero frequency by 2
    freq = seq(0, fs/2, length.out=length(psdx)) # Up until Nyquist frequency = fs / 2
    logf = log10(freq)
    logp = log10(psdx)

    # Select section of frequency range
    if (method == "Wijnants2013"){
      idx_min = 2 # Skip freq = 0
      idx_max = 51
      if (idx_max > length(freq)){
        message("Timeseries is too short to fit the slope of the spectral exponent over the first 50 frequencies (method = Wijnants2013)!")
        return()
      }
    } else if (method == "Prettyman2020"){
      f_min = 10**(-2)
      f_max = 10**(-1)
      if (f_max > (fs/2)){
        message("Error: The highest frequency f_max = %.4f you can estimate with a sampling frequency of %.4f is fs/2 = %.4f.", f_max, fs, fs/2)
        return()
      } else {
        idx_min = which.min(abs(freq - f_min))
        idx_max = which.min(abs(freq - f_max))
      }
    }

    pfit = stats::lm(logp[idx_min:idx_max] ~ logf[idx_min:idx_max], # Fit linear regression
              na.action=stats::na.omit)
    pse_value = -unname(pfit$coefficients[2])
  }
  return(pse_value)
}


#' Compute Spectral Ratio (Biggs, 2009)
#'
#' @inheritParams get_Smax
#' @inheritParams get_spectral_exp
#' @param f_min_to_f_max Low and high frequency to compare
#' @param n.freq Number of points to estimate frequency at
#'
#' @return Spectral ratio of spectral density estimated at specified low frequency to spectral density at specified high frequency
#' @export
#'
#' @examples
get_spectral_ratio <- function(x, fs, nr_timesteps,
                               f_min_to_f_max =  c(.05, .5),
                               n.freq = 500){

  if (f_min_to_f_max[2] > (fs/2)){
    message("Error: The highest frequency f_max = %.4f you can estimate with a sampling frequency of %.4f is fs/2 = %.4f.", f_min_to_f_max[2], fs, fs/2)
    return()
  }

  col_names  = paste0("LF", f_min_to_f_max[1], "_HF", f_min_to_f_max[2])
  if (identical(all.equal(stats::sd(x), 0), TRUE)){
    return(matrix(0, nrow = 1, ncol = length(f_min_to_f_max)) %>% magrittr::set_colnames(col_names) %>% as.data.frame() )

  } else {

    ARSPEC=stats::spec.ar(
      stats::ts(x, frequency = fs, end = nr_timesteps),
      n.freq=n.freq,
      # order=1, # Don't specify order to get order as decided by AIC
      plot=FALSE
    )
    freq = ARSPEC$freq[-1]
    spec = ARSPEC$spec[-1]

    # Spectral ratio dataframe
      f_min = f_min_to_f_max[1]
      f_max = f_min_to_f_max[2]
      idx_min = which.min(abs(freq - f_min))
      idx_max = which.min(abs(freq - f_max))
      spectral_ratio = spec[idx_min] / spec[idx_max]
      spectral_ratio_df =  matrix(spectral_ratio) %>% magrittr::set_colnames(col_names) %>% as.data.frame()

    return(spectral_ratio_df)

  }
}
