#' Compute early warning signals (EWS)
#'
#' @param df Dataframe
#' @param uni_metrics List of univariate EWS
#' @param multi_metrics List of multivariate EWS
#' @param EWS_args List of parameters to pass to EWS functions, empty list if no parameters need to be passed
#'
#' @return Dataframe with early warning signals
#' @export
#'
#' @examples
run_EWS <- function(df, uni_metrics, multi_metrics, EWS_args = list()){

  # print(names(list(...)))

  uni_EWS = lapply(names(uni_metrics), function(j){
    unlist(lapply(1:ncol(df), function(i){
      do.call(uni_metrics[[j]], utils::modifyList(list(x = df[,i]), EWS_args[j]))
    }))
    # apply(df, 2, uni_metric)
  }) %>%
    do.call(rbind, .) %>%
    magrittr::set_rownames(NULL) %>%
    as.data.frame() %>% cbind(metric = names(uni_metrics))  %>% tidyr::pivot_longer(!metric) %>%
    dplyr::mutate(metric = paste0(metric, "_", name)) %>% select(-name)

  multi_EWS = lapply(names(multi_metrics), function(j){
    do.call(multi_metrics[[j]], utils::modifyList(list(x = df), EWS_args[j]))
    # multi_metric(df)
  }) %>%
    do.call(rbind, .) %>% magrittr::set_rownames(NULL) %>%
    as.data.frame() %>% cbind(metric = names(multi_metrics)) %>%
    dplyr::rename(value = V1)

  return(rbind(uni_EWS, multi_EWS))
}

#' Run EWS for all bifurcation parameter values
#'
#' @param df Dataframe
#' @param X_names Names of variables in model
#' @param uni_metrics List of univariate EWS
#' @param multi_metrics List of multivariate EWS
#' @param EWS_args List of parameters to pass to EWS functions, empty list if no parameters need to be passed
#'
#' @return Dataframe with EWS
#' @export
#'
#' @examples
run_bifEWS <- function(df, X_names, uni_metrics, multi_metrics,
                       EWS_args = list("RQA" = list(emDim = 1, emLag = 1, theiler = 1, distNorm = "max", targetValue = .05))){

  # Split dataframe per bifurcation parameter
  split_df = df %>% as.data.frame() %>% dplyr::group_by(bifpar_idx) %>% dplyr::group_split()
  split_df_EWS = split_df %>%
    lapply(function(df_){run_EWS(df_ %>% dplyr::arrange(time_idx) %>% select(all_of(X_names)) %>% as.matrix(), uni_metrics, multi_metrics, EWS_args = EWS_args)} %>% dplyr::mutate(bifpar_idx = unique(df_$bifpar_idx))) %>%
    do.call(rbind, .) %>% as.data.frame()
  head(split_df_EWS)

  return(split_df_EWS)
}


#' Get warnings per critical value of sigma
#'
#' @param y Vector with EWS
#' @param bifpar_idx Bifurcation parameter index
#' @param z_score Z-score
#' @param crit_values Sequence of critical values of sigma
#' @param nr_consecutive_warnings Number of consecutive warnings to look for
#'
#' @return Dataframe with warnings
#' @export
#'
#' @examples
get_warnings_per_sigma <- function(y, bifpar_idx, z_score, crit_values, nr_consecutive_warnings = 1){
  lapply(crit_values, function(crit_value){
    idx_warnings = which(abs(z_score) >= crit_value)

    list_conseq_seq = split(idx_warnings, cumsum(c(1, diff(idx_warnings) != 1)))
    conseq_seq = list_conseq_seq[unlist(purrr::map(list_conseq_seq, length)) >= nr_consecutive_warnings]
    nr_patches = length(conseq_seq)
    nr_warnings = purrr::map(conseq_seq, length) %>% unlist() %>% sum()

    if (rlang::is_empty(conseq_seq)){
      first_warning_bifpar_idx = NA
      score = NA
    } else {
      first_warning_bifpar_idx = bifpar_idx[conseq_seq[[1]][1]]
      score = z_score[conseq_seq[[1]][1]]
    }
    return(data.frame(crit_value = crit_value, first_warning_bifpar_idx = first_warning_bifpar_idx, score=score,
                      nr_warnings = nr_warnings, nr_patches=nr_patches))
  }) %>% do.call(rbind, .) %>% as.data.frame() %>% return()

}


#' Get warnings per EWS metric
#'
#' @param split_df_EWS Dataframe of EWS split into lists, with one entry per value of the bifurcation parameter
#' @param baseline_steps Number of baseline steps in bifurcation parameter
#' @param transition_steps Number of transition steps in bifurcation parameter
#' @param sigmas_crit Sequence of critical values of sigma
#' @param nr_consecutive_warnings Number of consecutive warnings to look for
#'
#' @return Dataframe with warnings per EWS metric
#' @export
#'
#' @examples
get_warnings <- function(split_df_EWS, baseline_steps, transition_steps, sigmas_crit = seq(.25, 6, by = .25), nr_consecutive_warnings = 1){
  # Compute baseline mean and standard deviation
  EWS_df_CI = split_df_EWS %>% dplyr::arrange(bifpar_idx) %>% dplyr::group_by(metric) %>%
    dplyr::filter(bifpar_idx <= baseline_steps) %>%
    dplyr::summarise(mean_w0 = mean(value), sd_w0 = sd(value),
                     quantile_000 = as.numeric(quantile(value, 0, na.rm=TRUE)),
                     quantile_100 = as.numeric(quantile(value, 1, na.rm=TRUE)),
                     .groups = 'drop')

  # Compute z-scores - we're interested in absolute deviations, so either below or above mu + alpha_crit*sd
  winEWS_df = merge(split_df_EWS, EWS_df_CI) %>% dplyr::arrange(bifpar_idx) %>%
    dplyr::group_by(metric) %>%
    # Compute z-scores (normalized): Convert z-score so that it can be compared to any sigma (i.e. you don't have to rerun the procedure for every sigma) using the following formula:
    # current > mu_baseline + (sigma_baseline * sigma_crit)
    # (current - mu_baseline) > (sigma_baseline * sigma_crit)
    # ((current - mu_baseline) / sigma_baseline) > (sigma_crit)
    # z = (current - mu_baseline) / sigma_baseline
    dplyr::mutate(z_score_sd = ((value - mean_w0) / sd_w0)) %>% ungroup() %>% apply(2, unlist) %>% as.data.frame %>%
    dplyr::mutate_at(c("bifpar_idx", "value", "z_score_sd", "mean_w0", "sd_w0", "quantile_000", "quantile_100"), ~ as.numeric(as.character(.x))) %>%
    dplyr::group_by(metric) %>%
    dplyr::mutate(bifpar_idx = round(as.numeric(as.character(bifpar_idx)))) %>%
    ungroup()

  # Get warnings per critical sigma
  warning_df = winEWS_df %>%
    dplyr::filter(bifpar_idx > baseline_steps, bifpar_idx <= (baseline_steps + transition_steps)) %>%
    dplyr::group_by(metric) %>%
    group_modify(~ get_warnings_per_sigma(y = .y, bifpar_idx = .x$bifpar_idx, z_score = .x$z_score_sd, crit_values= sigmas_crit, nr_consecutive_warnings = nr_consecutive_warnings)) %>% ungroup

  return(list(winEWS_df = winEWS_df,
              warning_df = warning_df))
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
    distinct() %>% arrange(fpr, tpr)
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
  eigen(cov(x))$values[1]
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
  mean(abs(cor(x)[upper.tri(cor(x))]))
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


#' Recurrence Quantification Analysis
#'
#' @param x
#' @param emDim
#' @param emLag
#' @param theiler
#' @param distNorm
#' @param targetValue
#'
#' @return RQA measures
#' @export
#'
#' @examples
runRQA <- function(x, emDim = 1, emLag = 1, theiler = 1, distNorm = "max", targetValue = .05){
  RM <- casnet::rp(
    x,
    emDim   = emDim,
    emLag   = emLag,
    emRad   = NA, #pars$emRad,
    theiler = theiler,
    method = distNorm,
    targetValue    = targetValue
  )
  return(casnet::rp_measures(RM))
}


## Univariate EWS Metrics

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
  acf(x, plot = FALSE, lag = 1)$acf[[2]]
}

