#' Generate timeseries with changing bifurcation parameter
#'
#' @param model ODE in the style of deSolve
#' @param model_pars Fixed model parameters
#' @param bifpar_list List with changing bifurcation parameter values
#' @param X0 Initial state
#' @param X_names Names of variables in model
#' @param seed_nr Seed number
#' @param timestep Step size in time
#' @param nr_timesteps Number of time steps
#' @param deSolve_method Method of generating ODE passed to deSolve::ode
#' @param stopifregime End generating timeseries if this function is satisfied
#' @param do_downsample Reduce dataframe size by downsampling?
#' @param downsample_pars List of parameters for downsampling
#' @param silent Don't output progress if TRUE
#' @param max_iter Maximum number of tries to restart from a new initial condition and not hit undesirable regime
#' @param seed_change Function to update seed number when restarting after hitting undesirable regime
#'
#' @return List with timeseries matrix and parameters used to generate the timeseries
#' @export
#'
#' @examples
bifurcation_ts <- function(model, model_pars, bifpar_list,
                           X0 = c(),
                           X_names = names(X0),
                           seed_nr = 123,
                           timestep = .01, nr_timesteps = 100,
                           deSolve_method = c("lsoda", "euler", "rk4")[1],
                           # stopifregime = function(out){FALSE},
                           stopifregime = function(out){any(apply(out, 2, is.infinite)) | any(apply(out, 2, is.nan))},

                           do_downsample = TRUE,
                           downsample_pars = list(type = c("average", "one_sample")[1],
                                                  win_size = 50,
                                                  which_X = c(50, "first", "middle", "last", "random")[1],
                                                  seed = 123),
                           silent = FALSE,
                           max_iter = 1,
                           seed_change = function(seed_nr){seed_nr + 1}){
  if (rlang::is_empty(X0) & rlang::is_empty(X_names)){
    message("Specify either a named initial condition X0 (e.g. X0 = (X1 = .1, X2 = .3)) or the names of your variables to set a random initial condition using X_names (e.g. X_names = c('X1', 'X2')).")
    return()
  }

  # Initialize
  success = FALSE; iter = 1
  times = seq(0, nr_timesteps, by = timestep)
  min_t = times[1]
  bifpar_idxs = seq.int(length(bifpar_list))
  X0s <- matrix(nrow = length(bifpar_list), ncol = length(X_names)+1) %>% magrittr::set_colnames(c("bifpar_idx", X_names))

  while (!success & (iter <= max_iter)){
    # Initialize
    if (rlang::is_empty(X0)){
      set.seed(seed_nr)
      X0      <- stats::runif(length(X_names)) %>% stats::setNames(X_names)
      if (!silent){
        print(X0)
      }
    }
    tmps <- list() #as.list(rep("", length(bifpar_idxs))) #list()

    # Generate data for each bifurcation parameter
    for (bifpar_idx in bifpar_idxs){
      if (!silent){
        print(sprintf("%s Generating data for bifurcation index %d/%d", Sys.time(), bifpar_idx, length(bifpar_idxs)))
      }
      # Adjust bifurcation parameter
      bifpar = bifpar_list[[bifpar_idx]]
      model_pars = utils::modifyList(model_pars, bifpar)

      # Save initial condition
      X0s[bifpar_idx,] <- c(bifpar_idx, X0)

      # Generate data
      out <- deSolve::ode(y = X0, times = times + min_t, func = model,
                          parms = model_pars,
                          method = deSolve_method) %>% cbind(bifpar_idx=bifpar_idx)
      # head(out)
      # tail(out)
      # nrow(out)

      # Overwrite initial state
      X0 <- out[nrow(out),names(X0)]
      min_t = out[nrow(out),c("time")]

      # Remove last time point
      out <- out[-nrow(out),]
      # print(X0)

      # Stop if condition is met
      if (stopifregime(out) & (bifpar_idx > bifpar_idxs[1])){
        message("Hit undesirable regime! Deleting generated timeseries and starting from new initial condition")
        break
      }

      if (do_downsample){
        out = downsample(out, X_names, type = downsample_pars$type,
                         win_size = downsample_pars$win_size,
                         which_X = downsample_pars$which_X,
                         seed_nr = downsample_pars$seed_nr)
      }
      # Save intermediate result for efficiency
      tmp <- tempfile(fileext = ".RDS")
      saveRDS(out, tmp)
      tmps[bifpar_idx] <- tmp
      rm(out)
    }

    if (length(tmps) == length(bifpar_idxs)){
      # Collect all generate data and remove temporary files
      OUT = do.call(rbind, lapply(bifpar_idxs,
                                  function(bifpar_idx){
                                    tmp <- tmps[[bifpar_idx]]
                                    if (file.exists(tmp)){
                                    out <- readRDS(tmp)
                                    unlink(tmp)
                                    return(out)

                                    # return(cbind(out, bifpar_idx))
                                    } else {
                                      return(NULL)
                                    }
                                  }))
      success = TRUE
    } else {
      # Start again
      X0 = c()
      seed_nr = seed_change(seed_nr)
      iter = iter + 1
      print(sprintf("Iteration %d, seed number: %.8f", iter, seed_nr))
    }
  }

  # OUT$time_idx = 1:nrow(OUT)
  # OUT %>% head

  return(list(df = as.data.frame(cbind(OUT, time_idx=1:nrow(OUT))),
              X_names = X_names,
              X0s = X0s,
              timestep = timestep,
              nr_timesteps = nr_timesteps,
              model_pars = model_pars,
              seed_nr = seed_nr,
              stopifregime = stopifregime,
              deSolve_method = deSolve_method,
              bifpar_list = bifpar_list,
              do_downsample = do_downsample,
              downsample_pars = downsample_pars))
}


#' Find both local peaks and troughs of timeseries
#'
#' @param df Dataframe
#' @param X_names Column names of df to find peaks and troughs for
#'
#' @return Dataframe with peaks and troughs per variable, including index of the beginning, height, and end of each peak
#' @export
#' @importFrom dplyr select mutate filter slice group_by ungroup all_of rename .data
#'
#' @examples
peaks_bifdiag <- function(df, X_names){

  peaks = do.call(rbind, lapply(X_names, function(i){
    # Find maximum and minimum peaks
    maxpeaks = pracma::findpeaks(df[,i])
    # Negative to find minimum peaks
    minpeaks = pracma::findpeaks(-df[,i])

    if (!rlang::is_empty(maxpeaks)){
      maxpeaks_df = cbind(df[maxpeaks[,2], c(i, setdiff(colnames(df), X_names))],
                          peak_idx = maxpeaks[,2],
                          begin_peak_idx = maxpeaks[,3], end_peak_idx = maxpeaks[,4]
      ) %>% as.data.frame() %>% mutate(minmax = "maxpeak")
    }
    if (!rlang::is_empty(minpeaks)){
      minpeaks_df = cbind(df[minpeaks[,2], c(i, setdiff(colnames(df), X_names))],
                          peak_idx = minpeaks[,2],
                          begin_peak_idx = minpeaks[,3], end_peak_idx = minpeaks[,4]
      ) %>% as.data.frame() %>% mutate(minmax = "minpeak")
    }
    # Make sure each bifurcation parameter has an entry
    peaks_df = rbind(maxpeaks_df, minpeaks_df)
    missing_bifpar_idx = setdiff(unique(df[,c("bifpar_idx")]), peaks_df$bifpar_idx)
    node_df = df[,c(i, setdiff(colnames(df), X_names))] %>% as.data.frame() %>%
      slice(1, .by = .data$bifpar_idx) %>% ungroup() %>%
      filter(.data$bifpar_idx %in% missing_bifpar_idx) %>%
      mutate(minmax="node",
                    peak_idx = .data$time_idx,
                    begin_peak_idx = .data$time_idx, end_peak_idx = .data$time_idx)
    peaks_df_complete = rbind(peaks_df, node_df) %>% mutate(variable = i) %>% rename("X" = all_of(i))

    return(peaks_df_complete)
  }))

  return(peaks)
}



#' Generate bifurcation parameter sequence
#'
#' @param bifpar_start Bifurcation parameter value at start
#' @param bifpar_end Bifurcation parameter value at end; NA for null models
#' @param pre_steps Number of pre-steps, useful if transition did not start when planned
#' @param baseline_steps Number of baseline steps where the control parameter doesn't change
#' @param transition_steps Number of transition steps where the control parameter changes for every index in a step-wise manner
#' @param post_steps Number of post-transition steps
#'
#' @return List of bifurcation parameter values
#' @export
#'
#' @examples
get_bifurcation_range <- function(bifpar_start, bifpar_end, pre_steps = 0,
                                  baseline_steps = 50, transition_steps = 50, post_steps = 0){

  if (is.na(bifpar_end)){
    # Null model
    s_seq = rep(bifpar_start, pre_steps + baseline_steps + transition_steps + post_steps)
  } else {
    s_seq = c(rep(bifpar_start, pre_steps + baseline_steps),
              seq(bifpar_start, bifpar_end, length.out = transition_steps),
              rep(bifpar_end, post_steps))
  }
  return(purrr::transpose(expand.grid(s = s_seq)))
}


#' Find consecutive values in sequence
#'
#' @param bifpar_idx_ Sequence of indices
#'
#' @importFrom dplyr .data arrange mutate rowwise
#' @return Start and end of each consecutive value sequence
#'
#' @examples
find_consec_seq = function(bifpar_idx_) {
  bifpar_idx = sort(unique(bifpar_idx_))
  list_conseq_seq = split(bifpar_idx, cumsum(c(1, diff(bifpar_idx) != 1)))
  df_conseq_seq = rbind(
    purrr::map(list_conseq_seq, min) %>% as.data.frame(),
    purrr::map(list_conseq_seq, max) %>% as.data.frame(),
    purrr::map(list_conseq_seq, length) %>% as.data.frame()
  ) %>% t() %>%
    magrittr::set_colnames(c("start_bifpar_idx", "end_bifpar_idx", "length_region")) %>% as.data.frame() %>%
    mutate(region_nr = 1:nrow(.), nr_regions = nrow(.)) %>% rowwise()

  return(df_conseq_seq %>% arrange(.data$start_bifpar_idx))
}



#' Find whether and when basin boundaries are hit
#'
#' @param peaks_df Dataframe with peaks
#' @param variable_name Column name in dataframe to assess
#' @param min_edge Minimum value reached
#' @param max_edge Maximum value reached
#'
#' @return Dataframe with bifurcation parameter value for which the system touches the basin boundaries as defined by min_edge and max_edge, if it hits these. Otherwise, dataframe filled with NA.
#' @importFrom dplyr select mutate filter slice group_by ungroup all_of rename row_number summarise pull arrange n .data
#'
#' @examples
find_basin_boundary <- function(peaks_df, variable_name = "X1", min_edge = 0, max_edge = 1){

  minmax_peaks_df = peaks_df %>% select(.data$bifpar_idx, .data$variable, .data$minmax, .data$time_idx, .data$X) %>%
    filter(.data$variable == variable_name) %>%
    group_by(.data$bifpar_idx, .data$variable, .data$minmax) %>%
    summarise(max = round(max(.data$X), 2), min = round(min(.data$X), 2), .groups='drop') %>%
    tidyr::pivot_wider(names_from = "minmax", values_from = c("max", "min")) %>%
    select(.data$bifpar_idx, .data$variable, .data$max_maxpeak,.data$min_minpeak) %>% ungroup %>%
    filter(.data$max_maxpeak == max_edge & .data$min_minpeak == min_edge) %>% arrange(.data$bifpar_idx) %>%
    filter(row_number()==1 | row_number()==n())
  # If the edges were not touched
  if (nrow(minmax_peaks_df) == 0){
    # minmax_peaks_df = data.frame(bifpar_idx = NA)
    return(data.frame())
  }

  basin_bound = data.frame(start_bifpar_idx = minmax_peaks_df %>% slice(1) %>% pull(.data$bifpar_idx),
                           end_bifpar_idx = minmax_peaks_df %>% slice(nrow(.)) %>% pull(.data$bifpar_idx),
                           regime = "Basin-Boundary"
  ) %>%
    mutate(length_region = .data$end_bifpar_idx - .data$start_bifpar_idx + 1,
           region_nr = ifelse(is.na(.data$start_bifpar_idx), NA, 1),
           nr_regions = ifelse(is.na(.data$start_bifpar_idx), NA, 1))

  return(basin_bound)
}


#' Find type of regime switch
#'
#' @param from_regime The regime from which the transition occurs
#' @param to_regime The regime which is transitioned to
#' @param X_names Names of variables
#'
#' @return Regime switch type
#' @importFrom dplyr select mutate mutate_at filter slice group_by ungroup all_of rename row_number summarise pull bind_rows .data
#'
#' @examples
get_regime_switch_type <- function(from_regime, to_regime, X_names){

  regime_df = bind_rows(from_regime, to_regime)
  # Regime switches involving only periodic regimes
  # if (!grepl("Chaotic", from_regime$regime, fixed = TRUE) & !grepl("Chaotic", to_regime$regime, fixed = TRUE)){
  if (all(!grepl("Chaotic", regime_df$regime)) & all(!grepl("None", regime_df$regime)) & !any(regime_df %>% select(all_of(X_names)) %>%is.na()) ){
  period_switch = regime_df %>%
    mutate_at(X_names, ~ as.numeric(stringr::str_replace(., "Period-", ""))) %>%
    select(all_of(X_names)) %>% as.matrix()

  period_switch_factor = period_switch[2,] / period_switch[1,]
  nr_doubling = sum(period_switch_factor == 2)
  nr_halving = sum(period_switch_factor == 0.5)
  nr_same = sum(period_switch_factor == 1)
  nr_increasing = sum(period_switch_factor > 1 & period_switch_factor != 2)
  nr_decreasing = sum(period_switch_factor < 1 & period_switch_factor != 0.5)
  N = length(X_names)

  if (nr_doubling == N){
   return("Period-Doubling")
  } else if (nr_halving == N){
    return("Period-Halving")
  } else if (nr_increasing == N){
    return("Period-Increasing")
  } else if (nr_decreasing == N){
    return("Period-Decreasing")
  } else if (nr_same == N){
    return("Same-Period")
  } else if ((nr_increasing > 0 | nr_doubling > 0) & (nr_decreasing > 0 | nr_halving > 0)){
    return("Both-Increasing-Decreasing-Period")
  } else if (nr_doubling == 1 & nr_same==(N-1)){
    return("OneVar-Period-Doubling")
  } else if (nr_halving == 1 & nr_same==(N-1)){
    return("OneVar-Period-Halving")
  } else if (nr_increasing == 1 & nr_same==(N-1)){
    return("OneVar-Period-Increasing")
  } else if (nr_decreasing == 1 & nr_same==(N-1)){
    return("OneVar-Period-Decreasing")
  } else if (nr_doubling > 1 & nr_same==(N-sum(nr_doubling))){
    return("MultVar-Period-Doubling")
  } else if (nr_halving > 1 & nr_same==(N-sum(nr_halving))){
    return("MultVar-Period-Halving")
  } else if (nr_increasing > 1 & nr_same==(N-sum(nr_increasing))){
    return("MultVar-Period-Increasing")
  } else if (nr_decreasing > 1 & nr_same==(N-sum(nr_decreasing))){
    return("MultVar-Period-Decreasing")
    } else if ((nr_decreasing + nr_halving) == N){
      return("MultVar-Period-Decreasing-or-Halving")
    } else if ((nr_increasing + nr_doubling) == N){
      return("MultVar-Period-Increasing-or-Doubling")
    } else {
    return("Unknown Periodic")
  }

  } else {
    # Regime switches involving chaotic behaviour
    broad_period_switch = regime_df %>%
      tidyr::gather(X, value, -setdiff(colnames(.), X_names)) %>%
      group_by(.data$regime) %>%
      summarise(nr_periods = length(unique(.data$value)), period = paste0(unique(.data$value), collapse = ",")) %>%
      mutate(broad_regime = ifelse(grepl("Chaotic or Transitioning", .data$regime), "Chaotic or Transitioning",
                                   ifelse(grepl("Basin-Boundary", .data$regime), "Basin-Boundary",
                                          ifelse(grepl("None", .data$regime), "None",
                           ifelse(.data$nr_periods == 1, .data$period, "Mixed-Periodic")))))
    return(sprintf("%s to %s",
            broad_period_switch[from_regime$regime == broad_period_switch$regime, ] %>% pull(.data$broad_regime),
                broad_period_switch[to_regime$regime == broad_period_switch$regime,] %>% pull(.data$broad_regime)))
  }
}


#' Find regime boundaries
#'
#' @param regimes Dataframe with periodicity per value of the bifurcation parameter
#' @inheritParams find_regimes
#' @inheritParams bifurcation_ts
#'
#' @importFrom dplyr arrange ungroup mutate .data group_by_at
#' @return Dataframe with regime boundaries
#'
#' @examples
find_regime_bounds <- function(regimes, min_length_regime, X_names){
  # Find regimes that satisfy a certain size
  regimes = regimes %>% ungroup() %>% arrange(.data$start_bifpar_idx)
  regime_idx = which(regimes$length_region >= min_length_regime)

  if (rlang::is_empty(regime_idx)){
    return(data.frame(regime1 = NA,
                      regime2 = NA,
                      regime1_start_idx = NA,
                      regime1_halfway_idx = NA,
                      regime1_end_idx = NA,
                      regime1_length = NA,
                      regime2_start_idx = NA,
                      regime2_end_idx = NA,
                      regime2_length = NA
    ))
  } else {
  regime_bounds_df = lapply(1:(length(regime_idx)), function(i){

    from_regime = regimes[regime_idx[i],] %>% mutate(regime = ifelse(is.na(.data$regime), "None", .data$regime))
    to_regime = regimes[regime_idx[i+1],] %>% mutate(regime = ifelse(is.na(.data$regime), "None", .data$regime))
    regime_switch_type = get_regime_switch_type(from_regime, to_regime, X_names)

    return(data.frame(regime_switch_type = regime_switch_type,
                      regime1 = from_regime$regime,
                      regime2 = to_regime$regime,
                      regime1_start_idx = from_regime$start_bifpar_idx,
                      regime1_halfway_idx = from_regime$start_bifpar_idx + ceiling((from_regime$end_bifpar_idx - from_regime$start_bifpar_idx) / 2),
                      regime1_end_idx = from_regime$end_bifpar_idx,
                      regime1_length = from_regime$length_region,
                      regime2_start_idx = to_regime$start_bifpar_idx,
                      regime2_end_idx = to_regime$end_bifpar_idx,
                      regime2_length = to_regime$length_region
    ))
  }) %>% do.call(rbind, .) %>% as.data.frame()

  return(regime_bounds_df)
  }
}


#' Smooth over exceptions in otherwise stable regime
#'
#' @param periods Dataframe with periodicity per bifurcation parameter
#' @inheritParams find_regimes
#'
#' @return Dataframe with smoothed periodicity
#' @importFrom dplyr arrange mutate pull lag rowwise .data group_by ungroup rename
#'
#' @examples
smooth_periods <- function(periods, nr_smooth, min_length_regime){

   # If a consistent regime shows one exception, smooth over.
  find_nr_surrounding_regimes <- function(x, min_length_regime, nr_smooth){
    # Split vector
    m <- zoo::rollapply(x, min_length_regime*2 + nr_smooth, by = 1, FUN = c)
    # Find number of unique periods surrounding nr_smooth, excluding nr_smooth. Look back min_length_regime + nr_smooth - 1 steps (ignoring nr_smooth steps) and forward min_length_regime steps
    c(rep(NA, min_length_regime + nr_smooth - 1), split(m, row(m)) %>% purrr::map(
      function(x){length(unique(x[-seq(min_length_regime + 1, min_length_regime + nr_smooth)]))}) %>%
        unname %>% unlist, rep(NA, min_length_regime) ) %>% return()
  }
  # If a consistent regime shows one exception, smooth over.
  periods_smooth = periods %>% group_by(.data$variable) %>% arrange(.data$bifpar_idx, .by_group = TRUE) %>%
    mutate(lag_p = lag(.data$period, n = nr_smooth),
           # Apply rolling function to check for each row whether the previous min_length_regime values were the same and the future min_length_regime values were the same. If they were, fill in leading value.
           embedded_in_same_regime =
             # zoo::rollapply(.data$period, min_length_regime*2 + nr_smooth, function(x){length(unique(x[-seq(min_length_regime + 1, min_length_regime + nr_smooth)]))}, by = 1, fill=NA, align = 'center')
             find_nr_surrounding_regimes(.data$period, min_length_regime, nr_smooth) == 1 ) %>%
    ungroup() %>% rowwise() %>%
    rename(period_unsmooth = .data$period) %>%
    mutate(
      # check_ = (.data$embedded_in_same_regime == TRUE) & (.data$period != .data$lag_p),
      # check1 = (.data$embedded_in_same_regime == TRUE),
      # check2 =  (.data$period != .data$lag_p),
      period = ifelse(is.na(.data$lag_p) | is.na(.data$embedded_in_same_regime), .data$period_unsmooth,
                      ifelse((.data$embedded_in_same_regime == TRUE) & (.data$period_unsmooth != .data$lag_p), .data$lag_p, .data$period_unsmooth))) %>% ungroup() %>%
    select(-c(.data$lag_p, .data$embedded_in_same_regime)) # %>% pull(.data$period_smooth) %>% as.logical() %>% which()

  return(periods_smooth)

}

#' Find regimes of dynamical system
#'
#' @param GLV List object generated by bifurcation_ts()
#' @param thresh_node Threshold under which timeseries is classified as node
#' @param thresh_coord_spread Threshold for distance in peak and trough coordinates that determines whether a cluster (i.e. period length) fits; if exceeded, denoted as chaotic or transitioning
#' @param thresh_peak_idx_spread Same as thresh_coord_spread but for peak and trough indices
#' @param min_length_regime Minimum number of consecutive steps in the bifurcation parameter that have the same periodicity to qualify as a regime
#' @param max_k Maximum cluster size to look for
#' @param nr_smooth Number of exceptions in a stable periodicity window to smooth over; nr_smooth = 0 means no smoothing
#' @param factor_k Weighting of period length k; heavier weight means shorter k is preferred; factor_k = 0 means the optimal period length is chosen based solely on minimum spread
#' @param method_best_k Method of choosing best period length
#'
#' @return List of dataframes with periodicity per variable, periodicity per bifurcation parameter value, regimes, and regime boundaries
#' @importFrom dplyr arrange group_modify ungroup filter select group_by mutate rename bind_rows all_of slice_tail .data
#' @export
#'
#' @examples
find_regimes <- function(GLV,
                         thresh_node = .1,
                         thresh_coord_spread = .025,
                         thresh_peak_idx_spread=2,
                         min_length_regime = 5,
                         nr_smooth = 0,
                         factor_k = 1,
                         method_best_k = c("log", "min_scaled")[1],
                         max_k = NULL){

  # max_k = NULL
  #  thresh_node = .1
  # thresh_coord_spread = .025
  # thresh_peak_idx_spread=2
  # nr_smooth=0
  # factor_k = .1
  # min_length_regime = 5

  # Get dataframe with peaks
  peaks_df = peaks_bifdiag(GLV$df, GLV$X_names)

  # For each value of the bifurcation parameter, find the period length k which has a minimum spread
  print("Finding best fitting period length for all bifurcation parameter values")
  k_spread = peaks_df %>%
    # filter(bifpar_idx > 10, bifpar_idx < 30) %>%
    group_by(.data$variable, .data$bifpar_idx) %>%
    arrange(.data$time_idx, .by_group=TRUE) %>%
    group_modify(~ find_spread_per_k(coord = .x$X, peak_idx = .x$peak_idx,
                                     max_k = max_k)) %>%
    ungroup()

  # Select best k
  best_k = k_spread %>%
    group_by(.data$variable, .data$bifpar_idx) %>%
    group_modify(~ choose_best_k(spread_df = .x, thresh_node = thresh_node,
                                 factor_k = factor_k,
                                 method_best_k = method_best_k)) %>%
    ungroup()

  period_per_var = best_k %>%
    # tidyr::unnest(names(.)) %>%
    # Our algorithm will *always* find a period length k with a minimum spread, also for chaotic data. Set a threshold that decides which spread is too large to be classified as a neat periodic sequence.
    mutate(period = ifelse(.data$max_spread_coord > thresh_coord_spread & .data$max_spread_peak_idx > thresh_peak_idx_spread,
                           "Chaotic or Transitioning", paste0("Period-", .data$k) )) %>%
    # mutate(period =
    #          # ifelse(.data$max_spread_coord > thresh_coord_spread & .data$max_spread_peak_idx > thresh_peak_idx_spread,
    #                        # "Chaotic or Transitioning",
    #                        paste0("Period-", .data$k)) %>%
    # Smooth over period abnormalities that form nr_smooth exceptions
    smooth_periods(nr_smooth = nr_smooth, min_length_regime=min_length_regime)

  periods = period_per_var %>% select(.data$bifpar_idx, .data$period, .data$variable) %>%
    group_by(.data$bifpar_idx, .data$period) %>%
    mutate(period_group = paste0(
      unique(.data$period),
      " (",
      paste0(sort(unique(.data$variable)), collapse = ","),
      ")"
    )
    ) %>% group_by(.data$bifpar_idx) %>%
    mutate(period_bifpar = paste0(sort(unique(.data$period_group)), collapse = ' AND ')
    ) %>% select(-.data$period_group) %>% ungroup() %>%
    tidyr::pivot_wider(names_from = .data$variable, values_from = .data$period) %>%
    arrange(.data$bifpar_idx)

  # Get basin boundaries (when system hits edges of basin)
  basin_bound = find_basin_boundary(peaks_df, variable_name = "X1", min_edge = 0, max_edge = 1)

  # Find regimes with any chaotic behaviour and regimes with periodic behaviour
  broad_regimes = periods %>%
    mutate(period_bifpar2 = ifelse(grepl("Chaotic or Transitioning", .data$period_bifpar, fixed = TRUE),
                                   "Chaotic or Transitioning", "Periodic")) %>%
    group_by(.data$period_bifpar2) %>%
    group_modify( ~ find_consec_seq(.x$bifpar_idx)) %>%
    arrange(.data$start_bifpar_idx) %>%
    ungroup() %>%
    rename(regime = .data$period_bifpar2)

  # Compile regimes
  regimes = periods %>%
    filter(!grepl("Chaotic or Transitioning", .data$period_bifpar, fixed = TRUE)) %>%
    # group_by_at(.data$period_bifpar, all_of(X_names)) %>%
    group_by_at(setdiff(colnames(.), "bifpar_idx")) %>%
    group_modify( ~ find_consec_seq(.x$bifpar_idx)) %>%
    arrange(.data$start_bifpar_idx) %>%
    ungroup() %>% rename(regime = .data$period_bifpar) %>%
    bind_rows(basin_bound) %>%
    bind_rows(broad_regimes %>% filter(.data$regime != "Periodic")) %>%
    arrange(.data$start_bifpar_idx)

  # Find regime boundaries
  regime_bounds_ = find_regime_bounds(regimes, min_length_regime = min_length_regime, X_names = GLV$X_names)
  # Add corresponding initial conditions - the timepoint right before the bifurcation parameter changed to the starting value of the first regime
  regime_bounds =  merge(regime_bounds_,
                         GLV$X0s  %>% as.data.frame() %>%
                           filter(.data$bifpar_idx %in% c(regime_bounds_$regime1_start_idx)) %>%
                           rename(regime1_start_idx = .data$bifpar_idx),
                         #              df %>%
                         # filter(.data$bifpar_idx %in% c(regime_bounds_$regime1_start_idx - 1)) %>%
                         #   slice_tail(n=1, by = .data$bifpar_idx) %>%
                         # mutate(regime1_start_idx= .data$bifpar_idx + 1) %>%
                         # select(.data$regime1_start_idx, all_of(X_names)),
                         all.x = TRUE)

  regime_list = utils::modifyList(GLV, list(df = NULL,
                                            peaks_df = peaks_df,
                                            period_per_var = period_per_var,
                                            periods = periods,
                                            regimes = regimes,
                                            broad_regimes = broad_regimes,
                                            regime_bounds = regime_bounds,
                                            thresh_coord_spread = thresh_coord_spread,
                                            thresh_peak_idx_spread=thresh_peak_idx_spread,
                                            min_length_regime = min_length_regime,
                                            max_k = max_k,
                                            nr_smooth = nr_smooth,
                                            factor_k = factor_k))

  return(regime_list)
}




#' Find distance between points in chosen temporal sequence
#'
#' @param vec Vector
#' @param cluster_idx Vector of cluster indices
#'
#' @return Mean, minimum and maximum distance between vector points when partitioned using the indices in ks
#'
#' @examples
max_dist <- function(vec, cluster_idx){
  # From the package sarafrr/basicClEval
  vec = as.matrix(vec)
  sizeCl <- summary(as.factor(cluster_idx))
  # Compute maximal distance within each cluster
  max_dist_per_k = unlist(lapply(split(vec, cluster_idx), function(x){max(stats::dist(x))}))

  # return(mean(max_dist_per_k))
  return(c(
    mean_spread = mean(max_dist_per_k),
    median_spread = stats::median(max_dist_per_k),
    # min_spread = min(max_dist_per_k),
    max_spread = max(max_dist_per_k)
  ))
}


#' Compute mean maximum distance per cluster
#'
#'
#' @param ks Vector of possible period lengths
#' @inheritParams find_spread_per_k
#'
#' @return Distance in peak and trough coordinates and indices per cluster partitioning
#'
#' @examples
find_dist_per_k <- function(ks, coord, peak_idx){
  # Compute mean maximum distance per cluster - essentially the spread of each cluster averaged per partitioning; for peak coordinates and peak indices
  coord_spread = lapply(ks, function(k){max_dist(vec = coord,
                                                 cluster_idx = rep(1:k, length(coord))[1:length(coord)])}) %>%
    do.call(rbind, .) %>%
    magrittr::set_colnames(paste0(colnames(.), "_coord")) %>%
    as.data.frame()


  peak_idx_spread = lapply(ks, function(k){max_dist(vec = diff(peak_idx),
                                                    cluster_idx = rep(1:k, length(coord))[2:length(coord)])}) %>%
    do.call(rbind, .) %>%
    magrittr::set_colnames(paste0(colnames(.), "_peak_idx")) %>%
    as.data.frame()
  spread_df = cbind(k = ks, cbind(coord_spread, peak_idx_spread))

  # plot(cluster_idx, coord)
  # plot(cluster_idx, peak_idx)
  return(as.data.frame(spread_df))
}

#' Scale vector between range (a,b)
#'
#' @param x Vector
#' @param a Minimum value
#' @param b Maximum value
#'
#' @return Scaled vector
#' @export
#'
#' @examples
scale_range = function(x, a = 0, b = 1){
  return((b-a) * ((x - min(x)) / (max(x) - min(x))) + a)
}

#' Reverse scaling of vector between range (a,b)
#'
#' @inheritParams scale_range
#' @param min_x Previous minimum of x
#' @param max_x Previous maximum of x
#'
#' @return Scaled vector
#' @export
#'
#' @examples
rev_scale_range <- function(x,a,b,min_x,max_x){return((x - a) / (b-a) * (max_x - min_x) + min_x) }


#' Find best fitting period length k
#'
#' @param coord Peak and trough coordinates
#' @param peak_idx Peak and trough indices
#' @param max_k Maximum cluster size to look for
#'
#' @return Dataframe with best fitting period length k and the corresponding minimum within-cluster distance and between-cluster distance for peak and trough coordinates and indices
#' @export
#'
#' @examples
find_spread_per_k <- function(coord, peak_idx, max_k = NULL){
  # Finding the best period is tricky: The global minimum might be a subharmonic, but the first local minimum might not be optimal for a more complex oscillation. If the timeseries is truly periodic, local minima in distance will occur for every subharmonic (e.g. period = 4, minima at period = c(4,8,12,...)).

  # Only consider partitionings with at least two repetitions
  if (is.null(max_k)){
    max_k = (floor(length(coord)/2)-1)
  } else {
    # Check whether max_k is possible to look for
    if (max_k > (floor(length(coord)/2)-1)){
      max_k = (floor(length(coord)/2)-1)
    }
  }
  ks = seq(1, max_k)

  spread_df = find_dist_per_k(ks, coord, peak_idx)
  return(cbind(as.data.frame(spread_df), max_k=max_k))
}


#' Choose best period length
#'
#' @param spread_df Output of find_spread_per_k()
#' @inheritParams find_regimes
#'
#' @return Dataframe containing one row with best period length
#' @export
#' @importFrom dplyr select mutate_at mutate .data
#'
#' @examples
choose_best_k <- function(spread_df, thresh_node, factor_k, method_best_k = c("log", "min_scaled")[2]){

  # Detect nodes for which k = 1 should be a good fit
  if (spread_df[spread_df$k == 1, c("max_spread_coord")] < thresh_node){
    idx_min=1
  } else {
  # Method 1: log
  if (method_best_k == "log"){
    idx_min = spread_df %>% ungroup() %>%
      select(c("k", "max_spread_coord", "median_spread_coord",
                      # "max_spread_peak_idx"
      )) %>%
      # Penalize for period length k twice
      # dplyr::mutate(k = k * factor_k) %>%
      # dplyr::mutate_all(~ (.+.00001)) %>%
      mutate_at(c("max_spread_coord", "median_spread_coord"), ~ scale_range(., a = 0.00001, b = 0.00001 + 1)) %>%
      # dplyr::mutate_at(c("max_spread_coord", "median_spread_coord"), ~ ifelse(. == 0, 0, scale_range(log(.), a = 0.00001, b = 1))) %>%
      # dplyr::mutate_at(c("max_spread_coord", "median_spread_coord", "max_spread_peak_idx"), ~ ifelse(. == 0, 0, log(. * k))) %>%
      # scale(center = F) %>%
      # apply(2,scale_range, a = 1, b = 2) %>%
      # as.data.frame %>%
      mutate(k = scale_range(.data$k, a = 0.00001, b = 0.00001+factor_k)) %>%
      # as.data.frame %>%
    apply(1, mean) %>% which.min()

  } else if (method_best_k == "min_scaled"){
    # Method 2
    # Divide by minimum spread - how much more spread does each k have as compared to the minimum?
    max_spread_coord_scaled=spread_df$max_spread_coord / min(spread_df$max_spread_coord)
    # Divide by the period length k - how much longer or shorter is each k compared to the k corresponding to the minimum spread?
    k_scaled=spread_df$k / which.min(spread_df$max_spread_coord)
    # Choose k that balances minimum spread with short period length. If factor_k is very high, it's near impossible for a high period length to be the optimal period - it needs to have an excellent fit compared to other period lengths.
    idx_min = which.min(scale_range(max_spread_coord_scaled) + scale_range(k_scaled, a = 0, b = factor_k))
  }
  }
  return(as.data.frame(spread_df[idx_min,]))
}
