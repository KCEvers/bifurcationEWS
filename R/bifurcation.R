#' Generate timeseries with changing bifurcation parameter
#'
#' @param model ODE in the style of deSolve
#' @param model_pars Fixed model parameters
#' @param bifpar_list List with changing bifurcation parameter values
#' @param bifpar_pars List with parameters to generate bifpar_list, needs to include list(bifpar_start = .96, bifpar_end = .96), optionally specify list(pre_steps = 0, baseline_steps = 100, transition_steps = 0, post_steps = 0)
#' @param X0 Initial state
#' @param X_names Names of variables in model
#' @param seed_nr Seed number
#' @param timestep Step size in time
#' @param nr_timesteps Number of time steps
#' @param deSolve_method Method of generating ODE passed to deSolve::ode
#' @param stopifregime End generating timeseries if this function is satisfied
#' @param do_downsample Reduce dataframe size by downsampling?
#' @param fs Sampling frequency for downsampling if do_downsample = TRUE
#' @param which_X Time point(s) in downsampling window to select if do_downsample = TRUE
#' @param silent Don't output progress if TRUE
#' @param max_iter Maximum number of tries to restart from a new initial condition and not hit undesirable regime
#' @param seed_change Function to update seed number when restarting after hitting undesirable regime
#'
#' @return List with timeseries matrix and parameters used to generate the timeseries
#' @export
#'
#' @examples
bifurcation_ts <- function(model, model_pars, bifpar_list = NULL, bifpar_pars = NULL,
                           X0 = c(),
                           X_names = names(X0),
                           seed_nr = 123,
                           timestep = .01, nr_timesteps = 100,
                           deSolve_method = c("lsoda", "euler", "rk4")[1],
                           # stopifregime = function(out){FALSE},
                           stopifregime = function(out){any(apply(out, 2, is.infinite)) | any(apply(out, 2, is.nan))},

                           do_downsample = TRUE,
                           fs = 50,
                           which_X = "all",
                           silent = FALSE,
                           max_iter = 1,
                           seed_change = function(seed_nr){seed_nr + 1}){
  if (rlang::is_empty(X0) & rlang::is_empty(X_names)){
    message("Specify either a named initial condition X0 (e.g. X0 = (X1 = .1, X2 = .3)) or the names of your variables to set a random initial condition using X_names (e.g. X_names = c('X1', 'X2')).")
    return()
  }

  if (is.null(bifpar_list) & is.null(bifpar_pars)){
    message("Specify either a named list with bifurcation parameter values using bifpar_list or specify a list in bifpar_pars to construct such a list, e.g. bifpar_list = list(list('s' = .96)) or bifpar_pars = list(bifpar_start = .96, bifpar_end = .97, transition_steps = 100).")
    return()
  }  else if (is.null(bifpar_list)){
    bifpar_list = utils::modifyList(list(pre_steps = 0, baseline_steps = 100, transition_steps = 0, post_steps = 0), bifpar_pars) %>%
      do.call(get_bifurcation_range, .)
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
                          # rtol = 1e-6*10, atol = 1e-6*10,
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
        print(as.data.frame(out)[seq(1, nr_timesteps, length.out = 10),])
        break
      }

      if (do_downsample){
        # out = downsample(out, X_names,
        #                  win_size = downsample_pars$win_size,
        #                  which_X = downsample_pars$which_X,
        #                  seed_nr = downsample_pars$seed_nr)

        # out = do.call(downsample, utils::modifyList(downsample_pars,
                                                    # list(df = out, X_names = X_names)))
        out = downsample(out, X_names,
                         win_size = solve_sampling_par(fs = fs, timestep = timestep)$sample_interval,
                         which_X = which_X)
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

  if (success){
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
              bifpar_pars = bifpar_pars,
              do_downsample = do_downsample,
              # downsample_pars = downsample_pars
              fs = fs,
              which_X = which_X
              ))
  } else {
    message("Simulation NOT successful!")
    return()
  }
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
    } else {
      maxpeaks_df = data.frame()
    }
    if (!rlang::is_empty(minpeaks)){
      minpeaks_df = cbind(df[minpeaks[,2], c(i, setdiff(colnames(df), X_names))],
                          peak_idx = minpeaks[,2],
                          begin_peak_idx = minpeaks[,3], end_peak_idx = minpeaks[,4]
      ) %>% as.data.frame() %>% mutate(minmax = "minpeak")
    } else {
      minpeaks_df = data.frame()
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



#' Find type of regime switch
#'
#' @param from_regime The regime from which the transition occurs
#' @param to_regime The regime which is transitioned to
#' @param X_names Names of variables
#'
#' @return Regime switch type
#' @importFrom dplyr select mutate mutate_at filter slice group_by ungroup all_of rename row_number summarise pull bind_rows across .data
#'
#' @examples
get_regime_switch_type <- function(from_regime, to_regime, X_names){

  regime_df = bind_rows(from_regime, to_regime)

    # Check for hitting basin-boundary in the from_regime -> boundary-crisis
  if (grepl("Basin-Boundary", regime_df$regime[1]) & !grepl("Basin-Boundary", regime_df$regime[2])){
    bound_crisis = "Boundary-Crisis"
  } else {
    bound_crisis = ""
  }

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

  if (bound_crisis != ""){
    return(bound_crisis)
    # return(sprintf("%s: %s to %s",  bound_crisis, regime_df$regime[1], regime_df$regime[2]))
  } else if (nr_doubling == N){
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

  } else {    # Regime switches involving chaotic behaviour

    # Check for mixed periodicity
  regime_df_ = regime_df %>%
      tidyr::unite("z", X_names, sep = ";") %>%
    rowwise() %>%
      dplyr::mutate(check1 = paste0(unlist(stringr::str_split(.data$z, ";")), collapse = ","),
                    check2 = paste0(grepl("Period", unlist(stringr::str_split(.data$z, ";")), fixed = T), collapse = ","),
        nr_periods = ifelse( !all(grepl("Period", unlist(stringr::str_split(.data$z, ";")), fixed = T)), NA, length(unique(unlist(stringr::str_split(.data$z, ";")))))) %>%
      mutate(broad_regime = ifelse(is.na(.data$nr_periods), .data$regime, ifelse(.data$nr_periods == 1, .data$regime, "Mixed-Periodic"))) %>% ungroup()

     # Check for chaos expansion or reduction
    both_chaotic =  grepl("Chaotic", regime_df$regime[1], fixed = TRUE)& grepl("Chaotic", regime_df$regime[2], fixed = TRUE)
    type_chaos = grepl("Basin-Boundary", regime_df$regime, fixed = TRUE) | grepl("Merged-Band", regime_df$regime, fixed = TRUE)
    chaos_exp = !type_chaos[1] & type_chaos[2] & both_chaotic
    chaos_red = type_chaos[1] & !type_chaos[2] & both_chaotic
    # chaos_exp = grepl("Basin-Boundary", regime_df$regime[2], fixed = TRUE) & !grepl("Basin-Boundary", regime_df$regime[1], fixed = TRUE) & both_chaotic
    # chaos_red = grepl("Basin-Boundary", regime_df$regime[1], fixed = TRUE) & !grepl("Basin-Boundary", regime_df$regime[2], fixed = TRUE) & both_chaotic

    if (length(unique(regime_df_$regime)) == 1){
      return("Same-Chaos")
    } else if (any(grepl("None", regime_df$regime))){
      return(sprintf("%s to %s",
                     regime_df$regime[1], regime_df$regime[2]))
    } else if (chaos_red){
      return(trimws(sprintf("Chaos-Reduction")))
    } else if (chaos_exp){
      return(trimws(sprintf("Chaos-Expansion")))
    } else if (bound_crisis != ""){
      return(sprintf("%s: %s to %s",
                     bound_crisis, regime_df_$broad_regime[1], regime_df_$broad_regime[2]))
    } else {
      return(sprintf("%s to %s",
             regime_df_$broad_regime[1], regime_df_$broad_regime[2]))
    }
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
    return(data.frame(regime_switch_type = NA,
                      regime1 = NA,
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
#' @param rough_df Dataframe with periodicity per bifurcation parameter
#' @param col_to_smooth Column to smooth over
#' @param smooth_along_col Columns to smooth along too
#' @inheritParams find_regimes
#'
#' @return Dataframe with smoothed periodicity
#' @importFrom dplyr arrange mutate pull lag rowwise .data group_by ungroup rename
#'
#' @examples
smooth_df <- function(rough_df, col_to_smooth = "period", smooth_along_col = c(), nr_smooth = 0, min_length_regime = 10){

  min_length = rough_df %>% nrow() #dplyr::summarise(n = n(), .groups = 'drop') %>% pull(n) %>% min()
  if (min_length < min_length_regime*2){
    print("Warning: in smooth_df(), min_length_regime is not long enough. Skipping smoothing...")
    return(rough_df)
  } else {
   # If a consistent regime shows one exception, smooth over.
  find_nr_surrounding_regimes <- function(x, min_length_regime, nr_smooth){
    # Split vector
    m <- zoo::rollapply(x, min_length_regime*2 + nr_smooth, by = 1, FUN = c)
    # Find number of unique periods surrounding nr_smooth, excluding nr_smooth. Look back min_length_regime + nr_smooth - 1 steps (ignoring nr_smooth steps) and forward min_length_regime steps
    c(rep(NA, min_length_regime + ceiling(nr_smooth/2) - 1), split(m, row(m)) %>% purrr::map(
      function(x){length(unique(x[-seq(min_length_regime + 1, min_length_regime + nr_smooth)]))}) %>%
        unname %>% unlist, rep(NA, min_length_regime + (nr_smooth - ceiling(nr_smooth/2))) ) %>% return()
  }
  smooth_along <- function(x, to_smooth, nr_smooth){
    x[to_smooth] = lag(x, n = nr_smooth)[to_smooth]
    return(x)
  }
  # If a consistent regime shows one exception, smooth over.
  smooth_df = rough_df %>%
    mutate(lag_p = lag(.data[[col_to_smooth]], n = nr_smooth),
           # Apply rolling function to check for each row whether the previous min_length_regime values were the same and the future min_length_regime values were the same. If they were, fill in leading value.
           embedded_in_same_regime =
             # zoo::rollapply(.data$period, min_length_regime*2 + nr_smooth, function(x){length(unique(x[-seq(min_length_regime + 1, min_length_regime + nr_smooth)]))}, by = 1, fill=NA, align = 'center')
             find_nr_surrounding_regimes(.data[[col_to_smooth]], min_length_regime, nr_smooth) == 1 ) %>%
    ungroup() %>% rowwise() %>%
    rename(rough = !!col_to_smooth) %>% #select(-col_to_smooth) %>%
    mutate(to_smooth = ifelse(is.na(.data$lag_p) | is.na(.data$embedded_in_same_regime), F,
           ifelse((.data$embedded_in_same_regime == TRUE) & (.data$rough != .data$lag_p), T, F)),
      smooth = ifelse(.data$to_smooth, .data$lag_p, .data$rough)) %>% ungroup() %>%
    mutate_at(smooth_along_col, ~smooth_along(., .data$to_smooth, nr_smooth)) %>%
    select(-c(.data$lag_p, .data$embedded_in_same_regime, .data$to_smooth)) %>%
    rename(!!col_to_smooth := .data$smooth)


  return(smooth_df)
  }
}


#' Get bifurcation diagram band properties
#'
#' @param x Peak vector
#' @inheritParams periods_to_regimes
#' @param step_size Step size in histogram
#'
#' @return Band properties
#' @export
#'
#' @examples
get_bands = function(x, min_edge = 0, max_edge = 1, step_size = .05){
  # hist(x, breaks = seq(min_x, max_x, by = step_size))
  HIST = graphics::hist(x, breaks = seq(min_edge, max_edge, by = step_size), plot = F)

  dist_bins = diff(which(HIST$density > 0))
  dist_bins = dist_bins[dist_bins > 1]
  if (length(dist_bins) == 0){
    dist_bins = 0
  }
  HIST_df = data.frame(nr_bins = sum(HIST$density > 0),
                       max_dist_bands = max(dist_bins) * step_size,
                       total_dist_separate_bands = sum(dist_bins * step_size),
                       occupied_bins_in_band = sum(HIST$density > 0) / (diff(range(which(HIST$density > 0))) + 1),
                       min_band = HIST$mids[dplyr::first(which(HIST$density > 0))],
                       max_band = HIST$mids[dplyr::last(which(HIST$density > 0))]
  )
  return(HIST_df)
}


#' Update periodicities with basin boundaries
#'
#' @inheritParams find_regimes
#' @param peaks_df Dataframe with peaks
#' @param periods Dataframe with periodicities
#' @param X_names Names of variables
#'
#' @return Updated dataframe
#' @export
#' @importFrom dplyr select group_by group_by_at group_modify slice summarise arrange mutate filter ungroup all_of rename bind_rows n .data
#'
#' @examples
periods_to_regimes <- function(peaks_df, periods,
                               X_names,
                               min_length_regime,
                               variable_name = "X1", min_edge = 0, max_edge = 1,
                               thresh_full_band = .8, nr_smooth = 0
){

  # If there are no minima and maxima (i.e. only nodes), no chaotic regimes can be found
  if (!all(c("minpeak", "maxpeak") %in% peaks_df$minmax)){
    updated_periods=periods
  } else {

    minmax_peaks_df = peaks_df %>%
      select(.data$bifpar_idx, .data$variable, .data$minmax, .data$X) %>%
      group_by(.data$bifpar_idx, .data$variable, .data$minmax) %>%
      summarise(max = round(max(.data$X), 2), min = round(min(.data$X), 2), .groups='drop') %>%
      tidyr::pivot_wider(names_from = "minmax", values_from = c("max", "min")) %>%
      select(.data$bifpar_idx, .data$variable, .data$max_maxpeak,.data$min_minpeak)

  # Basin-Boundary
  basin_bound =  minmax_peaks_df %>%
    filter(.data$variable == !!variable_name) %>% select(-.data$variable) %>%
    filter(.data$max_maxpeak == !!max_edge & .data$min_minpeak == !!min_edge) %>%
    arrange(.data$bifpar_idx)

  # Merged bands
  bifpar_idx_chaotic = periods %>% filter(grepl("Chaotic", .data$period_bifpar, fixed = TRUE)) %>% pull(.data$bifpar_idx) %>% unique()
  merged_band_df = peaks_df %>%
    select(.data$bifpar_idx, .data$variable, .data$minmax, .data$X) %>%
    dplyr::filter(.data$variable == !!variable_name) %>%
    group_by(.data$bifpar_idx, .data$variable, .data$minmax) %>%
    group_modify(~ get_bands(x = .x$X, min_edge = min_edge, max_edge = max_edge,
                             # Step size dependent on how many data points there are
                             step_size = (max_edge - min_edge)/length(.x$X)  )) %>% ungroup() %>%
    dplyr::filter(.data$occupied_bins_in_band >= thresh_full_band) %>%
    filter(.data$bifpar_idx %in% bifpar_idx_chaotic) %>%
    ungroup()

  # Update which periods touch the basin-boundary
  updated_periods = periods %>% ungroup() %>% rowwise() %>%
    mutate(period_bifpar = ifelse(.data$bifpar_idx %in% basin_bound$bifpar_idx,
                                  paste0(.data$period_bifpar, " (Touching Basin-Boundary)"),
                                  .data$period_bifpar)
    ) %>%
    mutate(period_bifpar = ifelse(.data$bifpar_idx %in% merged_band_df$bifpar_idx,
                                  paste0(.data$period_bifpar, " (Merged-Band)"),
                                  .data$period_bifpar)
    ) %>% ungroup()

  }

  regimes_A = updated_periods %>%
    group_by_at(setdiff(colnames(.), "bifpar_idx")) %>%
    group_modify( ~ find_consec_seq(.x$bifpar_idx)) %>%
    arrange(.data$start_bifpar_idx) %>%
    ungroup() %>% rename(regime = .data$period_bifpar)

  regimes1 = combine_mixed_regimes(regimes_A, X_names, min_length_regime,
                                   name_mixed_regime = "Mixture: Periodic and Chaotic",
                                   grepl_mixed_regime = function(regime){(grepl("Chaotic", regime, fixed = TRUE) & grepl("Period", regime, fixed = TRUE) & !grepl("Merged-Band", regime, fixed = TRUE) & !grepl("Basin-Boundary", regime, fixed = TRUE))  }
  )

  regimes2 = combine_mixed_regimes(regimes1, X_names, min_length_regime,
                                   name_mixed_regime = "Mixture: Periodic and Chaotic (Merged-Band)",
                                   grepl_mixed_regime = function(regime){(((grepl("Chaotic", regime, fixed = TRUE) & grepl("Period", regime, fixed = TRUE)) | grepl("Chaotic or Transitioning (X1,X2,X3,X4)",regime,fixed=T)) & !grepl("Basin-Boundary", regime, fixed = TRUE))  }
  )


  regimes = regimes2 %>%
  tibble::rownames_to_column() %>%
    slice(rep(1:n(), each = .data$length_region)) %>%
    group_by(.data$rowname) %>%
    mutate(bifpar_idx = seq(unique(.data$start_bifpar_idx),
                            unique(.data$end_bifpar_idx))) %>% ungroup() %>%
    smooth_df(col_to_smooth = "regime", smooth_along_col = X_names, nr_smooth = nr_smooth, min_length_regime=min_length_regime) %>%
    select(all_of(X_names), .data$bifpar_idx, .data$regime) %>%
    # Apply regular regime-finding
    group_by_at(setdiff(colnames(.), "bifpar_idx")) %>%
    group_modify( ~ find_consec_seq(.x$bifpar_idx)) %>%
    arrange(.data$start_bifpar_idx) %>%
    ungroup()

  return(regimes)

}




#' Combine regimes that are too short to constitute their own regime
#'
#' @inheritParams periods_to_regimes
#' @param regimes_A Dataframe with regime
#' @param name_mixed_regime Name of new regime
#' @param grepl_mixed_regime Function to detect regimes that should be grouped together
#'
#' @return Dataframe with updated regimes
#' @export
#'
#' @examples
combine_mixed_regimes <- function(regimes_A, X_names, min_length_regime,

                                  name_mixed_regime = "Mixture: Periodic, Chaotic or Transitioning (X1,X2,X3,X4)",
                                  grepl_mixed_regime = function(regime){grepl("Chaotic or Transitioning (X1,X2,X3,X4)", regime, fixed = TRUE) | (grepl("Chaotic or Transitioning", regime, fixed = TRUE) & grepl("Period", regime, fixed = TRUE))}

){

  # Sometimes, the chaotic regime does not consistently touch the basin boundary, but switches about every bifurcation index. Merge regimes that are too short to form their own regime but have chaos in them.
  regimes_B = regimes_A %>% rowwise() %>%
    # By excluding those regimes that already satisfy the minimum length, we group together 'stray' regimes of short length
    mutate(regime_mixed = ifelse(.data$length_region < min_length_regime &
                                   # !grepl("Merged-Band", .data$regime, fixed = TRUE) &
                                   grepl_mixed_regime(.data$regime),
                                 name_mixed_regime, NA))


  if (all(is.na(regimes_B$regime_mixed))){
    regimes_C = data.frame()
  } else {
    regimes_C_ = regimes_B %>%
      filter(!is.na(.data$regime_mixed)) %>%
      select(-.data$regime) %>%
      # Turn regime of start to end_bifpar_idx into separate rows per bifpar_idx
      tibble::rownames_to_column() %>%
      slice(rep(1:n(), each = .data$length_region)) %>%
      group_by(.data$rowname) %>%
      mutate(bifpar_idx = seq(unique(.data$start_bifpar_idx),
                              unique(.data$end_bifpar_idx))) %>% ungroup() %>%
      dplyr::mutate_at(all_of(X_names), ~ "Mixture") %>%
      select(all_of(X_names), .data$bifpar_idx, .data$regime_mixed) %>%
      # Apply regular regime-finding
      group_by_at(setdiff(colnames(.), "bifpar_idx")) %>%
      group_modify( ~ find_consec_seq(.x$bifpar_idx)) %>%
      arrange(.data$start_bifpar_idx) %>%
      ungroup() %>% rename(regime = .data$regime_mixed) %>%
      filter(.data$length_region >= min_length_regime)

    # Add back rows whose regimes were not long enough
    if (nrow(regimes_C_) > 0){
      added_bifpar_idx = plyr::llply(1:nrow(regimes_C_), function(idx){seq(regimes_C_$start_bifpar_idx[idx], regimes_C_$end_bifpar_idx[idx] )}) %>% unlist() %>% unique()
      regimes_C = dplyr::bind_rows(regimes_C_,
                                   regimes_B %>%
                                     filter(!is.na(.data$regime_mixed)) %>%
                                     select(-.data$regime_mixed) %>% rowwise() %>%
                                     dplyr::filter(any(!(seq(.data$start_bifpar_idx, .data$end_bifpar_idx) %in% added_bifpar_idx))) %>% ungroup()
        )
    } else {
      regimes_C = regimes_B %>%
        filter(!is.na(.data$regime_mixed)) %>%
        select(-.data$regime_mixed)
    }

    # # If there's any member of the regime that touches the basin boundary, update
    # regimes_C = regimes_C_ %>% rowwise() %>%
    #   mutate(regime = ifelse(any(seq(.data$start_bifpar_idx, .data$end_bifpar_idx) %in% basin_bound$bifpar_idx),
    #                          paste0(.data$regime, " (Not) Touching Basin-Boundary)"),
    #                          .data$regime)
    #   )

  }
  # Merge
  regimes = bind_rows(
    regimes_B %>%
      filter(is.na(.data$regime_mixed)) %>% select(-.data$regime_mixed),
    regimes_C
  ) %>% arrange(.data$start_bifpar_idx)
  return(regimes)
}



#' Find regimes of dynamical system
#'
#' @param GLV List object generated by bifurcation_ts()
#' @param thresh_node Threshold under which timeseries is classified as node
#' @param thresh_coord_spread Threshold for distance in peak and trough coordinates that determines whether a cluster (i.e. period length) fits; if exceeded, denoted as chaotic or transitioning
#' @param thresh_peak_idx_spread Same as thresh_coord_spread but for peak and trough indices
#' @param thresh_full_band Percentage of bins that need to be occupied to qualify it as a fully filled banned
#' @param min_length_regime Minimum number of consecutive steps in the bifurcation parameter that have the same periodicity to qualify as a regime
#' @param max_k Maximum cluster size to look for
#' @param nr_smooth Number of exceptions in a stable periodicity window to smooth over; nr_smooth = 0 means no smoothing
#' @param factor_k Weighting of period length k; heavier weight means shorter k is preferred; factor_k = 0 means the optimal period length is chosen based solely on minimum spread
#' @param variable_name Column name in dataframe to assess for hitting basin boundaries
#' @param min_edge Minimum basin boundary
#' @param max_edge Maximum basin boundary
#' @param keep_nr_timesteps Number of timesteps to keep for each bifurcation parameter value; number of timesteps to retain after discarding the transient. If "all", none are removed
#'
#' @return List of dataframes with periodicity per variable, periodicity per bifurcation parameter value, regimes, and regime boundaries
#' @importFrom dplyr arrange group_modify ungroup filter select group_by mutate rename bind_rows all_of slice_tail .data
#' @export
#'
#' @examples
find_regimes <- function(GLV,
                         thresh_node = .001,
                         thresh_coord_spread = .025,
                         thresh_peak_idx_spread=2,
                         thresh_full_band = .9,
                         min_length_regime = 5,
                         nr_smooth = 0,
                         factor_k = .1,
                         variable_name = "X1", min_edge = 0, max_edge = 1,
                         max_k = NULL,
                         keep_nr_timesteps = c("all", 1000)){
#
# max_k = NULL
#  thresh_node = .001
# thresh_coord_spread = .025
# thresh_peak_idx_spread=2
# nr_smooth=0
# factor_k = .1
# min_length_regime = 10
# variable_name = "X1"
# min_edge = 0
# max_edge = 1
# thresh_expansion = .1
# thresh_full_band = .9
# keep_nr_timesteps = "all"
# peaks_df=regime_list$peaks_df
# k_spread = regime_list$k_spread
# period_per_var=regime_list$period_per_var
# GLV = regime_list
# X_names = GLV$X_names

  # Get dataframe with peaks
  if (keep_nr_timesteps == "all"){
    keep_nr_timesteps = GLV$df %>% dplyr::filter(bifpar_idx == 1) %>% nrow()
  } else if (is.numeric(keep_nr_timesteps)){
    keep_nr_timesteps = keep_nr_timesteps * GLV$fs
  }
  peaks_df = peaks_bifdiag(GLV$df %>%
                             # Discard transient if desired
                             dplyr::slice_tail(n = keep_nr_timesteps,
                                                        by = "bifpar_idx"),
                           GLV$X_names) %>%
    filter(.data$bifpar_idx >= 5) # Skip initial settling in points

  # For each value of the bifurcation parameter, find the period length k which has a minimum spread
  print("Finding best fitting period length for all bifurcation parameter values")
  k_spread = peaks_df %>%
    # filter(bifpar_idx > 560, bifpar_idx < 570) %>%
    group_by(.data$variable, .data$bifpar_idx) %>%
    arrange(.data$time_idx, .by_group=TRUE) %>%
    group_modify(~ find_spread_per_k(coord = .x$X, peak_idx = .x$peak_idx,
                                     max_k = max_k)) %>%
    ungroup()

  # Select best k
  best_k = k_spread %>%
    group_by(.data$variable, .data$bifpar_idx) %>%
    group_modify(~ choose_best_k(spread_df = .x, thresh_node = thresh_node,
                                 factor_k = factor_k)) %>%
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
    group_by(.data$variable) %>%
    arrange(.data$bifpar_idx, .by_group = TRUE) %>%
    smooth_df(col_to_smooth = "period", smooth_along_col = c(), nr_smooth = nr_smooth, min_length_regime=min_length_regime) %>% ungroup()

  periods = period_per_var %>%
    select(.data$bifpar_idx, .data$period, .data$variable) %>%
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
    # group_by(.data$bifpar_idx, .data$period_bifpar) %>%
    # group_modify(~ function(x = .x, y = .y){
    #   if (grepl("Chaotic or Transitioning",
    #             y$period_bifpar, fixed = TRUE)){
    #     return(x %>% mutate_at(X_names,
    #                            function(x){"Chaotic or Transitioning"}))
    #   } else {
    #     return(x)
    #   }
    # })
    # mutate_at(X_names, function(x){
    #                         ifelse(grepl("Chaotic or Transitioning",
    #                                      .data$period_bifpar, fixed = TRUE), "Chaotic or Transitioning", .)}) %>% ungroup() %>%
    arrange(.data$bifpar_idx)

  # Any time point containing at least one variable displaying chaotic behaviour is labelled as chaotic overall
  periods[grepl("Chaotic or Transitioning", periods$period_bifpar, fixed = TRUE), GLV$X_names] = "Chaotic or Transitioning"

  regimes = periods_to_regimes(peaks_df, periods,
                                 X_names = GLV$X_names,
                               min_length_regime=min_length_regime,
                                 variable_name = variable_name,
                               min_edge = min_edge,
                               max_edge = max_edge,
                               thresh_full_band = thresh_full_band, nr_smooth = nr_smooth)

  # Find regime boundaries
  regime_bounds_ = find_regime_bounds(regimes, min_length_regime = min_length_regime,
                                      X_names = GLV$X_names)

  # Add corresponding initial conditions - the timepoint right before the bifurcation parameter changed to the starting value of the first regime
  regime_bounds =  merge(regime_bounds_,
                         GLV$X0s  %>% as.data.frame() %>%
                           filter(.data$bifpar_idx %in% c(regime_bounds_$regime1_start_idx)) %>%
                           rename(regime1_start_idx = .data$bifpar_idx),
                         all.x = TRUE)

  regime_list = utils::modifyList(GLV, list(df = NULL,
                                            peaks_df = peaks_df,
                                            k_spread = k_spread,
                                            period_per_var = period_per_var,
                                            periods = periods,
                                            regimes = regimes,
                                            regime_bounds = regime_bounds,
                                            thresh_coord_spread = thresh_coord_spread,
                                            thresh_peak_idx_spread=thresh_peak_idx_spread,
                                            thresh_full_band=thresh_full_band,
                                            min_length_regime = min_length_regime,
                                            max_k = max_k,
                                            nr_smooth = nr_smooth,
                                            variable_name = variable_name,
                                            min_edge=min_edge,
                                            max_edge=max_edge,
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
  if (length(vec) <= 1){
    max_dist_per_k = 0
  } else {
  # From the package sarafrr/basicClEval
  vec = as.matrix(vec)
  sizeCl <- summary(as.factor(cluster_idx))
  # Compute maximal distance within each cluster
  max_dist_per_k = unlist(lapply(split(vec, cluster_idx), function(x){max(stats::dist(x))}))
}

  # return(mean(max_dist_per_k))
  return(c(
    mean_spread = mean(max_dist_per_k),
    median_spread = stats::median(max_dist_per_k),
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
                                                 cluster_idx = rep(1:k, length.out = length(coord)))}) %>%
    do.call(rbind, .) %>%
    magrittr::set_colnames(paste0(colnames(.), "_coord")) %>%
    as.data.frame()


  peak_idx_spread = lapply(ks, function(k){max_dist(vec = diff(peak_idx),
                                                    cluster_idx = rep(1:k, length.out = length(coord)))}) %>%
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
    max_k = floor(length(coord)/2)-1
  } else {
    # Check whether max_k is possible to look for
    if (max_k > (floor(length(coord)/2)-1)){
      max_k = (floor(length(coord)/2)-1)
    }
  }
  max_k = ifelse(max_k < 1, 1, max_k)
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
choose_best_k <- function(spread_df, thresh_node, factor_k){

  # Detect nodes for which k = 1 should be a good fit
  # if (spread_df[spread_df$k == 1, c("max_spread_coord")] < thresh_node){
  if (length(spread_df$max_spread_coord) == 1){
    idx_min = 1
  } else if (max(stats::dist(spread_df$max_spread_coord)) < thresh_node){ # If all fits are the same
    idx_min=1
  } else {
  # Method 1: log
  spread_df = spread_df[spread_df$k != 1,] # If it's not a node, k = 1 is not a good partitioning
    idx_min = spread_df %>% ungroup() %>%
      select(c("k", "max_spread_coord", "median_spread_coord"
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

  }
  return(as.data.frame(spread_df[idx_min,]))
}
