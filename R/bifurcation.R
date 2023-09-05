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
#'
#' @return List with timeseries matrix and parameters used to generate the timeseries
#' @importFrom magrittr `%>%`
#' @export
#'
#' @examples
bifurcation_ts <- function(model, model_pars, bifpar_list,
                           X0 = c(),
                           X_names = names(X0),
                           seed_nr = 123,
                           timestep = .01, nr_timesteps = 100,
                           deSolve_method = c("lsoda", "euler", "rk4")[1],
                           stopifregime = function(out){FALSE},
                           do_downsample = TRUE,
                           downsample_pars = list(type = c("average", "one_sample")[1],
                                                  win_size = 10,
                                                  which_X = c(10, "first", "middle", "last", "random")[1],
                                                  seed = 123),
                           silent = FALSE){
  if (rlang::is_empty(X0) & rlang::is_empty(X_names)){
    message("Specify either a named initial condition X0 (e.g. X0 = (X1 = .1, X2 = .3)) or the names of your variables to set a random initial condition using X_names (e.g. X_names = c('X1', 'X2')).")
    return()
  }
  # Initialize
  if (rlang::is_empty(X0)){
    set.seed(seed_nr)
    X0      <- stats::runif(length(X_names)) %>% stats::setNames(X_names)
  }
  times = seq(0, nr_timesteps, by = timestep)
  min_t = times[1]
  bifpar_idxs = seq.int(length(bifpar_list))
  X0s <- matrix(nrow = length(bifpar_list), ncol = length(X_names)+1)
  tmps <- list()

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
                        method = deSolve_method)
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
    if (stopifregime(out)){
      message("Hit undesirable regime! Deleting generated timeseries")
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

  # Collect all generate data and remove temporary files
  OUT = do.call(rbind, lapply(bifpar_idxs,
                              function(bifpar_idx){
                                tmp <- tmps[[bifpar_idx]]
                                out <- readRDS(tmp)
                                unlink(tmp)
                                return(cbind(out, bifpar_idx))
                              }))
  # OUT$time_idx = 1:nrow(OUT)
  # OUT %>% head

  return(list(df = as.data.frame(cbind(OUT, time_idx=1:nrow(OUT))),
              X_names = X_names,
              X0s = X0s,
              times = times,
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
#' @importFrom dplyr select mutate filter slice group_by ungroup all_of rename
#' @importFrom magrittr `%>%`
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
      slice(1, .by = bifpar_idx) %>% ungroup() %>%
      filter(bifpar_idx %in% missing_bifpar_idx) %>%
      mutate(minmax="node",
                    peak_idx = time_idx,
                    begin_peak_idx = time_idx, end_peak_idx = time_idx)
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
#' @importFrom magrittr `%>%`
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
#' @return Start and end of each consecutive value sequence
#' @importFrom magrittr `%>%`
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
    dplyr::mutate(region_nr = 1:nrow(.), nr_regions = nrow(.)) %>%  dplyr::rowwise()

  return(df_conseq_seq %>% dplyr::arrange(start_bifpar_idx))
}



#' Find whether and when basin boundaries are hit
#'
#' @param peaks_df Dataframe with peaks
#' @param variable_name Column name in dataframe to assess
#' @param min_edge Minimum value reached
#' @param max_edge Maximum value reached
#'
#' @return Dataframe with bifurcation parameter value for which the system touches the basin boundaries as defined by min_edge and max_edge, if it hits these. Otherwise, dataframe filled with NA.
#' @importFrom dplyr select mutate filter slice group_by ungroup all_of rename row_number summarise pull arrange n
#' @importFrom magrittr `%>%`
#'
#' @examples
find_basin_boundary <- function(peaks_df, variable_name = "X1", min_edge = 0, max_edge = 1){

  minmax_peaks_df = peaks_df %>% select(bifpar_idx, variable, minmax, time_idx, X) %>%
    filter(variable == variable_name) %>%
    group_by(bifpar_idx, variable, minmax) %>%
    summarise(max = round(max(X), 2), min = round(min(X), 2), .groups='drop') %>%
    tidyr::pivot_wider(names_from = "minmax", values_from = c("max", "min")) %>%
    select(bifpar_idx, variable, max_maxpeak, min_minpeak) %>% ungroup %>%
    filter(max_maxpeak == max_edge & min_minpeak == min_edge) %>% arrange(bifpar_idx) %>%
    filter(row_number()==1 | row_number()==n())
  # If the edges were not touched
  if (nrow(minmax_peaks_df) == 0){
    minmax_peaks_df = data.frame(bifpar_idx = NA)
  }

  basin_bound = data.frame(start_bifpar_idx = minmax_peaks_df %>% slice(1) %>% pull(bifpar_idx),
                           end_bifpar_idx = minmax_peaks_df %>% slice(nrow(.)) %>% pull(bifpar_idx),
                           regime = "Basin-Boundary"
  ) %>%
    mutate(length_region = end_bifpar_idx- start_bifpar_idx + 1, region_nr = ifelse(is.na(start_bifpar_idx), NA, 1), nr_regions = ifelse(is.na(start_bifpar_idx), NA, 1))

  return(basin_bound)
}


#' Find type of regime switch
#'
#' @param from_regime The regime from which the transition occurs
#' @param to_regime The regime which is transitioned to
#' @param X_names Names of variables
#'
#' @return Regime switch type
#' @importFrom dplyr select mutate mutate_at filter slice group_by ungroup all_of rename row_number summarise pull bind_rows
#' @importFrom magrittr `%>%`
#'
#' @examples
get_regime_switch_type <- function(from_regime, to_regime, X_names){

  regime_df = bind_rows(from_regime, to_regime)
  # Regime switches involving only periodic regimes
  # if (!grepl("Chaotic", from_regime$regime, fixed = TRUE) & !grepl("Chaotic", to_regime$regime, fixed = TRUE)){
  if (all(!grepl("Chaotic", regime_df$regime)) & all(!grepl("None", regime_df$regime))){
  period_switch = regime_df %>%
    mutate_at(X_names, ~ as.numeric(stringr::str_replace(., "Period-", ""))) %>%
    select(all_of(X_names)) %>% as.matrix()

  period_switch_factor = period_switch[2,] / period_switch[1,]

  if (all(period_switch_factor == 2)){
   return("Period-Doubling")
  } else if (all(period_switch_factor == 0.5)){
    return("Period-Halving")
  } else if (all(period_switch_factor > 1)){
    return("Period-Increasing")
  } else if (all(period_switch_factor < 1)){
    return("Period-Decreasing")
  } else if (all(period_switch_factor == 1)){
    return("Same-Period")
  }

  } else {
    # Regime switches involving chaotic behaviour
    broad_period_switch = regime_df %>%
      tidyr::gather(X, value, -setdiff(colnames(.), X_names)) %>%
      group_by(regime) %>%
      summarise(nr_periods = length(unique(value)), period = paste0(unique(value), collapse = ",")) %>%
      mutate(broad_regime = ifelse(grepl("Chaotic or Transitioning", regime), "Chaotic or Transitioning",
                                          ifelse(grepl("None", regime), "None",
                           ifelse(nr_periods == 1, period, "Mixed-Periodic"))))
    return(sprintf("%s to %s",
            broad_period_switch[from_regime$regime == broad_period_switch$regime, ] %>% pull(broad_regime),
                broad_period_switch[to_regime$regime == broad_period_switch$regime,] %>% pull(broad_regime)))
  }
}


#' Find regime boundaries
#'
#' @param regimes Dataframe with periodicity per value of the bifurcation parameter
#' @param min_length_regime Minimum number of consecutive steps in the bifurcation parameter that have the same periodicity to qualify as a regime
#' @param X_names Names of variables in model
#'
#' @return Dataframe with regime boundaries
#' @importFrom magrittr `%>%`
#'
#' @examples
find_regime_bounds <- function(regimes, min_length_regime, X_names){
  # Find regimes that satisfy a certain size
  regimes = regimes %>% dplyr::ungroup() %>% dplyr::arrange(start_bifpar_idx)
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
    from_regime = regimes[regime_idx[i],] %>% dplyr::mutate(regime = ifelse(is.na(regime), "None", regime))
    to_regime = regimes[regime_idx[i+1],] %>% dplyr::mutate(regime = ifelse(is.na(regime), "None", regime))
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
#' @param min_length_regime Minimum length of regime
#'
#' @return Dataframe with smoothed periodicity
#' @importFrom dplyr arrange mutate pull lag rowwise
#' @importFrom magrittr `%>%`
#'
#' @examples
smooth_periods <- function(periods, min_length_regime){
  # If a consistent regime shows one exception, smooth over.
  smooth_idxs = periods %>% arrange(bifpar_idx) %>%
    mutate(lag_p = lag(period_bifpar),
                  # Apply rolling function to check for each row whether the previous min_length_regime values were the same and the future min_length_regime values were the same. If they were, fill in leading value.
  embedded_in_same_regime = zoo::rollapply(period_bifpar, min_length_regime*2 + 1, function(x){length(unique(x[-(min_length_regime + 1)]))}, by = 1, fill=NA, align = 'center') ) %>% rowwise() %>%
    mutate(period_bifpar_smooth = ifelse(embedded_in_same_regime == 1 & period_bifpar != lag_p, 1, 0)) %>% pull(period_bifpar_smooth) %>% as.logical() %>% which()
  if (!rlang::is_empty(smooth_idxs)){
  periods[smooth_idxs, setdiff(colnames(periods), "bifpar_idx")] = periods[smooth_idxs - 1, setdiff(colnames(periods), "bifpar_idx")]
  }

  return(periods)
}

#' Find regimes of dynamical system
#'
#' @param df Dataframe
#' @param X_names Names of variables in model
#' @param thresh_coord_spread Threshold for within cluster sum of squares in peak and trough coordinates that determines a cluster (i.e. period length) fits; if exceeded, denoted as chaotic or transitioning
#' @param thresh_peak_idx_spread Same as thresh_coord_WCSS but for peak and trough indices
#' @param min_length_regime Minimum number of consecutive steps in the bifurcation parameter that have the same periodicity to qualify as a regime
#' @param max_k Maximum cluster size to look for
#'
#' @return List of dataframes with periodicity per variable, periodicity per bifurcation parameter value, regimes, and regime boundaries
#' @importFrom dplyr arrange group_modify ungroup filter select group_by mutate rename bind_rows all_of slice_tail
#' @importFrom magrittr `%>%`
#' @export
#'
#' @examples
find_regimes <- function(df,
                         X_names,
                         thresh_coord_spread = .025,
                         thresh_peak_idx_spread=2,
                         min_length_regime = 5, max_k = NULL){

  # max_k = NULL
  # thresh_coord_spread = .025
  # thresh_peak_idx_spread=2
  # min_length_regime = 5

  # Get dataframe with peaks
  peaks_df = peaks_bifdiag(df, X_names)

  # For each value of the bifurcation parameter, find the period length k which has a minimum WCSS.
  print("Finding best fitting period length for all bifurcation parameter values")
  period_per_var_ = peaks_df %>% group_by(variable, bifpar_idx) %>%
    arrange(time_idx, .by_group=TRUE) %>%
    group_modify(~ find_best_k(coord = .x$X, peak_idx = .x$peak_idx, max_k = max_k)) %>%
    ungroup()

  period_per_var = period_per_var_ %>%
    # tidyr::unnest(names(.)) %>%
    # Our algorithm will *always* find a period length k with a minimum WCSS, also for chaotic data. Set a threshold that decides which WCSS is too large to be classified as a neat periodic sequence.
    # mutate(period = ifelse(coord_WCSS > thresh_coord_WCSS & peak_idx_WCSS > thresh_peak_idx_WCSS,
    #                               "Chaotic or Transitioning", paste0("Period-", k) ))
    mutate(period = ifelse(max_spread_coord > thresh_coord_spread & max_spread_peak_idx > thresh_peak_idx_spread,
                                  "Chaotic or Transitioning", paste0("Period-", k) ))

  periods_ = period_per_var %>% select(bifpar_idx, period, variable) %>%
    group_by(bifpar_idx, period) %>%
    mutate(period_group = paste0(
      unique(period),
      " (",
      paste0(sort(unique(variable)), collapse = ","),
      ")"
    )
    ) %>% group_by(bifpar_idx) %>%
    mutate(period_bifpar = paste0(sort(unique(period_group)), collapse = ' AND ')
    ) %>% select(-period_group) %>% ungroup() %>%
    tidyr::pivot_wider(names_from = variable, values_from = period) %>%
    arrange(bifpar_idx)
  periods = smooth_periods(periods_, min_length_regime)

  # Get basin boundaries (when system hits edges of basin)
  basin_bound = find_basin_boundary(peaks_df, variable_name = "X1", min_edge = 0, max_edge = 1)

  # Find regimes with any chaotic behaviour and regimes with periodic behaviour
  broad_regimes = periods %>%
    mutate(period_bifpar2 = ifelse(grepl("Chaotic or Transitioning", period_bifpar, fixed = TRUE),
                                          "Chaotic or Transitioning", "Periodic")) %>%
    group_by(period_bifpar2) %>%
    group_modify( ~ find_consec_seq(.x$bifpar_idx)) %>% arrange(start_bifpar_idx) %>%
    ungroup() %>%
    rename(regime = period_bifpar2)

   # Compile regimes
  regimes = periods %>%
    filter(!grepl("Chaotic or Transitioning", period_bifpar, fixed = TRUE)) %>%
    group_by(period_bifpar, X1, X2, X3, X4) %>%
    group_modify( ~ find_consec_seq(.x$bifpar_idx)) %>% arrange(start_bifpar_idx) %>%
    ungroup() %>% rename(regime = period_bifpar) %>%
    bind_rows(basin_bound) %>%
    bind_rows(broad_regimes %>% filter(regime != "Periodic")) %>% arrange(start_bifpar_idx)

  # Find regime boundaries
  regime_bounds_ = find_regime_bounds(regimes, min_length_regime = min_length_regime, X_names = X_names)
    # Add corresponding initial conditions - the timepoint right before the bifurcation parameter changed to the starting value of the first regime
  regime_bounds =  merge(regime_bounds_,
                         df %>%
            filter(bifpar_idx %in% c(regime_bounds_$regime1_start_idx - 1)) %>%
              slice_tail(n=1, by = bifpar_idx) %>%
            mutate(regime1_start_idx= bifpar_idx + 1) %>%
            select(regime1_start_idx, all_of(X_names)), all.x = TRUE)

  return(list(peaks_df = peaks_df,
    period_per_var = period_per_var,
              periods = periods,
              regimes = regimes,
              broad_regimes = broad_regimes,
              regime_bounds = regime_bounds,
              thresh_coord_spread = thresh_coord_spread,
              thresh_peak_idx_spread=thresh_peak_idx_spread,
              min_length_regime = min_length_regime, max_k = max_k))
}




#' Find distance between points in chosen temporal sequence
#'
#' @param vec Vector
#' @param cluster_idx Vector of cluster indices
#'
#' @importFrom magrittr `%>%`
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
    min_spread = min(max_dist_per_k),
    max_spread = max(max_dist_per_k)
  ))
}


#' Compute mean maximum distance per cluster
#'
#' @param ks Vector of cluster sizes
#' @param coord Peak and trough coordinates
#' @param peak_idx Peak and trough indices
#'
#' @importFrom magrittr `%>%`
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


#' Find best fitting period length k
#'
#' @param coord Peak and trough coordinates
#' @param peak_idx Peak and trough indices
#' @param max_k Maximum cluster size to look for
#'
#' @return Dataframe with best fitting period length k and the corresponding minimum within-cluster distance and between-cluster distance for peak and trough coordinates and indices
#' @importFrom magrittr `%>%`
#' @export
#'
#' @examples
find_best_k <- function(coord, peak_idx, max_k = NULL){
  # Finding the best period is tricky: The global minimum might be a subharmonic, but the first local minimum might not be optimal for a more complex oscillation. If the timeseries is truly periodic, local minima in WCSS will occur for every subharmonic (e.g. period = 4, minima at period = c(4,8,12,...)).
  if (is.null(max_k)){
    max_k = (floor(length(coord)/2)-1)
  } else {
    # Check whether max_k is possible to look for
    if (max_k > (floor(length(coord)/2)-1)){
      max_k = (floor(length(coord)/2)-1)
    }
  }
  ks = seq(2, max_k) # Only consider partitionings with at least two repetitions

  spread_df = find_dist_per_k(ks, coord, peak_idx)

  # plot(log(spread_df[,c("mean_spread_coord")] * spread_df[,c("k")]) )
  # plot(log(spread_df[,c("max_spread_coord")] * spread_df[,c("k")]) )
  # spread_df[,c("max_spread_coord")][c(9,19)]
  # log(spread_df[,c("max_spread_coord")] * spread_df[,c("k")] )[c(9,19)]
  # log(spread_df[,c("max_spread_coord")] * log(spread_df[,c("k")]) )[c(9,19)]
  # log(spread_df[,c("max_spread_coord")] * cumsum(spread_df[,c("k")]) )[c(9,19)]
  # (spread_df[,c("max_spread_coord")] )[c(3,7)]
  #
  # dev.new()
  # plot(log(spread_df[,c("mean_spread_coord")] ) )
  # plot(log(spread_df[,c("max_spread_coord")] ) )
  #
  # plot(log(spread_df[,c("mean_spread_peak_idx")] * spread_df[,c("k")]) )
  # plot(log(spread_df[,c("max_spread_peak_idx")] * spread_df[,c("k")]) )
  #
  # plot(log(spread_df[,c("mean_spread_peak_idx")] ) )
  # plot(log(spread_df[,c("min_spread_peak_idx")] ) )
  # plot(log(spread_df[,c("max_spread_peak_idx")] ) )

  # Find best fitting period length k. First apply a transformation where the spread is multiplied by k in order to compensate for better fits simply due to more splitting. Apply a log-transformation to better see small differences, and choose the minimum value across the mean and maximum spread in peak coordinates and peak indices
  # idx_min = log(spread_df[,c("mean_spread_coord", "max_spread_coord",
                             # "mean_spread_peak_idx", "max_spread_peak_idx")] * spread_df[,c("k")]) %>% scale(center = F) %>% apply(1, mean) %>% which.min()


  # idx_min = log(spread_df[,c("max_spread_coord", "max_spread_peak_idx")] * spread_df[,c("k")]) %>% scale(center = F) %>% apply(1, mean) %>% which.min()
  # idx_min = spread_df %>% mutate_at(c("mean_spread_coord", "max_spread_coord",
                                      # "mean_spread_peak_idx", "max_spread_peak_idx"), ~ ifelse(. == 0, 0, log(. * k))) %>% scale(center = F) %>% apply(1, mean) %>% which.min()
  idx_min = spread_df %>% dplyr::mutate_at(c("max_spread_coord", "max_spread_peak_idx"), ~ ifelse(. == 0, 0, log(. * k))) %>% scale(center = F) %>% apply(1, mean) %>% which.min()

  # idx_min = log(spread_df[,c("max_spread_coord")] * spread_df[,c("k")]) %>%
  #   # scale(center = F) %>% apply(1, mean)  %>%
  #    which.min()

  # # Verifying the best period length k
  # cluster_idx = rep(1:(spread_df[idx_min,]$k), length(coord))[1:length(coord)]
  # plot(cluster_idx, coord, cex = .5, pch = 16, type = 'p')

  return(as.data.frame(spread_df[idx_min,]))
}
