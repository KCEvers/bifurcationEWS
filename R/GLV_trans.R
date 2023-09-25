#' @title Example transition model
#' @description Data set of a period-4 to period-8 transition generated using bifurcation_ts()
#' @format A named list:
#' \describe{
#'   \item{df}{Dataframe with timeseries}
#'   \item{X_names}{Names of variables in model}
#'   \item{X0s}{Dataframe with initial condition per bifurcation parameter value}
#'   \item{timestep}{Step size in time}
#'   \item{nr_timesteps}{Number of time steps}
#'   \item{model_pars}{Fixed model parameters}
#'   \item{seed_nr}{Seed number}
#'   \item{stopifregime}{End generating timeseries if this function is satisfied}
#'   \item{deSolve_method}{ODE in the style of deSolve}
#'   \item{bifpar_list}{List with changing bifurcation parameter values}
#'   \item{bifpar_pars}{List with parameters to generate bifpar_list, needs to include list(bifpar_start = .96, bifpar_end = .96), optionally specify list(pre_steps = 0, baseline_steps = 100, transition_steps = 0, post_steps = 0)}
#'   \item{do_downsample}{Reduce dataframe size by downsampling?}
#'   \item{downsample_pars}{List of parameters for downsampling}
#'}
"GLV_trans"
