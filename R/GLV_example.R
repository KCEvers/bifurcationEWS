#' Example data set of Generalized Lotka Volterra Model
#'
#' Description of the dataset.
#'
#' @format A list with 11 entries:
#' \describe{
#'   \item{"df"}{Dataframe}
#'   \item{"X_names"}{Names of variablesin model}
#'   \item{"X0s"}{Initial state}
#'   \item{"timestep"}{Step size in time}
#'   \item{"nr_timesteps"}{Number of time steps}
#'   \item{"model_pars"}{Model parameters}
#'   \item{"seed_nr"}{Seed number}
#'   \item{"stopifregime"}{End generating timeseries if this function is satisfied}
#'   \item{"deSolve_method"}{Method of generating ODE passed to deSolve::ode}
#'   \item{"bifpar_list"}{List with changing bifurcation parameter values}
#'   \item{"do_downsample"}{Reduce dataframe size by downsampling?}
#'   \item{"downsample_pars"}{List of parameters for downsampling}
#' }
#' @source \url{data_URL}
"GLV_example"
