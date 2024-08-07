#' Lorenz model
#'
#' @param t Time
#' @param state Named vector of variable state
#' @param parms List of model parameters
#'
#' @return Updated state
#' @export
Lorenz_model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dX <-  a * X + Y * Z
    dY <-  b * (Y - Z)
    dZ <- -X * Y + c * Y - Z
    list(c(dX, dY, dZ))
  })
}

#' Generalized Lotka-Volterra Model
#'
#' @param t Time
#' @param state Named vector of variable state
#' @param parms List of model parameters
#'
#' @return Updated state
#' @export
GLV_model <- function(t, state, parms) {
  with(as.list(c(parms)), {
    C = C0 * s
    diag(C) = diag(C0)
    dX <- r * state * (1 - C %*% state) + mu
    list(dX)
  })
}
