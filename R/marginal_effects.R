#' Compute Marginal Effects for a BetaBayesSpatial Model
#'
#' Computes posterior marginal effects for one or more mean model predictors
#' from a fitted `BetaBayesSpatial` model. For continuous variables, the
#' marginal effect is the average change in the predicted mean when the
#' predictor increases by one unit. For binary/categorical variables, it is
#' the average change when the predictor switches from 0 to 1. Effects are
#' summarised across individuals using the mean, median, or mode of the
#' posterior draws.
#'
#' @param fit A fitted model object of class `BetaBayesSpatial`, as returned
#'   by \code{bayes_badger()}.
#' @param variables A character vector of predictor names for which marginal
#'   effects should be computed. Must be a subset of the mean model terms in
#'   `fit`. Defaults to all mean model predictors.
#' @param type A named character vector indicating the type of each variable
#'   in `variables`. Each element should be either `"continuous"` (default)
#'   or `"binary"`. Names must match entries in `variables`. Any variable not
#'   listed defaults to `"continuous"`.
#' @param method A character string specifying how to summarize individual-level
#'   marginal effects across observations within each MCMC draw. One of
#'   `"mean"` (default), `"median"`, or `"mode"`.
#' @param burnin Integer number of initial MCMC draws to discard when
#'   `method = "mode"`. Default is `0`.
#' @param thin Integer thinning interval applied when `method = "mode"`.
#'   Default is `1`.
#'
#' @return If a single variable is requested, a numeric vector of length
#'   `n_draws` containing the posterior marginal effect for that variable.
#'   If multiple variables are requested, a named list of such vectors, one
#'   per variable.
#'
#' @importFrom rstan extract
#' @importFrom stats density median
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- bayes_badger(
#'   formula          = y ~ x1 + x2 | c1,
#'   individual_data  = my_data,
#'   spatial_id       = "region",
#'   spatial_data     = region_data,
#'   adjacency_matrix = adj_mat
#' )
#'
#' # Marginal effect of x1 (continuous) and x2 (binary)
#' me <- marginal_effects(
#'   fit       = fit,
#'   variables = c("x1", "x2"),
#'   type      = c(x1 = "continuous", x2 = "binary"),
#'   method    = "mean"
#' )
#' }
marginal_effects <- function(fit,
                             variables = NULL,
                             type      = NULL,
                             method    = c("mean", "median", "mode"),
                             burnin    = 0,
                             thin      = 1) {
  
  method <- match.arg(method)
  
  expit <- function(x) 1 / (1 + exp(-x))
  
  mapEst <- function(chain, burnin = 0, thin = 5) {
    s   <- seq(burnin + 1, length(chain), by = thin)
    den <- density(chain[s])
    den$x[which.max(den$y)]
  }
  
  # Summarise a vector of individual-level deltas (one MCMC draw) to a scalar
  summarise_draw <- switch(method,
                           mean   = mean,
                           median = median,
                           mode   = function(col) mapEst(col, burnin = burnin, thin = thin)
  )
  
  # --- extract draws --------------------------------------------------------
  draws       <- rstan::extract(fit$stanfit)
  beta_draws  <- draws$beta    # [n_draws x p]
  alpha_draws <- draws$alpha   # [n_draws x S]
  n_draws     <- nrow(beta_draws)
  
  X          <- fit$X          # [n x p], no intercept
  Z          <- fit$Z          # [n x S]
  mean_terms <- fit$mean_terms
  p          <- ncol(X)
  n          <- nrow(X)
  
  # alpha contribution per draw: [n x n_draws]
  Zalpha <- Z %*% t(alpha_draws)
  
  # --- defaults -------------------------------------------------------------
  if (is.null(variables)) variables <- mean_terms
  
  bad <- setdiff(variables, mean_terms)
  if (length(bad) > 0)
    stop(sprintf("Variable(s) not in mean model: %s", paste(bad, collapse = ", ")))
  
  if (is.null(type)) {
    type <- setNames(rep("continuous", length(variables)), variables)
  } else {
    for (v in variables)
      if (!v %in% names(type)) type[v] <- "continuous"
  }
  
  # --- loop over variables --------------------------------------------------
  result <- vector("list", length(variables))
  names(result) <- variables
  
  for (v in variables) {
    j     <- which(mean_terms == v)
    vtype <- type[v]
    
    # linear predictor without variable j: [n x n_draws]
    if (p == 1) {
      eta_rest <- Zalpha
    } else {
      eta_rest <- X[, -j, drop = FALSE] %*% t(beta_draws[, -j, drop = FALSE]) + Zalpha
    }
    
    beta_j <- matrix(beta_draws[, j], nrow = 1)  # [1 x n_draws]
    
    if (vtype == "continuous") {
      xj   <- X[, j, drop = FALSE]
      eta1 <- eta_rest + (xj + 1) %*% beta_j     # [n x n_draws]
      eta0 <- eta_rest + xj       %*% beta_j
    } else {
      # binary / categorical: x_j = 1 for all vs x_j = 0 for all
      eta1 <- eta_rest + matrix(1, n, 1) %*% beta_j
      eta0 <- eta_rest + matrix(0, n, 1) %*% beta_j
    }
    
    # individual-level deltas: [n x n_draws]
    delta_mat <- (expit(eta1) - expit(eta0))
    
    # Per draw: collapse over individuals using the chosen method -> [n_draws]
    result[[v]] <- apply(delta_mat, 2, summarise_draw)
  }
  
  # Return a named vector if a single variable, otherwise a named list of vectors
  if (length(variables) == 1) result[[1]] else result
}