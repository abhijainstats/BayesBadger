#' Compute Spatial Marginal Effects for a BetaBayesSpatial Model
#'
#' Computes posterior marginal effects for each spatial unit from a fitted
#' `BetaBayesSpatial` model. For each spatial unit, the marginal effect is
#' defined as the average change in the predicted mean when that spatial
#' unit's random effect (`alpha_s`) replaces the cross-unit mean random
#' effect (`alpha_0`). This captures how much the predicted outcome shifts
#' for a spatial unit relative to the average spatial effect.
#'
#' @param fit A fitted model object of class `BetaBayesSpatial`, as returned
#'   by \code{bayes_badger()}.
#' @param method A character string specifying how to summarise individual-level
#'   marginal effects across observations within each MCMC draw. One of
#'   `"mean"` (default), `"median"`, or `"mode"`.
#' @param burnin Integer number of initial MCMC draws to discard when
#'   `method = "mode"`. Default is `0`.
#' @param thin Integer thinning interval applied when `method = "mode"`.
#'   Default is `1`.
#'
#' @return A numeric matrix of dimensions `n_draws x S`, where `n_draws` is
#'   the number of posterior draws and `S` is the number of spatial units.
#'   Column names correspond to spatial unit identifiers. Each entry
#'   `[d, s]` gives the summarized marginal effect for spatial unit `s`
#'   in draw `d`.
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
#' # Compute spatial marginal effects using posterior mean
#' sme <- marginal_effects_spatial(fit, method = "mean")
#'
#' # Posterior mean AME per spatial unit
#' colMeans(sme)
#' }
marginal_effects_spatial <- function(fit,
                                     method = c("mean", "median", "mode"),
                                     burnin = 0,
                                     thin   = 1) {
  
  method <- match.arg(method)
  
  expit <- function(x) 1 / (1 + exp(-x))
  
  mapEst <- function(chain, burnin = 0, thin = 5) {
    s   <- seq(burnin + 1, length(chain), by = thin)
    den <- density(chain[s])
    den$x[which.max(den$y)]
  }
  
  # Collapse n individuals to a scalar for a single MCMC draw
  summarise_draw <- switch(method,
                           mean   = mean,
                           median = median,
                           mode   = function(x) mapEst(x, burnin = burnin, thin = thin)
  )
  
  # --- extract draws --------------------------------------------------------
  draws       <- rstan::extract(fit$stanfit)
  beta_draws  <- draws$beta    # [n_draws x p]
  alpha_draws <- draws$alpha   # [n_draws x S]
  n_draws     <- nrow(beta_draws)
  
  X             <- fit$X              # [n x p]
  spatial_units <- fit$spatial_units  # length S, ordered as in alpha columns
  
  S <- length(spatial_units)
  n <- nrow(X)
  
  # linear predictor from covariates only (no alpha): [n x n_draws]
  Xbeta <- X %*% t(beta_draws)
  
  # cross-spatial-unit mean of alpha per draw: [n_draws] vector
  alpha0_per_draw <- rowMeans(alpha_draws)
  
  # baseline eta with alpha0: [n x n_draws]
  eta0 <- Xbeta + matrix(alpha0_per_draw,
                         nrow  = n,
                         ncol  = n_draws,
                         byrow = TRUE)
  
  # --- loop over spatial units ----------------------------------------------
  # ame_draws: [n_draws x S] — per-draw AME for each spatial unit
  ame_draws <- matrix(NA_real_, nrow = n_draws, ncol = S)
  colnames(ame_draws) <- spatial_units
  
  for (s in seq_len(S)) {
    alpha_s <- matrix(alpha_draws[, s],
                      nrow  = n,
                      ncol  = n_draws,
                      byrow = TRUE)
    eta1 <- Xbeta + alpha_s
    
    delta_mat      <- (expit(eta1) - expit(eta0))
    ame_draws[, s] <- apply(delta_mat, 2, summarise_draw)
  }
  
  ame_draws
}