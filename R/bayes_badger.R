#' Bayesian Double Generalized Beta Regression with Spatial and Temporal Borrowing
#'
#' Fits a Bayesian double generalized beta regression model with spatial effects
#' using a graph Laplacian prior and temporal borrowing incorporated through the
#' mean of the prior of alpha (spatial effects). The mean model is specified via a
#' standard R formula and the dispersion model is optionally specified after a
#' `|` separator. Sampling is performed via Stan.
#'
#' @param formula A two-sided formula specifying the mean model, and optionally
#'   the dispersion model separated by `|`. For example,
#'   `y ~ x1 + x2 | c1 + c2`. If no dispersion model is provided, an
#'   intercept-only dispersion model is used.
#' @param individual_data A data frame of individual-level observations containing
#' the response variable, mean model predictors, and the spatial ID column.
#' @param spatial_id A character string naming the column in `individual_data`
#'   (and optionally `spatial_data`) that identifies the spatial unit each
#'   observation belongs to.
#' @param spatial_data An optional data frame of spatial-unit-level covariates
#'   for the dispersion model. Required if the dispersion model includes
#'   predictors. Default is `NULL`.
#' @param adjacency_matrix A symmetric, non-negative adjacency matrix with row and column
#'   names identifying spatial units. Used to construct the graph Laplacian.
#' @param alpha_prior_mean A numeric vector of length S (number of spatial units)
#'   specifying the prior mean for the spatial random effects. Defaults to a
#'   zero vector.
#' @param beta_prior_mean Prior mean for the mean model coefficients. Default
#'   is `0`.
#' @param beta_prior_var Prior variance for the mean model coefficients. Default
#'   is `10`.
#' @param omega_prior_mean Prior mean for the dispersion model coefficients.
#'   Default is `0`.
#' @param omega_prior_var Prior variance for the dispersion model coefficients.
#'   Default is `10`.
#' @param lambda_prior_shape Shape parameter for the Gamma prior on lambda.
#'   Default is `10`.
#' @param lambda_prior_rate Rate parameter for the Gamma prior on lambda.
#'   Default is `1`.
#' @param gamma_prior_shape Shape parameter for the Gamma prior on gamma.
#'   Default is `1`.
#' @param gamma_prior_rate Rate parameter for the Gamma prior on gamma.
#'   Default is `10`.
#' @param chains Number of Markov chains. Default is `3`.
#' @param iter Total number of iterations per chain (including warmup). Default
#'   is `2000`.
#' @param warmup Number of warmup iterations per chain. Default is
#'   `floor(iter / 2)`.
#' @param thin Thinning interval for saving samples. Default is `1`.
#' @param seed Optional integer seed for reproducibility. Default is `NULL`.
#' @param ... Additional arguments passed to [rstan::stan()].
#'
#' @return An object of class `BetaBayesSpatial`, which is a list containing:
#'   \describe{
#'     \item{stanfit}{The fitted `stanfit` object returned by `rstan::stan()`.}
#'     \item{mean_terms}{Character vector of mean model predictor names.}
#'     \item{dispersion_terms}{Character vector of dispersion model predictor names.}
#'     \item{spatial_units}{Character vector of spatial unit identifiers.}
#'     \item{n_obs}{Number of observations.}
#'     \item{n_spatial_units}{Number of spatial units (S).}
#'     \item{laplacian}{The graph Laplacian matrix L.}
#'     \item{X}{The mean model design matrix.}
#'     \item{Z}{The spatial assignment matrix.}
#'   }
#'
#' @importFrom rstan stan
#' @importFrom stats model.frame model.response model.matrix terms as.formula na.fail
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
#' print(fit)
#' }
bayes_badger <- function(formula,
                         individual_data,
                         spatial_id,
                         spatial_data       = NULL,
                         adjacency_matrix,
                         alpha_prior_mean   = NULL,
                         beta_prior_mean    = 0,
                         beta_prior_var     = 10,
                         omega_prior_mean   = 0,
                         omega_prior_var    = 10,
                         lambda_prior_shape = 10,
                         lambda_prior_rate  = 1,
                         gamma_prior_shape  = 1,
                         gamma_prior_rate   = 10,
                         chains             = 3,
                         iter               = 2000,
                         warmup             = floor(iter / 2),
                         thin               = 1,
                         seed               = NULL,
                         ...) {
  library(rstan)
  .beta_spatial_stan <- "
data {
  int<lower=1> n;
  int<lower=1> p;
  int<lower=1> q;
  int<lower=1> S;
  vector<lower=0,upper=1>[n] y;
  matrix[n,p] X;
  matrix[S,q] C;
  matrix[n,S] Z;
  matrix[S,S] L;
  vector[S] alphaPriorMean;
  matrix[S,S] I;
  real betaPriorMean;
  real<lower=0> betaPriorVar;
  real omegaPriorMean;
  real<lower=0> omegaPriorVar;
  real<lower=0> lambdaPrior1;
  real<lower=0> lambdaPrior2;
  real<lower=0> gammaPrior1;
  real<lower=0> gammaPrior2;
}
parameters {
  vector[p] beta;
  vector[S] alpha;
  vector[q] omega;
  real<lower=0> lambda;
  real<lower=0> gamma;
}
transformed parameters {
  matrix[S,S] M;
  M = lambda * L + lambda * gamma * I;
  matrix[S,S] Minv;
  Minv = inverse(M);
  vector[n] eta;
  vector[S] zeta;
  eta = X * beta + Z * alpha;
  zeta = C * omega;
  vector[n] mu;
  vector[S] phi;
  vector[n] phi2;
  vector[n] A;
  vector[n] B;
  mu   = inv_logit(eta);
  phi  = exp(zeta);
  phi2 = Z * phi;
  A    = mu .* phi2;
  B    = (1.0 - mu) .* phi2;
}
model {
  beta  ~ normal(betaPriorMean,  betaPriorVar);
  omega ~ normal(omegaPriorMean, omegaPriorVar);
  alpha ~ multi_normal(alphaPriorMean, Minv);
  lambda ~ gamma(lambdaPrior1, lambdaPrior2);
  gamma  ~ gamma(gammaPrior1,  gammaPrior2);
  y ~ beta(A, B);
}
"

# 1. Parse the formula
# Accepts either:
#   y ~ x1 + x2 | c1 + c2   (mean model | dispersion model)
#   y ~ x1 + x2              (mean model only; dispersion defaults to ~1)
fchar <- paste(deparse(formula), collapse = " ")
parts <- strsplit(fchar, "\\|")[[1]]
if (length(parts) > 2)
  stop("Formula must have at most one '|', e.g.:  y ~ x1 + x2 | c1 + c2")

mean_formula       <- as.formula(trimws(parts[1]))
dispersion_formula <- if (length(parts) == 2)
  as.formula(paste("~", trimws(parts[2])))
else
  as.formula("~ 1")

# 2. Build individual-level design matrix X and response y
if (!spatial_id %in% names(individual_data))
  stop(sprintf("Column '%s' not found in `individual_data`.", spatial_id))

mf <- model.frame(mean_formula, data = individual_data, na.action = na.fail)
y  <- model.response(mf)
X  <- model.matrix(mean_formula, data = mf)[, -1, drop = FALSE]  # drop intercept (alpha absorbs it)

if (!is.numeric(y) || any(y <= 0) || any(y >= 1))
  stop("Response must be numeric with all values strictly in (0, 1).")

n <- nrow(X)

# 3. Identify spatial units and build Z matrix
A_mat         <- as.matrix(adjacency_matrix)
spatial_units <- rownames(A_mat)

if (is.null(spatial_units))
  stop("`adjacency_matrix` must have row/column names identifying spatial units.")

unit_col <- as.character(individual_data[[spatial_id]])
unknown  <- setdiff(unique(unit_col), spatial_units)
if (length(unknown) > 0)
  stop(sprintf(
    "Units in `individual_data` not found in `adjacency_matrix`: %s",
    paste(unknown, collapse = ", ")
  ))

S <- length(spatial_units)
Z <- matrix(0.0, nrow = n, ncol = S)
colnames(Z) <- spatial_units
for (i in seq_len(n)) {
  Z[i, which(spatial_units == unit_col[i])] <- 1.0
}

# 4. Build spatial-unit-level dispersion matrix C
rhs_terms <- attr(terms(dispersion_formula), "term.labels")

if (length(rhs_terms) == 0) {
  # Intercept-only dispersion model
  C <- matrix(
    1.0,
    nrow = S,
    ncol = 1,
    dimnames = list(spatial_units, "(Intercept)")
  )
} else {
  if (is.null(spatial_data))
    stop("`spatial_data` must be supplied when the dispersion model has predictors.")
  
  # Identify unit labels in spatial_data
  if (spatial_id %in% names(spatial_data)) {
    uid_col <- as.character(spatial_data[[spatial_id]])
  } else {
    uid_col <- rownames(spatial_data)
  }
  if (is.null(uid_col))
    stop(
      "Cannot identify spatial units in `spatial_data`. ",
      "Add a column named '",
      spatial_id,
      "' or set row names."
    )
  
  missing_units <- setdiff(spatial_units, uid_col)
  if (length(missing_units) > 0)
    stop(sprintf(
      "Units missing from `spatial_data`: %s",
      paste(missing_units, collapse = ", ")
    ))
  
  sud <- spatial_data[match(spatial_units, uid_col), , drop = FALSE]
  rownames(sud) <- spatial_units
  C <- model.matrix(dispersion_formula, data = sud)
}

# 5. Compute weighted graph Laplacian L = D_w - A
if (!isSymmetric(A_mat, tol = .Machine$double.eps^0.5))
  stop("`adjacency_matrix` must be symmetric.")
if (any(A_mat < 0))
  stop("`adjacency_matrix` must have non-negative entries.")
if (any(diag(A_mat) != 0))
  stop("`adjacency_matrix` must have zeros on the diagonal (no self-loops).")

D <- diag(rowSums(A_mat))
L <- D - A_mat

# 6. Alpha prior
if (is.null(alpha_prior_mean))
  alpha_prior_mean <- rep(0.0, S)
if (length(alpha_prior_mean) != S)
  stop(sprintf("`alpha_prior_mean` must have length S = %d.", S))

# 7. Assemble Stan data list
stan_data <- list(
  n              = n,
  p              = ncol(X),
  q              = ncol(C),
  S              = S,
  y              = as.numeric(y),
  X              = X,
  C              = C,
  Z              = Z,
  L              = L,
  alphaPriorMean = as.numeric(alpha_prior_mean),
  I              = diag(S),
  betaPriorMean  = beta_prior_mean,
  betaPriorVar   = beta_prior_var,
  omegaPriorMean = omega_prior_mean,
  omegaPriorVar  = omega_prior_var,
  lambdaPrior1   = lambda_prior_shape,
  lambdaPrior2   = lambda_prior_rate,
  gammaPrior1    = gamma_prior_shape,
  gammaPrior2    = gamma_prior_rate
)

# 8. Compile and sample
stan_args <- c(
  list(
    model_code = .beta_spatial_stan,
    data       = stan_data,
    pars       = c("beta", "omega", "alpha", "lambda", "gamma"),
    chains     = chains,
    iter       = iter,
    warmup     = warmup,
    thin       = thin
  ),
  list(...)
)
if (!is.null(seed))
  stan_args$seed <- seed

stanfit <- do.call(stan, stan_args)

# 9. Return results
out <- list(
  stanfit          = stanfit,
  mean_terms       = colnames(X),
  dispersion_terms = colnames(C),
  spatial_units    = spatial_units,
  n_obs            = n,
  n_spatial_units  = S,
  laplacian        = L,
  X                = X,
  # stored for marginal effects computation
  Z                = Z    # stored for marginal effects computation
)
class(out) <- "BetaBayesSpatial"
out
}

