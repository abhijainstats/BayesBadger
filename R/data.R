#' Simulated Individual-Level Data
#'
#' A simulated dataset of individual-level observations nested within spatial
#' clusters, intended for use with [bayes_badger()].
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{unit}{Cluster identifier (character).}
#'   \item{y}{Response variable, in (0, 1).}
#'   \item{female}{Binary indicator (0/1).}
#'   \item{inc}{Standardized continuous covariate.}
#'   \item{educ}{Factor with levels "low", "mid", "high".}
#' }
"indiv_df"

#' Simulated Cluster-Level Data
#'
#' Cluster-level covariates corresponding to the spatial units in
#' \code{\link{indiv_df}}, used for the dispersion model in [bayes_badger()].
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{unit}{Cluster identifier (character).}
#'   \item{diversity}{Standardized diversity index.}
#'   \item{ses}{Standardized socioeconomic status.}
#' }
"spatial_df"

#' Adjacency Matrix for Simulated 10x10 Grid
#'
#' A 100x100 symmetric adjacency matrix encoding the neighbour structure of
#' a 10x10 grid of spatial units. Each unit is connected to its immediate
#' horizontal and vertical neighbours.
#'
#' @format A 100x100 numeric matrix with row and column names matching the
#'   `unit` column in \code{\link{indiv_df}} and \code{\link{spatial_df}}.
"A_mat"