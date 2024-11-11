
#' Generic summary for a `thurstem` object
#'
#' @param model_obj A fitted `thurstem` object.
#'
#' @return `lavaan`-esque summary, including sections for:
#'        \itemize{
#'          \item A summary of model fitting procedures and model fit.
#'          \item A summary of the latent variables and the factor loadings.
#'          \item A summary of the regression component of the model, modeling
#'                the intercepts of the latent variables.
#'          \item A summary of the threshold parameters, for which \eqn{1 - \Phi(\tau_\ell)},
#'                where \eqn{\Phi()} is the standard normal CDF, gives the estimated probability
#'                of comparison \eqn{\ell}, that is, item \eqn{i} ranked ahead of item \eqn{j}.
#'          \item A summary of the variances and covariances of the latent variables. Note that
#'                these will be fixed in most cases.
#'          \item Scales of \eqn{y^*}. Note that the user will not often find use for this section.
#'          \item A summary of the defined parameters. This is the section in which users will find
#'                information on marginal means of the latent utility distributions.
#'        }
#'
#' @importFrom lavaan summary
#'
#' @exportS3Method thurStEM::summary
#'
summary.thurstem <- function(model_obj){
  lavaan::summary(model_obj$mfit)
}


residuals.thurstem <- function(model_obj, type = "cor"){
  lavaan::resid(model_obj$mfit, type = type)
}




#' Get different measures of SEM fit
#'
#' @param model_obj Fitted `thurstem` model object.
#' @param ... Extra arguments to be pass on to \code{\link[lavaan]{fitMeasures}} from
#' \pkg{lavaan}.
#'
#' @return A `lavaan`-esque summary of the various SEM fit measures.
#' See \code{\link[lavaan]{fitMeasures}}.
#'
#' @importFrom lavaan fitMeasures
#'
#' @export
#'
fit_measures <- function(model_obj, ...){
  all_args <- c(
    list(object = model_obj$mfit),
    ...
  )
  do.call(lavaan::fitMeasures, all_args)
}

