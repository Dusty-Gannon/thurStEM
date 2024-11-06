
#' Generic summary for `thurstem` object
#'
#' @param model_obj Fitted `thurstem` object
#'
#' @return `lavaan`-esque summary
#' @exportS3Method thurStEM::summary
#'
summary.thurstem <- function(model_obj){
  lavaan::summary(model_obj$mfit)
}
