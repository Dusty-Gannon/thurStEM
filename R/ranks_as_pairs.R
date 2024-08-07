
#' Convert Rankings to Pairwise Comparisons
#'
#' This function converts a matrix of rankings into pairwise comparisons.
#' It can either aggregate the comparisons across all rows, giving the number
#' of times item \eqn{i} was ranked before item \eqn{j} across all rankings,
#' or provide them for each individual row as a set of binary variables.
#'
#' @param R A matrix or data frame of rankings. Each row represents a set of rankings,
#'          and each column represents an item being ranked. If R is a matrix,
#'          column names should be provided for the items.
#' @param cols A vector of column indices specifying which columns provide the rankings.
#'             If NULL (default), all columns are included.
#' @param agg A logical value indicating whether to aggregate the comparisons across
#'            all rows. If FALSE (default), the comparisons are provided for each
#'            individual row.
#'
#' @return If \code{agg} is TRUE, the function returns a matrix where each element
#'         \eqn{(i, j)} represents the number of times item \eqn{i} is ranked before
#'         item \eqn{j}. If \code{agg} is FALSE, the function returns a data frame
#'         with pairwise comparisons for each row and any remaining columns from \code{R}.
#'
#' @examples
#' # Example data
#' R <- data.frame(A = c(1, 2, 3), B = c(2, 1, 1), C = c(3, 3, 2))
#'
#' # Pairwise comparisons without aggregation
#' ranks_as_pairs(R)
#'
#' # Pairwise comparisons with aggregation
#' ranks_as_pairs(R, agg = TRUE)
#'
#' @export
ranks_as_pairs <- function(R, cols = NULL, agg = F){
  if(is.null(cols)){
    cols <- 1:ncol(R)
  }
  m <- length(cols)
  n <- nrow(R)
  items <- colnames(R)[cols]
  if(is.null(items)){
    stop("Please name the items being ranked by giving R column names.\n")
  }
  if(agg){
    R_sub <- R[, cols]
    P <- matrix(nrow = m, ncol = m)
    for(i in 1:m){
      for(j in 1:m){
        P[i, j] <- sum(R_sub[, i] < R_sub[, j])
      }
    }
    return(P)
  } else{
    P <- matrix(data = 0, nrow = n, ncol = choose(m, 2))
    pairnames <- vector(mode = "character")
    indexes <- vector(mode = "list")
    for(i in 1:(m - 1)){
      for(j in (i + 1):m){
        pairnames <- c(pairnames, paste(items[i], "before", items[j], sep = "_"))
        indexes <- c(indexes, list(c(i, j)))
      }
    }
    colnames(P) <- pairnames
    R_sub <- R[, cols]
    for(i in 1:n){
      for(j in 1:choose(m, 2)){
        P[i, j] <- as.numeric(R_sub[i, indexes[[j]][1]] < R_sub[i, indexes[[j]][2]])
      }
    }
    return(
      as.data.frame(cbind(P, R[-cols]))
    )
  }
}

