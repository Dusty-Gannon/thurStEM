% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranks_as_pairs.R
\name{ranks_as_pairs}
\alias{ranks_as_pairs}
\title{Convert Rankings to Pairwise Comparisons}
\usage{
ranks_as_pairs(R, cols = NULL, agg = F, scores = F)
}
\arguments{
\item{R}{A matrix or data frame of rankings. Each row represents a set of rankings,
and each column represents an item being ranked. If R is a matrix,
column names should be provided for the items.}

\item{cols}{A vector of column indices specifying which columns provide the rankings.
If \code{NULL} (default), all columns are included.}

\item{agg}{A logical value indicating whether to aggregate the comparisons across
all rows. If \code{FALSE} (default), the comparisons are provided for each
individual row.}

\item{scores}{If \code{FALSE} (default), the function assumes the data are as ranks such that
a rank of 1 is considered most favored. If \code{TRUE}, then a 1 is considered
least favorable since it got the lowest \emph{score}.}
}
\value{
If \code{agg} is TRUE, the function returns a matrix where each element
\eqn{(i, j)} represents the number of times item \eqn{i} is ranked before
item \eqn{j}. If \code{agg} is FALSE, the function returns a data frame
with pairwise comparisons for each row and any remaining columns from \code{R}.
}
\description{
This function converts a matrix of rankings into pairwise comparisons.
It can either aggregate the comparisons across all rows, giving the number
of times item \eqn{i} was ranked before item \eqn{j} across all rankings,
or provide them for each individual row as a set of binary variables.
}
\examples{
# Example data
R <- data.frame(A = c(1, 2, 3), B = c(2, 1, 1), C = c(3, 3, 2))

# Pairwise comparisons without aggregation
ranks_as_pairs(R)

# Pairwise comparisons with aggregation
ranks_as_pairs(R, agg = TRUE)

}
