
#' Fit a Thurstone Rank Model
#'
#' This function fits a Thurstone rank model to a dataset.
#' The function allows for optional inclusion of covariates and handles the conversion of pairwise
#' ranking data into a form suitable for fitting using the `lavaan` package.
#'
#' @param data A data frame containing the pairs data and any covariates.
#' @param pairs_cols A vector of column indices or names specifying the pairwise comparison columns in the data.
#'                   If \code{NULL}, all columns are used.
#' @param form An optional formula specifying the covariates to include in the model. If \code{NULL}, then the
#'             *intercept only* model will be fitted.
#' @param mclass A character string specifying the Thurstone model class to fit. Currently, only "V" is implemented.
#'
#' @return A list of class \code{"thurstem"} containing:
#'   \item{mfit}{The fitted lavaan model object.}
#'   \item{df_pvals}{A data frame of chi-square test statistics, degrees of freedom, and p-values for each covariate.}
#'   \item{mclass}{The specified Thurstone model class.}
#'
#'
#' @importFrom lavaan lavaan coef
#' @export
fit_rank_mod <- function(data, pairs_cols = NULL, form = NULL, mclass = "V"){

  if(!(mclass == "V")){
    stop("Sorry! Still working on implementing the other Thurstone models.
    Only Class V available right now.\n")
  }
  if(is.null(pairs_cols)){
    P <- data
  } else{
    P <- data[, pairs_cols]
  }

  m <- ncol(P)
  # solve for K based on the number of pairs
  K <- (1 + c(-1, 1) * (sqrt(1 + 8 * m))) / 2
  K <- K[which(K > 0)]

  # define model matrix, if necessary
  if(!is.null(form)){
    X <- stats::model.matrix(form, data = data)
  }

  # building model code
  indic <- colnames(P)

  # make the "A" matrix of contrasts
  A <- matrix(data = 0, nrow = 1, ncol = K)
  for(i in 1:(K - 1)){
    x_i <- matrix(0, nrow = (K - i), ncol = K)
    x_i[, i] <- 1
    for(j in (i+1):K){
      x_i[j - i, j] <- -1
    }
    A <- rbind(A, x_i)
  }
  # remove the initialization row
  A <- A[-1, ]

  # define loadings block of model code
  loadings <- vector("character", length = K)
  for(i in 1:K){
    rhs <- paste(A[,i], indic, sep = "*")
    loadings[i] <- paste0("factor_", i, " =~ ", paste(rhs, collapse = " + "))
  }
  loadings <- paste(loadings, collapse = "\n")

  # define latent model
  if(is.null(form)){
    lat_mod <- vector("character", length = K)
    for(i in 1:(K-1)){
      lat_mod[i] <- paste0("factor_", i, " ~ 1")
    }
    lat_mod[K] <- paste0("factor_", K, " ~ 0 * 1")
    lat_mod <- paste(lat_mod, collapse = "\n")
  } else{
    X_names <- colnames(X)[-1]
    lat_mod <- vector("character", length = K)
    for(i in 1:(K-1)){
      lat_mod[i] <- paste0("factor_", i, " ~ 1 + ", paste(X_names, collapse = " + "))
    }
    lat_mod[K] <- paste0("factor_", K, " ~ 0 * 1")
    lat_mod <- paste(lat_mod, collapse = "\n")
  }


  # now fix factor variances and covariances
  factor_var <- vector("character", length = (K * (K - 1) / 2))
  counter <- 1
  for(i in 1:K){
    for(j in i:K){
      if(i == j){
        factor_var[counter] <- paste0("factor_", i, " ~~ 1*factor_", j)
      } else{
        factor_var[counter] <- paste0("factor_", i, " ~~ 0*factor_", j)
      }
      counter <- counter + 1
    }
  }
  factor_var <- paste(factor_var, collapse = "\n")

  P <- as.data.frame(P)
  P[,1:K] <- apply(P[, 1:K], 2, ordered)
  if(is.null(form)){
    lav_dat <- P
  } else{
    lav_dat <- cbind(P, X[,-1])
    names(lav_dat) <- c(colnames(P), X_names)
  }
  # fit the model
  mfit <- lavaan::lavaan(
    model = c(loadings, lat_mod, factor_var),
    data = lav_dat,
    int.lv.free = T,
    ordered = colnames(P),
    parameterization = "theta",
    meanstructure = T
  )

  # p-values for each variable in X
  coefs <- coef(mfit)
  V <- mfit@vcov$vcov
  L <- rep(0, length(coefs))
  L[grep("~1", names(coefs))] <- 1
  chisq_int <- as.double(
    t(L * coefs) %*% solve(V) %*% (L * coefs)
  )
  df_pvals <- data.frame(
    chisq = chisq_int,
    df = sum(L),
    p_val = 1 - stats::pchisq(chisq_int, sum(L))
  )
  if(!is.null(form)){
    for(j in 2:ncol(X)){
      L <- rep(0, length(coefs))
      L[grep(colnames(X)[j], names(coefs))] <- 1
      chisq <- as.double(
        t(L * coefs) %*% solve(V) %*% (L * coefs)
      )
      df_pvals <- rbind(
        df_pvals,
        c(chisq, sum(L), 1 - stats::pchisq(chisq, sum(L)))
      )
    }
  }
  rownames(df_pvals) <- colnames(X)

  ret <- list(
    mfit = mfit,
    df_pvals = df_pvals,
    mclass = mclass
  )

  class(ret) <- "thurstem"

  return(ret)

}
