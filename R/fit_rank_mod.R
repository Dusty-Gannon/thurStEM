
#' Fit a Thurstone Rank Model
#'
#' This function fits a Thurstone rank model to a dataset.
#' The function allows for optional inclusion of covariates and handles the conversion of pairwise
#' ranking data into a form suitable for fitting using the `lavaan` package.
#'
#' @param data A data frame containing the pairs data and any covariates.
#' @param pairs_cols A vector of column indices or names specifying the pairwise comparison columns in the data.
#'                   If \code{NULL}, all columns are assumed to be ranks or scores.
#' @param form An optional formula specifying the covariates to include in the model. If \code{NULL}, then the
#'             *intercept only* model will be fitted.
#' @param mclass A character string specifying the Thurstone model class to fit. Currently, only "V" is implemented.
#'
#' @return A list of class \code{"thurstem"} containing:
#'   \item{mfit}{The fitted lavaan model object.}
#'   \item{emm_grid}{The model matrix for constructing the marginal means.}
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

  # define the regressions block of model code
  if(is.null(form)){
    reg_mod <- vector("character", length = K)
    for(i in 1:(K-1)){
      reg_mod[i] <- paste0("i", i, " ~ mu", i, " * 1")
    }
    reg_mod[K] <- paste0("i", K, " ~ 0 * 1 + mu", K, " * 1")
    reg_mod <- paste(reg_mod, collapse = "\n")
  } else{
    X_names <- colnames(X)[-1]
    reg_mod <- vector("character", length = K + length(X_names) - 1)
    counter <- 1
    for(i in 1:(K - 1)){
      intercept_i <- paste0("i", i, " ~ mu", i, " * 1")
        reg_ij <- ""
        for(j in 1:(ncol(X) - 1)){
          reg_ij <- paste0(reg_ij, " + b", i, j, " * ", X_names[j])
        }
        reg_mod[counter] <- paste0(intercept_i, reg_ij)
        counter <- counter + 1
    }
    reg_mod[length(reg_mod)] <- paste0("i", K, " ~ 0 * 1 + mu", K, " * 1")
  }
  reg_mod <- paste(reg_mod, collapse = "\n")


  # define loadings block of model code
  meas <- vector("character", length = K)
  for(i in 1:K){
    rhs <- paste(A[,i], indic, sep = "*")
    meas[i] <- paste0("i", i, " =~ ", paste(rhs, collapse = " + "))
  }
  meas <- paste(meas, collapse = "\n")


  # now fix factor variances and covariances
  Sigma_t <- diag(1, K)
  Sigma_d <- A %*% Sigma_t %*% t(A)
  covar <- vector("character", length = K + m)
  for(i in 1:K){
    lhs_i <- paste0("i", i, " ~~ ", Sigma_t[i, 1], " * i1")
    rhs_i <- ""
    if(i > 1){
      for(j in 2:i){
        rhs_i <- paste(rhs_i, paste0(Sigma_t[i,j], " * i", j), sep = " + ")
      }
    }
    covar[i] <- paste0(lhs_i, rhs_i)
  }
  for(i in (K + 1):length(covar)){
    lhs_i <- paste0(indic[i - K], " ~~ ", Sigma_d[(i - K), 1], " * ", indic[1])
    rhs_i <- ""
    if(i > (K + 1)){
      for(j in 2:(i - K)){
        rhs_i <- paste(rhs_i, paste0(Sigma_d[(i - K), j], " * ", indic[j]), sep = " + ")
      }
    }
    covar[i] <- paste0(lhs_i, rhs_i)
  }
  covar <- paste(covar, collapse = "\n")


  # now create the derived portion of the model
  # get the names of the factors
  fctrs_dat <- lapply(data, is.factor) |> unlist()
  fctrs_dat <- names(fctrs_dat)[fctrs_dat]
  if(is.null(form)){
    X_new <- matrix(1, nrow = 1, ncol = 1)
    derived <- vector("character", length = K)
    for(i in 1:(K - 1)){
      derived[i] <- paste0("mu", i, " := sqrt(2) * d", i, K)
    }
    derived[K] <- paste0("mu", K, " := 0")
  } else if(length(fctrs_dat) > 0){
    # first, find which columns of X are factors
    fctrs <- c(
      1,
      apply(as.matrix(X[,-1]), 2, function(x){!any(x != 1 & x != 0)}) |>
        which() + 1
    )

    # now find unique combinations of those along with medians of the continuous
    # variables
    X_new <- X

    if(any(!(1:ncol(X) %in% fctrs))){
      contins <- which(!(1:ncol(X) %in% fctrs))
      X_new[, contins] <- apply(as.matrix(X[, contins]), 2, median)
    }
    X_new <- unique(X_new)

    # make some rownames to track things
    data_grps <- as.data.frame(unique(data[, fctrs_dat]))
    names(data_grps) <- fctrs_dat
    rnames <- vector("character", nrow(data_grps))
    for(i in 1:nrow(data_grps)){
      rname_i <- vector("character", length = length(fctrs_dat))
      for(j in 1:length(fctrs_dat)){
        rname_i[j] <- paste0(fctrs_dat[j], data_grps[i, j])
      }
      rnames[i] <- paste(rname_i, collapse = "_")
    }

    # now construct derived variables
    derived <- vector("character", length = (K - 1) * nrow(X_new))

    for(i in 1:nrow(X_new)){
      for(k in 1:(K - 1)){
        intercept_ik <- paste0(
          rnames[i], "_mu", k, " := sqrt(2) * mu", k
        )
        rhs_ik <- ""
        for(l in 1:(ncol(X) - 1)){
          rhs_ik <- paste0(rhs_ik, " + ", X_new[i, (l + 1)], " * b", k, l)
        }
        derived[(K-1)*(i - 1) + k] <- paste0(intercept_ik, rhs_ik)
      }
    }
  } else {
    X_new <- NULL
    derived <- ""
  }
  derived <- paste(derived, collapse = "\n")


  P <- as.data.frame(P)
  if(is.null(form)){
    lav_dat <- P
  } else{
    lav_dat <- cbind(P, X[,-1])
    names(lav_dat) <- c(colnames(P), X_names)
  }
  # fit the model
  mfit <- lavaan::lavaan(
    model = c(reg_mod, meas, covar, derived),
    data = lav_dat,
    ordered = names(lav_dat)[1:ncol(P)],
    parameterization = "theta",
    meanstructure = TRUE
  )

  ret <- list(
    mfit = mfit,
    emm_grid = X_new,
    mclass = mclass,
    model = list(
      model_code = c(reg_mod, meas, covar, derived),
      Sigma_t = Sigma_t,
      Sigma_d = Sigma_d,
      A = A
    ),
    K = K, m = m
  )

  class(ret) <- "thurstem"

  return(ret)

}


