
#' Transform the parameter vector for use in pmwg
#'
#' See also the description for ll_single. This is a helper function that
#' transforms the parameters. To transform for pmwg the `r` and `s` values
#' are exponentiated in order to move to positive only and `w` is pnormed so
#' that values are between 0 and 1.
#'
#' These operation are reversed when `fwd=TRUE`
#'
#' @param x The named vector of parameter estimates
#' @param fwd Move certain parameter to real number line or back
#'
#' @return transformed parameter vector x
#' @importFrom stats qnorm pnorm
#' @export
transform_pars <- function(x, fwd=TRUE) {
  if(fwd) {
    x[c("r","s")] <- log(x[c("r","s")])
    x["w"] <- qnorm(x["w"])
  }  else {
    x[c("r","s")] <- exp(x[c("r","s")])
    x["r"] <- max(x["r"], 1e-10)  # Make sure not 0
    x["s"] <- max(x["s"], 1e-10)  # Make sure not 0
    x["w"] <- pnorm(x["w"])
  }
  x
}

#' Transform the parameter vector (dbl option, dbl attr) for use in pmwg
#'
#' See also the description for ll_double. This is a helper function that
#' transforms the parameters. To transform for pmwg the `r` and `gamma` values
#' are exponentiated in order to move to positive only and `w` is pnormed so
#' that values are between 0 and 1.
#'
#' These operation are reversed when `fwd=TRUE`
#'
#' @param x The named vector of parameter estimates
#' @param fwd Move certain parameter to real number line or back
#'
#' @return transformed parameter vector x
#' @importFrom stats qnorm pnorm
#' @export
transform_pars_dbl <- function(x, fwd=TRUE) {
  if(fwd) {
    x[c("r","gamma")] <- log(x[c("r","gamma")])
    x["w"] <- qnorm(x["w"])
  }  else {
    x[c("r","gamma")] <- exp(x[c("r","gamma")])
    x["r"] <- max(x["r"], 1e-10)  # Make sure not 0
    x["gamma"] <- max(x["gamma"], 1e-10)  # Make sure not 0
    x["w"] <- pnorm(x["w"])
  }
  x
}

#' Transform the parameter vector (2 x N) cut method for use in pmwg
#'
#' See also the description for rp_cut This is a helper function that
#' transforms the parameters. To transform for pmwg the `r` and `gamma` values
#' are exponentiated in order to move to positive only and cut points `c1`..`cn`
#' are pnormed so that values are between 0 and 1.
#'
#' These operation are reversed when `fwd=TRUE`
#'
#' @param x The named vector of parameter estimates
#' @param fwd Move certain parameter to real number line or back
#'
#' @return transformed parameter vector x
#' @importFrom stats qnorm pnorm
#' @export
transform_pars_cut <- function(x, fwd=TRUE) {
  cuts <- grep("^c", names(x))
  if(fwd) {
    x[c("r","gamma")] <- log(x[c("r","gamma")])
    x[cuts] <- qnorm(x[cuts])
    return(x)
  }
  x[c("r","gamma")] <- exp(x[c("r","gamma")])
  x["r"] <- max(x["r"], 1e-10)  # Make sure not 0
  x["gamma"] <- max(x["gamma"], 1e-10)  # Make sure not 0
  x[cuts] <- pnorm(x[cuts])

	cuts <- x[grep("^c", names(x))]

	weights <- c()
	filled_proportion <- 0
	for (cut_pt in cuts) {
		nxt_weight <- cut_pt * (1-filled_proportion)
		weights <- c(weights, nxt_weight)
		filled_proportion <- filled_proportion + nxt_weight
	}
	weights <- c(weights, 1 - sum(weights))
  names(weights) <- paste0("w", seq_along(weights))
  c(x, weights)
}

#' Transform the parameter vector (2 x N) dirichlet method for use in pmwg
#'
#' See also the description for rp_dir This is a helper function that
#' transforms the parameters. To transform for pmwg all parameter values
#' are exponentiated in order to move to positive only.
#'
#' These operation are reversed when `fwd=TRUE`
#'
#' @param x The named vector of parameter estimates
#' @param fwd Move certain parameter to real number line or back
#'
#' @return transformed parameter vector x
#' @export
transform_pars_dir <- function(x, fwd=TRUE) {
  if(fwd) {
    x <- log(x)
    return(x)
  }
  x <- exp(x)
  x <- pmax(x, 1e-10)  # Make sure not 0
  # Enforce alphas to be greater than 0.01 and less than 5000
	alpha_indices <- grep("^alpha", names(x))
  if (any(x[alpha_indices] > 5e3) | any(x[alpha_indices] < 0.01)) {
    stop("UNLIKELY")
  }
  x
}

#' Transform the parameter vector for use in pmwg
#'
#' See also the description for ll_single. This is a helper function that
#' transforms the parameters. To transform for pmwg the `r` and `s` values
#' and any weight values are exponentiated in order to move to positive only.
#'
#' These operation are reversed when `fwd=TRUE`
#'
#' @param x The named vector of parameter estimates
#' @param fwd Move certain parameter to real number line or back
#'
#' @return transformed parameter vector x
#' @export
transform_pars_weights <- function(x, fwd=TRUE) {
  weights <- grep("^w", names(x))
  others <- na.omit(match(c("r", "s", "delta", "gamma"), names(x)))
  if(fwd) {
    x[c(weights, others)] <- log(x[c(weights, others)])
  }  else {
    x[c(weights, others)] <- exp(x[c(weights, others)])
    x["r"] <- max(x["r"], 1e-10)  # Make sure not 0
  }
  x
}
