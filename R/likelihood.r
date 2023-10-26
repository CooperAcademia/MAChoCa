#' Calculate likelihood of x given the data and return it
#'
#' This function transforms parameters from the real number line to
#' the appropriate regions using `transform_pars`. The using a variation
#' of the generalised context model for categorisation it calculates
#' the likelihood of the parameter vector x at leading to the accept/reject
#' responses to the prices and ratings in the data.
#'
#' If `sample = TRUE` it will instead replace the accept column of the data with
#' values sampled from the model using the provided parameter values x
#'
#' The `x` parameter vector is expected to be a named vector containing values
#' for `w` weight on the price attribute, `r` form of the distance metric,
#' `s` sensitivity to distances in attreibute space and `delta` a choice offset
#' relative to the anchor.
#'
#' The data object is assumed to be a data.frame compatible object with at
#' least three columns, a `price_n` column with values for the price attribute
#' normalised to 0-1 and reversed so that 0 is the "worst" ie highest price and
#' 1 is the best price. It should also have a rating_n column with ratings
#' normalised to 0-1 (worst-best) and finally an `accept` column with T/F for
#' whether the option was accepted.
#'
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @param sample A boolean which specifies whether to return the likelihood, or
#'   sample new responses.
#'
#' @return A single value (log-likelihood) or a copy of `data` with new
#'   responses
#' @export
ll_pmwg <- function(x, data, sample=FALSE) {
  names(x) <- par_names
  x <- transform_pars(x, fwd=FALSE)
  d <- ((x["w"] * (data$price_n - a[1])^x["r"]) +
        ((1 - x["w"]) * (data$rating_n - a[2])^x["r"]))^(1/x["r"])
  resp_p <- pnorm(x["delta"] + x["s"] * d)
  # Make sure not too close to 1
  resp_p <- pmin(resp_p, 1 - 1e-10)
  # Make sure not too close to 0
  resp_p <- pmax(resp_p, 1e-10)

  if(!sample) {
    ll <- sum(log(resp_p[data$accept]), log((1 - resp_p)[!data$accept]))
    if (is.nan(ll)) {
      return(-Inf)
    }
    return(ll)
  } else {
    data$resp_p <- resp_p
    data$accept <- rbinom(n=nrow(data), size=1, prob=resp_p) %>% as.logical()
    return(data)
  }
}

#' Likelihood for single option, but minimises.
#'
#' See description for ll_pmwg. This is a thin wrapper that negates the output
#' for use with the R built-in optim.
#' 
#' @inheritParams ll_pmwg
#'
#' @return negated likelihood
#' @export
ll_optim <- function(x, data) {
  return -ll_pmwg(x, data, sample=FALSE)
}

#' Transform the parameter vector for use in pmwg
#'
#' See also the description for ll_pmwg. This is a helper function that
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
#' @export
transform_pars <- function(x, fwd=TRUE) {
  if(fwd) {
    x[c("r","s")] <- log(x[c("r","s")])
    x["w"] <- qnorm(x["w"])
  }  else {
    x[c("r","s")] <- exp(x[c("r","s")])
    x["w"] <- pnorm(x["w"])
  }
  x
}
