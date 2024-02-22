#' Calculate likelihood of x given single option data
#'
#' This function transforms parameters from the real number line to
#' the appropriate regions using `transform_pars`. Then using a variation
#' of the generalised context model for categorisation it calculates
#' the likelihood of the parameter vector x at leading to the accept/reject
#' responses to the prices and ratings in the single-option data.
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
#' @param rp_func A function for calculating response probabilities for choosing
#'   the left presented item.
#' @param ... Other parameters passed to `rp_func`
#'
#' @return A single value (log-likelihood)
#' @export
ll_single <- function(x, data, rp_func = rp_single, ...) {
  resp_p <- tryCatch(rp_func(x, data, ...),
           error = function(e) {
             if (e$message == "UNLIKELY") {
               return(-Inf)
             }
             stop(e)
           })
  if (any(is.infinite(resp_p))) {
    return(-1e10)
  }

  ll <- sum(log(resp_p[data$accept]), log((1 - resp_p)[!data$accept]))
  if (is.nan(ll)) {
    return(-Inf)
  }
  ll
}

#' Calculate likelihood of x given double option data (N attrs)
#'
#' This function transforms parameters from the real number line to
#' the appropriate regions using `transform_pars_dbl`. It then uses a
#' sensitivity adjusted Luce Choice Axiom rule and distances from the
#' exemplar as from the single option models to calculate the likelihood
#' of the parameter vector x given the data.
#'
#' The data object is assumed to be a data.frame compatible object with at
#' least columns for each attribute used in the `rp_func` and additionally
#' have a `response` column with T/F for whether the option selected was on the
#' left.
#'
#' @param x A named vector of parameters. Expects certain values as in
#'   the description for the `rp_func`
#' @param data A data.frame compatible object with specific rows as described.
#' @param rp_func A function for calculating response probabilities for choosing
#'   the left presented item.
#' @param ... Other parameters passed to `rp_func`
#'
#' @return A single value (log-likelihood)
#' @export
ll_double <- function(x, data, rp_func = rp_double, ...) {
  resp_p <- tryCatch(rp_func(x, data, ...),
           error = function(e) {
             if (e$message == "UNLIKELY") {
               return(-Inf)
             }
             stop(e)
           })
  if (any(is.infinite(resp_p))) {
    return(-1e10)
  }
  ll <- sum(log(resp_p[data$response]), log((1 - resp_p)[!data$response]))
  if (is.nan(ll)) {
    return(-1e10)
  }
  ll
}
