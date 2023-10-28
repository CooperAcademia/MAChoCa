#' Sample new responses from the parameter vector and single option data
#'
#' Given a parameter vector `x` and a dtaset `data` generate response
#' probabilities for the rows in the data for single option trials.
#'
#' The function \link{ll_single} explains the data format in more detail
#'
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @importFrom stats rbinom
#' @export
sample_single <- function(x, data) {
  resp_p <- rp_single(x, data)
  data$resp_p <- resp_p
  data$accept <- rbinom(n=nrow(data), size=1, prob=resp_p) |> as.logical()
  return(data)
}

