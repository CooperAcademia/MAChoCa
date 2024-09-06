
#' Calculate response probabilities for single option double attribute data
#'
#' Used for both calculating likelihoods or sampling from a parameter estimate
#' This function calculates the probabilities of responses to a dataset
#' based on the provided parameters (in the `x` vector).
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @param attr1 The name of the column containing normalised values for the
#'   first attribute (which has the weight value applied to it)
#' @param attr2 The name of the column containing normalised values for the
#'   second attribute(which has the (1 - weight) values applied to it
#' @importFrom stats pnorm
#' @export
rp_single <- function(x, data, attr1 = "price_n", attr2 = "rating_n") {
  x <- transform_pars(x, fwd=FALSE)
  d <- distance(c(x["w"], 1 - x["w"]),
                cbind(data[attr1], data[attr2]),
                c(0, 0),
                x["r"])
  resp_p <- pnorm(x["delta"] + x["s"] * d)
  # Make sure not too close to 1
  resp_p <- pmin(resp_p, 1 - 1e-10)
  # Make sure not too close to 0
  resp_p <- pmax(resp_p, 1e-10)
  resp_p
}

#' Calculate response probabilities for single option double attribute data
#'
#' Used for both calculating likelihoods or sampling from a parameter estimate
#' This function calculates the probabilities of responses to a dataset
#' based on the provided parameters (in the `x` vector).
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @param attrs The column names containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 1
#' @importFrom stats pnorm
#' @export
rp_single_weight <- function(x, data, attrs = c("price_n", "rating_n")) {
  x <- transform_pars_weights(x, fwd=FALSE)
  x["w2"] <- pmax(pmin(x["w2"], 1e3), 1e-3)
  w <- c(1, x['w2'])
  w <- w/sum(w)
  n_attr <- length(attrs)
  d <- distance(w,
                data.matrix(data[attrs]),
                rep(0, n_attr),
                x["r"])
  resp_p <- pnorm(x["delta"] + x["s"] * d)
  # Make sure not too close to 1
  resp_p <- pmin(resp_p, 1 - 1e-10)
  # Make sure not too close to 0
  resp_p <- pmax(resp_p, 1e-10)
  resp_p
}

#' Calculate response probabilities using Dirichlet for single option
#'
#' Used for both calculating likelihoods or sampling from a parameter estimate
#' This function calculates the probabilities of responses to a dataset
#' based on the provided parameters (in the `x` vector).
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @param attrs The column names containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 1
#' @importFrom stats pnorm
#' @export
rp_single_dir <- function(x, data, attrs = c("price_n", "rating_n")) {
  x <- transform_pars_dir(x, fwd=FALSE)
  w <- rdirichlet(nrow(data), x[grep("^alpha", names(x))])
  n_attr <- length(attrs)
  d <- distance(w,
                data.matrix(data[attrs]),
                rep(0, n_attr),
                x["r"])
  resp_p <- pnorm(x["delta"] + x["s"] * d)
  # Make sure not too close to 1
  resp_p <- pmin(resp_p, 1 - 1e-10)
  # Make sure not too close to 0
  resp_p <- pmax(resp_p, 1e-10)
  resp_p
}


#' Calculate response probabilities for double option double attribute data
#'
#' Used for both calculating likelihoods or sampling from a parameter estimate
#' This function calculates the probabilities of responses to a dataset
#' based on the provided parameters (in the `x` vector).
#'
#' The `x` parameter vector is expected to be a named vector containing values
#' for `w` weight on the price attribute, `r` form of the distance metric and
#' `gamma` sensitivity to differences in attribute space.
#'
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @param loc1_attr1 The name of the column containing normalised values for the
#'   first attribute (which has the weight value applied to it) in location 1
#' @param loc1_attr2 The name of the column containing normalised values for the
#'   second attribute(which has the (1 - weight) values applied to it in
#'   location 1
#' @param loc2_attr1 The name of the column containing normalised values for the
#'   first attribute (which has the weight value applied to it) in location 2
#' @param loc2_attr2 The name of the column containing normalised values for the
#'   second attribute(which has the (1 - weight) values applied to it in
#'   location 2
#' @importFrom stats pnorm
#' @export
rp_double <- function(x, data,
                      loc1_attr1 = "left_price", loc1_attr2 = "left_memory",
                      loc2_attr1 = "right_price", loc2_attr2 = "right_memory") {
  x <- transform_pars_dbl(x, fwd=FALSE)
  da <- distance(c(x["w"], 1 - x["w"]),
                 cbind(data[loc1_attr1], data[loc1_attr2]),
                 c(0, 0),
                 x["r"])
  db <- distance(c(x["w"], 1 - x["w"]),
                 cbind(data[loc2_attr1], data[loc2_attr2]),
                 c(0, 0),
                 x["r"])
  resp_p <- da^x["gamma"] / (da^x["gamma"] + db^x["gamma"])
#  resp_p <- pnorm(x["delta"] + x["s"] * d)
  # Make sure not too close to 1
  resp_p <- pmin(resp_p, 1 - 1e-10)
  # Make sure not too close to 0
  resp_p <- pmax(resp_p, 1e-10)
  resp_p
}

#' Calculate response probabilities for 2xN using the cut_remain method
#'
#' Used for both calculating likelihoods or sampling from a parameter estimate
#' This function calculates the probabilities of responses to a dataset
#' based on the provided parameters (in the `x` vector).
#'
#' The `x` parameter vector is expected to be a named vector containing values
#' for `\eqn{c_1..c_{N-1}}` where each \eqn{c_i} corresponds to the cut point
#' for each attribute \eqn{1..(N-1)}. The cut points are transformed to
#' attribute weights by cutting the remaining proportion of the 0-1 scale for
#' each new cut point, and the last weight being the remaining proportion.
#' The parameter vector `x` also contains `r` form of the distance metric and
#' `gamma` sensitivity to differences in attribute space.
#'
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @param loc1_attrs The column names containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 1
#' @param loc2_attrs The name of the column containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 2
#' @importFrom stats pnorm
#' @export
rp_cut <- function(x, data, loc1_attrs = c("left_price", "left_memory"),
                      loc2_attrs = c("right_price", "right_memory")) {
  x <- transform_pars_cut(x, fwd=FALSE)
  stopifnot(length(loc1_attrs) == length(loc2_attrs))
  n_attr <- length(loc1_attrs)
  da <- distance(x[grep("^w", names(x))],
                 data.matrix(data[loc1_attrs]),
                 rep(0, n_attr),
                 x["r"])
  db <- distance(x[grep("^w", names(x))],
                 data.matrix(data[loc2_attrs]),
                 rep(0, n_attr),
                 x["r"])
  resp_p <- da^x["gamma"] / (da^x["gamma"] + db^x["gamma"])
#  resp_p <- pnorm(x["delta"] + x["s"] * d)
  # Make sure not too close to 1
  resp_p <- pmin(resp_p, 1 - 1e-10)
  # Make sure not too close to 0
  resp_p <- pmax(resp_p, 1e-10)
  resp_p
}

#' Calculate response probabilities for 2xN using the dirichlet method
#'
#' Used for both calculating likelihoods or sampling from a parameter estimate
#' This function calculates the probabilities of responses to a dataset
#' based on the provided parameters (in the `x` vector).
#'
#' The `x` parameter vector is expected to be a named vector containing values
#' for `\eqn{alpha_1..alpha_{N-1}}` where each \eqn{alpha_i} corresponds to the
#' shape parameter for a dirichlet distribution for each attribute
#' \eqn{1..(N-1)}. The dirichlet shape parameters are used to draw weights
#' for each trial from the data.
#' The parameter vector `x` also contains `r` form of the distance metric and
#' `gamma` sensitivity to differences in attribute space.
#'
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @param loc1_attrs The column names containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 1
#' @param loc2_attrs The name of the column containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 2
#' @importFrom stats pnorm
#' @export
rp_dir <- function(x, data, loc1_attrs = c("left_price", "left_memory"),
                      loc2_attrs = c("right_price", "right_memory")) {
  x <- transform_pars_dir(x, fwd=FALSE)
  w <- rdirichlet(nrow(data), x[grep("^alpha", names(x))])
  stopifnot(length(loc1_attrs) == length(loc2_attrs))
  n_attr <- length(loc1_attrs)
  da <- distance(w,
                 data.matrix(data[loc1_attrs]),
                 rep(0, n_attr),
                 x["r"])
  db <- distance(w,
                 data.matrix(data[loc2_attrs]),
                 rep(0, n_attr),
                 x["r"])
  resp_p <- da^x["gamma"] / (da^x["gamma"] + db^x["gamma"])
#  resp_p <- pnorm(x["delta"] + x["s"] * d)
  # Make sure not too close to 1
  resp_p <- pmin(resp_p, 1 - 1e-10)
  # Make sure not too close to 0
  resp_p <- pmax(resp_p, 1e-10)
  resp_p
}

#' Calculate response probabilities for 2xN using the normalised weights method
#'
#' Used for both calculating likelihoods or sampling from a parameter estimate
#' This function calculates the probabilities of responses to a dataset
#' based on the provided parameters (in the `x` vector).
#'
#' The `x` parameter vector is expected to be a named vector containing values
#' for `\eqn{w_2..w_{N}}` where each \eqn{w_i} corresponds to the
#' unnormalised weight for an attribute, with the first weight (typically price)
#' fixed at 1.
#' The parameter vector `x` also contains `r` form of the distance metric and
#' `gamma` sensitivity to differences in attribute space.
#'
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @param loc1_attrs The column names containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 1
#' @param loc2_attrs The name of the column containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 2
#' @importFrom stats pnorm
#' @export
rp_weight <- function(x, data, loc1_attrs = c("left_price", "left_memory"),
                      loc2_attrs = c("right_price", "right_memory")) {
  x <- transform_pars_weights(x, fwd=FALSE)
  weights <- grep("^w", names(x))
  x[weights] <- pmax(pmin(x[weights], 1e3), 1e-3)
  w <- c(1, x[weights])
  w <- w/sum(w)
  stopifnot(length(loc1_attrs) == length(loc2_attrs))
  n_attr <- length(loc1_attrs)
  da <- distance(w,
                 data.matrix(data[loc1_attrs]),
                 rep(0, n_attr),
                 x["r"])
  db <- distance(w,
                 data.matrix(data[loc2_attrs]),
                 rep(0, n_attr),
                 x["r"])
  resp_p <- da^x["gamma"] / (da^x["gamma"] + db^x["gamma"])
  # Make sure not too close to 1
  resp_p <- pmin(resp_p, 1 - 1e-10)
  # Make sure not too close to 0
  resp_p <- pmax(resp_p, 1e-10)
  resp_p
}


#' Calculate response probabilities for heatmap dirichlet MDS
#'
#' Used for both calculating likelihoods or sampling from a parameter estimate
#' This function calculates the probabilities of responses to a dataset
#' based on the provided parameters (in the `x` vector).
#'
#' The `x` parameter vector is expected to be a named vector containing values
#' for `w` weight on all attributes, `r` form of the distance metric and
#' `gamma` sensitivity to differences in attribute space.
#' All attribute values are expected to be already transformed.
#'
#' @param x A named vector of parameters. Expects certain values as in
#'   description
#' @param data A data.frame compatible object with specific rows as described.
#' @param loc1_attrs The column names containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 1
#' @param loc2_attrs The name of the column containing normalised values for the
#'   attributes (which has the weight value applied to it) in location 2
#' @importFrom stats pnorm
#' @export
rp_heatdir <- function(x, data, loc1_attrs = c("left_price", "left_memory"),
                       loc2_attrs = c("right_price", "right_memory")) {
  stopifnot(length(loc1_attrs) == length(loc2_attrs))
  n_attr <- length(loc1_attrs)
  w <- x[grep("^w", names(x))]
  da <- distance(w,
                 data.matrix(data[loc1_attrs]),
                 rep(0, n_attr),
                 x["r"])
  db <- distance(w,
                 data.matrix(data[loc2_attrs]),
                 rep(0, n_attr),
                 x["r"])
  resp_p <- da^x["gamma"] / (da^x["gamma"] + db^x["gamma"])
  # Make sure not too close to 1
  resp_p <- pmin(resp_p, 1 - 1e-10)
  # Make sure not too close to 0
  resp_p <- pmax(resp_p, 1e-10)
  resp_p
}
