#' Calculate distances in psychological space
#'
#' Take a collection of weights and an anchor point, and a set of attribute
#' values and return the distance between each point represented by the
#' attribute values and the anchor, adjusted by the metric form `r`.
#'
#' @param w A vector or matrix of attribute weights.
#' @param ... Unused, for extensibility
#' @param attrs A matrix with ncol(attrs) = length(weights)
#' @param anchor An "exemplar" or anchor point to which attrs will be compared
#' @param r A positive value between that changes the metric form
#'
#' @return A vector of distances for each row of attrs
#' @export
distance <- function (w, ...) {
	UseMethod("distance", w)
}


#' @export
#' @rdname distance
distance.numeric <- function(w, attrs, anchor, r) {
	stopifnot(length(w) == ncol(attrs),
						length(w) == length(anchor),
						r > 0,
						isTRUE(all.equal(sum(w), 1)))

	attr_dist <- sapply(seq_along(w), FUN = function(i) {
												w[i] * abs(attrs[,i] - anchor[i])^r
						})
	option_dist <- apply(attr_dist, 1, sum)^(1/r)

	option_dist
}


#' @export
#' @rdname distance
distance.matrix <- function(w, attrs, anchor, r) {
	stopifnot(ncol(w) == ncol(attrs),
						ncol(w) == length(anchor),
						r > 0)
  summed_weights <- apply(w, 1, sum)

  if (!isTRUE(all.equal(summed_weights, rep(1, nrow(w))))) {
    # Check rows where sum is not close enough to 1. Borrowed from approxeq {cgwtools}
    too_far <- abs(summed_weights - 1) >= .Machine$double.eps ^ 0.5
    stop(paste("ERROR: found weights where sum is not close enough to 1\n",
              w[too_far], "\n"))
  }

	attr_dist <- sapply(seq_len(ncol(w)), FUN = function(i) {
												w[,i] * abs(attrs[,i] - anchor[i])^r
						})
	option_dist <- apply(attr_dist, 1, sum)^(1/r)

	option_dist
}
