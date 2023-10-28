#' Calculate distance in psychological space
#'
#' Take a vector of weights and an anchor point, and a set of attribute values
#' and return the distance between each point represented by the attribute
#' values and the anchor, adjusted by the metric form `r`.
#'
#' @param w A vector of attribute weights.
#' @param attrs A matrix with ncol(attrs) = length(weights)
#' @param anchor An "exemplar" or anchor point to which attrs will be compared
#' @param r A positive value between that changes the metric form
#'
#' @return A vector of distances for each row of attrs
#' @export
distance <- function(w, attrs, anchor, r) {
  stopifnot(length(w) == ncol(attrs),
            length(w) == length(anchor),
            r > 0,
            sum(w) == 1)
  
  attr_dist <- sapply(seq_along(w), FUN = function(i) {
           w[i] * abs(attrs[,i] - anchor[i])^r
            })
  option_dist <- apply(attr_dist, 1, sum)^(1/r)

  option_dist
}
