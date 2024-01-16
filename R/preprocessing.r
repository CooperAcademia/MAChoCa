#' Normalises a column of values to the 0-1 range for estimation
#'
#' Take a numeric vector of raw attribute values and a vector containing the
#' full range of possible values and normalise the values such that the
#' smallest real value is tranformed to 0 and the largest to 1. It is also
#' possible to invert the returned values for attributes such as price where
#' higher raw values should be lower.
#'
#' @param val_col A vector of raw attribute values.
#' @param val_rng A vector containning the range of possible values
#' @param invert If true return 1 - normalised value
#'
#' @return A vector of normalised values
#' @export
normalise <- function(val_col, val_rng = val_col, invert = FALSE) {
  shifted <- val_col - min(val_rng)
  val_rng <- val_rng - min(val_rng)
  normalised <- shifted / max(val_rng)
  if (invert) {
    normalised <- 1 - normalised
  }
  normalised
}
