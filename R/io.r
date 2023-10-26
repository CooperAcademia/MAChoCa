#' saveRDS and return .data
#'
#' This function uses the built-in saveRDS to save the data to a file,
#' but then returns the data so that it can be used for dplyr pipeline.
#'
#' @param .data object to save
#' @param file filename to save the object to.
#' @param ... passed onto saveRDS
#'
#' @return whatever object passed in.
#' @export
saveRDS_ <- function(.data, file, ...){
  saveRDS(object = .data, file = file, ...)
  return(.data)
}
