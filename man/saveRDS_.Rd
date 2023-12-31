% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.r
\name{saveRDS_}
\alias{saveRDS_}
\title{saveRDS and return .data}
\usage{
saveRDS_(.data, file, ...)
}
\arguments{
\item{.data}{object to save}

\item{file}{filename to save the object to.}

\item{...}{passed onto saveRDS}
}
\value{
whatever object passed in.
}
\description{
This function uses the built-in saveRDS to save the data to a file,
but then returns the data so that it can be used for dplyr pipeline.
}
