% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformers.r
\name{transform_pars_cut}
\alias{transform_pars_cut}
\title{Transform the parameter vector (2 x N) cut method for use in pmwg}
\usage{
transform_pars_cut(x, fwd = TRUE)
}
\arguments{
\item{x}{The named vector of parameter estimates}

\item{fwd}{Move certain parameter to real number line or back}
}
\value{
transformed parameter vector x
}
\description{
See also the description for rp_cut This is a helper function that
transforms the parameters. To transform for pmwg the \code{r} and \code{gamma} values
are exponentiated in order to move to positive only and cut points \code{c1}..\code{cn}
are pnormed so that values are between 0 and 1.
}
\details{
These operation are reversed when \code{fwd=TRUE}
}
