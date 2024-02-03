##
## Dirichlet (Multivariate Beta)
##

#' The Dirichlet Distribution
#'
#' Density function and random generation from the Dirichlet distribution.
#'
#' The Dirichlet distribution is the multidimensional generalization of the
#' beta distribution.
#'
#' @aliases ddirichlet rdirichlet
#'
#' @name Dirichlet
#'
#' @param x A vector containing a single deviate or matrix containing one
#' random deviate per row.
#'
#' @param n Number of random vectors to generate.
#'
#' @param alpha Vector of shape parameters, or matrix of shape parameters
#' corresponding to the number of draw.
#'
#' @return \code{ddirichlet} gives the density. \code{rdirichlet} returns a
#' matrix with \code{n} rows, each containing a single Dirichlet random
#' deviate.
#'
#' @author Code is taken from MCMCpack: Markov Chain Monte Carlo (MCMC) Package
#' (MCMCpack). His code was based on code posted by Ben Bolker to R-News on 15
#' Dec 2000.
#'
#' @seealso \code{\link[stats]{Beta}}
#'
#' @keywords distribution
#'
#' @examples
#'
#'   density <- ddirichlet(c(.1,.2,.7), c(1,1,1))
#'   draws <- rdirichlet(20, c(1,1,1) )
#'
NULL

#' @rdname Dirichlet
#' @export
ddirichlet <- function(x, alpha) {

	dirichlet1 <- function(x, alpha) {
		logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
		s <- sum((alpha-1)*log(x))
		exp(sum(s)-logD)
	}

	# make sure x is a matrix
	if(!is.matrix(x))
		if(is.data.frame(x))
			x <- as.matrix(x)
		else
			x <- t(x)
	if(!is.matrix(alpha))
		alpha <- matrix( alpha, ncol=length(alpha), nrow=nrow(x), byrow=TRUE)

	if( any(dim(x) != dim(alpha)) )
		stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")

	pd <- vector(length=nrow(x))
	for(i in 1:nrow(x))
		pd[i] <- dirichlet1(x[i,],alpha[i,])

	# Enforce 0 <= x[i,j] <= 1, sum(x[i,]) = 1
	pd[ apply( x, 1, function(z) any( z <0 | z > 1)) ] <- 0
	pd[ apply( x, 1, function(z) all.equal(sum( z ),1) !=TRUE) ] <- 0
	return(pd)
}


#' @rdname Dirichlet
#' @export
rdirichlet <- function(n, alpha) {
	l <- length(alpha)
	x <- matrix(rgamma(l*n,alpha),ncol=l,byrow=TRUE)
	sm <- x%*%rep(1,l)
	return(x/as.vector(sm))
}
