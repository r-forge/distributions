\name{probdistr}
\alias{example.continuous}
\alias{example.discrete}
\title{Auxiliary functions for 'probdistr' package.}
\description{
    Stirling numbers of the second kind and permutation of positive integers.
}
\usage{
example.continuous(distr="unif")
example.discrete(distr="binom")
}
\arguments{
  \item{distr}{a name of the probability distribution.}
}
\details{
        For continuous distributions, 
	\code{distr} shoud be one of the following character string \code{unif}, \code{triangle},
	\code{beta}, \code{kumar},\code{genbeta}, \code{exp}, \code{shiftexp},
	\code{invexp}, \code{shiftlnorm}, \code{gamma},\code{lgamma}, \code{invgamma},
	\code{invtrgamma}, \code{transgamma}, \code{unif}, \code{gig}, \code{generlang},
	\code{laplace}, \code{chisq},\code{ncchisq}, \code{chi}, \code{weibull}, 
	\code{weibull}, \code{invweibull}, \code{inchisq}, \code{student}, \code{cauchy},
	\code{pareto1}, \code{pareto2}, \code{pareto3}, \code{pareto4}, \code{invpareto},
	\code{burr}, \code{gpd}, \code{invburr}, \code{gumbel}.
	
	For discrete distributions,
	\code{distr} shoud be one of the following character string \code{binom}, \code{nbinom},
	\code{pois}, \code{geom}.
}
\value{
  a matrix with the abcisse and the plotted values.
}
\author{
 Christophe Dutang.
}
\examples{
example.discrete()
example.continuous()
}
\keyword{distribution}

