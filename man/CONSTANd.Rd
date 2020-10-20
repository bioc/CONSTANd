\name{CONSTANd}
\alias{CONSTANd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data normalization by matrix raking}
\description{
    Normalizes the data matrix by raking the Nrows by Ncols matrix such that the row means and column means equal Ncols and Nrows, respectively.
}
\usage{
CONSTANd(data, precision=1e-5, maxIterations=50, target=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{Nrows by Ncols matrix.}
    \item{precision}{Combined allowed deviation (residual error) of col and row means from target value.}
    \item{maxIterations}{Maximum amount of iterations (1x row and 1x col per iteration).}
    \item{target}{The mean value of quantifications in each row and column after normalization.}
}

\details{
  Normalizes the data matrix <data> by raking (using the RAS method by Bacharach, see references) the Nrows by Ncols matrix such that the row means and column means equal 1. The result is a normalized data matrix K=RAS, a product of row mulipliers R and column multipliers S with the original matrix A. Missing information needs to be presented as nan values and not as zero values, because CONSTANd is able to ignore nan-values when calculating the mean. The variable <maxIterations> is an integer value that denotes the number of raking cycles. The variable <precision> defines the stopping criteria based on the L1-norm as defined by Friedrich Pukelsheim, Bruno Simeone in "On the Iterative Proportional Fitting Procedure: Structure of Accumulation Points and L1-Error Analysis".
}

\value{
    \item{normalized_data}{Normalized data matrix `K=RAS` in the RAS-formulation of the problem.}
    \item{convergence_trail}{Precision acquired after each raking iteration (last value is the final precision).}
    \item{R}{Row multipliers in the `K=RAS` formulation of the problem.}
    \item{S}{Column multipliers in the `K=RAS` formulation of the problem.}
}

\author{Joris Van Houtven (joris.vanhoutven@uhasselt.be), Geert Jan Bex <geertjan.bex@uhasselt.be>, Dirk Valkenborg <dirk.valkenborg@uhasselt.be>}

\examples{
# generic use (mock data)
data_matrix <- matrix(runif(20), c(5,4))
normalized_matrix <- CONSTANd(data_matrix)$normalized_data

# customize parameters
result <- CONSTANd(data_matrix, precision=1e-3, maxIterations=30)$normalized_data

# explore parts of the result object
normalized_matrix <- result$normalized_matrix
num_iterations_performed <- length(result$convergence_trail)
attained_precision <- result$convergence_trail[num_iterations_performed]
}

\references{
    Maes, Evelyne, et al. "CONSTANd: A normalization method for isobaric labeled spectra by constrained optimization." Molecular & Cellular Proteomics 15.8 (2016): 2779-2790. \url{https://doi.org/10.1074/mcp.M115.056911}. Accessed 18 Oct. 2020.

    Bacharach, Michael. "Estimating Nonnegative Matrices from Marginal Data." International Economic Review, vol. 6, no. 3, 1965, pp. 294–310. JSTOR, \url{https://doi.org/10.2307\%2F2525582}. Accessed 18 Oct. 2020.
}

\keyword{manip}% at least one, from doc/KEYWORDS
