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
  Normalizes the data matrix <data> by raking the Nrows by Ncols matrix such that the row means and column means equal 1. Missing information needs to be presented as nan values and not as zero values, because CONSTANd is able to ignore nan-values when calculating the mean. The variable <maxIterations> is an integer value that denotes the number of raking cycles. The variable <precision> defines the stopping criteria based on the L1-norm as defined by Friedrich Pukelsheim, Bruno Simeone in "On the Iterative Proportional Fitting Procedure: Structure of Accumulation Points and L1-Error Analysis"
}

\value{
  An \code{ExpressionSet} object. The \code{assayDataElement(eset, "exprs")}
  will be populated with normalized signals, \code{assayDataElement(eset, "snDetect")}
  will be populated with S/N ratio values, and the \code{phenoData} slot will
  be populated with information from \code{designFile}. Further analysis can
  be performed on the \code{ExpressionSet} object with various R and Bioconductor
  packages.
}

\author{Joris Van Houtven (joris.vanhoutven@uhasselt.be), Jef Hooyberghs <jef.hooyberghs@uhasselt.be>, Dirk Valkenborg <dirk.valkenborg@uhasselt.be>}

\examples{
#- normalized_matrix <- CONSTANd(data_matrix)
#- normalized_matrix <- CONSTANd(data_matrix, precision=1e-6, maxIterations=30)
}

\keyword{manip}% at least one, from doc/KEYWORDS
