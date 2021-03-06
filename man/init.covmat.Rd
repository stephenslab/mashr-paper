% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sarah_mixem.R
\name{init.covmat}
\alias{init.covmat}
\title{init.covmat}
\usage{
init.covmat(t.stat = t.stat, factor.mat = factors, lambda.mat = lambda, K,
  P)
}
\arguments{
\item{t.stat}{matrix of strong t statistics (MxR) from which to derive covariance matrices}

\item{factor.mat}{KxR matrix of factors from SFA}

\item{lambda.mat}{KxR matric of loadings, also from SFA}

\item{K}{number of components to fit (i.e., the first dimension of the array)}

\item{P}{rank of SVD approximation}
}
\value{
KxRxR matric of prior covariance matrices to initialize the EM
}
