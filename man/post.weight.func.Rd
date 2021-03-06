% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{post.weight.func}
\alias{post.weight.func}
\title{post.weight.func}
\usage{
post.weight.func(pis, lik.mat)
}
\arguments{
\item{pis}{= object from EM output with prior weight P(Z=K) as computed from}

\item{lik.mat}{= a JxK matrix of likelihoods (may be training set) for the P(D|Z=K)}
}
\value{
a 1xK vector of posterior weight for each gene snp pait
}
\details{
converts the matrix of likelihood for each gene snp pairs to matrix of posterior weights
}
