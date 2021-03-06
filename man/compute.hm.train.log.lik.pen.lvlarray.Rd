% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{compute.hm.train.log.lik.pen.lvlarray}
\alias{compute.hm.train.log.lik.pen.lvlarray}
\title{compute.hm.train.log.lik.pen.lvlarray}
\usage{
compute.hm.train.log.lik.pen.lvlarray(train.b, lvlarray, covmat, A, pen)
}
\arguments{
\item{train.b}{=  JxR matrix of training beta hats}

\item{covmat}{= LxK dimenstional (unlisted list) of prior covariance matrices}

\item{A}{output file name}

\item{se.train}{= JxR matrix of training standard errors}

\item{Pen}{likelihood penalty, default1}
}
\value{
An object containing pis and model fit from the EM, and a pdf barplot
}
\description{
= takes a matrix of training sumamry statistics and their standard errors and computes the likelihood matrix according to a list of covariance matrices, using the exponent of the log likelihood - max (llj)
}
