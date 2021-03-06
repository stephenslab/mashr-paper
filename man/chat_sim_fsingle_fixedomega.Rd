% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulationswithmean.R
\name{chat_sim_fsingle_fixedomega}
\alias{chat_sim_fsingle_fixedomega}
\title{chat_sim_fsingle_fixedomega}
\usage{
chat_sim_fsingle_fixedomega(n = 1000, d = 8, omega = 2, esd = 0.1,
  mu = 0, reals = 0.1)
}
\arguments{
\item{n}{number of SNPs}

\item{d}{number of subgroups}

\item{esd}{the standard error of E in chat=mu+c+E, i.e., E~N(0,diag(esd^2))}

\item{betasd}{the vriance of true effect}
}
\value{
omega the size of the effects

f the factors from which they were simulated
}
