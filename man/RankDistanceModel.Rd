% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/rankdist.R
\name{RankDistanceModel}
\alias{RankDistanceModel}
\title{Fit A Mixture of Distance-based Models}
\usage{
RankDistanceModel(dat, init, ctrl)
}
\arguments{
\item{dat}{A \linkS4class{RankData} object.}

\item{init}{A \linkS4class{RankInit} object.}

\item{ctrl}{A \linkS4class{RankControl} object.}
}
\value{
A list containing the following components:
\describe{
\item{\code{modal_ranking.est}}{the estimated pi0 for each cluster.}
\item{\code{p}}{the probability of each cluster.}
\item{\code{w.est}}{the estimated weights of each cluster.}
\item{\code{param.est}}{the param parametrisation of weights of each cluster.}
\item{\code{SSR}}{the sum of squares of Pearson residuals}
\item{\code{log_likelihood}}{the fitted log_likelihood}
\item{\code{BIC}}{the fitted Bayesian Information Criterion value}
\item{\code{free_params}}{the number of free parameters in the model}
\item{\code{expectation}}{the expected value of each observation given by the model}
\item{\code{iteration}}{the number of EM iteration}
\item{\code{model.call}}{the function call}
}
}
\description{
\code{RankDistModel} fits a mixture of ranking models based on weighted Kendall distance.
}
\details{
The procedure will estimate central rankings, the probability of each cluster and weights.
}

