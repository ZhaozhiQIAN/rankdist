\name{rankdist-package}
\alias{rankdist-package}
\alias{rankdist}
\docType{package}
\title{
A package for fitting distance based ranking models
}
\description{
The rankdist package implements distance based probability models for ranking data. Mixture models are also supported. Ranking data are stored as S4 objects to avoid confusions about representations.
}
\details{
\tabular{ll}{
Package: \tab rankdist\cr
Type: \tab Package\cr
Version: \tab 0.3.1\cr
Date: \tab 2015-05-01\cr
License: \tab GPL (>= 2)\cr
}
Distance based models are effective ways to model ranking data. This package supports models based on Kendall distance and weighted Kendall distance. Mixture models can be easily fitted as well. The package includes a well-studied ranking data set, the APA Election data.
}
\author{
Zhaozhi Qian

Maintainer: Zhaozhi Qian <qianzhaozhi@connect.hku.hk>
}
\references{
Alvo, M., & Yu, Philip L.H. (2014), Statistical Methods for Ranking Data
}

\keyword{ package }

\examples{
\dontrun{
fitted_model <- RankDistanceModel(rankdata,rankinit,rankctrl)
}
}
