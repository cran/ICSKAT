% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ACAT.R
\name{ACAT}
\alias{ACAT}
\title{Aggregated Cauchy Association Test}
\usage{
ACAT(Pvals, Weights = NULL)
}
\arguments{
\item{Pvals}{a numeric vector of p-values to be combined by ACAT.}

\item{Weights}{a numeric vector of non-negative weights for the combined p-values. When it is NULL, the equal weights are used.}
}
\value{
p-value of ACAT.
}
\description{
A p-value combination method using the Cauchy distribution.
Code provided by Dr. Yaowu Liu.
}
\examples{
p.values<-c(2e-02,4e-04,0.2,0.1,0.8)
ACAT(Pvals=p.values)

}
\author{
Yaowu Liu
}
