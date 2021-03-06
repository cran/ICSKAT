% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chiSqMatchFast.R
\name{chiSqMatchFast}
\alias{chiSqMatchFast}
\title{chiSqMatchFast.R}
\usage{
chiSqMatchFast(lambdaVec, alwaysCentral = FALSE)
}
\arguments{
\item{lambdaVec}{Numeric vector holding the eigenvalues of the A term, where we are interested in x^TAX
and x is multivariate normal.}

\item{alwaysCentral}{Boolean determining whether to always set the noncentrality parameter to 0, as in SKAT package.}
}
\value{
A list with the elements:
\item{sigmaQrho}{Standard deviation of the mixture distribution}
\item{muQrho}{Mean of the mixture distribution}
\item{delta}{Noncentrality parameter of the matched distribution}
\item{l}{Degrees of freedom of the matched distribution}
}
\description{
Match the moments of a mixture of scaled chi-square random variables to a single non-central chi-square,
assumes the quadratic form case where the mean of the multivariate normal V=RV is 0.
}
\examples{
set.seed(2)
gMat <- matrix(data=rbinom(n=2000, size=2, prob=0.3), nrow=100)
xMat <- matrix(data=rnorm(200), nrow=100)
bhFunInv <- function(x) {x}
obsTimes <- 1:5
etaVec <- rep(0, 100)
outcomeDat <- gen_IC_data(bhFunInv = bhFunInv, obsTimes = obsTimes, windowHalf = 0.1,
probMiss = 0.1, etaVec = etaVec)
lt <- outcomeDat$leftTimes
rt <- outcomeDat$rightTimes
tpos_ind <- as.numeric(lt > 0)
obs_ind <- as.numeric(rt != Inf)
dmats <- make_IC_dmat(xMat, lt, rt, obs_ind, tpos_ind)
nullFit <- ICSKAT_fit_null(init_beta = rep(0, 5), left_dmat = dmats$left_dmat,
right_dmat=dmats$right_dmat, obs_ind = obs_ind, tpos_ind = tpos_ind, lt = lt, rt = rt)
icskatOut <- ICskat(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat,
lt = lt, rt = rt, obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat,
null_beta = nullFit$beta_fit, Itt = nullFit$Itt)
Rrho <- matrix(data=0.5, nrow=20, ncol=20)
diag(Rrho) <- 1
toDecomp <- Rrho \%*\% icskatOut$sig_mat
tempEvals <- eigen(toDecomp, symmetric = TRUE, only.values = TRUE)$values
idx1 <- which(tempEvals >= 0)
idx2 <- which(tempEvals > mean(tempEvals[idx1])/100000)
tempEvals <- tempEvals[idx2]
chiSqMatchFast(lambdaVec = tempEvals)

}
