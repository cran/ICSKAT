% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ICSKATO_bootstrap.R
\name{ICSKATO_bootstrap}
\alias{ICSKATO_bootstrap}
\title{ICSKATO_bootstrap.R}
\usage{
ICSKATO_bootstrap(
  icskatOut,
  B,
  intervalProbs,
  allVisits,
  quant_r,
  seed = NULL,
  null_fit,
  gMat,
  xMat,
  fitAgain,
  checkpoint = FALSE,
  downsample = 1,
  rhoVec
)
}
\arguments{
\item{icskatOut}{The output list from ICSKAT().}

\item{B}{Number of bootstrap replications.}

\item{intervalProbs}{n*(s+1) matrix where n is number of subjects and s is the number of visits possible, probability of falling in each interval.}

\item{allVisits}{n*s matrix with all the visit times for each subject.}

\item{quant_r}{Quantiles of time from make_IC_dmats, to keep them constant through bootstrapping.}

\item{seed}{Seed to start the bootstrapping.}

\item{null_fit}{The null fit output from ICSKAT_fit_null.}

\item{gMat}{Genotype matrix used in original test.}

\item{xMat}{n*p matrix of non-genetic covariates.}

\item{fitAgain}{Boolean, whether to fit the null model again in each bootstrap.}

\item{checkpoint}{Boolean, whether to print every time 100 bootstraps finish.}

\item{downsample}{A number in (0, 1], will use this fraction of the bootstrap iterations to try running the test with fewer bootstraps.}

\item{rhoVec}{Vector of rhos to search over in SKATO.}
}
\value{
A list with the elements:
\item{kurtQvec}{Vector of bootstrapped excess kurtosis of each Qrho.}
\item{varQvec}{Vector of bootstrapped variance of each Qrho.}
\item{meanQvec}{Vector of bootstrapped mean of each Qrho.}
\item{kurtKappa}{Bootstrapped kurtosis of kappa term without zeta.}
\item{kurtKappaAll}{Bootstrapped kurtosis of full kappa term with zeta.}
\item{varKappaAll}{Bootstrapped variance of full kappa term with zeta.}
\item{meanKappaAll}{Bootstrapped mean of full kappa term with zeta.}
\item{bootDF}{Matrix with B rows containing all the bootstrapped quantities over all iterations.}
\item{QrhoBoot}{Matrix with B rows containing all the bootstrapped Qrho values, one column for each rho.}
\item{listDS}{A list containing all of the other elements in this return list, except using the downsampled iterations.}
\item{nonNA}{Number of bootstraps that did not result in NA (and thus were not removed).}
}
\description{
The version of ICSKATO to run when bootstrapping to match kurtosis.
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
right_dmat=dmats$right_dmat, obs_ind = obs_ind, tpos_ind = tpos_ind,
lt = lt, rt = rt)
icskatOut <- ICskat(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat,
lt = lt, rt = rt, obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat,
null_beta = nullFit$beta_fit, Itt = nullFit$Itt)
intervalProbOutput <- construct_interval_probs(allTimes = outcomeDat$allVisits,
dmats = dmats, nullBeta = nullFit$beta_fit, p = ncol(xMat), nKnots=1)
ICSKATO_bootstrap(icskatOut = icSkatOut, B = 100, intervalProbs = intervalProbOutput$probMat,
allVisits = intervalProbOutput$allTimesFilled, quant_r = dmats$quant_r, seed = 0,
null_fit = nullFit, gMat = gMat, xMat, fitAgain = TRUE,
rhoVec=c(0, 0.01, 0.04, 0.09, 0.25, 0.5, 1))
}
