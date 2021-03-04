## -----------------------------------------------------------------------------
library(ICSKAT)
set.seed(0)
n <- 10^4
q <- 50

# generate data
# all SNPs have minor allele frequency 0.3 in this toy example
gMat <- matrix(data=rbinom(n=n*q, size=2, prob=0.3), nrow=n)
# gender and whether they take daily vitamins
xMat <- matrix(data=rbinom(n=n*2, size=1, prob=0.5), nrow=n)
# the baseline cumulative hazard function
bhFunInv <- function(x) {x}
# observation times
obsTimes <- 1:4
# no effect of either gender or daily vitamins
etaVec <- rep(0, n)

# generate data
outcomeDat <- gen_IC_data(bhFunInv = bhFunInv, obsTimes = obsTimes, windowHalf = 0.25,
                          probMiss = 0.1, etaVec = etaVec)
lt <- outcomeDat$leftTimes
rt <- outcomeDat$rightTimes
# indicators of left- and right-censoring
tpos_ind <- as.numeric(lt > 0)
obs_ind <- as.numeric(rt != Inf)

# make design matrix with cubic spline terms
dmats <- make_IC_dmat(xMat, lt, rt, obs_ind, tpos_ind)
# fit null model - only need to do this once for each genetic set (note there is no information
# on the SNPs used here)
nullFit <- ICSKAT_fit_null(init_beta = rep(0, 5), left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, 
                           obs_ind = obs_ind, tpos_ind = tpos_ind, lt = lt, rt = rt)
# perform the ICSKAT and Burden tests
icskatOut <- ICskat(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, lt = lt, rt = rt,
       obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat, null_beta = as.numeric(nullFit$beta_fit), Itt = nullFit$Itt)
icskatOut$p_SKAT
icskatOut$p_burden
# perform the ICSKATO test
ICSKATO(icskatOut = icskatOut)$pval

## -----------------------------------------------------------------------------

# another gene with 100 SNPs in it
gMat2 <- matrix(data=rbinom(n=n*100, size=2, prob=0.3), nrow=n)

# we don't need to run the make_IC_dmat or ICSKAT_fit_null functions again as long
# as the event times and non-genetic covariates haven't changed.
icskatOut2 <- ICskat(left_dmat = dmats$left_dmat, right_dmat=dmats$right_dmat, lt = lt, rt = rt,
       obs_ind = obs_ind, tpos_ind = tpos_ind, gMat = gMat2, null_beta = as.numeric(nullFit$beta_fit), Itt = nullFit$Itt)
icskatOut2$p_SKAT
icskatOut2$p_burden
# perform the ICSKATO test
ICSKATO(icskatOut = icskatOut2)$pval


