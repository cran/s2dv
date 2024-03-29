% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UltimateBrier.R
\name{UltimateBrier}
\alias{UltimateBrier}
\title{Compute Brier scores}
\usage{
UltimateBrier(
  exp,
  obs,
  dat_dim = NULL,
  memb_dim = "member",
  time_dim = "sdate",
  quantile = TRUE,
  thr = c(5/100, 95/100),
  type = "BS",
  decomposition = TRUE,
  ncores = NULL
)
}
\arguments{
\item{exp}{A numeric array of forecast anomalies with named dimensions that
at least include 'memb_dim', and 'time_dim'. It can be provided
by \code{Ano()}.}

\item{obs}{A numeric array of observational reference anomalies with named
dimensions that at least include 'time_dim'. If it has 
'memb_dim', the length  must be 1. The dimensions should be consistent with
'exp' except 'dat_dim' and 'memb_dim'. It can be provided by \code{Ano()}.}

\item{dat_dim}{A character string indicating the name of the dataset 
dimension in 'exp' and 'obs'. The default value is NULL (no dataset).
dimension, set NULL.}

\item{memb_dim}{A character string indicating the name of the member 
dimension in 'exp' (and 'obs') for ensemble mean calculation. The default
value is 'member'.}

\item{time_dim}{A character string indicating the dimension along which to 
compute the probabilistic scores. The default value is 'sdate'.}

\item{quantile}{A logical value to decide whether a quantile (TRUE) or a 
threshold (FALSE) is used to estimate the forecast and observed 
probabilities. If 'type' is 'FairEnsembleBS' or 'FairEnsembleBSS', it must
be TRUE. The default value is TRUE.}

\item{thr}{A numeric vector to be used in probability calculation (for 'BS', 
'FairStartDatesBS', 'BSS', and 'FairStartDatesBSS') and binary event 
judgement (for 'FairEnsembleBS' and 'FairEnsembleBSS'). It is as
quantiles if 'quantile' is TRUE or as thresholds if 'quantile' is FALSE. 
The default value is \code{c(0.05, 0.95)} for 'quantile = TRUE'.}

\item{type}{A character string of the desired score type. It can be the 
 following values:
\itemize{
 \item{'BS': Simple Brier Score. Use SpecsVerification::BrierDecomp inside.}
 \item{'FairEnsembleBS': Corrected Brier Score computed across ensemble 
   members. Use SpecsVerification::FairBrier inside.}
 \item{'FairStartDatesBS': Corrected Brier Score computed across starting 
   dates. Use s2dv:::.BrierScore inside.}
 \item{'BSS': Simple Brier Skill Score. Use s2dv:::.BrierScore inside.}
 \item{'FairEnsembleBSS': Corrected Brier Skill Score computed across 
   ensemble members. Use SpecsVerification::FairBrierSs inside.}
 \item{'FairStartDatesBSS': Corrected Brier Skill Score computed across 
   starting dates. Use s2dv:::.BrierScore inside.}
}
 The default value is 'BS'.}

\item{decomposition}{A logical value to determine whether the decomposition 
of the Brier Score should be provided (TRUE) or not (FALSE). It is only 
used when 'type' is 'BS' or 'FairStartDatesBS'. The default value is TRUE.}

\item{ncores}{An integer indicating the number of cores to use for parallel 
computation. The default value is NULL.}
}
\value{
If 'type' is 'BS' or 'FairStartDatesBS' and 'decomposition' is TRUE, the 
output is a list of 4 arrays (see details below.) In other cases, the output 
is an array of Brier scores or Brier skill scores. All the arrays have the 
same dimensions:
c(nexp, nobs, no. of bins, the rest dimensions of 'exp' except 'time_dim' and
'memb_dim'). 'nexp' and 'nobs' is the length of dataset dimension in 'exp'
and 'obs' respectively. If dat_dim is NULL, nexp and nobs are omitted.\cr
The list of 4 includes: 
 \itemize{
   \item{$bs: Brier Score}
   \item{$rel: Reliability component}
   \item{$res: Resolution component}
   \item{$unc: Uncertainty component}
 }
}
\description{
Interface to compute probabilistic scores (Brier Score, Brier Skill Score) 
from the forecast and observational data anomalies. It provides six types
to choose.
}
\examples{
 \dontshow{
startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
sampleData <- s2dv:::.LoadSampleData('tos', c('experiment'),
                                    c('observation'), startDates,
                                    leadtimemin = 1,
                                    leadtimemax = 4,
                                    output = 'lonlat',
                                    latmin = 27, latmax = 48,
                                    lonmin = -12, lonmax = 40)
 }
sampleData$mod <- Season(sampleData$mod, monini = 11, moninf = 12, monsup = 2)
sampleData$obs <- Season(sampleData$obs, monini = 11, moninf = 12, monsup = 2)
clim <- Clim(sampleData$mod, sampleData$obs)
exp <- Ano(sampleData$mod, clim$clim_exp)
obs <- Ano(sampleData$obs, clim$clim_obs)
bs <- UltimateBrier(exp, obs, dat_dim = 'dataset')
bss <- UltimateBrier(exp, obs, type = 'BSS', dat_dim = 'dataset')

}
