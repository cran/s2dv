% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSE.R
\name{MSE}
\alias{MSE}
\title{Compute mean square error}
\usage{
MSE(
  exp,
  obs,
  time_dim = "sdate",
  dat_dim = NULL,
  memb_dim = NULL,
  comp_dim = NULL,
  limits = NULL,
  conf = TRUE,
  alpha = 0.05,
  ncores = NULL
)
}
\arguments{
\item{exp}{A named numeric array of experimental data, with at least 
'time_dim' dimension. It can also be a vector with the same length as 'obs'.}

\item{obs}{A named numeric array of observational data, same dimensions as
parameter 'exp' except along 'dat_dim' and 'memb_dim'. It can also be a 
vector with the same length as 'exp'.}

\item{time_dim}{A character string indicating the name of dimension along
which the correlations are computed. The default value is 'sdate'.}

\item{dat_dim}{A character string indicating the name of dataset or member 
(nobs/nexp) dimension. The datasets of exp and obs will be paired and 
computed MSE for each pair. The default value is NULL.}

\item{memb_dim}{A character string indicating the name of the member dimension
to compute the ensemble mean; it should be set to NULL if the input data are
already the ensemble mean. The default value is NULL.}

\item{comp_dim}{A character string indicating the name of dimension along which
obs is taken into account only if it is complete. The default value
is NULL.}

\item{limits}{A vector of two integers indicating the range along comp_dim to 
be completed. The default value is c(1, length(comp_dim dimension)).}

\item{conf}{A logical value indicating whether to retrieve the confidence 
intervals or not. The default value is TRUE.}

\item{alpha}{A numeric indicating the significance level for the statistical
significance test. The default value is 0.05.}

\item{ncores}{An integer indicating the number of cores to use for parallel 
computation. The default value is NULL.}
}
\value{
A list containing the numeric arrays with dimension:\cr 
 c(nexp, nobs, all other dimensions of exp except time_dim).\cr
nexp is the number of experiment (i.e., dat_dim in exp), and nobs is the 
number of observation (i.e., dat_dim in obs).\cr
\item{$mse}{
 The mean square error. 
}
\item{$conf.lower}{
 The lower confidence interval. Only present if \code{conf = TRUE}.
}
\item{$conf.upper}{
 The upper confidence interval. Only present if \code{conf = TRUE}.
}
}
\description{
Compute the mean square error for an array of forecasts and an array of
observations. The MSEs are computed along time_dim, the dimension which 
corresponds to the start date dimension. If comp_dim is given, the MSEs are 
computed only if obs along the comp_dim dimension are complete between 
limits[1] and limits[2], i.e. there are no NAs between limits[1] and 
limits[2]. This option can be activated if the user wants to account only 
for the forecasts for which the corresponding observations are available at 
all leadtimes.\cr
The confidence interval is computed by the chi2 distribution.\cr
}
\examples{
# Load sample data as in Load() example:
example(Load)
clim <- Clim(sampleData$mod, sampleData$obs)
ano_exp <- Ano(sampleData$mod, clim$clim_exp)
ano_obs <- Ano(sampleData$obs, clim$clim_obs)
smooth_ano_exp <- Smoothing(ano_exp, runmeanlen = 12, time_dim = 'ftime')
smooth_ano_obs <- Smoothing(ano_obs, runmeanlen = 12, time_dim = 'ftime')
res <- MSE(smooth_ano_exp, smooth_ano_obs, memb_dim = 'member', 
          comp_dim = 'ftime', limits = c(7, 54))

# Synthetic data:
exp1 <- array(rnorm(120), dim = c(dat = 3, sdate = 10, ftime = 4))
obs1 <- array(rnorm(80),  dim = c(dat = 2, sdate = 10, ftime = 4))
res1 <- MSE(exp1, obs1, comp_dim = 'ftime', dat_dim = 'dat')
 
exp2 <- array(rnorm(20), dim = c(sdate = 5, member = 4))
obs2 <- array(rnorm(10),  dim = c(sdate = 5, member = 2))
res2 <- MSE(exp2, obs2, memb_dim = 'member')

}
