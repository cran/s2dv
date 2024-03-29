% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Corr.R
\name{Corr}
\alias{Corr}
\title{Compute the correlation coefficient between an array of forecast and their corresponding observation}
\usage{
Corr(
  exp,
  obs,
  time_dim = "sdate",
  dat_dim = NULL,
  comp_dim = NULL,
  limits = NULL,
  method = "pearson",
  memb_dim = NULL,
  memb = TRUE,
  pval = TRUE,
  conf = TRUE,
  sign = FALSE,
  alpha = 0.05,
  ncores = NULL
)
}
\arguments{
\item{exp}{A named numeric array of experimental data, with at least dimension
'time_dim'.}

\item{obs}{A named numeric array of observational data, same dimensions as
parameter 'exp' except along 'dat_dim' and 'memb_dim'.}

\item{time_dim}{A character string indicating the name of dimension along
which the correlations are computed. The default value is 'sdate'.}

\item{dat_dim}{A character string indicating the name of dataset (nobs/nexp) 
dimension. The default value is NULL (no dataset).}

\item{comp_dim}{A character string indicating the name of dimension along which
obs is taken into account only if it is complete. The default value
is NULL.}

\item{limits}{A vector of two integers indicating the range along comp_dim to 
be completed. The default is c(1, length(comp_dim dimension)).}

\item{method}{A character string indicating the type of correlation: 
'pearson', 'spearman', or 'kendall'. The default value is 'pearson'.}

\item{memb_dim}{A character string indicating the name of the member 
dimension. It must be one dimension in 'exp' and 'obs'. If there is no 
member dimension, set NULL. The default value is NULL.}

\item{memb}{A logical value indicating whether to remain 'memb_dim' dimension
(TRUE) or do ensemble mean over 'memb_dim' (FALSE). Only functional when 
'memb_dim' is not NULL. The default value is TRUE.}

\item{pval}{A logical value indicating whether to return or not the p-value 
of the test Ho: Corr = 0. The default value is TRUE.}

\item{conf}{A logical value indicating whether to return or not the confidence 
intervals. The default value is TRUE.}

\item{sign}{A logical value indicating whether to retrieve the statistical
significance of the test Ho: Corr = 0 based on 'alpha'. The default value is
FALSE.}

\item{alpha}{A numeric indicating the significance level for the statistical
significance test. The default value is 0.05.}

\item{ncores}{An integer indicating the number of cores to use for parallel 
computation. The default value is NULL.}
}
\value{
A list containing the numeric arrays with dimension:\cr 
 c(nexp, nobs, exp_memb, obs_memb, all other dimensions of exp except 
 time_dim and memb_dim).\cr
nexp is the number of experiment (i.e., 'dat_dim' in exp), and nobs is the 
number of observation (i.e., 'dat_dim' in obs). If dat_dim is NULL, nexp and 
nobs are omitted. exp_memb is the number of member in experiment (i.e., 
'memb_dim' in exp) and obs_memb is the number of member in observation (i.e.,
'memb_dim' in obs). If memb = F, exp_memb and obs_memb are omitted.\cr\cr
\item{$corr}{
 The correlation coefficient. 
}
\item{$p.val}{
 The p-value. Only present if \code{pval = TRUE}.
}
\item{$conf.lower}{
 The lower confidence interval. Only present if \code{conf = TRUE}.
}
\item{$conf.upper}{
 The upper confidence interval. Only present if \code{conf = TRUE}.
}
\item{$sign}{
 The statistical significance. Only present if \code{sign = TRUE}.
}
}
\description{
Calculate the correlation coefficient (Pearson, Kendall or Spearman) for 
an array of forecast and an array of observation. The correlations are 
computed along 'time_dim' that usually refers to the start date dimension. If
'comp_dim' is given, the correlations are computed only if obs along comp_dim
dimension are complete between limits[1] and limits[2], i.e., there is no NA 
between limits[1] and limits[2]. This option can be activated if the user 
wants to account only for the forecasts which the corresponding observations 
are available at all leadtimes.\cr 
The confidence interval is computed by the Fisher transformation and the 
significance level relies on an one-sided student-T distribution.\cr 
The function can calculate ensemble mean before correlation by 'memb_dim' 
specified and 'memb = F'. If ensemble mean is not calculated, correlation will
be calculated for each member.
If there is only one dataset for exp and obs, you can simply use cor() to 
compute the correlation.
}
\examples{
# Case 1: Load sample data as in Load() example: 
example(Load) 
clim <- Clim(sampleData$mod, sampleData$obs) 
ano_exp <- Ano(sampleData$mod, clim$clim_exp) 
ano_obs <- Ano(sampleData$obs, clim$clim_obs) 
runmean_months <- 12 

# Smooth along lead-times   
smooth_ano_exp <- Smoothing(ano_exp, runmeanlen = runmean_months) 
smooth_ano_obs <- Smoothing(ano_obs, runmeanlen = runmean_months) 
required_complete_row <- 3  # Discard start dates which contain any NA lead-times 
leadtimes_per_startdate <- 60 
corr <- Corr(MeanDims(smooth_ano_exp, 'member'),              
            MeanDims(smooth_ano_obs, 'member'),              
            comp_dim = 'ftime', dat_dim = 'dataset', 
            limits = c(ceiling((runmean_months + 1) / 2),                         
            leadtimes_per_startdate - floor(runmean_months / 2))) 

# Case 2: Keep member dimension
corr <- Corr(smooth_ano_exp, smooth_ano_obs, memb_dim = 'member', dat_dim = 'dataset')
# ensemble mean
corr <- Corr(smooth_ano_exp, smooth_ano_obs, memb_dim = 'member', memb = FALSE,
            dat_dim = 'dataset')

}
