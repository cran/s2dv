% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NAO.R
\name{NAO}
\alias{NAO}
\title{Compute the North Atlantic Oscillation (NAO) Index}
\usage{
NAO(
  exp = NULL,
  obs = NULL,
  lat,
  lon,
  time_dim = "sdate",
  memb_dim = "member",
  space_dim = c("lat", "lon"),
  ftime_dim = "ftime",
  ftime_avg = 2:4,
  obsproj = TRUE,
  ncores = NULL
)
}
\arguments{
\item{exp}{A named numeric array of North Atlantic SLP (20N-80N, 80W-40E) 
forecast anomalies from \code{Ano()} or \code{Ano_CrossValid()} with 
dimensions 'time_dim', 'memb_dim', 'ftime_dim', and 'space_dim' at least.
If only NAO of observational data needs to be computed, this parameter can
be left to NULL. The default value is NULL.}

\item{obs}{A named numeric array of North Atlantic SLP (20N-80N, 80W-40E) 
observed anomalies from \code{Ano()} or \code{Ano_CrossValid()} with 
dimensions 'time_dim', 'ftime_dim', and 'space_dim' at least.
If only NAO of experimental data needs to be computed, this parameter can 
be left to NULL. The default value is NULL.}

\item{lat}{A vector of the latitudes of 'exp' and 'obs'.}

\item{lon}{A vector of the longitudes of 'exp' and 'obs'.}

\item{time_dim}{A character string indicating the name of the time dimension
of 'exp' and 'obs'. The default value is 'sdate'.}

\item{memb_dim}{A character string indicating the name of the member 
dimension of 'exp' (and 'obs', optional). If 'obs' has memb_dim, the length
must be 1. The default value is 'member'.}

\item{space_dim}{A vector of two character strings. The first is the dimension
name of latitude of 'ano' and the second is the dimension name of longitude
of 'ano'. The default value is c('lat', 'lon').}

\item{ftime_dim}{A character string indicating the name of the forecast time 
dimension of 'exp' and 'obs'. The default value is 'ftime'.}

\item{ftime_avg}{A numeric vector of the forecast time steps to average
across the target period. If average is not needed, set NULL. The default 
value is 2:4, i.e., from 2nd to 4th forecast time steps.}

\item{obsproj}{A logical value indicating whether to compute the NAO index by
projecting the forecast anomalies onto the leading EOF of observational 
reference (TRUE) or compute the NAO by first computing the leading 
EOF of the forecast anomalies (in cross-validation mode, i.e. leaving the 
year you are evaluating out), and then projecting forecast anomalies onto 
this EOF (FALSE). The default value is TRUE.}

\item{ncores}{An integer indicating the number of cores to use for parallel 
computation. The default value is NULL.}
}
\value{
A list which contains:
\item{exp}{
 A numeric array of forecast NAO index in verification format with the same 
 dimensions as 'exp' except space_dim and ftime_dim. If ftime_avg is NULL, 
 ftime_dim remains.
 }
\item{obs}{
 A numeric array of observed NAO index in verification format with the same
 dimensions as 'obs' except space_dim and ftime_dim. If ftime_avg is NULL,
 ftime_dim remains.
}
}
\description{
Compute the North Atlantic Oscillation (NAO) index based on the leading EOF 
of the sea level pressure (SLP) anomalies over the north Atlantic region 
(20N-80N, 80W-40E). The PCs are obtained by projecting the forecast and 
observed anomalies onto the observed EOF pattern or the forecast 
anomalies onto the EOF pattern of the other years of the forecast. 
By default (ftime_avg = 2:4), NAO() computes the NAO index for 1-month 
lead seasonal forecasts that can be plotted with PlotBoxWhisker(). It returns
cross-validated PCs of the NAO index for forecast (exp) and observations 
(obs) based on the leading EOF pattern.
}
\examples{
# Make up synthetic data
set.seed(1)
exp <- array(rnorm(1620), dim = c(member = 2, sdate = 3, ftime = 5, lat = 6, lon = 9))
set.seed(2)
obs <- array(rnorm(1620), dim = c(member = 1, sdate = 3, ftime = 5, lat = 6, lon = 9))
lat <- seq(20, 80, length.out = 6)
lon <- seq(-80, 40, length.out = 9) 
nao <- NAO(exp = exp, obs = obs, lat = lat, lon = lon)

# plot the NAO index
 \dontrun{
nao$exp <- Reorder(nao$exp, c(2, 1))
nao$obs <- Reorder(nao$obs, c(2, 1))
PlotBoxWhisker(nao$exp, nao$obs, "NAO index, DJF", "NAO index (PC1) TOS",
       monini = 12, yearini = 1985, freq = 1, "Exp. A", "Obs. X")
 }

}
\references{
Doblas-Reyes, F.J., Pavan, V. and Stephenson, D. (2003). The skill of 
 multi-model seasonal forecasts of the wintertime North Atlantic 
 Oscillation. Climate Dynamics, 21, 501-514. 
 DOI: 10.1007/s00382-003-0350-4
}
