% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Trend.R
\name{Trend}
\alias{Trend}
\title{Compute the trend}
\usage{
Trend(
  data,
  time_dim = "ftime",
  interval = 1,
  polydeg = 1,
  conf = TRUE,
  conf.lev = 0.95,
  pval = TRUE,
  ncores = NULL
)
}
\arguments{
\item{data}{An numeric array including the dimension along which the trend 
is computed.}

\item{time_dim}{A character string indicating the dimension along which to 
compute the trend. The default value is 'ftime'.}

\item{interval}{A positive numeric indicating the unit length between two 
points along 'time_dim' dimension. The default value is 1.}

\item{polydeg}{A positive integer indicating the degree of polynomial 
regression. The default value is 1.}

\item{conf}{A logical value indicating whether to retrieve the confidence 
intervals or not. The default value is TRUE.}

\item{conf.lev}{A numeric indicating the confidence level for the 
regression computation. The default value is 0.95.}

\item{pval}{A logical value indicating whether to compute the p-value or not. 
The default value is TRUE.}

\item{ncores}{An integer indicating the number of cores to use for parallel 
computation. The default value is NULL.}
}
\value{
A list containing:
\item{$trend}{
 A numeric array with the first dimension 'stats', followed by the same 
 dimensions as parameter 'data' except the 'time_dim' dimension. The length
 of the 'stats' dimension should be \code{polydeg + 1}, containing the 
 regression coefficients from the lowest order (i.e., intercept) to the 
 highest degree.   
}
\item{$conf.lower}{
 A numeric array with the first dimension 'stats', followed by the same 
 dimensions as parameter 'data' except the 'time_dim' dimension. The length
 of the 'stats' dimension should be \code{polydeg + 1}, containing the 
 lower limit of the \code{conf.lev}\% confidence interval for all the 
 regression coefficients with the same order as \code{$trend}. Only present 
 \code{conf = TRUE}.
}
\item{$conf.upper}{
 A numeric array with the first dimension 'stats', followed by the same 
 dimensions as parameter 'data' except the 'time_dim' dimension. The length
 of the 'stats' dimension should be \code{polydeg + 1}, containing the 
 upper limit of the \code{conf.lev}\% confidence interval for all the 
 regression coefficients with the same order as \code{$trend}. Only present 
 \code{conf = TRUE}.
}
\item{$p.val}{
 The p-value calculated by anova(). Only present if \code{pval = TRUE}.
}
\item{$detrended}{
 A numeric array with the same dimensions as paramter 'data', containing the 
 detrended values along the 'time_dim' dimension.
}
}
\description{
Compute the linear trend or any degree of polynomial regression along the 
forecast time. It returns the regression coefficients (including the intercept)
and the detrended array. The confidence intervals and p-value are also 
provided if needed.\cr
The confidence interval relies on the student-T distribution, and the p-value 
is calculated by ANOVA.
}
\examples{
# Load sample data as in Load() example:
example(Load)
months_between_startdates <- 60
trend <- Trend(sampleData$obs, polydeg = 2, interval = months_between_startdates)

}
