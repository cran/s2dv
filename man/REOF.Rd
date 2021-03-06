% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/REOF.R
\name{REOF}
\alias{REOF}
\title{Area-weighted empirical orthogonal function analysis with varimax rotation using SVD}
\usage{
REOF(
  ano,
  lat,
  lon,
  ntrunc = 15,
  time_dim = "sdate",
  space_dim = c("lat", "lon"),
  corr = FALSE,
  ncores = NULL
)
}
\arguments{
\item{ano}{A numerical array of anomalies with named dimensions to calculate
REOF. The dimensions must have at least 'time_dim' and 'space_dim'.}

\item{lat}{A vector of the latitudes of 'ano'.}

\item{lon}{A vector of the longitudes of 'ano'.}

\item{ntrunc}{A positive integer of the number of eofs to be kept for varimax 
rotation. This function uses this value as 'neof' too, which is the number
of eofs to return by \code{.EOF()}. The default value is 15. If time length
or the product of latitude length and longitude length is less than 
'ntrunc', 'ntrunc' is equal to the minimum of the three values.}

\item{time_dim}{A character string indicating the name of the time dimension
of 'ano'. The default value is 'sdate'.}

\item{space_dim}{A vector of two character strings. The first is the dimension
name of latitude of 'ano' and the second is the dimension name of longitude
of 'ano'. The default value is c('lat', 'lon').}

\item{corr}{A logical value indicating whether to base on a correlation (TRUE)
or on a covariance matrix (FALSE). The default value is FALSE.}

\item{ncores}{An integer indicating the number of cores to use for parallel 
computation. The default value is NULL.}
}
\value{
A list containing:
\item{REOFs}{
 An array of REOF patterns normalized to 1 (unitless) with dimensions 
 (number of modes, the rest of the dimensions of 'ano' except 
 'time_dim'). Multiplying 'REOFs' by 'RPCs' gives the original 
 reconstructed field.
}
\item{RPCs}{
 An array of principal components with the units of the original field to 
 the power of 2, with dimensions (time_dim, number of modes, the rest of the
 dimensions of 'ano' except 'space_dim'). 
} 
\item{var}{
 An array of the percentage (%) of variance fraction of total variance 
 explained by each mode. The dimensions are (number of modes, the rest of 
 the dimension except 'time_dim' and 'space_dim').
}
\item{wght}{
 An array of the area weighting with dimensions 'space_dim'. It is calculated
 by the square root of cosine of 'lat' and used to compute the fraction of 
 variance explained by each REOFs.
}
}
\description{
Perform an area-weighted EOF analysis with varimax rotation using single 
value decomposition (SVD) based on a covariance matrix or a correlation matrix if 
parameter 'corr' is set to TRUE. The internal s2dv function \code{.EOF()} is used 
internally.
}
\examples{
# This example computes the REOFs along forecast horizons and plots the one 
# that explains the greatest amount of variability. The example data has low  
# resolution so the result may not be explanatory, but it displays how to 
# use this function.
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
ano <- Ano_CrossValid(sampleData$mod, sampleData$obs)
ano <- MeanDims(ano$exp, c('dataset', 'member'))
res <- REOF(ano, lat = sampleData$lat, lon = sampleData$lon, ntrunc = 5)
\dontrun{
PlotEquiMap(eof$EOFs[1, , , 1], sampleData$lat, sampleData$lon)
}

}
\seealso{
EOF
}
