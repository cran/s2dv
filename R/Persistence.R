#'Compute persistence 
#'
#'Compute a persistence forecast based on a lagged autoregression of 
#'observational data along the time dimension, with a measure of forecast 
#'uncertainty (prediction interval) based on Coelho et al., 2004.\cr\cr
#'
#'@param data A numeric array corresponding to the observational data
#'  including the time dimension along which the autoregression is computed.
#'  The data should start at least 40 time steps (years or days) before 
#'  'start'.
#'@param dates A sequence of 4-digit integers (YYYY) or dates (YYYY-MM-DD) 
#'  indicating the dates available in the observations.
#'@param time_dim A character string indicating the dimension along which to 
#'  compute the autoregression. The default value is 'time'.
#'@param start A 4-digit integer (YYYY) or a date in the ISOdate format 
#'  (YYYY-MM-DD) indicating the first start date of the persistence forecast.
#'@param end A 4-digit integer (YYYY) or a date in the ISOdate format
#'  (YYYY-MM-DD) indicating the last start date of the persistence forecast.
#'@param ft_start An integer indicating the forecast time for which the 
#'  persistence forecast should be calculated, or the first forecast time of
#'  the average forecast times for which persistence should be calculated.
#'@param ft_end An (optional) integer indicating the last forecast time of the
#'  average forecast times for which persistence should be calculated in the
#'  case of a multi-timestep average persistence. The default value is 
#'  'ft_start'.
#'@param max_ft An integer indicating the maximum forecast time possible for 
#'  'data'. For example, for decadal prediction 'max_ft' would correspond to 10
#'  (years). The default value is 10.
#'@param nmemb An integer indicating the number of ensemble members to generate 
#'  for the persistence forecast. The default value is 1.
#'@param na.action A function or an integer. A function (e.g., na.omit, 
#'  na.exclude, na.fail, na.pass) indicates what should happen when the data 
#'  contain NAs. A numeric indicates the maximum number of NA position allowed 
#'  to compute regression. The default value is 10.
#'@param ncores An integer indicating the number of cores to use for parallel
#'  computation. The default value is NULL.
#'
#'@return 
#'A list containing:
#'\item{$persistence}{
#'  A numeric array with dimensions 'memb', time (start dates), latitudes and longitudes
#'  containing the persistence forecast.   
#'}
#'\item{$persistence.mean}{
#'  A numeric array with same dimensions as 'persistence', except the 'memb' dimension
#'  which is of length 1, containing the ensemble mean persistence forecast.
#'}
#'\item{$persistence.predint}{
#'  A numeric array with same dimensions as 'persistence', except the 'memb' dimension
#'  which is of length 1, containing the prediction interval of the persistence forecast.
#'}
#'\item{$AR.slope}{
#'  A numeric array with same dimensions as 'persistence', except the 'memb' dimension
#'  which is of length 1, containing the slope coefficient of the autoregression.
#'}
#'\item{$AR.intercept}{
#'  A numeric array with same dimensions as 'persistence', except the 'memb' dimension
#'  which is of length 1, containing the intercept coefficient of the autoregression.
#'}
#'\item{$AR.lowCI}{
#'  A numeric array with same dimensions as 'persistence', except the 'memb' dimension
#'  which is of length 1, containing the lower value of the confidence interval of the
#'  autoregression.
#'}
#'\item{$AR.highCI}{
#'  A numeric array with same dimensions as 'persistence', except the 'memb' dimension
#'  which is of length 1, containing the upper value of the confidence interval of the
#'  autoregression.
#'}
#'
#'@examples
#'#Building an example dataset with yearly start dates from 1920 to 2009
#'set.seed(1)
#'obs1 <- rnorm(1 * 70 * 6 * 7)
#'dim(obs1) <- c(member = 1, time = 70, lat = 6, lon = 7)
#'dates <- seq(1940, 2009, 1)
#'persist <- Persistence(obs1, dates = dates, start = 1961, end = 2005, ft_start = 1,
#'                       nmemb = 40)
#'
#'@import multiApply
#'@export
Persistence <- function(data, dates, time_dim = 'time', start, end, ft_start,
                        ft_end = ft_start, max_ft = 10, nmemb = 1, na.action = 10,
                        ncores = NULL) {

  # Check inputs 
  ## data
  if (is.null(data)) {
    stop("Parameter 'data' cannot be NULL.")
  }
  if (!is.numeric(data)) {
    stop("Parameter 'data' must be a numeric array.")
  }
  if (is.null(dim(data))) {  #is vector
    dim(data) <- c(length(data))
    names(dim(data)) <- time_dim
  }
  if(any(is.null(names(dim(data))))| any(nchar(names(dim(data))) == 0)) {
    stop("Parameter 'data' must have dimension names.")
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) > 1) {
    stop("Parameter 'time_dim' must be a character string.")
  }
  if (!time_dim %in% names(dim(data))) {
    stop("Parameter 'time_dim' is not found in 'data' dimension.")
  }
  ## dates
  if (length(dates) != dim(data)[time_dim]) {
    stop("Parameter 'dates' must have the same length as in 'time_dim'.")
  }
  ## start
#  if (!is.numeric(start) | start %% 1 != 0 | start < 0 |
#      length(start) > 1 | start < 1850 | start > 2020) {
#        stop("Parameter 'start' must be an integer between 1850 and 2020.")
#  }
#  if (start < dates[1] + 40) {
#        stop("Parameter 'start' must start at least 40 time steps after the 
#             first start date of 'data'.")
#  }
  ## end
#  if (!is.numeric(end) | end %% 1 != 0 | end < 0 |
#      length(end) > 1 | end < 1850 | end > 2020) {
#        stop("Parameter 'end' must be an integer between 1850 and 2020.")
#  }
#  if (end > dates[length(dates)] + 1) {
#        stop("Parameter 'end' must end at most 1 time step after the 
#             last start date of 'data'.")
#  }
  ## ft_start
  if (!is.numeric(ft_start) | ft_start %% 1 != 0 | ft_start < 0 |
      length(ft_start) > 1) {
        stop("Parameter 'ft_start' must be a positive integer.")
  }
  ## ft_end
  if (!is.numeric(ft_end) | ft_end %% 1 != 0 | ft_end < 0 |
      length(ft_end) > 1 | ft_end > max_ft) {
        stop("Parameter 'ft_end' must be a positive integer below 'max_ft'.")
  }
  ## max_ft
  if (!is.numeric(max_ft) | max_ft %% 1 != 0 | max_ft < 0 |
      length(max_ft) > 1) {
        stop("Parameter 'max_ft' must be a positive integer.")
  }
  ## nmemb
  if (!is.numeric(nmemb) | nmemb %% 1 != 0 | nmemb <= 0 |
      length(nmemb) > 1) {
        stop("Parameter 'nmemb' must be a positive integer.")
  }
  ## na.action
  if (!is.function(na.action) & !is.numeric(na.action)) {
    stop(paste0("Parameter 'na.action' must be a function for NA values or ",
                "a numeric indicating the number of NA values allowed ",
                "before returning NA."))
  }
  if (is.numeric(na.action)) {
    if (any(na.action %% 1 != 0) | any(na.action < 0) | length(na.action) > 1) {
      stop(paste0("Parameter 'na.action' must be a function for NA values or ",
                  "a numeric indicating the number of NA values allowed ",
                  "before returning NA."))
    }
  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 |
      length(ncores) > 1) {
      stop("Parameter 'ncores' must be a positive integer.")
    }
  }

  ###############################
  # Calculate Persistence
  dim_names <- names(dim(data))


  output <- Apply(list(data),
                  target_dims = time_dim,
                  fun = .Persistence,
#                  output_dims = output_dims,
                  dates = dates,
                  time_dim = time_dim,
                  start = start,
                  end = end,
                  ft_start = ft_start,
                  ft_end = ft_end,
                  max_ft = max_ft,
                  nmemb = nmemb, 
                  na.action = na.action,
                  ncores = ncores)

  return(output)
}
# x could be a vector timeseries
# start/end is a year (4-digit numeric) or a date (ISOdate)
# ft_start/ft_end are indices
.Persistence <- function(x, dates, time_dim = 'time', start, end, ft_start = 1, 
           ft_end = 1, max_ft = 10, nmemb = 1, na.action = 10) {

  tm <- end - start + 1
  max_date <- match(start, dates)
  interval <- ft_end - ft_start
  persistence.mean <- persistence.predint <- NULL
  AR.slope <- AR.intercept <- AR.lowCI <- AR.highCI <- NULL
  persistence <- matrix(NA, nrow = nmemb, ncol = tm)
  names(dim(persistence)) <- c('realization', time_dim)

  for (sdate in tm:1){
    min_y = max_ft + ft_start
    max_y = max_date + sdate - 2
    min_x = max_ft        # for extreme case: ex. forecast years 1-10, interval = 9
    max_x = max_date + sdate - 2 - ft_start

    regdates = max_y - min_y + 1

    for (val_x in min_x:max_x) {
      tmp_x <- mean(x[(val_x - interval):val_x])
      if (val_x == min_x){
        obs_x <- tmp_x
      } else {
        obs_x <- c(obs_x, tmp_x)
      }
    }

    for (val_y in min_y:max_y) {
      tmp_y <- mean(x[val_y:(val_y + interval)])
      if (val_y == min_y) {
        obs_y <- tmp_y
      } else {
        obs_y <- c(obs_y, tmp_y)
      }
    }
    reg <- .Regression(obs_y, obs_x,  na.action = na.action)
    len <- length(obs_x)
    a <- reg$regression[2] # slope
    b <- reg$regression[1] # intercept
    CI <- abs(reg$conf.upper[2] - reg$conf.lower[2]) # confidence interval
    stdev_reg <- CI / 1.96
    n <- max_x
    X_sq <- (obs_x[len] - mean(obs_x)) ** 2
    S_sq <- sum((obs_x[1:len] - mean(obs_x)) ** 2)

    persistence.mean[sdate] <- a * mean(x[(max_y - interval):max_y]) + b
    persistence.predint[sdate] <- stdev_reg * sqrt(1 + 1 / n + X_sq / S_sq)
    AR.slope[sdate] <- a
    AR.intercept[sdate] <- b
    AR.lowCI[sdate] <- reg$regression[1]
    AR.highCI[sdate] <- reg$regression[3]
    persistence[ ,sdate] <- rnorm(n = nmemb, mean = persistence.mean[sdate],
                                  sd = persistence.predint[sdate])
  }

  return(list(persistence = persistence, persistence.mean = persistence.mean, 
              persistence.predint = persistence.predint, AR.slope = AR.slope,
              AR.intercept = AR.intercept, AR.lowCI = AR.lowCI,
              AR.highCI = AR.highCI))

}