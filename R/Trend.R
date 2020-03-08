#'Compute the trend 
#'
#'Compute the linear trend or any degree of polynomial regression along the 
#'forecast time. It returns the regression coefficients (including the intercept)
#'and the confidence intervals if needed. The detrended array is also provided.\cr
#'The confidence interval relies on the student-T distribution.\cr\cr
#'
#'@param data An numeric array including the dimension along which the trend 
#'  is computed.
#'@param time_dim A character string indicating the dimension along which to 
#'  compute the trend. The default value is 'sdate'.
#'@param interval A positive numeric indicating the unit length between two 
#' points along 'time_dim' dimension. The default value is 1.
#'@param polydeg A positive integer indicating the degree of polynomial 
#'  regression. The default value is 1.
#'@param conf A logical value indicating whether to retrieve the confidence 
#'  intervals or not. The default value is TRUE.
#'@param conf.lev A numeric indicating the confidence level for the 
#'  regression computation. The default value is 0.95.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return 
#'A list containing:
#'\item{$trend}{
#'  A numeric array with the first dimension 'stats', followed by the same 
#'  dimensions as parameter 'data' except the 'time_dim' dimension. The length
#'  of the 'stats' dimension should be \code{polydeg + 1}, containing the 
#'  regression coefficients from the lowest order (i.e., intercept) to the 
#'  highest degree.   
#'}
#'\item{$conf.lower}{
#'  A numeric array with the first dimension 'stats', followed by the same 
#'  dimensions as parameter 'data' except the 'time_dim' dimension. The length
#'  of the 'stats' dimension should be \code{polydeg + 1}, containing the 
#'  lower limit of the \code{conf.lev}\% confidence interval for all the 
#'  regression coefficients with the same order as \code{$trend}. Only present 
#'  \code{conf = TRUE}.
#'}
#'\item{$conf.upper}{
#'  A numeric array with the first dimension 'stats', followed by the same 
#'  dimensions as parameter 'data' except the 'time_dim' dimension. The length
#'  of the 'stats' dimension should be \code{polydeg + 1}, containing the 
#'  upper limit of the \code{conf.lev}\% confidence interval for all the 
#'  regression coefficients with the same order as \code{$trend}. Only present 
#'  \code{conf = TRUE}.
#'}
#'\item{$detrended}{
#'  A numeric array with the same dimensions as paramter 'data', containing the 
#'  detrended values along the 'time_dim' dimension.
#'}
#'
#'@examples
#'# Load sample data as in Load() example:
#'example(Load)
#'months_between_startdates <- 60
#'trend <- Trend(sampleData$obs, polydeg = 2)
#'
#'@rdname Trend
#'@import multiApply
#'@export
Trend <- function(data, time_dim = 'sdate', interval = 1, polydeg = 1,
                  conf = TRUE, conf.lev = 0.95, ncores = NULL) {

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
  ## interval
  if (any(!is.numeric(interval) | interval <= 0 | length(interval) > 1)) {
    stop("Parameter 'interval' must be a positive number.")
  }
  ## polydeg
  if (!is.numeric(polydeg) | polydeg %% 1 != 0 | polydeg < 0 |
      length(polydeg) > 1) {
    stop("Parameter 'polydeg' must be a positive integer.")
  }
  ## conf
  if (!is.logical(conf) | length(conf) > 1) {
    stop("Parameter 'conf' must be one logical value.")
  }
  ## conf.lev
  if (!is.numeric(conf.lev) | conf.lev < 0 | conf.lev > 1 | length(conf.lev) > 1) {
    stop("Parameter 'conf.lev' must be a numeric number between 0 and 1.")
  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores < 0 |
      length(ncores) > 1) {
      stop("Parameter 'ncores' must be a positive integer.")
    }
  }

  ###############################
  # Calculate Trend
  dim_names <- names(dim(data))

  if (conf) {
    output_dims <- list(trend = 'stats', conf.lower = 'stats',
                        conf.upper = 'stats', detrended = time_dim)
  } else if (!conf) {
    output_dims <- list(trend = 'stats', detrended = time_dim)
  }


  output <- Apply(list(data),
                  target_dims = time_dim,
                  fun = .Trend,
                  output_dims = output_dims,
                  time_dim = time_dim, interval = interval, 
                  polydeg = polydeg, conf = conf,
                  conf.lev = conf.lev,
                  ncores = ncores)

  #output <- lapply(output, .reorder, time_dim = time_dim, dim_names = dim_names)

  return(output)
}

.Trend <- function(x, time_dim = 'sdate', interval = 1, polydeg = 1,
                   conf = TRUE, conf.lev = 0.95) {

  mon <- seq(x) * interval

  # remove NAs for potential poly()
  NApos <- 1:length(x)
  NApos[which(is.na(x))] <- NA
  x2 <- x[!is.na(NApos)]
  mon2 <- mon[!is.na(NApos)]

  if (length(x2) > 0) {
#    lm.out <- lm(x ~ mon, na.action = na.omit)
    lm.out <- lm(x2 ~ poly(mon2, degree = polydeg, raw = TRUE), na.action = na.omit)
    trend <- lm.out$coefficients  #intercept, slope1, slope2,...

    if (conf) {
      conf.lower <- confint(lm.out, level = conf.lev)[, 1]
      conf.upper <- confint(lm.out, level = conf.lev)[, 2]
    }

    detrended <- c()
    detrended[is.na(x) == FALSE] <- x[is.na(x) == FALSE] - lm.out$fitted.values
  } else {
    trend <- rep(NA, polydeg + 1)
    detrend <- NA
    if (conf) {
      conf.lower <- rep(NA, polydeg + 1)
      conf.upper <- rep(NA, polydeg + 1)
    }
  }

  if (conf) {
    return(list(trend = trend, conf.lower = conf.lower, conf.upper = conf.upper, 
                detrended = detrended))
  } else {
    return(list(trend = trend, detrended = detrended))
  }

}
