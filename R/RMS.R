#'Compute root mean square error
#'
#'Compute the root mean square error for an array of forecasts and an array of
#'observations. The RMSEs are computed along time_dim, the dimension which 
#'corresponds to the startdate dimension. If comp_dim is given, the RMSEs are 
#'computed only if obs along the comp_dim dimension are complete between 
#'limits[1] and limits[2], i.e. there are no NAs between limits[1] and 
#'limits[2]. This option can be activated if the user wishes to account only 
#'for the forecasts for which the corresponding observations are available at 
#'all leadtimes.\cr
#'The confidence interval is computed by the chi2 distribution.\cr
#'
#'@param exp A named numeric array of experimental data, with at least two 
#'  dimensions 'time_dim' and 'memb_dim'.
#'@param obs A named numeric array of observational data, same dimensions as  
#'  parameter 'exp' except along memb_dim.
#'@param time_dim A character string indicating the name of dimension along  
#'  which the correlations are computed. The default value is 'sdate'.
#'@param memb_dim A character string indicating the name of member (nobs/nexp) 
#'  dimension. The default value is 'member'.
#'@param comp_dim A character string indicating the name of dimension along which
#'  obs is taken into account only if it is complete. The default value
#'  is NULL.
#'@param limits A vector of two integers indicating the range along comp_dim to 
#'  be completed. The default value is c(1, length(comp_dim dimension)).
#'@param conf A logical value indicating whether to retrieve the confidence 
#'  intervals or not. The default value is TRUE.
#'@param conf.lev A numeric indicating the confidence level for the 
#'  regression computation. The default value is 0.95.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return 
#'A list containing the numeric arrays with dimension:\cr 
#'  c(nexp, nobs, all other dimensions of exp except time_dim).\cr
#'nexp is the number of experiment (i.e., memb_dim in exp), and nobs is the 
#'number of observation (i.e., memb_dim in obs).\cr
#'\item{$rms}{
#'  The root mean square error. 
#'}
#'\item{$conf.lower}{
#'  The lower confidence interval. Only present if \code{conf = TRUE}.
#'}
#'\item{$conf.upper}{
#'  The upper confidence interval. Only present if \code{conf = TRUE}.
#'}
#'
#'@examples
#'  set.seed(1)
#'  exp1 <- array(rnorm(120), dim = c(member = 3, sdate = 5, ftime = 2, lon = 1, lat = 4))
#'  set.seed(2)
#'  obs1 <- array(rnorm(80),  dim = c(member = 2, sdate = 5, ftime = 2, lon = 1, lat = 4))
#'  set.seed(2)
#'  na <- floor(runif(10, min = 1, max = 80))
#'  obs1[na] <- NA
#'  res <- RMS(exp1, obs1, comp_dim = 'ftime')
#'
#'@rdname RMS
#'@import multiApply
#'@importFrom ClimProjDiags Subset
#'@importFrom stats qchisq
#'@export
RMS <- function(exp, obs, time_dim = 'sdate', memb_dim = 'member',
                comp_dim = NULL, limits = NULL, 
                conf = TRUE, conf.lev = 0.95, ncores = NULL) {       
  # Check inputs 
  ## exp and obs (1)
  if (is.null(exp) | is.null(obs)) {
    stop("Parameter 'exp' and 'obs' cannot be NULL.")
  }
  if (!is.numeric(exp) | !is.numeric(obs)) {
    stop("Parameter 'exp' and 'obs' must be a numeric array.")
  }
  if (is.null(dim(exp)) | is.null(dim(obs))) {
    stop(paste0("Parameter 'exp' and 'obs' must be at least two dimensions ",
                "containing time_dim and memb_dim."))
  }
  if(any(is.null(names(dim(exp))))| any(nchar(names(dim(exp))) == 0) |
     any(is.null(names(dim(obs))))| any(nchar(names(dim(obs))) == 0)) {
    stop("Parameter 'exp' and 'obs' must have dimension names.")
  }
  if(!all(names(dim(exp)) %in% names(dim(obs))) | 
     !all(names(dim(obs)) %in% names(dim(exp)))) {
    stop("Parameter 'exp' and 'obs' must have same dimension name")
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) > 1) {
    stop("Parameter 'time_dim' must be a character string.")
  }
  if (!time_dim %in% names(dim(exp)) | !time_dim %in% names(dim(obs))) {
    stop("Parameter 'time_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## memb_dim
  if (!is.character(memb_dim) | length(memb_dim) > 1) {
    stop("Parameter 'memb_dim' must be a character string.")
  }
  if (!memb_dim %in% names(dim(exp)) | !memb_dim %in% names(dim(obs))) {
    stop("Parameter 'memb_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## comp_dim
  if (!is.null(comp_dim)) {
    if (!is.character(comp_dim) | length(comp_dim) > 1) {
      stop("Parameter 'comp_dim' must be a character string.")
    }
    if (!comp_dim %in% names(dim(exp)) | !comp_dim %in% names(dim(obs))) {
      stop("Parameter 'comp_dim' is not found in 'exp' or 'obs' dimension.")
    }
  }
  ## limits
  if (!is.null(limits)) {
    if (is.null(comp_dim)) {
      stop("Paramter 'comp_dim' cannot be NULL if 'limits' is assigned.")
    }
    if (!is.numeric(limits) | any(limits %% 1 != 0) | any(limits < 0) | 
        length(limits) != 2 | any(limits > dim(exp)[comp_dim])) {
      stop(paste0("Parameter 'limits' must be a vector of two positive ",
                  "integers smaller than the length of paramter 'comp_dim'."))
    }
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
  ## exp and obs (2)
  name_exp <- sort(names(dim(exp)))
  name_obs <- sort(names(dim(obs)))
  name_exp <- name_exp[-which(name_exp == memb_dim)]
  name_obs <- name_obs[-which(name_obs == memb_dim)]
  if(!all(dim(exp)[name_exp] == dim(obs)[name_obs])) {
    stop(paste0("Parameter 'exp' and 'obs' must have same length of ",
                "all dimension expect 'memb_dim'."))
  }
  if (dim(exp)[time_dim] < 2) {
    stop("The length of time_dim must be at least 2 to compute RMS.")
  }


  ###############################
  # Sort dimension
  name_exp <- names(dim(exp))
  name_obs <- names(dim(obs))
  order_obs <- match(name_exp, name_obs)
  obs <- Reorder(obs, order_obs)


  ###############################
  # Calculate RMS
  
  #  Remove data along comp_dim dim if there is at least one NA between limits
  if (!is.null(comp_dim)) {
    if (is.null(limits)) {
      limits <- c(1, dim(obs)[comp_dim])
    }
    pos <- which(names(dim(obs)) == comp_dim)
    obs_sub <- Subset(obs, pos, list(limits[1]:limits[2]))
    outrows <- is.na(MeanDims(obs_sub, pos, na.rm = FALSE))
    outrows <- InsertDim(outrows, pos, dim(obs)[comp_dim])
    obs[which(outrows)] <- NA
  }
  
  res <- Apply(list(exp, obs), 
               target_dims = list(c(time_dim, memb_dim), 
                                  c(time_dim, memb_dim)),
               fun = .RMS, 
               time_dim = time_dim, memb_dim = memb_dim,
               conf = conf, conf.lev = conf.lev, ncores = ncores)
  return(res)
}

.RMS <- function(exp, obs, time_dim = 'sdate', memb_dim = 'member',
                 conf = TRUE, conf.lev = 0.95) { 

  # exp: [sdate, member_exp]
  # obs: [sdate, member_obs]
  n_exp <- as.numeric(dim(exp)[2])
  n_obs <- as.numeric(dim(obs)[2])
  n_sdate <- as.numeric(dim(exp)[1])

  dif <- array(dim = c(sdate = n_sdate, n_exp = n_exp, n_obs = n_obs))
  chi <- array(dim = c(nexp = n_exp, nobs = n_obs))
  if (conf) {
    conflow <- (1 - conf.lev) / 2
    confhigh <- 1 - conflow
    conf.lower <- array(dim = c(nexp = n_exp, nobs = n_obs))
    conf.upper <- array(dim = c(nexp = n_exp, nobs = n_obs))
  }

  # dif
  for (i in 1:n_obs) {
    dif[, , i] <- sapply(1:n_exp, function(x) {exp[, x] - obs[, i]})
  }
  rms <- apply(dif^2, c(2, 3), mean, na.rm = TRUE)^0.5  #array(dim = c(n_exp, n_obs))

  if (conf) {
    #eno <- Eno(dif, 1) #count effective sample along sdate. dim = c(n_exp, n_obs)
    eno <- Eno(dif, time_dim) #change to this line when Eno() is done

    # conf.lower
    chi <- sapply(1:n_obs, function(i) {
                             qchisq(confhigh, eno[, i] - 1)
                           })
    conf.lower <- (eno * rms ** 2 / chi) ** 0.5

    # conf.upper
    chi <- sapply(1:n_obs, function(i) {
                             qchisq(conflow, eno[, i] - 1)
                           })
    conf.upper <- (eno * rms ** 2 / chi) ** 0.5
  }

  if (conf) {
    res <- list(rms = rms, conf.lower = conf.lower, conf.upper = conf.upper)
  } else {
    res <- list(rms = rms)
  } 

  return(res)

}
