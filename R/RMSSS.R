#'Compute root mean square error skill score
#'
#'Compute the root mean square error skill score (RMSSS) between an array of 
#'forecast 'exp' and an array of observation 'obs'. The two arrays should 
#'have the same dimensions except along dat_dim, where the length can be 
#'different, with the number of experiments/models (nexp) and the number of 
#'observational datasets (nobs).\cr
#'RMSSS computes the root mean square error skill score of each jexp in 1:nexp 
#'against each jobs in 1:nobs which gives nexp * nobs RMSSS for each other 
#'grid point of the array.\cr
#'The RMSSS are computed along the time_dim dimension which should corresponds
#'to the startdate dimension.\cr
#'The p-value is optionally provided by an one-sided Fisher test.\cr
#'
#'@param exp A named numeric array of experimental data which contains at least
#'  two dimensions for dat_dim and time_dim. It can also be a vector with the 
#'  same length as 'obs', then the vector will automatically be 'time_dim' and 
#'  'dat_dim' will be 1.
#'@param obs A named numeric array of observational data which contains at least
#'  two dimensions for dat_dim and time_dim. The dimensions should be the same 
#'  as paramter 'exp' except the length of 'dat_dim' dimension. The order of 
#'  dimension can be different. It can also be a vector with the same length as
#'  'exp', then the vector will automatically be 'time_dim' and 'dat_dim' will
#'  be 1.
#'@param dat_dim A character string indicating the name of dataset (nobs/nexp) 
#'  dimension. The default value is 'dataset'.
#'@param time_dim A character string indicating the name of dimension along  
#'  which the RMSSS are computed. The default value is 'sdate'.
#'@param pval A logical value indicating whether to compute or not the p-value 
#'  of the test Ho: RMSSS = 0. If pval = TRUE, the insignificant RMSSS will 
#'  return NA. The default value is TRUE.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return
#'A list containing the numeric arrays with dimension:\cr 
#'  c(nexp, nobs, all other dimensions of exp except time_dim).\cr
#'nexp is the number of experiment (i.e., dat_dim in exp), and nobs is the 
#'number of observation (i.e., dat_dim in obs).\cr
#'\item{$rmsss}{
#'  The root mean square error skill score. 
#'}
#'\item{$p.val}{
#'  The p-value. Only present if \code{pval = TRUE}.
#'}
#'
#'@examples
#' set.seed(1)
#' exp <- array(rnorm(30), dim = c(dataset = 2, time = 3, memb = 5))
#' set.seed(2)
#' obs <- array(rnorm(15), dim = c(time = 3, memb = 5, dataset = 1))
#' res <- RMSSS(exp, obs, time_dim = 'time')
#'
#'@rdname RMSSS
#'@import multiApply
#'@importFrom stats pf
#'@export
RMSSS <- function(exp, obs, time_dim = 'sdate', dat_dim = 'dataset', 
                  pval = TRUE, ncores = NULL) {
  
  # Check inputs 
  ## exp and obs (1)
  if (is.null(exp) | is.null(obs)) {
    stop("Parameter 'exp' and 'obs' cannot be NULL.")
  }
  if (!is.numeric(exp) | !is.numeric(obs)) {
    stop("Parameter 'exp' and 'obs' must be a numeric array.")
  }
  if (is.null(dim(exp)) & is.null(dim(obs))) {  #is vector
    if (length(exp) == length(obs)) {
      exp <- array(exp, dim = c(length(exp), 1))
      names(dim(exp)) <- c(time_dim, dat_dim)
      obs <- array(obs, dim = c(length(obs), 1))
      names(dim(obs)) <- c(time_dim, dat_dim)
    } else {
    stop(paste0("Parameter 'exp' and 'obs' must be array with as least two ",
                "dimensions time_dim and dat_dim, or vector of same length."))
    }
  } else if (is.null(dim(exp)) | is.null(dim(obs))) {
    stop(paste0("Parameter 'exp' and 'obs' must be array with as least two ",
                "dimensions time_dim and dat_dim, or vector of same length."))
  }
  if(any(is.null(names(dim(exp))))| any(nchar(names(dim(exp))) == 0) |
     any(is.null(names(dim(obs))))| any(nchar(names(dim(obs))) == 0)) {
    stop("Parameter 'exp' and 'obs' must have dimension names.")
  }
  if(!all(names(dim(exp)) %in% names(dim(obs))) | 
     !all(names(dim(obs)) %in% names(dim(exp)))) {
    stop("Parameter 'exp' and 'obs' must have same dimension name.")
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) > 1) {
    stop("Parameter 'time_dim' must be a character string.")
  }
  if (!time_dim %in% names(dim(exp)) | !time_dim %in% names(dim(obs))) {
    stop("Parameter 'time_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## dat_dim
  if (!is.character(dat_dim) | length(dat_dim) > 1) {
    stop("Parameter 'dat_dim' must be a character string.")
  }
  if (!dat_dim %in% names(dim(exp)) | !dat_dim %in% names(dim(obs))) {
    stop("Parameter 'dat_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## pval
  if (!is.logical(pval) | length(pval) > 1) {
    stop("Parameter 'pval' must be one logical value.")
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
  name_exp <- name_exp[-which(name_exp == dat_dim)]
  name_obs <- name_obs[-which(name_obs == dat_dim)]
  if(!all(dim(exp)[name_exp] == dim(obs)[name_obs])) {
    stop(paste0("Parameter 'exp' and 'obs' must have same length of ",
                "all dimension expect 'dat_dim'."))
  }
  if (dim(exp)[time_dim] <= 2) {
    stop("The length of time_dim must be more than 2 to compute RMSSS.")
  }
 

  ###############################
  # Sort dimension
  name_exp <- names(dim(exp))
  name_obs <- names(dim(obs))
  order_obs <- match(name_exp, name_obs)
  obs <- Reorder(obs, order_obs)


  ###############################
  # Calculate RMSSS
  
  res <- Apply(list(exp, obs), 
               target_dims = list(c(time_dim, dat_dim), 
                                  c(time_dim, dat_dim)),
               fun = .RMSSS, 
               time_dim = time_dim, dat_dim = dat_dim,
               pval = pval, #conf = conf, conf.lev = conf.lev,
               ncores = ncores)
  
  return(res)
}

.RMSSS <- function(exp, obs, time_dim = 'sdate', dat_dim = 'dataset', pval = TRUE) {
  # exp: [sdate, dat_exp]
  # obs: [sdate, dat_obs]
  nexp <- as.numeric(dim(exp)[2])
  nobs <- as.numeric(dim(obs)[2])
  nsdate <- as.numeric(dim(exp)[1])
  
  p_val <- array(dim = c(nexp = nexp, nobs = nobs))
  dif1 <- array(dim = c(nsdate, nexp, nobs))
  names(dim(dif1)) <- c(time_dim, 'nexp', 'nobs')

#  if (conf) {
#    conflow <- (1 - conf.lev) / 2
#    confhigh <- 1 - conflow
#    conf_low <- array(dim = c(nexp = nexp, nobs = nobs))
#    conf_high <- array(dim = c(nexp = nexp, nobs = nobs))
#  }

  # dif1
  for (i in 1:nobs) {
    dif1[, , i] <- sapply(1:nexp, function(x) {exp[, x] - obs[, i]})
  }

  # rms1 and eno1
  rms1 <- apply(dif1^2, c(2, 3), mean, na.rm = TRUE)^0.5  #array(dim = c(nexp, nobs))
  # rms2 and eno2
  rms2 <- array(colMeans(obs^2, na.rm = TRUE)^0.5, dim = c(nobs = nobs))
  rms2[which(abs(rms2) <= (max(abs(rms2), na.rm = TRUE) / 1000))] <- max(abs(
    rms2), na.rm = TRUE) / 1000
  #rms2 above: [nobs]
  rms2 <- array(rms2, dim = c(nobs = nobs, nexp = nexp))
  #rms2 above: [nobs, nexp]
  rms2 <- Reorder(rms2, c(2, 1))
  #rms2 above: [nexp, nobs]

  # use rms1 and rms2 to calculate rmsss
  rmsss <- 1 - rms1/rms2

  ## pval and conf 
  if (pval) {
    eno1 <- Eno(dif1, time_dim)  
    eno2 <- Eno(obs, time_dim)  
    eno2 <- array(eno2, dim = c(nobs = nobs, nexp = nexp))
    eno2 <- Reorder(eno2, c(2, 1))
  }

  # pval
  if (pval) {

    F.stat <- (eno2 * rms2^2 / (eno2- 1)) / ((eno1 * rms1^2 / (eno1- 1)))
    tmp <- !is.na(eno1) & !is.na(eno2) & eno1 > 2 & eno2 > 2
    p_val <- 1 - pf(F.stat, eno1 - 1, eno2 - 1)
    p_val[which(!tmp)] <- NA
    
    # change not significant rmsss to NA
    rmsss[which(!tmp)] <- NA
  }

  # output  
  if (pval) {
    res <- list(rmsss = rmsss, p.val = p_val)

  } else {
    res <- list(rmsss = rmsss)
  }

  return(res)
}
