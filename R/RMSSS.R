#'Compute root mean square error skill score
#'
#'Compute the root mean square error skill score (RMSSS) between an array of 
#'forecast 'exp' and an array of observation 'obs'. The two arrays should 
#'have the same dimensions except along dat_dim, where the length can be 
#'different, with the number of experiments/models (nexp) and the number of 
#'observational datasets (nobs).\cr
#'RMSSS computes the root mean square error skill score of each jexp in 1:nexp 
#'against each job in 1:nobs which gives nexp * nobs RMSSS for each grid point 
#'of the array.\cr
#'The RMSSS are computed along the time_dim dimension which should correspond
#'to the start date dimension.\cr
#'The p-value and significance test are optionally provided by an one-sided 
#'Fisher test or Random Walk test.\cr
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
#'@param ref A named numerical array of the reference forecast data with at 
#'  least time dimension, or 0 (typical climatological forecast) or 1 
#'  (normalized climatological forecast). If it is an array, the dimensions must
#'  be the same as 'exp' except 'memb_dim' and 'dat_dim'. If there is only one
#'  reference dataset, it should not have dataset dimension. If there is 
#'  corresponding reference for each experiment, the dataset dimension must 
#'  have the same length as in 'exp'. If 'ref' is NULL, the typical 
#'  climatological forecast is used as reference forecast (equivelant to 0.)
#'  The default value is NULL.
#'@param dat_dim A character string indicating the name of dataset (nobs/nexp) 
#'  dimension. The default value is 'dataset'.
#'@param time_dim A character string indicating the name of dimension along  
#'  which the RMSSS are computed. The default value is 'sdate'.
#'@param memb_dim A character string indicating the name of the member dimension
#'  to compute the ensemble mean; it should be set to NULL if the parameter 'exp'
#'  and 'ref' are already the ensemble mean. The default value is NULL.
#'@param pval A logical value indicating whether to compute or not the p-value 
#'  of the test Ho: RMSSS = 0. The default value is TRUE.
#'@param sign A logical value indicating whether to compute or not the 
#'  statistical significance of the test Ho: RMSSS = 0. The default value is 
#'  FALSE.
#'@param alpha A numeric of the significance level to be used in the 
#'  statistical significance test. The default value is 0.05.
#'@param sig_method A character string indicating the significance method. The
#'  options are "one-sided Fisher" (default) and "Random Walk". 
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return
#'A list containing the numeric arrays with dimension:\cr 
#'  c(nexp, nobs, all other dimensions of exp except time_dim).\cr
#'nexp is the number of experiment (i.e., dat_dim in exp), and nobs is the 
#'number of observation (i.e., dat_dim in obs). If dat_dim is NULL, nexp and 
#'nobs are omitted.\cr
#'\item{$rmsss}{
#'  A numerical array of the root mean square error skill score. 
#'}
#'\item{$p.val}{
#'  A numerical array of the p-value with the same dimensions as $rmsss.
#'  Only present if \code{pval = TRUE}.
#'}
#'\item{sign}{
#'  A logical array of the statistical significance of the RMSSS with the same
#'  dimensions as $rmsss. Only present if \code{sign = TRUE}.
#'}
#'
#'@examples
#' set.seed(1)
#' exp <- array(rnorm(30), dim = c(dataset = 2, time = 3, memb = 5))
#' set.seed(2)
#' obs <- array(rnorm(15), dim = c(time = 3, memb = 5, dataset = 1))
#' res <- RMSSS(exp, obs, time_dim = 'time', dat_dim = 'dataset')
#'
#'@rdname RMSSS
#'@import multiApply
#'@importFrom stats pf
#'@export
RMSSS <- function(exp, obs, ref = NULL, time_dim = 'sdate', dat_dim = 'dataset',
                  memb_dim = NULL, pval = TRUE, sign = FALSE, alpha = 0.05, 
                  sig_method = 'one-sided Fisher', ncores = NULL) {
  
  # Check inputs 
  ## exp, obs, and ref (1)
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
  if (!is.null(ref)) {
    if (!is.numeric(ref)) {
      stop("Parameter 'ref' must be numeric.")
    }
    if (is.array(ref)) {
      if (any(is.null(names(dim(ref))))| any(nchar(names(dim(ref))) == 0)) {
        stop("Parameter 'ref' must have dimension names.")
      }
    } else if (length(ref) != 1 | any(!ref %in% c(0, 1))) {
      stop("Parameter 'ref' must be a numeric array or number 0 or 1.")
    }
  }

  ## time_dim
  if (!is.character(time_dim) | length(time_dim) > 1) {
    stop("Parameter 'time_dim' must be a character string.")
  }
  if (!time_dim %in% names(dim(exp)) | !time_dim %in% names(dim(obs))) {
    stop("Parameter 'time_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## dat_dim
  if (!is.null(dat_dim)) {
    if (!is.character(dat_dim) | length(dat_dim) > 1) {
      stop("Parameter 'dat_dim' must be a character string or NULL.")
    }
    if (!dat_dim %in% names(dim(exp)) | !dat_dim %in% names(dim(obs))) {
      stop("Parameter 'dat_dim' is not found in 'exp' or 'obs' dimension.",
           " Set it as NULL if there is no dataset dimension.")
    }
  }
  ## memb_dim
  if (!is.null(memb_dim)) {
    if (!is.character(memb_dim) | length(memb_dim) > 1) {
      stop("Parameter 'memb_dim' must be a character string.")
    }
    if (!memb_dim %in% names(dim(exp))) {
      stop("Parameter 'memb_dim' is not found in 'exp' dimension.")
    }
    if (memb_dim %in% names(dim(obs))) {
      if (identical(as.numeric(dim(obs)[memb_dim]), 1)) {
        obs <- ClimProjDiags::Subset(x = obs, along = memb_dim, indices = 1, drop = 'selected')
      } else {
        stop("Not implemented for observations with members ('obs' can have 'memb_dim', ",
             "but it should be of length = 1).")
      }
    }
  }
  ## pval
  if (!is.logical(pval) | length(pval) > 1) {
    stop("Parameter 'pval' must be one logical value.")
  }
  ## sign
  if (!is.logical(sign) | length(sign) > 1) {
    stop("Parameter 'sign' must be one logical value.")
  }
  ## alpha
  if (!is.numeric(alpha) | length(alpha) > 1) {
    stop("Parameter 'alpha' must be one numeric value.")
  }
  ## sig_method
  if (length(sig_method) != 1 | !any(sig_method %in% c('one-sided Fisher', 'Random Walk'))) {
    stop("Parameter 'sig_method' must be one of 'one-sided Fisher' or 'Random Walk'.")
  }
  if (sig_method == "Random Walk" & pval == T) {
    warning("p-value cannot be calculated by significance method 'Random Walk'.")
    pval <- FALSE
  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 |
      length(ncores) > 1) {
      stop("Parameter 'ncores' must be a positive integer.")
    }
  } 
  ## exp and obs (2)
  name_exp <- sort(names(dim(exp)))
  name_obs <- sort(names(dim(obs)))
  if (!is.null(memb_dim)) {
    name_exp <- name_exp[-which(name_exp == memb_dim)]
  }
  if (!is.null(dat_dim)) {
  name_exp <- name_exp[-which(name_exp == dat_dim)]
  name_obs <- name_obs[-which(name_obs == dat_dim)]
  }
  if(!all(dim(exp)[name_exp] == dim(obs)[name_obs])) {
    stop(paste0("Parameter 'exp' and 'obs' must have same length of ",
                "all dimension except 'memb_dim' and 'dat_dim'."))
  }
  if (!is.null(ref)) {
    name_ref <- sort(names(dim(ref)))
    if (!is.null(memb_dim) && memb_dim %in% name_ref) {
      name_ref <- name_ref[-which(name_ref == memb_dim)]
    }
    if (!is.null(dat_dim)) {
      if (dat_dim %in% name_ref) {
        if (!identical(dim(exp)[dat_dim], dim(ref)[dat_dim])) {
          stop(paste0("If parameter 'ref' has dataset dimension, it must be ", 
                      "equal to dataset dimension of 'exp'."))
        }
        name_ref <- name_ref[-which(name_ref == dat_dim)]
      }
    }
    if (!identical(length(name_exp), length(name_ref)) |
        !identical(dim(exp)[name_exp], dim(ref)[name_ref])) {
      stop(paste0("Parameter 'exp' and 'ref' must have the same length of ",
                  "all dimensions except 'memb_dim' and 'dat_dim' if there is ",
                  "only one reference dataset."))
    }
  }

  if (dim(exp)[time_dim] <= 2) {
    stop("The length of time_dim must be more than 2 to compute RMSSS.")
  }
 

  ###############################
#  # Sort dimension
#  name_exp <- names(dim(exp))
#  name_obs <- names(dim(obs))
#  order_obs <- match(name_exp, name_obs)
#  obs <- Reorder(obs, order_obs)


  ###############################
  # Create ref array if needed
  if (is.null(ref)) ref <- 0
  if (!is.array(ref)) {
    ref <- array(data = ref, dim = dim(exp))
  }

  ###############################
  ## Ensemble mean
  if (!is.null(memb_dim)) {
    exp <- MeanDims(exp, memb_dim, na.rm = T)
    if (!is.null(ref) & memb_dim %in% names(dim(ref))) {
      ref <- MeanDims(ref, memb_dim, na.rm = T)
    }
  }

  ###############################
  # Calculate RMSSS

#  if (!is.null(ref)) { # use "ref" as reference forecast
#    if (!is.null(dat_dim) && dat_dim %in% names(dim(ref))) {
#      target_dims_ref <- c(time_dim, dat_dim)
#    } else {
#      target_dims_ref <- c(time_dim)
#    }
#    data <- list(exp = exp, obs = obs, ref = ref)
#    target_dims = list(exp = c(time_dim, dat_dim),
#                       obs = c(time_dim, dat_dim),
#                       ref = target_dims_ref)
#  } else {
#    data <- list(exp = exp, obs = obs)
#    target_dims = list(exp = c(time_dim, dat_dim),
#                       obs = c(time_dim, dat_dim))
#  }
  data <- list(exp = exp, obs = obs, ref = ref)
  if (!is.null(dat_dim)) {
    if (dat_dim %in% names(dim(ref))) {
       target_dims <- list(exp = c(time_dim, dat_dim),
                           obs = c(time_dim, dat_dim),
                           ref = c(time_dim, dat_dim))
    } else {
       target_dims <- list(exp = c(time_dim, dat_dim),
                           obs = c(time_dim, dat_dim),
                           ref = c(time_dim))
    }
  } else {
    target_dims <- list(exp = time_dim, obs = time_dim, ref = time_dim)
  }
    
  res <- Apply(data, 
               target_dims = target_dims,
               fun = .RMSSS, 
               time_dim = time_dim, dat_dim = dat_dim,
               pval = pval, sign = sign, alpha = alpha,
               sig_method = sig_method,
               ncores = ncores)
  
  return(res)
}

.RMSSS <- function(exp, obs, ref = NULL, time_dim = 'sdate', dat_dim = 'dataset', pval = TRUE, 
                   sign = FALSE, alpha = 0.05, sig_method = 'one-sided Fisher') {
  # exp: [sdate, (dat)]
  # obs: [sdate, (dat)]
  # ref: [sdate, (dat)] or NULL

  if (is.null(ref)) {
    ref <- array(data = 0, dim = dim(obs))
  } else if (identical(ref, 0) | identical(ref, 1)) {
    ref <- array(ref, dim = dim(exp))
  }

  if (is.null(dat_dim)) {
    # exp: [sdate]
    # obs: [sdate]
    nexp <- 1
    nobs <- 1
    nref <- 1
    # Add dat dim back temporarily
    dim(exp) <- c(dim(exp), dat = 1)
    dim(obs) <- c(dim(obs), dat = 1)
    dim(ref) <- c(dim(ref), dat = 1)

  } else {
  # exp: [sdate, dat_exp]
  # obs: [sdate, dat_obs]
    nexp <- as.numeric(dim(exp)[2])
    nobs <- as.numeric(dim(obs)[2])
    if (dat_dim %in% names(dim(ref))) {
      nref <- as.numeric(dim(ref)[2])
    } else {
      dim(ref) <- c(dim(ref), dat = 1)
      nref <- 1
    }
  }

  nsdate <- as.numeric(dim(exp)[1])

  # RMS of forecast
  dif1 <- array(dim = c(nsdate, nexp, nobs))
  names(dim(dif1)) <- c(time_dim, 'nexp', 'nobs')

  for (i in 1:nobs) {
    dif1[, , i] <- sapply(1:nexp, function(x) {exp[, x] - obs[, i]})
  }

  rms_exp <- apply(dif1^2, c(2, 3), mean, na.rm = TRUE)^0.5  #array(dim = c(nexp, nobs))

  # RMS of reference
#  if (!is.null(ref)) {
    dif2 <- array(dim = c(nsdate, nref, nobs))
    names(dim(dif2)) <- c(time_dim, 'nexp', 'nobs')
    for (i in 1:nobs) {
      dif2[, , i] <- sapply(1:nref, function(x) {ref[, x] - obs[, i]})
    }
    rms_ref <- apply(dif2^2, c(2, 3), mean, na.rm = TRUE)^0.5  #array(dim = c(nref, nobs))
    if (nexp != nref) {
      # expand rms_ref to nexp (nref is 1)
      rms_ref <- array(rms_ref, dim = c(nobs = nobs, nexp = nexp))
      rms_ref <- Reorder(rms_ref, c(2, 1))
    }
#  } else {
#    rms_ref <- array(colMeans(obs^2, na.rm = TRUE)^0.5, dim = c(nobs = nobs, nexp = nexp))
##    rms_ref[which(abs(rms_ref) <= (max(abs(rms_ref), na.rm = TRUE) / 1000))] <- max(abs(
##    rms_ref), na.rm = TRUE) / 1000
#    rms_ref <- Reorder(rms_ref, c(2, 1))
#    #rms_ref above: [nexp, nobs]
#  }

  rmsss <- 1 - rms_exp / rms_ref

#################################################

#  if (conf) {
#    conflow <- (1 - conf.lev) / 2
#    confhigh <- 1 - conflow
#    conf_low <- array(dim = c(nexp = nexp, nobs = nobs))
#    conf_high <- array(dim = c(nexp = nexp, nobs = nobs))
#  }

  if (sig_method == 'one-sided Fisher') {
    p_val <- array(dim = c(nexp = nexp, nobs = nobs))
    ## pval and sign 
    if (pval || sign) {
      eno1 <- Eno(dif1, time_dim)
      if (is.null(ref)) {
        eno2 <- Eno(obs, time_dim) 
        eno2 <- array(eno2, dim = c(nobs = nobs, nexp = nexp))
        eno2 <- Reorder(eno2, c(2, 1))
      } else {
        eno2 <- Eno(dif2, time_dim)
        if (nref != nexp) {
          eno2 <- array(eno2, dim = c(nobs = nobs, nexp = nexp))
          eno2 <- Reorder(eno2, c(2, 1))
        }
      }

      F.stat <- (eno2 * rms_ref^2 / (eno2 - 1)) / ((eno1 * rms_exp^2 / (eno1- 1)))
      tmp <- !is.na(eno1) & !is.na(eno2) & eno1 > 2 & eno2 > 2
      p_val <- 1 - pf(F.stat, eno1 - 1, eno2 - 1)
      if (sign) signif <- p_val <= alpha 
      # If there isn't enough valid data, return NA
      p_val[which(!tmp)] <- NA
      if (sign) signif[which(!tmp)] <- NA
    
      # change not enough valid data rmsss to NA
      rmsss[which(!tmp)] <- NA
    }

  } else if (sig_method == "Random Walk") {
    signif <- array(dim = c(nexp = nexp, nobs = nobs))
    for (i in 1:nexp) {
      for (j in 1:nobs) {

        # Error
        error_exp <- array(data = abs(exp[, i] - obs[, j]), dim = c(time = nsdate))
        if (nref == nexp) {
          error_ref <- array(data = abs(ref[, i] - obs[, j]), dim = c(time = nsdate))
        } else {
          # nref = 1
          error_ref <- array(data = abs(ref - obs[, j]), dim = c(time = nsdate))
        }
        signif[i, j] <- .RandomWalkTest(skill_A = error_exp, skill_B = error_ref)$sign
      }
    }
  }

  ###################################
  # Remove extra dimensions if dat_dim = NULL
  if (is.null(dat_dim)) {
    dim(rmsss) <- NULL
    dim(p_val) <- NULL
    if (sign) dim(signif) <- NULL
  }
  ###################################

  # output  
  res <- list(rmsss = rmsss)
  if (pval) {
    p.val <- list(p.val = p_val)
    res <- c(res, p.val)
  }
  if (sign) {
    signif <- list(sign = signif)
    res <- c(res, signif)
  }

  return(res)
}
