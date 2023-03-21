#'Compute the Relative Operating Characteristic Skill Score
#'
#'The Relative Operating Characteristic Skill Score (ROCSS; Kharin and Zwiers, 
#'2003) is based on the ROC curve, which gives information about the hit rates 
#'against the false-alarm rates for a particular category or event. The ROC
#'curve can be summarized with the area under the ROC curve, known as the ROC
#'score, to provide a skill value for each category. The ROCSS ranges between 
#'minus infinite and 1. A positive ROCSS value indicates that the forecast has
#'higher skill than the reference forecasts, meaning the contrary otherwise.
#'@param exp A named numerical array of the forecast with at least time and
#'  member dimension.
#'@param obs A named numerical array of the observation with at least time 
#'  dimension. The dimensions must be the same as 'exp' except 'memb_dim' and
#'  'dat_dim'.
#'@param ref A named numerical array of the reference forecast data with at 
#'  least time and member dimension. The dimensions must be the same as 'exp' 
#'  except 'memb_dim' and 'dat_dim'. If there is only one reference dataset,
#'  it should not have dataset dimension. If there is corresponding reference 
#'  for each experiement, the dataset dimension must have the same length as in
#'  'exp'. If 'ref' is NULL, the random forecast is used as reference forecast. 
#'  The default value is NULL.
#'@param time_dim A character string indicating the name of the time dimension.
#'  The default value is 'sdate'.
#'@param memb_dim A character string indicating the name of the member dimension
#'  to compute the probabilities of the forecast and the reference forecast. The
#'  default value is 'member'.
#'@param dat_dim A character string indicating the name of dataset dimension. 
#'  The length of this dimension can be different between 'exp' and 'obs'. 
#'  The default value is NULL.
#'@param prob_thresholds A numeric vector of the relative thresholds (from 0 to
#'  1) between the categories. The default value is c(1/3, 2/3), which 
#'  corresponds to tercile equiprobable categories.
#'@param indices_for_clim A vector of the indices to be taken along 'time_dim' 
#'  for computing the thresholds between the probabilistic categories. If NULL,
#'  the whole period is used. The default value is NULL.
#'@param cross.val A logical indicating whether to compute the thresholds 
#'  between probabilistic categories in cross-validation. The default value is 
#'  FALSE.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return
#'A numerical array of ROCSS with the same dimensions as 'exp' excluding 
#''time_dim' and 'memb_dim' dimensions and including 'cat' dimension, which is
#'each category. The length if 'cat' dimension corresponds to the number of  
#'probabilistic categories, i.e., 1 + length(prob_thresholds). If there are 
#'multiple datasets, two additional dimensions 'nexp' and 'nobs' are added.
#'
#'@references 
#'Kharin, V. V. and Zwiers, F. W. (2003): 
#'  https://doi.org/10.1175/1520-0442(2003)016%3C4145:OTRSOP%3E2.0.CO;2
#'
#'@examples
#'exp <- array(rnorm(1000), dim = c(lon = 3, lat = 2, sdate = 60, member = 10))
#'ref <- array(rnorm(1000), dim = c(lon = 3, lat = 2, sdate = 60, member = 10))
#'obs <- array(rnorm(1000), dim = c(lon = 3, lat = 2, sdate = 60))
#'ROCSS(exp = exp, obs = obs) ## random forecast as reference forecast
#'ROCSS(exp = exp, obs = obs, ref = ref) ## ref as reference forecast
#'
#'@import multiApply
#'@importFrom easyVerification EnsRoca
#'@export
ROCSS <- function(exp, obs, ref = NULL, time_dim = 'sdate', memb_dim = 'member',
                  dat_dim = NULL, prob_thresholds = c(1/3, 2/3), indices_for_clim = NULL, 
                  cross.val = FALSE, ncores = NULL) {
  
  # Check inputs
  ## exp, obs, and ref (1)
  if (!is.array(exp) | !is.numeric(exp)) {
    stop("Parameter 'exp' must be a numeric array.")
  }
  if (!is.array(obs) | !is.numeric(obs)) {
    stop("Parameter 'obs' must be a numeric array.")
  }
  if (any(is.null(names(dim(exp))))| any(nchar(names(dim(exp))) == 0) |
      any(is.null(names(dim(obs))))| any(nchar(names(dim(obs))) == 0)) {
    stop("Parameter 'exp' and 'obs' must have dimension names.")
  }
  if (!is.null(ref)) {
    if (!is.array(ref) | !is.numeric(ref))
      stop("Parameter 'ref' must be a numeric array.")
    if (any(is.null(names(dim(ref))))| any(nchar(names(dim(ref))) == 0)) {
      stop("Parameter 'ref' must have dimension names.")
    }
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) != 1) {
    stop("Parameter 'time_dim' must be a character string.")
  }
  if (!time_dim %in% names(dim(exp)) | !time_dim %in% names(dim(obs))) {
    stop("Parameter 'time_dim' is not found in 'exp' or 'obs' dimension.")
  }
  if (!is.null(ref) & !time_dim %in% names(dim(ref))) {
    stop("Parameter 'time_dim' is not found in 'ref' dimension.")
  }
  ## memb_dim
  if (!is.character(memb_dim) | length(memb_dim) > 1) {
    stop("Parameter 'memb_dim' must be a character string.")
  }
  if (!memb_dim %in% names(dim(exp))) {
    stop("Parameter 'memb_dim' is not found in 'exp' dimension.")
  }
  if (!is.null(ref) & !memb_dim %in% names(dim(ref))) {
    stop("Parameter 'memb_dim' is not found in 'ref' dimension.")
  }
  ## dat_dim
  if (!is.null(dat_dim)) {
    if (!is.character(dat_dim) | length(dat_dim) > 1) {
      stop("Parameter 'dat_dim' must be a character string.")
    }
    if (!dat_dim %in% names(dim(exp)) | !dat_dim %in% names(dim(obs))) {
      stop("Parameter 'dat_dim' is not found in 'exp' or 'obs' dimension.",
           " Set it as NULL if there is no dataset dimension.")
    }
   }
  ## exp, obs, and ref (2)
  name_exp <- sort(names(dim(exp)))
  name_obs <- sort(names(dim(obs)))
  name_exp <- name_exp[-which(name_exp == memb_dim)]
  if (memb_dim %in% name_obs) {
    name_obs <- name_obs[-which(name_obs == memb_dim)]
  }
  if (!is.null(dat_dim)) {
    name_exp <- name_exp[-which(name_exp == dat_dim)]
    name_obs <- name_obs[-which(name_obs == dat_dim)]
  }
  if (!identical(length(name_exp), length(name_obs)) |
      !identical(dim(exp)[name_exp], dim(obs)[name_obs])) {
    stop(paste0("Parameter 'exp' and 'obs' must have same length of ",
                "all dimensions except 'memb_dim' and 'dat_dim'."))
  }
  if (!is.null(ref)) {
    name_ref <- sort(names(dim(ref)))
    name_ref <- name_ref[-which(name_ref == memb_dim)]
    if (!is.null(dat_dim)) {
      if (dat_dim %in% name_ref) {
        if (!identical(dim(exp)[dat_dim], dim(ref)[dat_dim])) {
          stop(paste0("If parameter 'ref' has dataset dimension, it must be", 
                      " equal to dataset dimension of 'exp'."))
        }
        name_ref <- name_ref[-which(name_ref == dat_dim)]
      }
    }
    if (!identical(length(name_exp), length(name_ref)) |
        !identical(dim(exp)[name_exp], dim(ref)[name_ref])) {
      stop(paste0("Parameter 'exp' and 'ref' must have the same length of ",
                  "all dimensions except 'memb_dim' and 'dat_dim'",
                  " if there is only one reference dataset."))
    }
  }
  ## prob_thresholds
  if (!is.numeric(prob_thresholds) | !is.vector(prob_thresholds) |
      any(prob_thresholds <= 0) | any(prob_thresholds >= 1)) {
    stop("Parameter 'prob_thresholds' must be a numeric vector between 0 and 1.")
  }
  ## indices_for_clim
  if (is.null(indices_for_clim)) {
    indices_for_clim <- 1:dim(obs)[time_dim]
  } else {
    if (!is.numeric(indices_for_clim) | !is.vector(indices_for_clim)) {
      stop("Parameter 'indices_for_clim' must be NULL or a numeric vector.")
    } else if (length(indices_for_clim) > dim(obs)[time_dim] |
               max(indices_for_clim) > dim(obs)[time_dim] |
               any(indices_for_clim) < 1) {
      stop("Parameter 'indices_for_clim' should be the indices of 'time_dim'.")
    }
  }
  ## cross.val
  if (!is.logical(cross.val)  | length(cross.val) > 1) {
    stop("Parameter 'cross.val' must be either TRUE or FALSE.")
  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 |
        length(ncores) > 1) {
      stop("Parameter 'ncores' must be either NULL or a positive integer.")
    }
  }
  
  ###############################
  
  # Compute ROCSS
  ## output_dims
  if (is.null(dat_dim)) {
    output_dims <- 'cat'
  } else {
    output_dims <- c('nexp', 'nobs', 'cat')
  }
  ## target_dims
  if (!memb_dim %in% names(dim(obs))) {
    target_dims_obs <- c(time_dim, dat_dim)
  } else {
    target_dims_obs <- c(time_dim, memb_dim, dat_dim)
  }
  ## If ref doesn't have & dat_dim is not NULL
  if (!is.null(ref) && !is.null(dat_dim) &&!dat_dim %in% names(dim(ref))) {
    target_dims_ref <- c(time_dim, memb_dim)
  } else {
    target_dims_ref <- c(time_dim, memb_dim, dat_dim)
  }

  if (!is.null(ref)) { ## reference forecast is provided
    res <- Apply(data = list(exp = exp, obs = obs, ref = ref),
                             target_dims = list(exp = c(time_dim, memb_dim, dat_dim),
                                                obs = target_dims_obs,
                                                ref = target_dims_ref),
                             output_dims = output_dims,
                             fun = .ROCSS,
                             prob_thresholds = prob_thresholds, 
                             indices_for_clim = indices_for_clim,
                             time_dim = time_dim, dat_dim = dat_dim,
                             cross.val = cross.val,
                             ncores = ncores)$output1
    
  } else { ## Random forecast as reference forecast
    res <- Apply(data = list(exp = exp, obs = obs),
                             target_dims = list(exp = c(time_dim, memb_dim, dat_dim),
                                                obs = target_dims_obs),
                             output_dims = output_dims,
                             fun = .ROCSS,
                             ref = ref,
                             prob_thresholds = prob_thresholds, 
                             indices_for_clim = indices_for_clim,
                             time_dim = time_dim, dat_dim = dat_dim,
                             cross.val = cross.val,
                             ncores = ncores)$output1
  }
  
  return(res)
}

.ROCSS <- function(exp, obs, ref = NULL, time_dim = 'sdate', dat_dim = NULL, prob_thresholds = c(1/3, 2/3), 
                   indices_for_clim = NULL, cross.val = FALSE) {
  
  # exp: [sdate, memb, (dat)]
  # obs: [sdate, (dat)]
  # ref: [sdate, memb, (dat)] or NULL
  
  if (is.null(dat_dim)) {
    nexp <- 1
    nobs <- 1
    dat_dim <- 'dat'
    dim(exp) <- c(dim(exp), dat = 1)
    dim(obs) <- c(dim(obs), dat = 1)
    if (!is.null(ref)) {
      dim(ref) <- c(dim(ref), dat = 1)
    }
    remove_dat_dim <- TRUE
  } else {
    nexp <- as.numeric(dim(exp)[dat_dim])
    nobs <- as.numeric(dim(obs)[dat_dim])
    if (!is.null(ref) && !dat_dim %in% names(dim(ref))) { # make ref have the same dat dim as exp
      ref <- array(ref, dim = dim(exp))
    }
    remove_dat_dim <- FALSE
  }

  ncats <- 1 + length(prob_thresholds)
  rocs_exp <- array(dim = c(nexp = nexp, nobs = nobs, cat = ncats))
  if (!is.null(ref)) rocs_ref <- array(dim = dim(rocs_exp))

  for (exp_i in 1:nexp) {
    for (obs_i in 1:nobs) {
 
      # Input dim for .get_probs
      ## if exp: [sdate, memb]
      ## if obs: [sdate, (memb)]
      exp_probs <- .get_probs(ClimProjDiags::Subset(exp, dat_dim, exp_i, drop = 'selected'), 
                              indices_for_quantiles = indices_for_clim, 
                              prob_thresholds = prob_thresholds, cross.val = cross.val)
      obs_probs <- .get_probs(data = ClimProjDiags::Subset(obs, dat_dim, obs_i, drop = 'selected'), 
                              indices_for_quantiles = indices_for_clim, 
                              prob_thresholds = prob_thresholds, cross.val = cross.val)
      ## exp_probs and obs_probs: [bin, sdate]
    
      ## ROCS (exp)
      rocs_exp[exp_i, obs_i, ] <- unlist(EnsRoca(ens = Reorder(exp_probs, c(time_dim, 'bin')),
                                                 obs = Reorder(obs_probs, c(time_dim, 'bin')))[1:ncats])
    
      if (!is.null(ref)) {
        ref_probs <- .get_probs(ClimProjDiags::Subset(ref, dat_dim, exp_i, drop = 'selected'),
                                indices_for_quantiles = indices_for_clim, 
                                prob_thresholds = prob_thresholds, cross.val = cross.val)
        rocs_ref[exp_i, obs_i, ] <- unlist(EnsRoca(ens = Reorder(ref_probs, c(time_dim, 'bin')),
                                                   obs = Reorder(obs_probs, c(time_dim, 'bin')))[1:ncats])
      }
    }
  }

  ## ROCSS
  if (is.null(ref)) { ## Random forecast as reference forecast
    rocss <- 2 * rocs_exp - 1
  } else { ## Reference forecast is provided
    rocss <- (rocs_exp - rocs_ref) / (1 - rocs_ref)
  }
  if (remove_dat_dim) {
    rocss <- array(rocss, dim = dim(rocss)['cat'])
  }

  return(rocss)
}
