#'Compute the Continuous Ranked Probability Skill Score
#'
#'The Continuous Ranked Probability Skill Score (CRPSS; Wilks, 2011) is the 
#'skill score based on the Continuous Ranked Probability Score (CRPS; Wilks, 
#'2011). It can be used to assess whether a forecast presents an improvement or 
#'worsening with respect to a reference forecast. The CRPSS ranges between minus 
#'infinite and 1. If the CRPSS is positive, it indicates that the forecast has 
#'higher skill than the reference forecast, while a negative value means that it
#'has a lower skill. Examples of reference forecasts are the climatological 
#'forecast (same probabilities for all categories for all time steps), 
#'persistence, a previous model version, or another model. It is computed as 
#'CRPSS = 1 - CRPS_exp / CRPS_ref. The statistical significance is obtained 
#'based on a Random Walk test at the 95% confidence level (DelSole and Tippett,
#'2016).
#'
#'@param exp A named numerical array of the forecast with at least time 
#'  dimension.
#'@param obs A named numerical array of the observation with at least time 
#'  dimension. The dimensions must be the same as 'exp' except 'memb_dim'
#'  and 'dat_dim'. 
#'@param ref A named numerical array of the reference forecast data with at 
#'  least time and member dimension. The dimensions must be the same as 'exp' 
#'  except 'memb_dim' and 'dat_dim'. If there is only one reference dataset,
#'  it should not have dataset dimension. If there is corresponding reference 
#'  for each experiement, the dataset dimension must have the same length as in
#'  'exp'. If 'ref' is NULL, the climatological forecast is used as reference 
#'  forecast. The default value is NULL.
#'@param time_dim A character string indicating the name of the time dimension.
#'  The default value is 'sdate'.
#'@param memb_dim A character string indicating the name of the member dimension
#'  to compute the probabilities of the forecast and the reference forecast. The
#'  default value is 'member'.
#'@param dat_dim A character string indicating the name of dataset dimension. 
#'  The length of this dimension can be different between 'exp' and 'obs'. 
#'  The default value is NULL.
#'@param Fair A logical indicating whether to compute the FairCRPSS (the 
#'  potential CRPSS that the forecast would have with an infinite ensemble 
#'  size). The default value is FALSE.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return
#'\item{$crpss}{
#'  A numerical array of the CRPSS with dimensions c(nexp, nobs, the rest 
#'  dimensions of 'exp' except 'time_dim' and 'memb_dim' dimensions). nexp is 
#'  the number of experiment (i.e., dat_dim in exp), and nobs is the number of 
#'  observation (i.e., dat_dim in obs). If 'dat_dim' is NULL, nexp and nobs are
#'  omitted.
#'}
#'\item{$sign}{
#'  A logical array of the statistical significance of the CRPSS with the same 
#'  dimensions as $crpss.
#'}
#'
#'@references 
#'Wilks, 2011; https://doi.org/10.1016/B978-0-12-385022-5.00008-7
#'DelSole and Tippett, 2016; https://doi.org/10.1175/MWR-D-15-0218.1
#'
#'@examples
#'exp <- array(rnorm(1000), dim = c(lat = 3, lon = 2, member = 10, sdate = 50))
#'obs <- array(rnorm(1000), dim = c(lat = 3, lon = 2, sdate = 50))
#'ref <- array(rnorm(1000), dim = c(lat = 3, lon = 2, member = 10, sdate = 50))
#'res <- CRPSS(exp = exp, obs = obs) ## climatology as reference forecast
#'res <- CRPSS(exp = exp, obs = obs, ref = ref) ## ref as reference forecast
#'
#'@import multiApply
#'@importFrom ClimProjDiags Subset
#'@export
CRPSS <- function(exp, obs, ref = NULL, time_dim = 'sdate', memb_dim = 'member', dat_dim = NULL,
                  Fair = FALSE, ncores = NULL) {
  
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
  if (memb_dim %in% names(dim(obs))) {
    if (identical(as.numeric(dim(obs)[memb_dim]), 1)) {
      obs <- ClimProjDiags::Subset(x = obs, along = memb_dim, indices = 1, drop = 'selected')
    } else {
      stop("Not implemented for observations with members ('obs' can have ",
           "'memb_dim', but it should be of length = 1).")
    }
  }
  name_exp <- sort(names(dim(exp)))
  name_obs <- sort(names(dim(obs)))
  name_exp <- name_exp[-which(name_exp == memb_dim)]
  if (!is.null(dat_dim)) {
    name_exp <- name_exp[-which(name_exp == dat_dim)]
    name_obs <- name_obs[-which(name_obs == dat_dim)]
  }
  if (!identical(length(name_exp), length(name_obs)) |
      !identical(dim(exp)[name_exp], dim(obs)[name_obs])) {
    stop(paste0("Parameter 'exp' and 'obs' must have same length of all dimensions", 
                " except 'memb_dim' and 'dat_dim'."))
  }
  if (!is.null(ref)) {
    name_ref <- sort(names(dim(ref)))
    name_ref <- name_ref[-which(name_ref == memb_dim)]
    if (!is.null(dat_dim)) {
      if (dat_dim %in% name_ref) {
        if (!identical(dim(exp)[dat_dim], dim(ref)[dat_dim])) {
          stop(paste0("If parameter 'ref' has dataset dimension it must be", 
                      " equal to dataset dimension of 'exp'."))
        }
        name_ref <- name_ref[-which(name_ref == dat_dim)]
      }
    }
    if (!identical(length(name_exp), length(name_ref)) |
        !identical(dim(exp)[name_exp], dim(ref)[name_ref])) {
      stop(paste0("Parameter 'exp' and 'ref' must have same length of ",
                  "all dimensions except 'memb_dim' and 'dat_dim' if there is ",
                  "only one reference dataset."))
    }
  }
  ## Fair
  if (!is.logical(Fair)  | length(Fair) > 1) {
    stop("Parameter 'Fair' must be either TRUE or FALSE.")
  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 | length(ncores) > 1) {
      stop("Parameter 'ncores' must be either NULL or a positive integer.")
    }
  }
  
  ###############################
  
  # Compute CRPSS
  if (!is.null(ref)) { # use "ref" as reference forecast
    if (!is.null(dat_dim) && (dat_dim %in% names(dim(ref)))) {
      target_dims_ref <- c(time_dim, memb_dim, dat_dim)
    } else {
      target_dims_ref <- c(time_dim, memb_dim)
    }
    data <- list(exp = exp, obs = obs, ref = ref)
    target_dims = list(exp = c(time_dim, memb_dim, dat_dim),
                       obs = c(time_dim, dat_dim),
                       ref = target_dims_ref)
  } else {
    data <- list(exp = exp, obs = obs)
    target_dims = list(exp = c(time_dim, memb_dim, dat_dim),
                       obs = c(time_dim, dat_dim))
  }
  output <- Apply(data,
                  target_dims = target_dims,
                  fun = .CRPSS,
                  time_dim = time_dim, memb_dim = memb_dim, 
                  dat_dim = dat_dim, 
                  Fair = Fair,
                  ncores = ncores)
  
  return(output)
}

.CRPSS <- function(exp, obs, ref = NULL, time_dim = 'sdate', memb_dim = 'member', dat_dim = NULL,
                   Fair = FALSE) {
  
  # exp: [sdate, memb, (dat)]
  # obs: [sdate, (dat)]
  # ref: [sdate, memb, (dat)] or NULL

  if (is.null(dat_dim)) {
    nexp <- 1
    nobs <- 1
  } else {
    nexp <- as.numeric(dim(exp)[dat_dim])
    nobs <- as.numeric(dim(obs)[dat_dim])
  }
  
  #----- CRPS of the forecast
  # [sdate, (nexp), (nobs)]
  crps_exp <- .CRPS(exp = exp, obs = obs, time_dim = time_dim, memb_dim = memb_dim, 
                      dat_dim = dat_dim, Fair = Fair)
  
  #----- CRPS of the reference forecast
  if (is.null(ref)) {
    ## using climatology as reference forecast
    ## all the time steps are used as if they were members
    ## then, ref dimensions are [sdate, memb]
    ## memb dimension has length(sdate) - 1 due to cross-validation

    obs_time_len <- dim(obs)[time_dim]
    if (is.null(dat_dim)) {
      
      ## Without cross-validation: 
      ## ref <- array(data = rep(obs, each = obs_time_len), dim = c(obs_time_len, obs_time_len))
      ## With cross-validation (excluding the value of that year to create ref for that year):
      ref <- array(data = NA, dim = c(obs_time_len, obs_time_len - 1))
      for (i in 1:obs_time_len) {
        ref[i, ] <- obs[-i]
      }
      
      names(dim(ref)) <- c(time_dim, memb_dim)
      # ref: [sdate, memb]; obs: [sdate]
      crps_ref <- .CRPS(exp = ref, obs = obs, time_dim = time_dim, memb_dim = memb_dim,
                        dat_dim = dat_dim, Fair = Fair)
      # crps_ref should be [sdate]

    } else {
      crps_ref <- array(dim = c(obs_time_len, nobs))
      names(dim(crps_ref)) <- c(time_dim, 'nobs')
      for (i_obs in 1:nobs) {
        
        ## Without cross-validation: 
        ## ref <- array(data = rep(obs[, i_obs], each = obs_time_len), dim = c(obs_time_len, obs_time_len))
        ## With cross-validation (excluding the value of that year to create ref for that year):
        ref <- array(data = NA, dim = c(obs_time_len, obs_time_len - 1))
        for (i in 1:obs_time_len) {
          ref[i, ] <- obs[-i, i_obs]
        }
        
        names(dim(ref)) <- c(time_dim, memb_dim)
        crps_ref[, i_obs] <- .CRPS(exp = ref, obs = ClimProjDiags::Subset(obs, dat_dim, i_obs, drop = 'selected'), 
                                   time_dim = time_dim, memb_dim = memb_dim, dat_dim = NULL, Fair = Fair)
      }
      # crps_ref should be [sdate, nobs]
    }

  } else { # ref is not NULL
    if (!is.null(dat_dim) && (!dat_dim %in% names(dim(ref)))) {
      remove_dat_dim <- TRUE
      ref <- InsertDim(data = ref, posdim = length(dim(ref)) + 1 , lendim = 1, name = dat_dim)
    } else {
      remove_dat_dim <- FALSE
    }
    crps_ref <- .CRPS(exp = ref, obs = obs, time_dim = time_dim, memb_dim = memb_dim,
                      dat_dim = dat_dim, Fair = Fair)
    # crps_ref should be [sdate, (nexp), (nobs)]

    if (!is.null(dat_dim)) {
      if (isTRUE(remove_dat_dim)) {
        dim(crps_ref) <- dim(crps_ref)[-2]
      }
    }
  }

  #----- CRPSS
  if (!is.null(dat_dim)) {
    # If ref != NULL & ref has dat_dim, crps_ref = [sdate, nexp, nobs]; else, crps_ref = [sdate, nobs]

    crps_exp_mean <- MeanDims(crps_exp, time_dim, na.rm = FALSE)
    crps_ref_mean <- MeanDims(crps_ref, time_dim, na.rm = FALSE)
    crpss <- array(dim = c(nexp = nexp, nobs = nobs))
    sign <- array(dim = c(nexp = nexp, nobs = nobs))

    if (length(dim(crps_ref_mean)) == 1) {
      for (i in 1:nexp) {
        for (j in 1:nobs) {
          crpss[i, j] <- 1 - crps_exp_mean[i, j] / crps_ref_mean[j]
          sign[i, j] <- .RandomWalkTest(skill_A = crps_exp_mean[i, j], skill_B = crps_ref_mean[j], sign = T, pval = F)$sign
        }
      }
    } else {
      for (i in 1:nexp) {
        for (j in 1:nobs) {
          crpss[i, j] <- 1 - crps_exp_mean[i, j] / crps_ref_mean[i, j]
          sign[i, j] <- .RandomWalkTest(skill_A = crps_exp_mean[i, j], skill_B = crps_ref_mean[i, j], sign = T, pval = F)$sign
        }
      }
    }

  } else {
    crpss <- 1 - mean(crps_exp) / mean(crps_ref)
    # Significance
    sign <- .RandomWalkTest(skill_A = crps_exp, skill_B = crps_ref, sign = T, pval = F)$sign
  }
  
  return(list(crpss = crpss, sign = sign))
}
