#'Compute the Ranked Probability Skill Score
#'
#'The Ranked Probability Skill Score (RPSS; Wilks, 2011) is the skill score 
#'based on the Ranked Probability Score (RPS; Wilks, 2011). It can be used to 
#'assess whether a forecast presents an improvement or worsening with respect to
#'a reference forecast. The RPSS ranges between minus infinite and 1. If the 
#'RPSS is positive, it indicates that the forecast has higher skill than the 
#'reference forecast, while a negative value means that it has a lower skill. 
#'Examples of reference forecasts are the climatological forecast (same 
#'probabilities for all categories for all time steps), persistence, a previous
#'model version, and another model. It is computed as 
#'\code{RPSS = 1 - RPS_exp / RPS_ref}. The statistical significance is obtained 
#'based on a Random Walk test at the 95% confidence level (DelSole and Tippett, 
#'2016). If there is more than one dataset, RPS will be computed for each pair 
#'of exp and obs data.
#'
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
#'@param prob_thresholds A numeric vector of the relative thresholds (from 0 to
#'  1) between the categories. The default value is c(1/3, 2/3), which 
#'  corresponds to tercile equiprobable categories.
#'@param indices_for_clim A vector of the indices to be taken along 'time_dim' 
#'  for computing the thresholds between the probabilistic categories. If NULL,
#'  the whole period is used. The default value is NULL.
#'@param Fair A logical indicating whether to compute the FairRPSS (the 
#'  potential RPSS that the forecast would have with an infinite ensemble size).
#'  The default value is FALSE.
#'@param weights Deprecated and will be removed in the next release. Please use
#'  'weights_exp' and 'weights_ref' instead.
#'@param weights_exp A named numerical array of the forecast ensemble weights.
#'  The dimension should include 'memb_dim', 'time_dim' and 'dat_dim' if there
#'  are multiple datasets. All dimension lengths must be equal to 'exp' 
#'  dimension lengths. The default value is NULL, which means no weighting is 
#'  applied. The ensemble should have at least 70 members or span at least 10 
#'  time steps and have more than 45 members if consistency between the weighted
#'   and unweighted methodologies is desired.
#'@param weights_ref Same as 'weights_exp' but for the reference forecast.
#'@param cross.val A logical indicating whether to compute the thresholds between 
#'  probabilistics categories in cross-validation.
#'  The default value is FALSE.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return
#'\item{$rpss}{
#'  A numerical array of RPSS with dimensions c(nexp, nobs, the rest dimensions 
#'  of 'exp' except 'time_dim' and 'memb_dim' dimensions). nexp is the number of 
#'  experiment (i.e., dat_dim in exp), and nobs is the number of observation 
#'  i.e., dat_dim in obs). If dat_dim is NULL, nexp and nobs are omitted.
#'}
#'\item{$sign}{
#'  A logical array of the statistical significance of the RPSS with the same 
#'  dimensions as $rpss.
#'}
#'
#'@references 
#'Wilks, 2011; https://doi.org/10.1016/B978-0-12-385022-5.00008-7
#'DelSole and Tippett, 2016; https://doi.org/10.1175/MWR-D-15-0218.1
#'
#'@examples
#'set.seed(1)
#'exp <- array(rnorm(3000), dim = c(lat = 3, lon = 2, member = 10, sdate = 50))
#'set.seed(2)
#'obs <- array(rnorm(300), dim = c(lat = 3, lon = 2, sdate = 50))
#'set.seed(3)
#'ref <- array(rnorm(3000), dim = c(lat = 3, lon = 2, member = 10, sdate = 50))
#'weights <- sapply(1:dim(exp)['sdate'], function(i) {
#'             n <- abs(rnorm(10))
#'             n/sum(n)
#'           })
#'dim(weights) <- c(member = 10, sdate = 50)
#'res <- RPSS(exp = exp, obs = obs) ## climatology as reference forecast
#'res <- RPSS(exp = exp, obs = obs, ref = ref) ## ref as reference forecast
#'res <- RPSS(exp = exp, obs = obs, ref = ref, weights_exp = weights, weights_ref = weights)
#'@import multiApply
#'@export
RPSS <- function(exp, obs, ref = NULL, time_dim = 'sdate', memb_dim = 'member',
                 dat_dim = NULL, prob_thresholds = c(1/3, 2/3), indices_for_clim = NULL,
                 Fair = FALSE, weights = NULL, weights_exp = NULL, weights_ref = NULL, 
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
                  "all dimensions except 'memb_dim' and 'dat_dim' if there is ",
                  "only one reference dataset."))
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
  ## Fair
  if (!is.logical(Fair) | length(Fair) > 1) {
    stop("Parameter 'Fair' must be either TRUE or FALSE.")
  }
  ## cross.val
  if (!is.logical(cross.val)  | length(cross.val) > 1) {
    stop("Parameter 'cross.val' must be either TRUE or FALSE.")
  }
  ## weights
  if (!is.null(weights)) {
    .warning(paste0("Parameter 'weights' is deprecated and will be removed in the next release. ",
                    "Use 'weights_exp' and 'weights_ref' instead. The value will be assigned ",
                    "to these two parameters now if they are NULL."), tag = '! Deprecation: ')
    if (is.null(weights_exp)) weights_exp <- weights
    if (is.null(weights_ref)) weights_ref <- weights
  }
  ## weights_exp
  if (!is.null(weights_exp)) {
    if (!is.array(weights_exp) | !is.numeric(weights_exp))
      stop("Parameter 'weights_exp' must be a named numeric array.")

    if (is.null(dat_dim)) {
      if (length(dim(weights_exp)) != 2 | any(!names(dim(weights_exp)) %in% c(memb_dim, time_dim)))
        stop("Parameter 'weights_exp' must have two dimensions with the names of 'memb_dim' and 'time_dim'.")
      if (dim(weights_exp)[memb_dim] != dim(exp)[memb_dim] |
          dim(weights_exp)[time_dim] != dim(exp)[time_dim]) {
        stop("Parameter 'weights_exp' must have the same dimension lengths as 'memb_dim' and 'time_dim' in 'exp'.")
      }
      weights_exp <- Reorder(weights_exp, c(time_dim, memb_dim))

    } else {
      if (length(dim(weights_exp)) != 3 | any(!names(dim(weights_exp)) %in% c(memb_dim, time_dim, dat_dim)))
        stop("Parameter 'weights_exp' must have three dimensions with the names of 'memb_dim', 'time_dim' and 'dat_dim'.")
      if (dim(weights_exp)[memb_dim] != dim(exp)[memb_dim] |
          dim(weights_exp)[time_dim] != dim(exp)[time_dim] |
          dim(weights_exp)[dat_dim] != dim(exp)[dat_dim]) {
        stop(paste0("Parameter 'weights_exp' must have the same dimension lengths ", 
                    "as 'memb_dim', 'time_dim' and 'dat_dim' in 'exp'."))
      }
      weights_exp <- Reorder(weights_exp, c(time_dim, memb_dim, dat_dim))
    }
    
  }
  ## weights_ref
  if (!is.null(weights_ref)) {
    if (!is.array(weights_ref) | !is.numeric(weights_ref))
      stop('Parameter "weights_ref" must be a named numeric array.')

    if (is.null(dat_dim) | ((!is.null(dat_dim)) && (!dat_dim %in% names(dim(ref))))) {
      if (length(dim(weights_ref)) != 2 | any(!names(dim(weights_ref)) %in% c(memb_dim, time_dim)))
        stop("Parameter 'weights_ref' must have two dimensions with the names of 'memb_dim' and 'time_dim'.")
      if (dim(weights_ref)[memb_dim] != dim(exp)[memb_dim] |
          dim(weights_ref)[time_dim] != dim(exp)[time_dim]) {
        stop("Parameter 'weights_ref' must have the same dimension lengths as 'memb_dim' and 'time_dim' in 'ref'.")
      }
      weights_ref <- Reorder(weights_ref, c(time_dim, memb_dim))

    } else {
      if (length(dim(weights_ref)) != 3 | any(!names(dim(weights_ref)) %in% c(memb_dim, time_dim, dat_dim)))
        stop("Parameter 'weights_ref' must have three dimensions with the names of 'memb_dim', 'time_dim' and 'dat_dim'.")
      if (dim(weights_ref)[memb_dim] != dim(ref)[memb_dim] |
          dim(weights_ref)[time_dim] != dim(ref)[time_dim] |
          dim(weights_ref)[dat_dim] != dim(ref)[dat_dim]) {
        stop(paste0("Parameter 'weights_ref' must have the same dimension lengths ", 
                    "as 'memb_dim', 'time_dim' and 'dat_dim' in 'ref'."))
      }
      weights_ref <- Reorder(weights_ref, c(time_dim, memb_dim, dat_dim))
    }

  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 |
      length(ncores) > 1) {
      stop("Parameter 'ncores' must be either NULL or a positive integer.")
    }
  }

  ###############################

  # Compute RPSS
  if (!memb_dim %in% names(dim(obs))) {
    target_dims_obs <- c(time_dim, dat_dim)
  } else {
    target_dims_obs <- c(time_dim, memb_dim, dat_dim)
  }

  if (!is.null(ref)) { # use "ref" as reference forecast
    if (!is.null(dat_dim) && (dat_dim %in% names(dim(ref)))) {
      target_dims_ref <- c(time_dim, memb_dim, dat_dim)
    } else {
        target_dims_ref <- c(time_dim, memb_dim)
    }
    data <- list(exp = exp, obs = obs, ref = ref)
    target_dims = list(exp = c(time_dim, memb_dim, dat_dim),
                       obs = target_dims_obs,
                       ref = target_dims_ref)
  } else {
    data <- list(exp = exp, obs = obs)
    target_dims = list(exp = c(time_dim, memb_dim, dat_dim),
                       obs = target_dims_obs)
  }
  output <- Apply(data,
                  target_dims = target_dims,
                  fun = .RPSS,
                  time_dim = time_dim, memb_dim = memb_dim, 
                  dat_dim = dat_dim, 
                  prob_thresholds = prob_thresholds,
                  indices_for_clim = indices_for_clim, Fair = Fair,
                  weights_exp = weights_exp,
                  weights_ref = weights_ref,
                  cross.val = cross.val, 
                  ncores = ncores)
  
  return(output)

}

.RPSS <- function(exp, obs, ref = NULL, time_dim = 'sdate', memb_dim = 'member', dat_dim = NULL,
                  prob_thresholds = c(1/3, 2/3), indices_for_clim = NULL, Fair = FALSE, 
                  weights_exp = NULL, weights_ref = NULL, cross.val = FALSE) {

  # exp: [sdate, memb, (dat)]
  # obs: [sdate, (memb), (dat)]
  # ref: [sdate, memb, (dat)] or NULL

  if (is.null(dat_dim)) {
    nexp <- 1
    nobs <- 1
  } else {
    nexp <- as.numeric(dim(exp)[dat_dim])
    nobs <- as.numeric(dim(obs)[dat_dim])
  }

  # RPS of the forecast
  rps_exp <- .RPS(exp = exp, obs = obs, time_dim = time_dim, memb_dim = memb_dim, dat_dim = dat_dim,
                  prob_thresholds = prob_thresholds, indices_for_clim = indices_for_clim,
                  Fair = Fair, weights = weights_exp, cross.val = cross.val)
  
  # RPS of the reference forecast
  if (is.null(ref)) { ## using climatology as reference forecast
    if (!memb_dim %in% names(dim(obs))) {
      obs <- InsertDim(obs, posdim = 2, lendim = 1, name = memb_dim)
    }
    if (is.null(dat_dim)) {
      dim(obs) <- c(dim(obs), nobs = nobs)
    }
    rps_ref <- array(dim = c(dim(obs)[time_dim], nobs = nobs))

    for (j in 1:nobs) {
      obs_data <- obs[ , , j]
      if (is.null(dim(obs_data))) dim(obs_data) <- c(dim(obs)[1:2])
      # obs_probs: [bin, sdate] 
      obs_probs <- .get_probs(data = obs_data, indices_for_quantiles = indices_for_clim, 
                              prob_thresholds = prob_thresholds, weights = NULL, cross.val = cross.val)
      # clim_probs: [bin, sdate]
      clim_probs <- c(prob_thresholds[1], diff(prob_thresholds), 1 - prob_thresholds[length(prob_thresholds)])
      clim_probs <- array(clim_probs, dim = dim(obs_probs))

      # Calculate RPS for each time step
      probs_clim_cumsum <- apply(clim_probs, 2, cumsum)
      probs_obs_cumsum <- apply(obs_probs, 2, cumsum)
      rps_ref[ , j] <- apply((probs_clim_cumsum - probs_obs_cumsum)^2, 2, sum)

  #   if (Fair) { # FairRPS
  #     ## adjustment <- rowSums(-1 * (1/R - 1/R.new) * ens.cum * (R - ens.cum)/R/(R - 1)) [formula taken from SpecsVerification::EnsRps]
  #     R <- dim(exp)[2]  #memb
  #     R_new <- Inf
  #     adjustment <- (-1) / (R - 1) * probs_clim_cumsum * (1 - probs_clim_cumsum)
  #     adjustment <- apply(adjustment, 2, sum)
  #     rps_ref <- rps_ref + adjustment
  #   }

    }
    if (is.null(dat_dim)) {
      dim(rps_ref) <- dim(exp)[time_dim]
    }

  } else { # use "ref" as reference forecast
    if (!is.null(dat_dim) && (!dat_dim %in% names(dim(ref)))) {
      remove_dat_dim <- TRUE
      ref <- InsertDim(ref, posdim = 3, lendim = 1, name = dat_dim)
        if (!is.null(weights_ref)) {
          weights_ref <- InsertDim(weights_ref, posdim = 3, lendim = 1, name = dat_dim)
        }
    } else {
      remove_dat_dim <- FALSE
    }

    rps_ref <- .RPS(exp = ref, obs = obs, time_dim = time_dim, memb_dim = memb_dim, dat_dim = dat_dim,
                    prob_thresholds = prob_thresholds, indices_for_clim = indices_for_clim,
                    Fair = Fair, weights = weights_ref, cross.val = cross.val)
    if (!is.null(dat_dim)) {
      if (isTRUE(remove_dat_dim)) {
        dim(rps_ref) <- dim(rps_ref)[-2]
      }
    }
  }
  
  if (!is.null(dat_dim)) {

    rps_exp_mean <- MeanDims(rps_exp, time_dim, na.rm = FALSE)
    rps_ref_mean <- MeanDims(rps_ref, time_dim, na.rm = FALSE)
    rpss <- array(dim = c(nexp = nexp, nobs = nobs))
    sign <- array(dim = c(nexp = nexp, nobs = nobs))

    if (length(dim(rps_ref_mean)) == 1) {
      for (i in 1:nexp) {
        for (j in 1:nobs) {
          rpss[i, j] <- 1 - rps_exp_mean[i, j] / rps_ref_mean[j]
          sign[i, j] <- .RandomWalkTest(skill_A = rps_exp[, i, j], skill_B = rps_ref[, j])$sign
        }
      }
    } else {
      for (i in 1:nexp) {
        for (j in 1:nobs) {
          rpss[i, j] <- 1 - rps_exp_mean[i, j] / rps_ref_mean[i, j]
          sign[i, j] <- .RandomWalkTest(skill_A = rps_exp[, i, j], skill_B = rps_ref[, i, j])$sign
        }
      }
    }
  } else {
    rpss <- 1 - mean(rps_exp) / mean(rps_ref)
    # Significance
    sign <- .RandomWalkTest(skill_A = rps_exp, skill_B = rps_ref, sign = T, pval = F)$sign
  }
  
  return(list(rpss = rpss, sign = sign))
}
