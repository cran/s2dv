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
#'model version, and another model. It is computed as RPSS = 1 - RPS_exp / RPS_ref.
#'The statistical significance is obtained based on a Random Walk test at the 
#'95% confidence level (DelSole and Tippett, 2016).
#'
#'@param exp A named numerical array of the forecast with at least time 
#'  dimension.
#'@param obs A named numerical array of the observation with at least time 
#'  dimension. The dimensions must be the same as 'exp' except 'memb_dim'. 
#'@param ref A named numerical array of the reference forecast data with at 
#'  least time dimension. The dimensions must be the same as 'exp' except 
#'  'memb_dim'. If it is NULL, the climatological forecast is used as reference
#'  forecast. The default value is NULL.
#'@param time_dim A character string indicating the name of the time dimension.
#'  The default value is 'sdate'.
#'@param memb_dim A character string indicating the name of the member dimension
#'  to compute the probabilities of the forecast and the reference forecast. The
#'  default value is 'member'.
#'@param prob_thresholds A numeric vector of the relative thresholds (from 0 to
#'  1) between the categories. The default value is c(1/3, 2/3), which 
#'  corresponds to tercile equiprobable categories.
#'@param indices_for_clim A vector of the indices to be taken along 'time_dim' 
#'  for computing the thresholds between the probabilistic categories. If NULL,
#'  the whole period is used. The default value is NULL.
#'@param Fair A logical indicating whether to compute the FairRPSS (the 
#'  potential RPSS that the forecast would have with an infinite ensemble size).
#'  The default value is FALSE.
#'@param weights A named two-dimensional numerical array of the weights for each
#'  member and time. The dimension names should include 'memb_dim' and 
#'  'time_dim'. The default value is NULL. The ensemble should have at least 70 
#'  members or span at least 10 time steps and have more than 45 members if 
#'  consistency between the weighted and unweighted methodologies is desired.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return
#'\item{$rpss}{
#'  A numerical array of the RPSS with the same dimensions as "exp" except the
#'  'time_dim' and 'memb_dim' dimensions.
#'}
#'\item{$sign}{
#'  A logical array of the statistical significance of the RPSS with the same
#'  dimensions as 'exp' except the 'time_dim' and 'memb_dim' dimensions.
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
#'res <- RPSS(exp = exp, obs = obs) ## climatology as reference forecast
#'res <- RPSS(exp = exp, obs = obs, ref = ref) ## ref as reference forecast
#'
#'@import multiApply
#'@export
RPSS <- function(exp, obs, ref = NULL, time_dim = 'sdate', memb_dim = 'member',
                 prob_thresholds = c(1/3, 2/3), indices_for_clim = NULL, Fair = FALSE, 
                 weights = NULL, ncores = NULL) {
 
  # Check inputs
  ## exp, obs, and ref (1)
  if (!is.array(exp) | !is.numeric(exp))
    stop('Parameter "exp" must be a numeric array.')
  if (!is.array(obs) | !is.numeric(obs))
    stop('Parameter "obs" must be a numeric array.')
  if (!is.null(ref)) {
    if (!is.array(ref) | !is.numeric(ref))
      stop('Parameter "ref" must be a numeric array.')
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) != 1)
    stop('Parameter "time_dim" must be a character string.')
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
  ## exp and obs (2)
  name_exp <- sort(names(dim(exp)))
  name_obs <- sort(names(dim(obs)))
  name_exp <- name_exp[-which(name_exp == memb_dim)]
  if (memb_dim %in% name_obs) {
    name_obs <- name_obs[-which(name_obs == memb_dim)]
  }
  if (!identical(length(name_exp), length(name_obs)) |
      !identical(dim(exp)[name_exp], dim(obs)[name_obs])) {
    stop(paste0("Parameter 'exp' and 'obs' must have same length of ",
                "all dimensions expect 'memb_dim'."))
  }
  if (!is.null(ref)) {
    name_ref <- sort(names(dim(ref)))
    name_ref <- name_ref[-which(name_ref == memb_dim)]
    if (!identical(length(name_exp), length(name_ref)) |
        !identical(dim(exp)[name_exp], dim(ref)[name_ref])) {
      stop(paste0("Parameter 'exp' and 'obs' must have same length of ",
                  "all dimensions expect 'memb_dim'."))
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
  if (!is.logical(Fair)  | length(Fair) > 1) {
    stop("Parameter 'Fair' must be either TRUE or FALSE.")
  }
  ## weights
  if (!is.null(weights)) {
    if (!is.array(weights) | !is.numeric(weights))
      stop('Parameter "weights" must be a two-dimensional numeric array.')
    if (length(dim(weights)) != 2 | any(!names(dim(weights)) %in% c(memb_dim, time_dim)))
      stop("Parameter 'weights' must have two dimensions with the names of memb_dim and time_dim.")
    if (dim(weights)[memb_dim] != dim(exp)[memb_dim] |
        dim(weights)[time_dim] != dim(exp)[time_dim]) {
      stop("Parameter 'weights' must have the same dimension lengths as memb_dim and time_dim in 'exp'.")
    }
    weights <- Reorder(weights, c(time_dim, memb_dim))
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
    target_dims_obs <- time_dim
  } else {
    target_dims_obs <- c(time_dim, memb_dim)
  }

  if (!is.null(ref)) { # use "ref" as reference forecast
    data <- list(exp = exp, obs = obs, ref = ref)
    target_dims = list(exp = c(time_dim, memb_dim),
                       obs = target_dims_obs,
                       ref = c(time_dim, memb_dim))
  } else {
    data <- list(exp = exp, obs = obs)
    target_dims = list(exp = c(time_dim, memb_dim),
                       obs = target_dims_obs)
  }
  output <- Apply(data,
                  target_dims = target_dims,
                  fun = .RPSS,
                  prob_thresholds = prob_thresholds,
                  indices_for_clim = indices_for_clim, Fair = Fair,
                  weights = weights,
                  ncores = ncores)
  
  return(output)
}

.RPSS <- function(exp, obs, ref = NULL, prob_thresholds = c(1/3, 2/3), 
                  indices_for_clim = NULL, Fair = FALSE, weights = NULL) {
  # exp: [sdate, memb]
  # obs: [sdate, (memb)]
  # ref: [sdate, memb] or NULL

  # RPS of the forecast
  rps_exp <- .RPS(exp = exp, obs = obs, prob_thresholds = prob_thresholds, 
                  indices_for_clim = indices_for_clim, Fair = Fair, weights = weights)
  
  # RPS of the reference forecast
  if (is.null(ref)) { ## using climatology as reference forecast
    obs_probs <- .get_probs(data = obs, indices_for_quantiles = indices_for_clim, 
                            prob_thresholds = prob_thresholds, weights = NULL)
    # obs_probs: [bin, sdate] 

    clim_probs <- c(prob_thresholds[1], diff(prob_thresholds), 1 - prob_thresholds[length(prob_thresholds)])
    clim_probs <- array(clim_probs, dim = dim(obs_probs))
    # clim_probs: [bin, sdate]  

    # Calculate RPS for each time step
    probs_clim_cumsum <- apply(clim_probs, 2, cumsum)
    probs_obs_cumsum <- apply(obs_probs, 2, cumsum)
    rps_ref <- apply((probs_clim_cumsum - probs_obs_cumsum)^2, 2, sum)
    # rps_ref: [sdate]

#    if (Fair) { # FairRPS
#      ## adjustment <- rowSums(-1 * (1/R - 1/R.new) * ens.cum * (R - ens.cum)/R/(R - 1)) [formula taken from SpecsVerification::EnsRps]
#      R <- dim(exp)[2]  #memb
#      R_new <- Inf
#      adjustment <- (-1) / (R - 1) * probs_clim_cumsum * (1 - probs_clim_cumsum)
#      adjustment <- apply(adjustment, 2, sum)
#      rps_ref <- rps_ref + adjustment
#    }

  } else { # use "ref" as reference forecast
    rps_ref <- .RPS(exp = ref, obs = obs, prob_thresholds = prob_thresholds, 
                    indices_for_clim = indices_for_clim, Fair = Fair, weights = weights)
  }
  
  # RPSS
  rpss <- 1 - mean(rps_exp) / mean(rps_ref)
  
  # Significance
  sign <- .RandomWalkTest(skill_A = rps_exp, skill_B = rps_ref)$signif
  
  return(list(rpss = rpss, sign = sign))
}

