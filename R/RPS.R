#'Compute the Ranked Probability Score
#'
#'The Ranked Probability Score (RPS; Wilks, 2011) is defined as the sum of the
#'squared differences between the cumulative forecast probabilities (computed 
#'from the ensemble members) and the observations (defined as 0% if the category
#'did not happen and 100% if it happened). It can be used to evaluate the skill 
#'of multi-categorical probabilistic forecasts. The RPS ranges between 0 
#'(perfect forecast) and n-1 (worst possible forecast), where n is the number of
#'categories. In the case of a forecast divided into two categories (the lowest
#'number of categories that a probabilistic forecast can have), the RPS 
#'corresponds to the Brier Score (BS; Wilks, 2011), therefore, ranges between 0 
#'and 1. If there is more than one dataset, RPS will be computed for each pair 
#'of exp and obs data.
#'
#'@param exp A named numerical array of the forecast with at least time and  
#'  member dimension.  
#'@param obs A named numerical array of the observation with at least time 
#'  dimension. The dimensions must be the same as 'exp' except 'memb_dim' and 
#'  'dat_dim'.
#'@param time_dim A character string indicating the name of the time dimension.
#'  The default value is 'sdate'.
#'@param dat_dim A character string indicating the name of dataset dimension. 
#'  The length of this dimension can be different between 'exp' and 'obs'. 
#'  The default value is NULL.
#'@param memb_dim A character string indicating the name of the member dimension
#'  to compute the probabilities of the forecast. The default value is 'member'.
#'@param prob_thresholds A numeric vector of the relative thresholds (from 0 to
#'  1) between the categories. The default value is c(1/3, 2/3), which 
#'  corresponds to tercile equiprobable categories.
#'@param indices_for_clim A vector of the indices to be taken along 'time_dim' 
#'  for computing the thresholds between the probabilistic categories. If NULL,
#'  the whole period is used. The default value is NULL.
#'@param Fair A logical indicating whether to compute the FairRPS (the 
#'  potential RPS that the forecast would have with an infinite ensemble size).
#'  The default value is FALSE.
#'@param weights A named numerical array of the weights for 'exp'. If 'dat_dim' 
#'  is NULL, the dimension should include 'memb_dim' and 'time_dim'. Else, the 
#'  dimension should also include 'dat_dim'. The default value is NULL. The 
#'  ensemble should have at least 70 members or span at least 10 time steps and
#'  have more than 45 members if consistency between the weighted and unweighted
#'  methodologies is desired.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return 
#'A numerical array of RPS with dimensions c(nexp, nobs, the rest dimensions of 
#''exp' except 'time_dim' and 'memb_dim' dimensions). nexp is the number of 
#'experiment (i.e., dat_dim in exp), and nobs is the number of observation 
#'(i.e., dat_dim in obs). If dat_dim is NULL, nexp and nobs are omitted.
#'
#'@references 
#'Wilks, 2011; https://doi.org/10.1016/B978-0-12-385022-5.00008-7
#'
#'@examples
#'exp <- array(rnorm(1000), dim = c(lat = 3, lon = 2, member = 10, sdate = 50))
#'obs <- array(rnorm(1000), dim = c(lat = 3, lon = 2, sdate = 50))
#'res <- RPS(exp = exp, obs = obs)
#'
#'@import multiApply
#'@importFrom easyVerification convert2prob
#'@export
RPS <- function(exp, obs, time_dim = 'sdate', memb_dim = 'member', dat_dim = NULL,
                prob_thresholds = c(1/3, 2/3), indices_for_clim = NULL, Fair = FALSE,
                weights = NULL, ncores = NULL) {
  
  # Check inputs
  ## exp and obs (1)
  if (!is.array(exp) | !is.numeric(exp))
    stop('Parameter "exp" must be a numeric array.')
  if (!is.array(obs) | !is.numeric(obs))
    stop('Parameter "obs" must be a numeric array.')
  if(any(is.null(names(dim(exp))))| any(nchar(names(dim(exp))) == 0) |
     any(is.null(names(dim(obs))))| any(nchar(names(dim(obs))) == 0)) {
    stop("Parameter 'exp' and 'obs' must have dimension names.")
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) != 1)
    stop('Parameter "time_dim" must be a character string.')
  if (!time_dim %in% names(dim(exp)) | !time_dim %in% names(dim(obs))) {
    stop("Parameter 'time_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## memb_dim
  if (!is.character(memb_dim) | length(memb_dim) > 1) {
    stop("Parameter 'memb_dim' must be a character string.")
  }
  if (!memb_dim %in% names(dim(exp))) {
    stop("Parameter 'memb_dim' is not found in 'exp' dimension.")
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
  ## exp and obs (2)
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
      stop("Parameter 'weights' must be a named numeric array.")
    if (is.null(dat_dim)) {
      if (length(dim(weights)) != 2 | any(!names(dim(weights)) %in% c(memb_dim, time_dim)))
        stop("Parameter 'weights' must have two dimensions with the names of 'memb_dim' and 'time_dim'.")
      if (dim(weights)[memb_dim] != dim(exp)[memb_dim] |
          dim(weights)[time_dim] != dim(exp)[time_dim]) {
        stop(paste0("Parameter 'weights' must have the same dimension lengths ", 
                    "as 'memb_dim' and 'time_dim' in 'exp'."))
      }
      weights <- Reorder(weights, c(time_dim, memb_dim))

    } else {
      if (length(dim(weights)) != 3 | any(!names(dim(weights)) %in% c(memb_dim, time_dim, dat_dim)))
        stop("Parameter 'weights' must have three dimensions with the names of 'memb_dim', 'time_dim' and 'dat_dim'.")
      if (dim(weights)[memb_dim] != dim(exp)[memb_dim] |
          dim(weights)[time_dim] != dim(exp)[time_dim] |
          dim(weights)[dat_dim] != dim(exp)[dat_dim]) {
        stop(paste0("Parameter 'weights' must have the same dimension lengths ", 
                    "as 'memb_dim', 'time_dim' and 'dat_dim' in 'exp'."))
      }
      weights <- Reorder(weights, c(time_dim, memb_dim, dat_dim))

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
 
  # Compute RPS
  if (!memb_dim %in% names(dim(obs))) {
    target_dims_obs <- c(time_dim, dat_dim)
  } else {
    target_dims_obs <- c(time_dim, memb_dim, dat_dim)
  }
  rps <- Apply(data = list(exp = exp, obs = obs), 
               target_dims = list(exp = c(time_dim, memb_dim, dat_dim), 
                                  obs = target_dims_obs),
               fun = .RPS,
               dat_dim = dat_dim, time_dim = time_dim,
               memb_dim = memb_dim,
               prob_thresholds = prob_thresholds, 
               indices_for_clim = indices_for_clim, Fair = Fair,
               weights = weights, ncores = ncores)$output1
 
  # Return only the mean RPS
  rps <- MeanDims(rps, time_dim, na.rm = FALSE)
 
  return(rps)
}


.RPS <- function(exp, obs, time_dim = 'sdate', memb_dim = 'member', dat_dim = NULL, 
                 prob_thresholds = c(1/3, 2/3), indices_for_clim = NULL, Fair = FALSE, weights = NULL) {

  # exp: [sdate, memb, (dat)]
  # obs: [sdate, (memb), (dat)]
  # weights: NULL or same as exp

  # Adjust dimensions to be [sdate, memb, dat] for both exp and obs
  if (!memb_dim %in% names(dim(obs))) obs <- InsertDim(obs, posdim = 2, lendim = 1, name = memb_dim)

  if (is.null(dat_dim)) {
    nexp <- 1
    nobs <- 1
    dim(exp) <- c(dim(exp), nexp = nexp)
    dim(obs) <- c(dim(obs), nobs = nobs)
    if (!is.null(weights)) dim(weights) <- c(dim(weights), nexp = nexp)
  } else {
    nexp <- as.numeric(dim(exp)[dat_dim])
    nobs <- as.numeric(dim(obs)[dat_dim])
  }

  rps <- array(dim = c(dim(exp)[time_dim], nexp = nexp, nobs = nobs))

  for (i in 1:nexp) {
    for (j in 1:nobs) {
      exp_data <- exp[ , , i]
      obs_data <- obs[ , , j]

      if (is.null(dim(exp_data))) dim(exp_data) <- c(dim(exp)[1:2])
      if (is.null(dim(obs_data))) dim(obs_data) <- c(dim(obs)[1:2])

      if (!is.null(weights)) {
        weights_data <- weights[ , , i]
        if (is.null(dim(weights_data))) dim(weights_data) <- c(dim(weights)[1:2])
      } else {
        weights_data <- weights
      }

      exp_probs <- .get_probs(data = exp_data, indices_for_quantiles = indices_for_clim, 
                              prob_thresholds = prob_thresholds, weights = weights_data)
      # exp_probs: [bin, sdate]
      obs_probs <- .get_probs(data = obs_data, indices_for_quantiles = indices_for_clim, 
                              prob_thresholds = prob_thresholds, weights = NULL)
      # obs_probs: [bin, sdate]
      probs_exp_cumsum <- apply(exp_probs, 2, cumsum)
      probs_obs_cumsum <- apply(obs_probs, 2, cumsum)

      # rps: [sdate, nexp, nobs]
      rps[ , i, j] <- apply((probs_exp_cumsum - probs_obs_cumsum)^2, 2, sum)
      
      if (Fair) { # FairRPS
        ## adjustment <- rowSums(-1 * (1/R - 1/R.new) * ens.cum * (R - ens.cum)/R/(R - 1)) [formula taken from SpecsVerification::EnsRps]
        R <- dim(exp)[2]  #memb
        R_new <- Inf
        adjustment <- (-1) / (R - 1) * probs_exp_cumsum * (1 - probs_exp_cumsum)
        adjustment <- apply(adjustment, 2, sum)
        rps[ , i, j] <- rps[ , i, j] + adjustment
      }
    }
  }

  if (is.null(dat_dim)) {
    dim(rps) <- dim(exp)[time_dim]
  }

  return(rps)
}

.get_probs <- function(data, indices_for_quantiles, prob_thresholds, weights = NULL) {
  # if exp: [sdate, memb]
  # if obs: [sdate, (memb)]

  # Add dim [memb = 1] to obs if it doesn't have memb_dim
  if (length(dim(data)) == 1) dim(data) <- c(dim(data), 1)

  # Absolute thresholds
  if (is.null(weights)) {
    quantiles <- quantile(x = as.vector(data[indices_for_quantiles, ]), probs = prob_thresholds, type = 8, na.rm = TRUE)
  } else {
    # weights: [sdate, memb]
    sorted_arrays <- .sorted_distributions(data[indices_for_quantiles, ], weights[indices_for_quantiles, ])
    sorted_data <- sorted_arrays$data
    cumulative_weights <- sorted_arrays$cumulative_weights
    quantiles <- approx(cumulative_weights, sorted_data, prob_thresholds, "linear")$y
  }

  # Probabilities
  probs <- array(dim = c(bin = length(quantiles) + 1, dim(data)[1])) # [bin, sdate]
  for (i_time in 1:dim(data)[1]) {
    if (anyNA(data[i_time, ])) {
      probs[, i_time] <- rep(NA, dim = length(quantiles) + 1)
    } else {
      if (is.null(weights)) {
        probs[, i_time] <- colMeans(easyVerification::convert2prob(data[i_time, ], threshold = quantiles))
      } else {
        sorted_arrays <- .sorted_distributions(data[i_time, ], weights[i_time, ])
        sorted_data <- sorted_arrays$data
        cumulative_weights <- sorted_arrays$cumulative_weights
        # find any quantiles that are outside the data range
        integrated_probs <- array(dim = c(bin = length(quantiles)))
        for (i_quant in 1:length(quantiles)) {
          # for thresholds falling under the distribution
          if (quantiles[i_quant] < min(sorted_data)) {
            integrated_probs[i_quant] <- 0 
          # for thresholds falling over the distribution
          } else if (max(sorted_data) < quantiles[i_quant]) {
            integrated_probs[i_quant] <- 1
          } else {
            integrated_probs[i_quant] <- approx(sorted_data, cumulative_weights, quantiles[i_quant], "linear")$y
          }
        }
        probs[, i_time] <- append(integrated_probs, 1) - append(0, integrated_probs)
        if (min(probs[, i_time]) < 0 | max(probs[, i_time]) > 1) {
          stop(paste0("Probability in i_time = ", i_time, " is out of [0, 1]."))
        } 
      }
    }
  }
  return(probs)
}

.sorted_distributions <- function(data_vector, weights_vector) {
  weights_vector <- as.vector(weights_vector)
  data_vector <- as.vector(data_vector)
  weights_vector <- weights_vector / sum(weights_vector) # normalize to 1
  sorter <- order(data_vector)
  sorted_weights <- weights_vector[sorter]
  cumulative_weights <- cumsum(sorted_weights) - 0.5 * sorted_weights
  cumulative_weights <- cumulative_weights - cumulative_weights[1] # fix the 0
  cumulative_weights <- cumulative_weights / cumulative_weights[length(cumulative_weights)] # fix the 1
  return(list("data" = data_vector[sorter], "cumulative_weights" = cumulative_weights))
}
