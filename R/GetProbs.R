#'Compute probabilistic forecasts or the corresponding observations
#'
#'Compute probabilistic forecasts from an ensemble based on the relative 
#'thresholds, or the probabilistic observations (i.e., which probabilistic 
#'category was observed). A reference period can be specified to calculate the 
#'absolute thresholds between each probabilistic category. The absolute 
#'thresholds can be computed in cross-validation mode. If data is an ensemble, 
#'the probabilities are calculated as the percentage of members that fall into 
#'each category. For observations (or forecast without member dimension), 1 
#'means that the event happened, while 0 indicates that the event did not 
#'happen. Weighted probabilities can be computed if the weights are provided for 
#'each ensemble member and time step.
#'
#'@param data A named numerical array of the forecasts or observations with, at 
#'  least, time dimension.
#'@param time_dim A character string indicating the name of the time dimension.
#'  The default value is 'sdate'.
#'@param memb_dim A character string indicating the name of the member dimension
#'  to compute the probabilities of the forecast, or NULL if there is no member
#'  dimension (e.g., for observations, or for forecast with only one ensemble
#'  member). The default value is 'member'.
#'@param prob_thresholds A numeric vector of the relative thresholds (from 0 to
#'  1) between the categories. The default value is c(1/3, 2/3), which 
#'  corresponds to tercile equiprobable categories.
#'@param indices_for_quantiles A vector of the indices to be taken along 
#'  'time_dim' for computing the absolute thresholds between the probabilistic
#'  categories. If NULL, the whole period is used. The default value is NULL.
#'@param weights A named numerical array of the weights for 'data' with 
#'  dimensions 'time_dim' and 'memb_dim' (if 'data' has them). The default value
#'  is NULL. The ensemble should have at least 70 members or span at least 10 
#'  time steps and have more than 45 members if consistency between the weighted 
#'  and unweighted methodologies is desired.
#'@param cross.val A logical indicating whether to compute the thresholds 
#'  between probabilistic categories in cross-validation mode. The default value
#'  is FALSE.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return 
#'A numerical array of probabilities with dimensions c(bin, the rest dimensions
#'of 'data' except 'memb_dim'). 'bin' dimension has the length of probabilistic 
#'categories, i.e., \code{length(prob_thresholds) + 1}.
#'
#'@examples
#'data <- array(rnorm(2000), dim = c(ensemble = 25, sdate = 20, time = 4))
#'res <- GetProbs(data = data, time_dim = 'sdate', memb_dim = 'ensemble', 
#'                indices_for_quantiles = 4:17)
#'
#'@import multiApply
#'@importFrom easyVerification convert2prob
#'@export
GetProbs <- function(data, time_dim = 'sdate', memb_dim = 'member', indices_for_quantiles = NULL, 
                     prob_thresholds = c(1/3, 2/3), weights = NULL, cross.val = FALSE, ncores = NULL) {

  # Check inputs
  ## data
  if (is.null(data)) {
    stop("Parameter 'data' cannot be NULL.")
  }
  if (!is.numeric(data)) {
    stop("Parameter 'data' must be a numeric array.")
  }
  if (any(is.null(names(dim(data)))) | any(nchar(names(dim(data))) == 0)) {
    stop("Parameter 'data' must have dimension names.")
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) != 1)
    stop('Parameter "time_dim" must be a character string.')
  if (!time_dim %in% names(dim(data))) {
    stop("Parameter 'time_dim' is not found in 'data' dimensions.")
  }
  ## memb_dim
  if (!is.null(memb_dim)) {
    if (!is.character(memb_dim) | length(memb_dim) > 1) {
      stop("Parameter 'memb_dim' must be a character string.")
    }
    if (!memb_dim %in% names(dim(data))) {
      stop("Parameter 'memb_dim' is not found in 'data' dimensions. If no member ",
           "dimension exists, set it as NULL.")
    }
  }
  ## prob_thresholds
  if (!is.numeric(prob_thresholds) | !is.vector(prob_thresholds) |
      any(prob_thresholds <= 0) | any(prob_thresholds >= 1)) {
    stop("Parameter 'prob_thresholds' must be a numeric vector between 0 and 1.")
  }
  ## indices_for_quantiles
  if (is.null(indices_for_quantiles)) {
    indices_for_quantiles <- 1:dim(data)[time_dim]
  } else {
    if (!is.numeric(indices_for_quantiles) | !is.vector(indices_for_quantiles)) {
      stop("Parameter 'indices_for_quantiles' must be NULL or a numeric vector.")
    } else if (length(indices_for_quantiles) > dim(data)[time_dim] |
               max(indices_for_quantiles) > dim(data)[time_dim] |
               any(indices_for_quantiles < 1)) {
      stop("Parameter 'indices_for_quantiles' should be the indices of 'time_dim'.")
    }
  }
  ## weights
  if (!is.null(weights)) {
    if (!is.array(weights) | !is.numeric(weights))
      stop("Parameter 'weights' must be a named numeric array.")

#    if (is.null(dat_dim)) {
      if (!is.null(memb_dim)) {
        lendim_weights <- 2
        namesdim_weights <- c(time_dim, memb_dim)
      } else {
        lendim_weights <- 1
        namesdim_weights <- c(time_dim)
      }
      if (length(dim(weights)) != lendim_weights | 
          any(!names(dim(weights)) %in% namesdim_weights)) {
        stop(paste0("Parameter 'weights' must have dimension ", 
                    paste0(namesdim_weights, collapse = ' and '), "."))
      }
      if (any(dim(weights)[namesdim_weights] != dim(data)[namesdim_weights])) {
        stop(paste0("Parameter 'weights' must have the same dimension length as ", 
                    paste0(namesdim_weights, collapse = ' and '), " dimension in 'data'."))
      }
      weights <- Reorder(weights, namesdim_weights)

#    } else {
#      if (length(dim(weights)) != 3 | any(!names(dim(weights)) %in% c(memb_dim, time_dim, dat_dim)))
#        stop("Parameter 'weights' must have three dimensions with the names of 'memb_dim', 'time_dim' and 'dat_dim'.")
#      if (dim(weights)[memb_dim] != dim(exp)[memb_dim] |
#          dim(weights)[time_dim] != dim(exp)[time_dim] |
#          dim(weights)[dat_dim] != dim(exp)[dat_dim]) {
#        stop(paste0("Parameter 'weights' must have the same dimension lengths ", 
#                    "as 'memb_dim', 'time_dim' and 'dat_dim' in 'exp'."))
#      }
#      weights <- Reorder(weights, c(time_dim, memb_dim, dat_dim))
#    }
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

  res <- Apply(data = list(data = data), 
               target_dims = c(time_dim, memb_dim), #, dat_dim), 
               output_dims = c("bin", time_dim),
               fun = .GetProbs,
#               dat_dim = dat_dim, 
               prob_thresholds = prob_thresholds, 
               indices_for_quantiles = indices_for_quantiles, 
               weights = weights, cross.val = cross.val, ncores = ncores)$output1
 
  return(res)
}

.GetProbs <- function(data, indices_for_quantiles,
                      prob_thresholds = c(1/3, 2/3), weights = NULL, cross.val = FALSE) {
  # .GetProbs() is used in RPS, RPSS, ROCSS
  # data
  ## if data is exp: [sdate, memb]
  ## if data is obs: [sdate, (memb)]
  # weights: [sdate, (memb)], same as data

  # Add dim [memb = 1] to data if it doesn't have memb_dim
  if (length(dim(data)) == 1) {
    dim(data) <- c(dim(data), 1)
    if (!is.null(weights)) dim(weights) <- c(dim(weights), 1)
  }
  # Absolute thresholds
  if (cross.val) {
    quantiles <- array(NA, dim = c(bin = length(prob_thresholds), sdate = dim(data)[1]))
    for (i_time in 1:dim(data)[1]) {
      if (is.null(weights)) {
        quantiles[, i_time] <- quantile(x = as.vector(data[indices_for_quantiles[which(indices_for_quantiles != i_time)], ]),
                                        probs = prob_thresholds, type = 8, na.rm = TRUE)
      } else {
        # weights: [sdate, memb]
        sorted_arrays <- .sorted_distributions(data[indices_for_quantiles[which(indices_for_quantiles != i_time)], ], 
                                               weights[indices_for_quantiles[which(indices_for_quantiles != i_time)], ])
        sorted_data <- sorted_arrays$data
        cumulative_weights <- sorted_arrays$cumulative_weights
        quantiles[, i_time] <- approx(cumulative_weights, sorted_data, prob_thresholds, "linear")$y
      }
    }

  } else {
    if (is.null(weights)) {
      quantiles <- quantile(x = as.vector(data[indices_for_quantiles, ]), 
                            probs = prob_thresholds, type = 8, na.rm = TRUE)
    } else {
      # weights: [sdate, memb]
      sorted_arrays <- .sorted_distributions(data[indices_for_quantiles, ], 
                                             weights[indices_for_quantiles, ])
      sorted_data <- sorted_arrays$data
      cumulative_weights <- sorted_arrays$cumulative_weights
      quantiles <- approx(cumulative_weights, sorted_data, prob_thresholds, "linear")$y
    }
    quantiles <- array(rep(quantiles, dim(data)[1]), dim = c(bin = length(quantiles), dim(data)[1]))
  }
  # quantiles: [bin-1, sdate]

  # Probabilities
  probs <- array(dim = c(dim(quantiles)[1] + 1, dim(data)[1])) # [bin, sdate]
  for (i_time in 1:dim(data)[1]) {
    if (anyNA(data[i_time, ])) {
      probs[, i_time] <- rep(NA, dim = dim(quantiles)[1] + 1)
    } else {
      if (is.null(weights)) {
        probs[, i_time] <- colMeans(easyVerification::convert2prob(data[i_time, ], 
                                                                   threshold = quantiles[, i_time]))
      } else {
        sorted_arrays <- .sorted_distributions(data[i_time, ], weights[i_time, ])
        sorted_data <- sorted_arrays$data
        cumulative_weights <- sorted_arrays$cumulative_weights
        # find any quantiles that are outside the data range
        integrated_probs <- array(dim = dim(quantiles))
        for (i_quant in 1:dim(quantiles)[1]) {
          # for thresholds falling under the distribution
          if (quantiles[i_quant, i_time] < min(sorted_data)) {
            integrated_probs[i_quant, i_time] <- 0 
          # for thresholds falling over the distribution
          } else if (max(sorted_data) < quantiles[i_quant, i_time]) {
            integrated_probs[i_quant, i_time] <- 1
          } else {
            integrated_probs[i_quant, i_time] <- approx(sorted_data, cumulative_weights, 
                                                        quantiles[i_quant, i_time], "linear")$y
          }
        }
        probs[, i_time] <- append(integrated_probs[, i_time], 1) - append(0, integrated_probs[, i_time])
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

