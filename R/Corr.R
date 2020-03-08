#'Compute the correlation coefficient between an array of forecast and their corresponding observation
#'
#'Calculate the correlation coefficient (Pearson, Kendall or Spearman) for 
#'an array of forecast and an array of observation. The correlations are 
#'computed along time_dim, the startdate dimension. If comp_dim is given, 
#'the correlations are computed only if obs along the comp_dim dimension are 
#'complete between limits[1] and limits[2], i.e., there is no NA between 
#'limits[1] and limits[2]. This option can be activated if the user wants to 
#'account only for the forecasts which the corresponding observations are 
#'available at all leadtimes.\cr 
#'The confidence interval is computed by the Fisher transformation and the 
#'significance level relies on an one-sided student-T distribution.\cr 
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
#'  be completed. The default is c(1, length(comp_dim dimension)).
#'@param method A character string indicating the type of correlation: 
#'  'pearson', 'spearman', or 'kendall'. The default value is 'pearson'.
#'@param pval A logical value indicating whether to compute or not the p-value 
#'  of the test Ho: Corr = 0. The default value is TRUE.
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
#'\item{$corr}{
#'  The correlation coefficient. 
#'}
#'\item{$p.val}{
#'  The p-value. Only present if \code{pval = TRUE}.
#'}
#'\item{$conf.lower}{
#'  The lower confidence interval. Only present if \code{conf = TRUE}.
#'}
#'\item{$conf.upper}{
#'  The upper confidence interval. Only present if \code{conf = TRUE}.
#'}
#'
#'@examples
#'# Load sample data as in Load() example: 
#'example(Load) 
#'clim <- Clim(sampleData$mod, sampleData$obs) 
#'corr <- Corr(clim$clim_exp, clim$clim_obs, time_dim = 'ftime')
#'
#'@rdname Corr
#'@import multiApply
#'@importFrom ClimProjDiags Subset
#'@importFrom stats cor pt qnorm 
#'@export
Corr <- function(exp, obs, time_dim = 'sdate', memb_dim = 'member', 
                 comp_dim = NULL, limits = NULL,
                 method = 'pearson', pval = TRUE, conf = TRUE,
                 conf.lev = 0.95, ncores = NULL) {

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
  ## method
  if (!(method %in% c("kendall", "spearman", "pearson"))) {
    stop("Parameter 'method' must be one of 'kendall', 'spearman' or 'pearson'.")
  }
  ## pval
  if (!is.logical(pval) | length(pval) > 1) {
    stop("Parameter 'pval' must be one logical value.")
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
  if (dim(exp)[time_dim] < 3) {
    stop("The length of time_dim must be at least 3 to compute correlation.")
  }


  ###############################
  # Sort dimension
  name_exp <- names(dim(exp))
  name_obs <- names(dim(obs))
  order_obs <- match(name_exp, name_obs)
  obs <- Reorder(obs, order_obs)


  ###############################
  # Calculate Corr

  # Remove data along comp_dim dim if there is at least one NA between limits
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
              fun = .Corr, 
              time_dim = time_dim, method = method,
              pval = pval, conf = conf, conf.lev = conf.lev, 
              ncores = ncores)
 return(res)
}

.Corr <- function(exp, obs, time_dim = 'sdate', method = 'pearson',
                  conf = TRUE, pval = TRUE, conf.lev = 0.95) {

  # exp: [sdate, member_exp]
  # obs: [sdate, member_obs]
  n_exp <- as.numeric(dim(exp)[2])
  n_obs <- as.numeric(dim(obs)[2])
  
  CORR <- array(dim = c(n_exp = n_exp, n_obs = n_obs))
  eno_expand <- array(dim = c(n_exp = n_exp, n_obs = n_obs))
  p.val <- array(dim = c(n_exp = n_exp, n_obs = n_obs))

  # ens_mean
  for (i in 1:n_obs) {

    CORR[, i] <- sapply(1:n_exp, 
                          function(x) {
    if (any(!is.na(exp[, x])) && sum(!is.na(obs[, i])) > 2) {
cor(exp[, x], obs[, i],
                                       use = "pairwise.complete.obs", 
                                       method = method)
} else {
    CORR[, i] <- NA
}
})
  }

#  if (pval) {
#    for (i in 1:n_obs) {
#      p.val[, i] <- try(sapply(1:n_exp,
#                           function(x) {(cor.test(exp[, x], obs[, i],
#                                         use = "pairwise.complete.obs",
#                                         method = method)$p.value)/2}), silent = TRUE)
#      if (class(p.val[, i]) == 'character') {
#        p.val[, i] <- NA
#      }
#    }
#  }

  if (pval | conf) {
    if (method == "kendall" | method == "spearman") {
      tmp <- apply(obs, 2, rank)
      names(dim(tmp))[1] <- time_dim
      eno <- Eno(tmp, time_dim)
    } else if (method == "pearson") {
      eno <- Eno(obs, time_dim)  
    }
    for (i in 1:n_exp) {
      eno_expand[i, ] <- eno
    }
  }
#############old#################
#This doesn't return error but it's diff from cor.test() when method is spearman and kendall
  if (pval) {
    t <-sqrt(CORR * CORR * (eno_expand - 2) / (1 - (CORR ^ 2)))
    p.val <- pt(t, eno_expand - 2, lower.tail = FALSE)
  } 
###################################
  if (conf) {
    conf.lower <- (1 - conf.lev) / 2
    conf.upper <- 1 - conf.lower
    conflow <- tanh(atanh(CORR) + qnorm(conf.lower) / sqrt(eno_expand - 3))
    confhigh <- tanh(atanh(CORR) + qnorm(conf.upper) / sqrt(eno_expand - 3))
  }

  if (pval & conf) {
    res <- list(corr = CORR, p.val = p.val, 
                conf.lower = conflow, conf.upper = confhigh)
  } else if (pval & !conf) {
    res <- list(corr = CORR, p.val = p.val)
  } else if (!pval & conf) {
    res <- list(corr = CORR,
                conf.lower = conflow, conf.upper = confhigh)
  } else {
    res <- list(corr = CORR)
  }

  return(res) 

}
