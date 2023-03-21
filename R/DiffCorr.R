#'Compute the correlation difference and its significance
#'
#'Compute the correlation difference between two deterministic forecasts. 
#'Positive values of the correlation difference indicate that the forecast is 
#'more skillful than the reference forecast, while negative values mean that 
#'the reference forecast is more skillful. The statistical significance of the 
#'correlation differences is computed with a one-sided or two-sided test for 
#'equality of dependent correlation coefficients (Steiger, 1980; Siegert et al.,
#'2017) using effective degrees of freedom to account for the autocorrelation of
#'the time series (Zwiers and von Storch, 1995).
#'
#'@param exp A named numerical array of the forecast data with at least time 
#'  dimension. 
#'@param obs A named numerical array with the observations with at least time
#'  dimension. The dimensions must be the same as "exp" except 'memb_dim'.
#'@param ref A named numerical array of the reference forecast data with at 
#'  least time dimension. The dimensions must be the same as "exp" except 
#'  'memb_dim'.
#'@param N.eff Effective sample size to be used in the statistical significance
#'  test. It can be NA (and it will be computed with the s2dv:::.Eno), a 
#'  numeric (which is used for all cases), or an array with the same dimensions
#'  as "obs" except "time_dim" (for a particular N.eff to be used for each case)
#'  . The default value is NA.
#'@param time_dim A character string indicating the name of the time dimension.
#'  The default value is 'sdate'.
#'@param memb_dim A character string indicating the name of the member 
#'  dimension to compute the ensemble mean of the forecast and reference 
#'  forecast. If it is NULL (default), the ensemble mean should be provided
#'  directly to the function.
#'@param method A character string indicating the correlation coefficient to be
#'  computed ("pearson" or "spearman"). The default value is "pearson".
#'@param alpha A numeric of the significance level to be used in the statistical
#'  significance test. If it is a numeric, "sign" will be returned. If NULL, the
#'  p-value will be returned instead. The default value is NULL.
#'@param handle.na A charcater string indicating how to handle missing values.
#'  If "return.na", NAs will be returned for the cases that contain at least one
#'  NA in "exp", "ref", or "obs". If "only.complete.triplets", only the time 
#'  steps with no missing values in all "exp", "ref", and "obs" will be used. If
#'  "na.fail", an error will arise if any of "exp", "ref", or "obs" contains any
#'  NA. The default value is "return.na".
#'@param test.type A character string indicating the type of significance test. 
#'  It can be "two-sided" (to assess whether the skill of "exp" and "ref" are 
#'  significantly different) or "one-sided" (to assess whether the skill of  
#'  "exp" is significantly higher than that of "ref") following Steiger (1980).
#'  The default value is "two-sided".
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return A list with:
#'\item{$diff.corr}{
#'  A numerical array of the correlation differences with the same dimensions as
#'  the input arrays except "time_dim" (and "memb_dim" if provided).
#'}
#'\item{$sign}{
#'  A logical array of the statistical significance of the correlation 
#'  differences with the same dimensions as the input arrays except "time_dim"
#'  (and "memb_dim" if provided). Returned only if "alpha" is a numeric.
#'}
#'\item{$p.val}{
#'  A numeric array of the p-values with the same dimensions as the input arrays
#'  except "time_dim" (and "memb_dim" if provided). Returned only if "alpha" is 
#'  NULL.
#'}
#'
#'@references 
#'Steiger, 1980; https://content.apa.org/doi/10.1037/0033-2909.87.2.245
#'Siegert et al., 2017; https://doi.org/10.1175/MWR-D-16-0037.1
#'Zwiers and von Storch, 1995; https://doi.org/10.1175/1520-0442(1995)008<0336:TSCIAI>2.0.CO;2
#'
#'@examples
#' exp <- array(rnorm(1000), dim = c(lat = 3, lon = 2, member = 10, sdate = 50))
#' obs <- array(rnorm(1000), dim = c(lat = 3, lon = 2, sdate = 50))
#' ref <- array(rnorm(1000), dim = c(lat = 3, lon = 2, member = 10, sdate = 50))
#' res_two.sided_sign <- DiffCorr(exp, obs, ref, memb_dim = 'member', 
#'                                test.type = 'two-sided', alpha = 0.05)
#' res_one.sided_pval <- DiffCorr(exp, obs, ref, memb_dim = 'member', 
#'                                test.type = 'one-sided', alpha = NULL)
#'
#'@import multiApply
#'@export
DiffCorr <- function(exp, obs, ref, N.eff = NA, time_dim = 'sdate', 
                     memb_dim = NULL, method = 'pearson', alpha = NULL, 
                     handle.na = 'return.na', test.type = "two-sided", ncores = NULL) {
  
  # Check inputs
  ## exp, ref, and obs (1)
  if (!is.array(exp) | !is.numeric(exp))
    stop('Parameter "exp" must be a numeric array.')
  if (!is.array(obs) | !is.numeric(obs))
    stop('Parameter "obs" must be a numeric array.')
  if (!is.array(ref) | !is.numeric(ref))
    stop('Parameter "ref" must be a numeric array.')
  ## N.eff
  if (is.array(N.eff)) {
    if (!is.numeric(N.eff)) stop("Parameter 'N.eff' must be numeric.")
    if (!all(names(dim(N.eff)) %in% names(dim(obs))) |
        any(dim(obs)[match(names(dim(N.eff)), names(dim(obs)))] != dim(N.eff))) {
      stop(paste0('If parameter "N.eff" is provided with an array, it must ',
                  'have the same dimensions as "obs" except "time_dim".'))
    }
  } else if (any((!is.na(N.eff) & !is.numeric(N.eff)) | length(N.eff) != 1)) {
    stop(paste0('Parameter "N.eff" must be NA, a numeric, or an array with ',
                'the same dimensions as "obs" except "time_dim".'))
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) != 1)
    stop('Parameter "time_dim" must be a character string.')
  if (!time_dim %in% names(dim(exp)) | !time_dim %in% names(dim(obs)) | 
      !time_dim %in% names(dim(ref))) {
    stop("Parameter 'time_dim' is not found in 'exp', 'obs', or 'ref' dimension.")
  }
  ## memb_dim
  if (!is.null(memb_dim)) {
    if (!is.character(memb_dim) | length(memb_dim) > 1) {
      stop("Parameter 'memb_dim' must be a character string.")
    }
    if (!memb_dim %in% names(dim(exp)) | !memb_dim %in% names(dim(ref))) {
      stop("Parameter 'memb_dim' is not found in 'exp' or 'ref' dimension.")
    }
  }
  ## exp, ref, and obs (2)
  name_exp <- sort(names(dim(exp)))
  name_obs <- sort(names(dim(obs)))
  name_ref <- sort(names(dim(ref)))
  if (!is.null(memb_dim)) {
    name_exp <- name_exp[-which(name_exp == memb_dim)]
    name_ref <- name_ref[-which(name_ref == memb_dim)]
  }
  if (length(name_exp) != length(name_obs) | length(name_exp) != length(name_ref) |
      !identical(dim(exp)[name_exp], dim(obs)[name_obs]) | !identical(dim(exp)[name_exp], dim(ref)[name_ref])) {
    stop(paste0("Parameter 'exp', 'obs', and 'ref' must have same length of ",
                "all dimensions except 'memb_dim'."))
  }
  ## method
  if (!method %in% c("pearson", "spearman")) {
    stop('Parameter "method" must be "pearson" or "spearman".')
  }
  if (method == "spearman") {
    .warning(paste0("The test used in this function is built on Pearson method. ", 
                    "To verify if Spearman method is reliable, you can run the ",
                    "Monte-Carlo simulations that are done in Siegert et al., 2017"))
  }
  ## alpha
  if (!is.null(alpha)) {
    if (any(!is.numeric(alpha) | alpha <= 0 | alpha >= 1 |  
        length(alpha) > 1)) {
      stop('Parameter "alpha" must be NULL or a number between 0 and 1.')
    }
  }
  ## handle.na
  if (!handle.na %in% c('return.na', 'only.complete.triplets', 'na.fail')) {
    stop('Parameter "handle.na" must be "return.na", "only.complete.triplets" or "na.fail".')
  }
  ## test.type
  if (!test.type %in% c('two-sided', 'one-sided')) {
    stop("Parameter 'test.type' must be 'two-sided' or 'one-sided'.")
  }

  ## ncores
  if (!is.null(ncores)) {
    if (any(!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 |
        length(ncores) > 1)) {
      stop('Parameter "ncores" must be either NULL or a positive integer.')
    }
  }
  
  ###############################

  # NA check: na.fail
  if (handle.na == "na.fail" & (anyNA(exp) | anyNA(ref) | anyNA(obs)))
    stop('The data contain NAs.')

  # Calculate ensemble means
  dim_exp <- dim(exp)
  dim_ref <- dim(ref)
  dim_obs <- dim(obs)

  if (!is.null(memb_dim)) {
    exp_memb_dim_ind <- which(names(dim_exp) == memb_dim)
    ref_memb_dim_ind <- which(names(dim_ref) == memb_dim)
    exp <- apply(exp, c(1:length(dim_exp))[-exp_memb_dim_ind], mean, na.rm = FALSE)
    ref <- apply(ref, c(1:length(dim_ref))[-ref_memb_dim_ind], mean, na.rm = FALSE)
    if (is.null(dim(exp))) exp <- array(exp, dim = c(dim_exp[time_dim]))
    if (is.null(dim(ref))) ref <- array(ref, dim = c(dim_ref[time_dim]))
  }

  # output_dims
  if (is.null(alpha)) {
    output_dims <- list(diff.corr = NULL, p.val = NULL)
  } else {
    output_dims <- list(diff.corr = NULL, sign = NULL)
  }
  # Correlation difference
  if (is.array(N.eff)) {
    output <- Apply(data = list(exp = exp, obs = obs, ref = ref,
                                N.eff = N.eff),
                    target_dims = list(exp = time_dim, obs = time_dim,
                                       ref = time_dim, N.eff = NULL),
                    output_dims = output_dims,
                    fun = .DiffCorr, method = method,
                    alpha = alpha, handle.na = handle.na, 
                    test.type = test.type, ncores = ncores)
  } else { 
    output <- Apply(data = list(exp = exp, obs = obs, ref = ref),
                    target_dims = list(exp = time_dim, obs = time_dim, 
                                       ref = time_dim), 
                    output_dims = output_dims, N.eff = N.eff,
                    fun = .DiffCorr, method = method, 
                    alpha = alpha, handle.na = handle.na, 
                    test.type = test.type, ncores = ncores)
  }

  return(output)
}

.DiffCorr <- function(exp, obs, ref, N.eff = NA, method = 'pearson', alpha = NULL, 
                      handle.na = 'return.na', test.type = 'two.sided') {

  .diff.corr <- function(exp, obs, ref, method = 'pearson', N.eff = NA, alpha = NULL, test.type = 'two.sided') {
    # Correlation difference
    cor.exp <- cor(x = exp, y = obs, method = method)
    cor.ref <- cor(x = ref, y = obs, method = method)
    output <- NULL
    output$diff.corr <- cor.exp - cor.ref
    if (is.na(N.eff)) {
      N.eff <- .Eno(x = obs, na.action = na.pass) ## effective degrees of freedom
    }

    # Significance with one-sided or two-sided test for equality of dependent correlation coefficients (Steiger, 1980)
    r12 <- cor.exp
    r13 <- cor.ref
    r23 <- cor(exp, ref)
    R <- (1 - r12 * r12 - r13 * r13 - r23 * r23) + 2 * r12 * r13 * r23
    t <- (r12 - r13) * sqrt((N.eff - 1) * (1 + r23) / (2 * ((N.eff - 1) / (N.eff - 3)) * R + 0.25 * (r12 + r13)^2 * (1 - r23)^3))
    
    if (test.type == 'one-sided') {
      
      ## H0: the skill of exp is not higher than that of ref
      ## H1: the skill of exp is higher than that of ref
      
      p.value <- pt(t, df = N.eff - 3, lower.tail = FALSE)
      
      if (is.null(alpha)) {
        output$p.val <- p.value
      } else {
        output$sign <- ifelse(!is.na(p.value) & p.value <= alpha & output$diff.corr > 0, TRUE, FALSE)
      }

    } else if (test.type == 'two-sided') {
      
      ## H0: the skill difference of exp and ref is zero
      ## H1: the skill difference of exp and ref is different from zero
      
      p.value <- pt(abs(t), df = N.eff - 3, lower.tail = FALSE)

      if (is.null(alpha)) {
        output$p.val <- p.value
      } else {
        output$sign <- ifelse(!is.na(p.value) & p.value <= alpha / 2, TRUE, FALSE)
      }
 
    } else {
      stop("Parameter 'test.type' is not supported.")
    }

    return(output)
  }
  
  #==================================================
  
  if (anyNA(exp) | anyNA(obs) | anyNA(ref)) { ## There are NAs 
    if (handle.na == 'only.complete.triplets') {
      nna <- is.na(exp) | is.na(obs) | is.na(ref) # A vector of T/F
      if (all(nna)) stop("There is no complete set of forecasts and observations.")
      # Remove the incomplete set
      exp <- exp[!nna]
      obs <- obs[!nna]
      ref <- ref[!nna]

      output <- .diff.corr(exp = exp, obs = obs, ref = ref, method = method, 
                           N.eff = N.eff, alpha = alpha, test.type = test.type)
      
    } else if (handle.na == 'return.na') {
      # Data contain NA, return NAs directly without passing to .diff.corr
      if (is.null(alpha)) {
        output <- list(diff.corr = NA, p.val = NA)
      } else {
        output <- list(diff.corr = NA, sign = NA)
      }
    }
    
  } else { ## There is no NA  
    output <- .diff.corr(exp = exp, obs = obs, ref = ref, method = method, 
                         N.eff = N.eff, alpha = alpha, test.type = test.type)
  }

  return(output)
}
