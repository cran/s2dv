#'Compute the anomaly correlation coefficient between the forecast and corresponding observation
#'
#'Calculate the anomaly correlation coefficient for the ensemble mean of each 
#'model and the corresponding references over a spatial domain. It can return a
#'forecast time series if the data contain forest time dimension, and also the 
#'start date mean if the data contain start date dimension. 
#'The domain of interest can be specified by providing the list 
#'of longitudes/latitudes (lon/lat) of the data together with the corners 
#'of the domain: lonlatbox = c(lonmin, lonmax, latmin, latmax).
#'
#'@param exp A numeric array of experimental anomalies with named dimensions.
#'  It must have at least 'dat_dim' and 'space_dim'.
#'@param obs A numeric array of observational anomalies with named dimensions.
#'  It must have the same dimensions as 'exp' except the length of 'dat_dim' 
#'  and 'memb_dim'.
#'@param dat_dim A character string indicating the name of dataset (nobs/nexp) 
#'  dimension. The default value is 'dataset'.
#'@param space_dim A character string vector of 2 indicating the name of the
#'  latitude and longitude dimensions (in order) along which ACC is computed. 
#'  The default value is c('lat', 'lon').
#'@param avg_dim A character string indicating the name of the dimension to be
#'  averaged. It must be one of 'time_dim'. The mean ACC is calculated along 
#'  averaged. If no need to calculate mean ACC, set as NULL. The default value 
#'  is 'sdate'.
#'@param memb_dim A character string indicating the name of the member 
#'  dimension. If the data are not ensemble ones, set as NULL. The default 
#'  value is 'member'.
#'@param lat A vector of the latitudes of the exp/obs grids. Only required when
#'  the domain of interested is specified. The default value is NULL.
#'@param lon A vector of the longitudes of the exp/obs grids. Only required when
#'  the domain of interested is specified. The default value is NULL.
#'@param lonlatbox A numeric vector of 4 indicating the corners of the domain of
#'  interested: c(lonmin, lonmax, latmin, latmax). Only required when the domain
#'  of interested is specified. The default value is NULL.
#'@param conf A logical value indicating whether to retrieve the confidence 
#'  intervals or not. The default value is TRUE.
#'@param conftype A charater string of "parametric" or "bootstrap". 
#'  "parametric" provides a confidence interval for the ACC computed by a 
#'  Fisher transformation and a significance level for the ACC from a one-sided
#'  student-T distribution. "bootstrap" provides a confidence interval for the
#'  ACC and MACC computed from bootstrapping on the members with 100 drawings 
#'  with replacement. To guarantee the statistical robustness of the result, 
#'  make sure that your experiment and observation always have the same number
#'  of members. "bootstrap" requires 'memb_dim' has value. The default value is
#'  'parametric'.
#'@param conf.lev A numeric indicating the confidence level for the 
#'  regression computation. The default value is 0.95.
#'@param pval A logical value indicating whether to compute the p-value or not.
#'  The default value is TRUE.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return 
#'A list containing the numeric arrays:\cr 
#'\item{acc}{
#'  The ACC with the dimensions c(nexp, nobs, the rest of the dimension except 
#'  space_dim and memb_dim). nexp is the number of experiment (i.e., dat_dim in
#'  exp), and nobs is the number of observation (i.e., dat_dim in obs).
#'}
#'\item{conf.lower (if conftype = "parametric") or acc_conf.lower (if 
#'      conftype = "bootstrap")}{
#'  The lower confidence interval of ACC with the same dimensions as ACC. Only
#'  present if \code{conf = TRUE}.
#'}
#'\item{conf.upper (if conftype = "parametric") or acc_conf.upper (if 
#'      conftype = "bootstrap")}{
#'  The upper confidence interval of ACC with the same dimensions as ACC. Only 
#'  present if \code{conf = TRUE}.
#'}
#'\item{p.val}{
#'  The p-value with the same dimensions as ACC. Only present if 
#'  \code{pval = TRUE} and code{conftype = "parametric"}.  
#'}
#'\item{macc}{
#'  The mean anomaly correlation coefficient with dimensions
#'  c(nexp, nobs, the rest of the dimension except space_dim, memb_dim, and 
#'  avg_dim). Only present if 'avg_dim' is not NULL.
#'}
#'\item{macc_conf.lower}{
#'  The lower confidence interval of MACC with the same dimensions as MACC. 
#'  Only present if \code{conftype = "bootstrap"}.
#'}
#'\item{macc_conf.upper}{
#'  The upper confidence interval of MACC with the same dimensions as MACC. 
#'  Only present if \code{conftype = "bootstrap"}.
#'}
#'
#'@examples
#'  \dontshow{
#'startDates <- c('19851101', '19901101', '19951101', '20001101', '20051101')
#'sampleData <- s2dv:::.LoadSampleData('tos', c('experiment'),
#'                                                c('observation'), startDates,
#'                                                leadtimemin = 1,
#'                                                leadtimemax = 4,
#'                                                output = 'lonlat',
#'                                                latmin = 27, latmax = 48,
#'                                                lonmin = -12, lonmax = 40)
#'  }
#'sampleData$mod <- Season(sampleData$mod, monini = 11, moninf = 12, monsup = 2)
#'sampleData$obs <- Season(sampleData$obs, monini = 11, moninf = 12, monsup = 2) 
#'clim <- Clim(sampleData$mod, sampleData$obs)
#'ano_exp <- Ano(sampleData$mod, clim$clim_exp)
#'ano_obs <- Ano(sampleData$obs, clim$clim_obs)
#'acc <- ACC(ano_exp, ano_obs)
#'acc_bootstrap <- ACC(ano_exp, ano_obs, conftype = 'bootstrap')
#'# Combine acc results for PlotACC
#'res <- array(c(acc$conf.lower, acc$acc, acc$conf.upper, acc$p.val), 
#'             dim = c(dim(acc$acc), 4))
#'res_bootstrap <- array(c(acc$acc_conf.lower, acc$acc, acc$acc_conf.upper, acc$p.val),
#'                       dim = c(dim(acc$acc), 4))
#'  \donttest{
#'PlotACC(res, startDates)
#'PlotACC(res_bootstrap, startDates)
#'  }
#'@references Joliffe and Stephenson (2012). Forecast Verification: A 
#'  Practitioner's Guide in Atmospheric Science. Wiley-Blackwell.
#'@import multiApply
#'@importFrom abind abind
#'@importFrom stats qt qnorm quantile
#'@importFrom ClimProjDiags Subset
#'@export
ACC <- function(exp, obs, dat_dim = 'dataset', space_dim = c('lat', 'lon'),
                avg_dim = 'sdate', memb_dim = 'member', 
                lat = NULL, lon = NULL, lonlatbox = NULL, 
                conf = TRUE, conftype = "parametric", conf.lev = 0.95, pval = TRUE,
                ncores = NULL) {

  # Check inputs 
  ## exp and obs (1)
  if (is.null(exp) | is.null(obs)) {
    stop("Parameter 'exp' and 'obs' cannot be NULL.")
  }
  if (!is.numeric(exp) | !is.numeric(obs)) {
    stop("Parameter 'exp' and 'obs' must be a numeric array.")
  }
  if (is.null(dim(exp)) | is.null(dim(obs))) {
    stop(paste0("Parameter 'exp' and 'obs' must have at least dimensions ",
                "dat_dim and space_dim."))
  }
  if(any(is.null(names(dim(exp))))| any(nchar(names(dim(exp))) == 0) |
     any(is.null(names(dim(obs))))| any(nchar(names(dim(obs))) == 0)) {
    stop("Parameter 'exp' and 'obs' must have dimension names.")
  }
  if(!all(names(dim(exp)) %in% names(dim(obs))) |
     !all(names(dim(obs)) %in% names(dim(exp)))) {
    stop("Parameter 'exp' and 'obs' must have same dimension names.")
  }
  ## dat_dim
  if (!is.character(dat_dim) | length(dat_dim) > 1) {
    stop("Parameter 'dat_dim' must be a character string.")
  }
  if (!dat_dim %in% names(dim(exp)) | !dat_dim %in% names(dim(obs))) {
    stop("Parameter 'dat_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## space_dim
  if (!is.character(space_dim) | length(space_dim) != 2) {
    stop("Parameter 'space_dim' must be a character vector of 2.")
  }
  if (any(!space_dim %in% names(dim(exp))) | any(!space_dim %in% names(dim(obs)))) {
    stop("Parameter 'space_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## avg_dim
  if (!is.null(avg_dim)) {
    if (!is.character(avg_dim) | length(avg_dim) > 1) {
      stop("Parameter 'avg_dim' must be a character string.")
    } 
    if (!avg_dim %in% names(dim(exp)) | !avg_dim %in% names(dim(obs))) {
      stop("Parameter 'avg_dim' is not found in 'exp' or 'obs' dimension.")
    }
  }
  ## memb_dim
  if (!is.null(memb_dim)) {
    if (!is.character(memb_dim) | length(memb_dim) > 1) {
      stop("Parameter 'memb_dim' must be a character string.")
    }
    if (!memb_dim %in% names(dim(exp)) | !memb_dim %in% names(dim(obs))) {
      stop("Parameter 'memb_dim' is not found in 'exp' or 'obs' dimension.")
    }
  }
  ## lat
  if (!is.null(lat)) {
    if (!is.numeric(lat) | length(lat) != dim(exp)[space_dim[1]]) {
      stop(paste0("Parameter 'lat' must be a numeric vector with the same ",
                  "length as the latitude dimension of 'exp' and 'obs'."))
    }
  }
  ## lon
  if (!is.null(lon)) {
    if (!is.numeric(lon) | length(lon) != dim(exp)[space_dim[2]]) {
      stop(paste0("Parameter 'lon' must be a numeric vector with the same ",
                  "length as the longitude dimension of 'exp' and 'obs'."))
    }
  }
  ## lonlatbox
  if (!is.null(lonlatbox)) {
    if (!is.numeric(lonlatbox) | length(lonlatbox) != 4) {
      stop("Parameter 'lonlatbox' must be a numeric vector of 4.")
    }
  }
  ## lat, lon, and lonlatbox
  if (!is.null(lon) & !is.null(lat) & !is.null(lonlatbox)) {
    select_lonlat <- TRUE
  } else if (is.null(lon) & is.null(lat) & is.null(lonlatbox)) {
    select_lonlat <- FALSE
  } else {
    stop(paste0("Parameters 'lon', 'lat', and 'lonlatbox' must be used or be ",
                "NULL at the same time."))
  }
  ## conf
  if (!is.logical(conf) | length(conf) > 1) {
    stop("Parameter 'conf' must be one logical value.")
  }
  if (conf) {
    ## conftype 
    if (!conftype %in% c('parametric', 'bootstrap')) {
      stop("Parameter 'conftype' must be either 'parametric' or 'bootstrap'.")
    }
    if (conftype == 'bootstrap' & is.null(memb_dim)) {
      stop("Parameter 'memb_dim' cannot be NULL when parameter 'conftype' is 'bootstrap'.")
    }
    ## conf.lev
    if (!is.numeric(conf.lev) | conf.lev < 0 | conf.lev > 1 | length(conf.lev) > 1) {
      stop("Parameter 'conf.lev' must be a numeric number between 0 and 1.")
    }
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
  if (!is.null(memb_dim)) {
    name_exp <- name_exp[-which(name_exp == memb_dim)]
    name_obs <- name_obs[-which(name_obs == memb_dim)]
  }
  if(!all(dim(exp)[name_exp] == dim(obs)[name_obs])) {
    stop(paste0("Parameter 'exp' and 'obs' must have same length of ",
                "all the dimensions expect 'dat_dim' and 'memb_dim'."))
  }

#-----------------------------------------------------------------


  ###############################
  # Sort dimension
  name_exp <- names(dim(exp))
  name_obs <- names(dim(obs))
  order_obs <- match(name_exp, name_obs)
  obs <- Reorder(obs, order_obs)
  ###############################

  # Select the domain 
  if (select_lonlat) {
    for (jind in 1:2) {
      while (lonlatbox[jind] < 0) {
        lonlatbox[jind] <- lonlatbox[jind] + 360
      }
      while (lonlatbox[jind] > 360) {
        lonlatbox[jind] <- lonlatbox[jind] - 360
      }
    }
    indlon <- which((lon >= lonlatbox[1] & lon <= lonlatbox[2]) | 
                    (lonlatbox[1] > lonlatbox[2] & (lon > lonlatbox[1] | lon < lonlatbox[2])))
    indlat <- which(lat >= lonlatbox[3] & lat <= lonlatbox[4])

    exp <- ClimProjDiags::Subset(exp, space_dim, list(indlat, indlon), drop = FALSE)
    obs <- ClimProjDiags::Subset(obs, space_dim, list(indlat, indlon), drop = FALSE)
  }

   # Ensemble mean
  if (!is.null(memb_dim)) {
    if (conftype == 'bootstrap') {
      exp_ori <- exp
      obs_ori <- obs
    }
    exp <- MeanDims(exp, memb_dim, na.rm = TRUE)
    obs <- MeanDims(obs, memb_dim, na.rm = TRUE)
  }

  if (is.null(avg_dim)) {
    res <- Apply(list(exp, obs),
                 target_dims = list(c(space_dim, dat_dim),
                                    c(space_dim, dat_dim)),
                 fun = .ACC,
                 dat_dim = dat_dim, avg_dim = avg_dim,
                 conftype = conftype, pval = pval, conf = conf, conf.lev = conf.lev,
                 ncores_input = ncores,
                 ncores = ncores)

    if (conftype == 'bootstrap') {
      res_conf <- Apply(list(exp_ori, obs_ori),
                 target_dims = list(c(memb_dim, dat_dim, space_dim),
                                    c(memb_dim, dat_dim, space_dim)),
                 fun = .ACC_bootstrap,
                 dat_dim = dat_dim, memb_dim = memb_dim, avg_dim = avg_dim,
                 conftype = conftype, pval = pval, conf = conf, conf.lev = conf.lev,
                 ncores_input = ncores,
                 ncores = ncores)
      #NOTE: pval?
      res <- list(acc = res$acc, 
                  acc_conf.lower = res_conf$acc_conf.lower, 
                  acc_conf.upper = res_conf$acc_conf.upper,
                  macc = res$macc,
                  macc_conf.lower = res_conf$macc_conf.lower, 
                  macc_conf.upper = res_conf$macc_conf.upper)
    }

  } else {
    res <- Apply(list(exp, obs),
                 target_dims = list(c(space_dim, avg_dim, dat_dim),
                                    c(space_dim, avg_dim, dat_dim)),
                 fun = .ACC,
                 dat_dim = dat_dim, avg_dim = avg_dim,
                 conftype = conftype, pval = pval, conf = conf, conf.lev = conf.lev,
                 ncores_input = ncores,
                 ncores = ncores)

    if (conftype == 'bootstrap') {
      res_conf <- Apply(list(exp_ori, obs_ori),
                        target_dims = list(c(memb_dim, dat_dim, avg_dim, space_dim),
                                           c(memb_dim, dat_dim, avg_dim, space_dim)),
                        fun = .ACC_bootstrap,
                        dat_dim = dat_dim, memb_dim = memb_dim, avg_dim = avg_dim,
                        conftype = conftype, pval = pval, conf = conf, conf.lev = conf.lev,
                        ncores_input = ncores,
                        ncores = ncores)
      res <- list(acc = res$acc, 
                  acc_conf.lower = res_conf$acc_conf.lower, 
                  acc_conf.upper = res_conf$acc_conf.upper,
                  macc = res$macc,
                  macc_conf.lower = res_conf$macc_conf.lower,
                  macc_conf.upper = res_conf$macc_conf.upper)

    }

  }

 return(res)
}

.ACC <- function(exp, obs,  dat_dim = 'dataset', #space_dim = c('lat', 'lon'),
                 avg_dim = 'sdate', #memb_dim = NULL,
                 lon = NULL, lat = NULL, lonlatbox = NULL,
                 conf = TRUE, conftype = "parametric", conf.lev = 0.95, pval = TRUE,
                 ncores_input = NULL) {

# if (is.null(avg_dim)) 
  # exp: [space_dim, dat_exp]
  # obs: [space_dim, dat_obs]
# if (!is.null(avg_dim)) 
  # exp: [space_dim, avg_dim, dat_exp]
  # obs: [space_dim, avg_dim, dat_obs]

  # .ACC() should use all the spatial points to calculate ACC. It returns [nexp, nobs].

  nexp <- as.numeric(dim(exp)[length(dim(exp))])
  nobs <- as.numeric(dim(obs)[length(dim(obs))])

  if (is.null(avg_dim)) {
    acc <- array(dim = c(nexp = nexp, nobs = nobs))
    if (pval) p.val <- array(dim = c(nexp = nexp, nobs = nobs))
    if (conf) {
      conf.upper <- array(dim = c(nexp = nexp, nobs = nobs))
      conf.lower <- array(dim = c(nexp = nexp, nobs = nobs))
      if (conftype == 'bootstrap') {
         ndraw <- 100
         acc_draw <- array(dim = c(nexp = nexp, nobs = nobs, ndraw)) 
      }
    }
    
  } else {
        acc <- array(dim = c(nexp = nexp, nobs = nobs, dim(exp)[length(dim(exp)) - 1]))
        names(dim(acc))[3] <- avg_dim
        macc <- array(dim = c(nexp = nexp, nobs = nobs))
        if (pval) p.val <- array(dim = c(nexp = nexp, nobs = nobs, dim(exp)[length(dim(exp)) - 1]))
    if (conf) {
        conf.upper <- array(dim = c(nexp = nexp, nobs = nobs, dim(exp)[length(dim(exp)) - 1]))
        conf.lower <- array(dim = c(nexp = nexp, nobs = nobs, dim(exp)[length(dim(exp)) - 1]))
      if (conftype == 'bootstrap') {
         ndraw <- 100
         acc_draw <- array(dim = c(nexp = nexp, nobs = nobs, dim(exp)[length(dim(exp)) - 1], ndraw))
         macc_draw <- array(dim = c(nexp = nexp, nobs = nobs, ndraw))
      }
    }
  }  
  
  # Per-paired exp and obs. NAs should be in the same position in both exp and obs
  for (iobs in 1:nobs) {
    for (iexp in 1:nexp) {
      exp_sub <- ClimProjDiags::Subset(exp, dat_dim, iexp, drop = 'selected')
      obs_sub <- ClimProjDiags::Subset(obs, dat_dim, iobs, drop = 'selected')
      # dim: [space_dim]

      # Variance(iexp) should not take into account any point 
      # that is not available in iobs and therefore not accounted for 
      # in covariance(iexp, iobs) and vice-versa 
      exp_sub[is.na(obs_sub)] <- NA
      obs_sub[is.na(exp_sub)] <- NA

      if (is.null(avg_dim)) {
        # ACC
        top <- sum(exp_sub*obs_sub, na.rm = TRUE)  #a number
        bottom <- sqrt(sum(exp_sub^2, na.rm = TRUE) * sum(obs_sub^2, na.rm = TRUE))
        acc[iexp, iobs] <- top/bottom #a number
        # handle bottom = 0
        if (is.infinite(acc[iexp, iobs])) acc[iexp, iobs] <- NA
        # pval and conf
        if (pval | conf) {
          if (conftype == "parametric") {
            # calculate effective sample size along space_dim
            # combine space_dim into one dim first
            obs_tmp <- array(obs_sub, dim = c(space = length(obs_sub)))
            eno <- Eno(obs_tmp, 'space', ncores = ncores_input)  # a number
            if (pval) {
              t <- qt(conf.lev, eno - 2)  # a number
              p.val[iexp, iobs] <- sqrt(t^2 / (t^2 + eno - 2))
            }
            if (conf) {
              conf.upper[iexp, iobs] <- tanh(atanh(acc[iexp, iobs]) + qnorm(1 - (1 - conf.lev) / 2) / sqrt(eno - 3))
              conf.lower[iexp, iobs] <- tanh(atanh(acc[iexp, iobs]) + qnorm((1 - conf.lev) / 2) / sqrt(eno - 3))
            }
          }
        }

      } else {  #avg_dim is not NULL
        # MACC
        top <- sum(exp_sub*obs_sub, na.rm = TRUE)  #a number
        bottom <- sqrt(sum(exp_sub^2, na.rm = TRUE) * sum(obs_sub^2, na.rm = TRUE))
        macc[iexp, iobs] <- top/bottom #a number
        # handle bottom = 0
        if (is.infinite(macc[iexp, iobs])) macc[iexp, iobs] <- NA
        # ACC
        for (i in 1:dim(acc)[3]) {   #NOTE: use sapply!!!
          exp_sub_i <- ClimProjDiags::Subset(exp_sub, avg_dim, i, drop = 'selected')
          obs_sub_i <- ClimProjDiags::Subset(obs_sub, avg_dim, i, drop = 'selected')
          #dim: [space_dim]
          top <- sum(exp_sub_i*obs_sub_i, na.rm = TRUE)  #a number
          bottom <- sqrt(sum(exp_sub_i^2, na.rm = TRUE) * sum(obs_sub_i^2, na.rm = TRUE))
          acc[iexp, iobs, i] <- top/bottom #a number
          # handle bottom = 0
          if (is.infinite(acc[iexp, iobs, i])) acc[iexp, iobs, i] <- NA
        }

        # pval and conf
        if (pval | conf) {
          if (conftype == "parametric") {
            # calculate effective sample size along space_dim 
            # combine space_dim into one dim first
            obs_tmp <- array(obs_sub, dim = c(space = prod(dim(obs_sub)[-length(dim(obs_sub))]), 
                                              dim(obs_sub)[length(dim(obs_sub))]))
            eno <- Eno(obs_tmp, 'space', ncores = ncores_input)  # a vector of avg_dim
            if (pval) {
              t <- qt(conf.lev, eno - 2)  # a vector of avg_dim
              p.val[iexp, iobs, ] <- sqrt(t^2 / (t^2 + eno - 2))
            }
            if (conf) {
              conf.upper[iexp, iobs, ] <- tanh(atanh(acc[iexp, iobs, ]) + qnorm(1 - (1 - conf.lev) / 2) / sqrt(eno - 3))
              conf.lower[iexp, iobs, ] <- tanh(atanh(acc[iexp, iobs, ]) + qnorm((1 - conf.lev) / 2) / sqrt(eno - 3))
            }
          }
        }

      }  # if avg_dim is not NULL

    }
  }

#------------------------------------------------
  


  # Return output
  if (is.null(avg_dim)) {
    if (conf & pval) {
      return(list(acc = acc, conf.lower = conf.lower, conf.upper = conf.upper,
                  p.val = p.val))
    } else if (conf & !pval) {
      return(list(acc = acc, conf.lower = conf.lower, conf.upper = conf.upper,
                  macc = macc))
    } else if (!conf & pval) {
      return(list(acc = acc, p.val = p.val))
    } else {
      return(list(acc = acc))
    }
  } else {
    if (conf & pval) {
      return(list(acc = acc, conf.lower = conf.lower, conf.upper = conf.upper, 
                  p.val = p.val, macc = macc))
    } else if (conf & !pval) {
      return(list(acc = acc, conf.lower = conf.lower, conf.upper = conf.upper,
                  macc = macc))
    } else if (!conf & pval) {
      return(list(acc = acc, p.val = p.val, macc = macc))
    } else {
      return(list(acc = acc, macc = macc))
    }
  }

}


.ACC_bootstrap <- function(exp, obs,  dat_dim = 'dataset', #space_dim = c('lat', 'lon'),
                           avg_dim = 'sdate', memb_dim = NULL,
                           lon = NULL, lat = NULL, lonlatbox = NULL,
                           conf = TRUE, conftype = "parametric", conf.lev = 0.95, pval = TRUE, 
                           ncores_input = NULL) {
# if (is.null(avg_dim)) 
  # exp: [memb_exp, dat_exp, space_dim]
  # obs: [memb_obs, dat_obs, space_dim]
# if (!is.null(avg_dim)) 
  # exp: [memb_exp, dat_exp, avg_dim, space_dim]
  # obs: [memb_obs, dat_obs, avg_dim, space_dim]

  nexp <- as.numeric(dim(exp)[2])
  nobs <- as.numeric(dim(obs)[2])
  nmembexp <- as.numeric(dim(exp)[1])
  nmembobs <- as.numeric(dim(obs)[1])

  ndraw <- 100
  if (is.null(avg_dim)) {
    acc_draw <- array(dim = c(nexp = nexp, nobs = nobs, ndraw))
  } else {
    acc_draw <- array(dim = c(nexp = nexp, nobs = nobs, dim(exp)[3], ndraw))
    macc_draw <- array(dim = c(nexp = nexp, nobs = nobs, ndraw))
  }

  for (jdraw in 1:ndraw) {
    #choose a randomly member index for each point of the matrix 
    indexp <- array(sample(nmembexp, size = prod(dim(exp)[-c(length(dim(exp)) - 1, length(dim(exp)))]),
                           replace = TRUE), 
                    dim = dim(exp))
    indobs <- array(sample(nmembobs, size = prod(dim(obs)[-c(length(dim(obs)) - 1, length(dim(obs)))]),               
                           replace = TRUE),
                     dim = dim(obs))

      #combine maxtrix of data and random index
      varindexp <- abind::abind(exp, indexp, along = length(dim(exp)) + 1)
      varindobs <- abind::abind(obs, indobs, along = length(dim(obs)) + 1)

    #select randomly the members for each point of the matrix
#    if (is.null(avg_dim)) {

    drawexp <- array( 
                    apply(varindexp, c(2:length(dim(exp))), function(x) x[,1][x[,2]] ),
                          dim = dim(exp)) 
    drawobs <- array(
                    apply(varindobs, c(2:length(dim(obs))), function(x) x[,1][x[,2]] ),
                          dim = dim(obs))

    # ensemble mean before .ACC
    drawexp <- MeanDims(drawexp, memb_dim, na.rm = TRUE)
    drawobs <- MeanDims(drawobs, memb_dim, na.rm = TRUE)
    # Reorder
    if (is.null(avg_dim)) {
      drawexp <- Reorder(drawexp, c(2, 3, 1))
      drawobs <- Reorder(drawobs, c(2, 3, 1))
    } else {
      drawexp <- Reorder(drawexp, c(3, 4, 2, 1))
      drawobs <- Reorder(drawobs, c(3, 4, 2, 1))
    } 
   
    #calculate the ACC of the randomized field
    tmpACC <- .ACC(drawexp, drawobs, conf = FALSE, pval = FALSE, avg_dim = avg_dim,
                   ncores_input = ncores_input)
    if (is.null(avg_dim)) {
      acc_draw[, , jdraw] <- tmpACC$acc
    } else {
      acc_draw[, , , jdraw] <- tmpACC$acc
      macc_draw[, , jdraw] <- tmpACC$macc
    }
  }

  #calculate the confidence interval
  if (is.null(avg_dim)) {
  acc_conf.upper <- apply(acc_draw, c(1, 2),
                          function (x) {
                            quantile(x, 1 - (1 - conf.lev) / 2, na.rm = TRUE)})
  acc_conf.lower <- apply(acc_draw, c(1, 2),
                          function (x) {
                            quantile(x, (1 - conf.lev) / 2, na.rm = TRUE)})

  } else {
  acc_conf.upper <- apply(acc_draw, c(1, 2, 3), 
                          function (x) {
                            quantile(x, 1 - (1 - conf.lev) / 2, na.rm = TRUE)})
  acc_conf.lower <- apply(acc_draw, c(1, 2, 3),
                          function (x) {
                            quantile(x, (1 - conf.lev) / 2, na.rm = TRUE)})
  macc_conf.upper <- apply(macc_draw, c(1, 2),
                          function (x) {
                            quantile(x, 1 - (1 - conf.lev) / 2, na.rm = TRUE)})
  macc_conf.lower <- apply(macc_draw, c(1, 2),
                          function (x) {
                            quantile(x, (1 - conf.lev) / 2, na.rm = TRUE)})
  }

  # Return output
  if (is.null(avg_dim)) {
      return(list(acc_conf.lower = acc_conf.lower,
                  acc_conf.upper = acc_conf.upper))
  } else {
      return(list(acc_conf.lower = acc_conf.lower,
                  acc_conf.upper = acc_conf.upper,
                  macc_conf.lower = macc_conf.lower,
                  macc_conf.upper = macc_conf.upper))
  }
                  
}
