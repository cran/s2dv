#'Compute the North Atlantic Oscillation (NAO) Index
#'
#'Compute the North Atlantic Oscillation (NAO) index based on the leading EOF 
#'of the sea level pressure (SLP) anomalies over the north Atlantic region 
#'(20N-80N, 80W-40E). The PCs are obtained by projecting the forecast and 
#'observed anomalies onto the observed EOF pattern or the forecast 
#'anomalies onto the EOF pattern of the other years of the forecast. 
#'By default (ftime_avg = 2:4), NAO() computes the NAO index for 1-month 
#'lead seasonal forecasts that can be plotted with PlotBoxWhisker(). It returns
#'cross-validated PCs of the NAO index for forecast (exp) and observations 
#'(obs) based on the leading EOF pattern.
#'
#'@param exp A named numeric array of North Atlantic SLP (20N-80N, 80W-40E) 
#'  forecast anomalies from \code{Ano()} or \code{Ano_CrossValid()} with 
#'  dimensions 'time_dim', 'memb_dim', 'ftime_dim', and 'space_dim' at least.
#'  If only NAO of observational data needs to be computed, this parameter can
#'  be left to NULL. The default value is NULL.
#'@param obs A named numeric array of North Atlantic SLP (20N-80N, 80W-40E) 
#'  observed anomalies from \code{Ano()} or \code{Ano_CrossValid()} with 
#'  dimensions 'time_dim', 'ftime_dim', and 'space_dim' at least.
#'  If only NAO of experimental data needs to be computed, this parameter can 
#'  be left to NULL. The default value is NULL.
#'@param lat A vector of the latitudes of 'exp' and 'obs'.
#'@param lon A vector of the longitudes of 'exp' and 'obs'.
#'@param time_dim A character string indicating the name of the time dimension
#' of 'exp' and 'obs'. The default value is 'sdate'. 
#'@param memb_dim A character string indicating the name of the member 
#'  dimension of 'exp' (and 'obs', optional). If 'obs' has memb_dim, the length
#'  must be 1. The default value is 'member'.
#'@param space_dim A vector of two character strings. The first is the dimension
#'  name of latitude of 'ano' and the second is the dimension name of longitude
#'  of 'ano'. The default value is c('lat', 'lon').
#'@param ftime_dim A character string indicating the name of the forecast time 
#'  dimension of 'exp' and 'obs'. The default value is 'ftime'.
#'@param ftime_avg A numeric vector of the forecast time steps to average
#'  across the target period. If average is not needed, set NULL. The default 
#'  value is 2:4, i.e., from 2nd to 4th forecast time steps.
#'@param obsproj A logical value indicating whether to compute the NAO index by
#'  projecting the forecast anomalies onto the leading EOF of observational 
#'  reference (TRUE) or compute the NAO by first computing the leading 
#'  EOF of the forecast anomalies (in cross-validation mode, i.e. leaving the 
#'  year you are evaluating out), and then projecting forecast anomalies onto 
#'  this EOF (FALSE). The default value is TRUE.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return 
#'A list which contains:
#'\item{exp}{
#'  A numeric array of forecast NAO index in verification format with the same 
#'  dimensions as 'exp' except space_dim and ftime_dim. If ftime_avg is NULL, 
#'  ftime_dim remains.
#'  }
#'\item{obs}{
#'  A numeric array of observed NAO index in verification format with the same
#'  dimensions as 'obs' except space_dim and ftime_dim. If ftime_avg is NULL,
#'  ftime_dim remains.
#'}
#'
#'@references
#'Doblas-Reyes, F.J., Pavan, V. and Stephenson, D. (2003). The skill of 
#'  multi-model seasonal forecasts of the wintertime North Atlantic 
#'  Oscillation. Climate Dynamics, 21, 501-514. 
#'  DOI: 10.1007/s00382-003-0350-4
#'
#'@examples
#'# Make up synthetic data
#'set.seed(1)
#'exp <- array(rnorm(1620), dim = c(member = 2, sdate = 3, ftime = 5, lat = 6, lon = 9))
#'set.seed(2)
#'obs <- array(rnorm(1620), dim = c(member = 1, sdate = 3, ftime = 5, lat = 6, lon = 9))
#'lat <- seq(20, 80, length.out = 6)
#'lon <- seq(-80, 40, length.out = 9) 
#'nao <- NAO(exp = exp, obs = obs, lat = lat, lon = lon)
#'
#'# plot the NAO index
#'  \dontrun{
#'nao$exp <- Reorder(nao$exp, c(2, 1))
#'nao$obs <- Reorder(nao$obs, c(2, 1))
#'PlotBoxWhisker(nao$exp, nao$obs, "NAO index, DJF", "NAO index (PC1) TOS",
#'        monini = 12, yearini = 1985, freq = 1, "Exp. A", "Obs. X")
#'  }
#'
#'@import multiApply
#'@importFrom ClimProjDiags Subset
#'@export
NAO <- function(exp = NULL, obs = NULL, lat, lon, time_dim = 'sdate',
                memb_dim = 'member', space_dim = c('lat', 'lon'),
                ftime_dim = 'ftime', ftime_avg = 2:4, 
                obsproj = TRUE, ncores = NULL) {

  # Check inputs 
  ## exp and obs (1)
  if (is.null(obs) & is.null(exp)) {
    stop("Parameter 'exp' and 'obs' cannot both be NULL.")
  }
  if (!is.null(exp)) {
    if (!is.numeric(exp)) {
      stop("Parameter 'exp' must be a numeric array.")
    }
    if (is.null(dim(exp))) {
      stop(paste0("Parameter 'exp' must have at least dimensions ",
                  "time_dim, memb_dim, space_dim, and ftime_dim."))
    }
    if(any(is.null(names(dim(exp)))) | any(nchar(names(dim(exp))) == 0)) {
      stop("Parameter 'exp' must have dimension names.")
    }
  }
  if (!is.null(obs)) {
    if (!is.numeric(obs)) {
      stop("Parameter 'obs' must be a numeric array.")
    }
    if (is.null(dim(obs))) {
      stop(paste0("Parameter 'obs' must have at least dimensions ",
                  "time_dim, space_dim, and ftime_dim."))
    }
    if(any(is.null(names(dim(obs)))) | any(nchar(names(dim(obs))) == 0)) {
      stop("Parameter 'obs' must have dimension names.")
    }
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) > 1) {
    stop("Parameter 'time_dim' must be a character string.")
  }
  if (!is.null(exp)) {
    if (!time_dim %in% names(dim(exp))) {
      stop("Parameter 'time_dim' is not found in 'exp' or 'obs' dimension.")
    }
  }
  if (!is.null(obs)) {
    if (!time_dim %in% names(dim(obs))) {
      stop("Parameter 'time_dim' is not found in 'exp' or 'obs' dimension.")
    }
  }
  ## memb_dim
  if (!is.character(memb_dim) | length(memb_dim) > 1) {
    stop("Parameter 'memb_dim' must be a character string.")
  }
  if (!is.null(exp)) {
    if (!memb_dim %in% names(dim(exp))) {
      stop("Parameter 'memb_dim' is not found in 'exp' dimension.")
    }
  }
  if (!is.null(obs)) {
    if (memb_dim %in% names(dim(obs))) {
      if (dim(obs)[memb_dim] != 1) {
        stop("The length of parameter 'memb_dim' in 'obs' must be 1.")
      } else {
        add_member_back <- TRUE
        obs <- ClimProjDiags::Subset(obs, memb_dim, 1, drop = 'selected')
      }
    } else {
      add_member_back <- FALSE
    }
  }
  ## space_dim
  if (!is.character(space_dim) | length(space_dim) != 2) {
    stop("Parameter 'space_dim' must be a character vector of 2.")
  }
  if (!is.null(exp)) {
    if (any(!space_dim %in% names(dim(exp)))) {
      stop("Parameter 'space_dim' is not found in 'exp' or 'obs' dimension.")
    }
  }
  if (!is.null(obs)) {
    if (any(!space_dim %in% names(dim(obs)))) {
      stop("Parameter 'space_dim' is not found in 'exp' or 'obs' dimension.")
    }
  }
  ## ftime_dim
  if (!is.character(ftime_dim) | length(ftime_dim) > 1) {
    stop("Parameter 'ftime_dim' must be a character string.")
  }
  if (!is.null(exp)) {
    if (!ftime_dim %in% names(dim(exp))) {
      stop("Parameter 'ftime_dim' is not found in 'exp' or 'obs' dimension.")
    }
  }
  if (!is.null(obs)) {
    if (!ftime_dim %in% names(dim(obs))) {
      stop("Parameter 'ftime_dim' is not found in 'exp' or 'obs' dimension.")
    }
  }
  ## exp and obs (2)
  if (!is.null(exp) & !is.null(obs)) {
    name_exp <- sort(names(dim(exp)))
    name_obs <- sort(names(dim(obs)))
    name_exp <- name_exp[-which(name_exp == memb_dim)]
    throw_error <- FALSE
    if (length(name_exp) != length(name_obs)) {
      throw_error <- TRUE
    } else if (any(name_exp != name_obs)) {
      throw_error <- TRUE
    } else if(!all(dim(exp)[name_exp] == dim(obs)[name_obs])) {
      throw_error <- TRUE
    }
    if (throw_error) {
      stop(paste0("Parameter 'exp' and 'obs' must have the same names and lengths ",
                  "of all the dimensions except 'memb_dim'."))
    }
  }
  ## ftime_avg
  if (!is.null(ftime_avg)) {
    if (!is.vector(ftime_avg) | !is.numeric(ftime_avg)) {
      stop("Parameter 'ftime_avg' must be an integer vector.")
    }
    if (!is.null(exp)) {
      if (max(ftime_avg) > dim(exp)[ftime_dim] | min(ftime_avg) < 1) {
        stop("Parameter 'ftime_avg' must be within the range of ftime_dim length.")
      }
    } else {
      if (max(ftime_avg) > dim(obs)[ftime_dim] | min(ftime_avg) < 1) {
        stop("Parameter 'ftime_avg' must be within the range of ftime_dim length.")
      }
    }
  }
  ## sdate >= 2
  if (!is.null(exp)) {
    if (dim(exp)[time_dim] < 2) {
      stop("The length of time_dim must be at least 2.")
    }
  } else {
    if (dim(obs)[time_dim] < 2) {
      stop("The length of time_dim must be at least 2.")
    }
  }
  ## lat and lon
  if (!is.null(exp)) {
    if (!is.numeric(lat) | length(lat) != dim(exp)[space_dim[1]]) {
      stop(paste0("Parameter 'lat' must be a numeric vector with the same ",
                  "length as the latitude dimension of 'exp' and 'obs'."))
    }
    if (!is.numeric(lon) | length(lon) != dim(exp)[space_dim[2]]) {
      stop(paste0("Parameter 'lon' must be a numeric vector with the same ",
                  "length as the longitude dimension of 'exp' and 'obs'."))
    }
  } else {
    if (!is.numeric(lat) | length(lat) != dim(obs)[space_dim[1]]) {
      stop(paste0("Parameter 'lat' must be a numeric vector with the same ",
                  "length as the latitude dimension of 'exp' and 'obs'."))
    }
    if (!is.numeric(lon) | length(lon) != dim(obs)[space_dim[2]]) {
      stop(paste0("Parameter 'lon' must be a numeric vector with the same ",
                  "length as the longitude dimension of 'exp' and 'obs'."))
    }
  } 
  stop_needed <- FALSE
  if (max(lat) > 80 | min(lat) < 20) {
    stop_needed <- TRUE
  }
  #NOTE: different from s2dverification
  # lon is not used in the calculation actually. EOF only uses lat to do the
  # weight. So we just need to ensure the data is in this region, regardless
  # the order. 
  if (any(lon < 0)) {  #[-180, 180] 
    if (!(min(lon) > -90 & min(lon) < -70 & max(lon) < 50 & max(lon) > 30)) {
      stop_needed <- TRUE
    }
  } else {  #[0, 360]
    if (any(lon >= 50 & lon <= 270)) {
      stop_needed <- TRUE
    } else {
      lon_E <- lon[which(lon < 50)]
      lon_W <- lon[-which(lon < 50)]
      if (max(lon_E) < 30 | min(lon_W) > 290) {
        stop_needed <- TRUE
      }
    }
  }
  if (stop_needed) {
    stop(paste0("The typical domain used to compute the NAO is 20N-80N, ",
                "80W-40E. 'lat' or 'lon' is out of range."))
  }
  ## obsproj
  if (!is.logical(obsproj)  | length(obsproj) > 1) {
    stop("Parameter 'obsproj' must be either TRUE or FALSE.")
  }
  if (obsproj) {
    if (is.null(obs)) {
      stop("Parameter 'obsproj' set to TRUE but no 'obs' provided.")
    }
    if (is.null(exp)) {
      .warning("parameter 'obsproj' set to TRUE but no 'exp' provided.")
    }
  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores == 0 |
      length(ncores) > 1) {
      stop("Parameter 'ncores' must be a positive integer.")
    }
  }

  # Average ftime 
  if (!is.null(ftime_avg)) {
    if (!is.null(exp)) {
      exp_sub <- ClimProjDiags::Subset(exp, ftime_dim, ftime_avg, drop = FALSE)
      exp <- MeanDims(exp_sub, ftime_dim, na.rm = TRUE)
      ## Cross-validated PCs. Fabian. This should be extended to
      ## nmod and nlt by simple loops. Virginie
    }
    if (!is.null(obs)) {
      obs_sub <- ClimProjDiags::Subset(obs, ftime_dim, ftime_avg, drop = FALSE)
      obs <- MeanDims(obs_sub, ftime_dim, na.rm = TRUE)
    }
  }

  # wght
  wght <- array(sqrt(cos(lat * pi/180)), dim = c(length(lat), length(lon)))

  if (!is.null(exp) & !is.null(obs)) {
    res <- Apply(list(exp, obs),
                 target_dims = list(exp = c(memb_dim, time_dim, space_dim),
                                    obs = c(time_dim, space_dim)),
                 fun = .NAO,
                 lat = lat, wght = wght, 
                 obsproj = obsproj, add_member_back = add_member_back,
                 ncores = ncores)
  } else if (!is.null(exp)) {
    res <- Apply(list(exp = exp),
                 target_dims = list(exp = c(memb_dim, time_dim, space_dim)),
                 fun = .NAO,
                 lat = lat, wght = wght, obs = NULL, 
                 obsproj = obsproj, add_member_back = FALSE,
                 ncores = ncores)
  } else if (!is.null(obs)) {
    if (add_member_back) {
      output_dims <- list(obs = c(time_dim, memb_dim))
    } else {
      output_dims <- list(obs = time_dim)
    }
    res <- Apply(list(obs = obs),
                 target_dims = list(obs = c(time_dim, space_dim)),
                 output_dims = output_dims,
                 fun = .NAO,
                 lat = lat, wght = wght, exp = NULL, 
                 obsproj = obsproj, add_member_back = add_member_back,
                 ncores = ncores)
  }
  return(res)
}

.NAO <- function(exp = NULL, obs = NULL, lat, wght, obsproj = TRUE, add_member_back = FALSE) {
  # exp: [memb_exp, sdate, lat, lon]
  # obs: [sdate, lat, lon]
  # wght: [lat, lon]

  if (!is.null(exp)) {
    ntime <- dim(exp)[2]
    nlat <- dim(exp)[3]
    nlon <- dim(exp)[4]
    nmemb_exp <- dim(exp)[1]
  } else {
    ntime <- dim(obs)[1]
    nlat <- dim(obs)[2]
    nlon <- dim(obs)[3]
  }

  if (!is.null(obs)) NAOO.ver <- array(NA, dim = ntime)
  if (!is.null(exp)) NAOF.ver <- array(NA, dim = c(ntime, nmemb_exp))

  for (tt in 1:ntime) {  #sdate

    if (!is.null(obs)) {
      ## Calculate observation EOF. Excluding one forecast start year.
      obs_sub <- obs[c(1:ntime)[-tt], , , drop = FALSE]
      obs_EOF <- .EOF(obs_sub, neofs = 1, wght = wght)  # $EOFs: [mode, lat, lon]

      ## Correct polarity of pattern.
      # dim(obs_EOF$EOFs): [mode, lat, lon]
      if (0 < mean(obs_EOF$EOFs[1, which.min(abs(lat - 65)), ], na.rm = T)) {
        obs_EOF$EOFs <- obs_EOF$EOFs * (-1)
#        obs_EOF$PCs <- obs_EOF$PCs * (-1)  # not used
      }
      ## Project observed anomalies.
      PF <- .ProjectField(obs, eof_mode = obs_EOF$EOFs[1, , ], wght = wght)  # [sdate]
      ## Keep PCs of excluded forecast start year. Fabian.
      NAOO.ver[tt] <- PF[tt]
    }

    if (!is.null(exp)) {
      if (!obsproj) {
        exp_sub <- exp[, c(1:ntime)[-tt], , , drop = FALSE]
        # Combine 'memb' and 'sdate' to calculate EOF
        dim(exp_sub) <- c(nmemb_exp * (ntime - 1), nlat, nlon)
        exp_EOF <- .EOF(exp_sub, neofs = 1, wght = wght)  # $EOFs: [mode, lat, lon]

        ## Correct polarity of pattern.
        ##NOTE: different from s2dverification, which doesn't use mean().
#        if (0 < exp_EOF$EOFs[1, which.min(abs(lat - 65)), ]) {
        if (0 < mean(exp_EOF$EOFs[1, which.min(abs(lat - 65)), ], na.rm = T)) {
          exp_EOF$EOFs <- exp_EOF$EOFs * (-1)
#          exp_EOF$PCs <- exp_EOF$PCs * sign  # not used
        }

        ### Lines below could be simplified further by computing
        ### ProjectField() only on the year of interest... (though this is
        ### not vital). Lauriane
        for (imemb in 1:nmemb_exp) {
          PF <- .ProjectField(exp[imemb, , , ], eof_mode = exp_EOF$EOFs[1, , ], wght = wght) # [sdate, memb]
          NAOF.ver[tt, imemb] <- PF[tt]
        }    
      } else {
        ## Project forecast anomalies on obs EOF
        for (imemb in 1:nmemb_exp) {
          PF <- .ProjectField(exp[imemb, , , ], eof_mode = obs_EOF$EOFs[1, , ], wght = wght)  # [sdate]
          NAOF.ver[tt, imemb] <- PF[tt]
        }
      }
    }

  }  # for loop sdate

  # add_member_back
  if (add_member_back) {
    suppressWarnings( 
    NAOO.ver <- InsertDim(NAOO.ver, 2, 1, name = names(dim(exp))[1])
    )
  }

  #NOTE: EOFs_obs is not returned because it's only the result of the last sdate
  #     (It is returned in s2dverification.)
  if (!is.null(exp) & !is.null(obs)) {
    return(list(exp = NAOF.ver, obs = NAOO.ver)) #, EOFs_obs = obs_EOF))
  } else if (!is.null(exp)) {
    return(list(exp = NAOF.ver))
  } else if (!is.null(obs)) {
    return(list(obs = NAOO.ver))
  }
}
