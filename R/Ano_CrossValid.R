#'Compute anomalies in cross-validation mode
#'
#'Compute the anomalies from the arrays of the experimental and observational 
#'data output by subtracting the climatologies computed with a leave-one-out 
#'cross validation technique and a per-pair method (Garcia-Serrano and 
#'Doblas-Reyes, CD, 2012).
#'Per-pair climatology means that only the start dates covered by the 
#'whole experiments/observational datasets will be used. In other words, the 
#'startdates which do not all have values along 'dat_dim' dimension of both
#'the 'exp' and 'obs' are excluded when computing the climatologies.
#'
#'@param exp A named numeric array of experimental data, with at least 
#'  dimensions 'time_dim' and 'dat_dim'.
#'@param obs A named numeric array of observational data, same dimensions as
#'  parameter 'exp' except along 'dat_dim'.
#'@param time_dim A character string indicating the name of the time dimension.  
#'  The default value is 'sdate'.
#'@param dat_dim A character vector indicating the name of the dataset and 
#'  member dimensions. When calculating the climatology, if data at one 
#'  startdate (i.e., 'time_dim') is not complete along 'dat_dim', this startdate
#'  along 'dat_dim' will be discarded. The default value is 
#'  "c('dataset', 'member')".
#'@param memb_dim A character string indicating the name of the member 
#'  dimension. Only used when parameter 'memb' is FALSE. It must be one element
#'  in 'dat_dim'. The default value is 'member'.
#'@param memb A logical value indicating whether to subtract the climatology 
#'  based on the individual members (TRUE) or the ensemble mean over all
#'  members (FALSE) when calculating the anomalies. The default value is TRUE.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return
#'A list of 2:
#'\item{$exp}{
#'  A numeric array with the same dimensions as 'exp'. The dimension order may
#'  change.
#'}
#'\item{$obs}{
#'  A numeric array with the same dimensions as 'obs'.The dimension order may
#'  change.
#'}
#'
#'@examples 
#'# Load sample data as in Load() example:
#'example(Load)
#'anomalies <- Ano_CrossValid(sampleData$mod, sampleData$obs)
#'\dontrun{
#'PlotAno(anomalies$exp, anomalies$obs, startDates, 
#'        toptitle = paste('anomalies'), ytitle = c('K', 'K', 'K'), 
#'        legends = 'ERSST', biglab = FALSE, fileout = 'tos_ano_crossvalid.eps')
#'}
#'@import multiApply
#'@importFrom ClimProjDiags Subset
#'@export
Ano_CrossValid <- function(exp, obs, time_dim = 'sdate', dat_dim = c('dataset', 'member'),
                           memb_dim = 'member', memb = TRUE, ncores = NULL) {

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
                "time_dim and dat_dim."))
  }
  if(any(is.null(names(dim(exp))))| any(nchar(names(dim(exp))) == 0) |
     any(is.null(names(dim(obs))))| any(nchar(names(dim(obs))) == 0)) {
    stop("Parameter 'exp' and 'obs' must have dimension names.")
  }
  if(!all(names(dim(exp)) %in% names(dim(obs))) |
     !all(names(dim(obs)) %in% names(dim(exp)))) {
    stop("Parameter 'exp' and 'obs' must have same dimension name.")
  }
  ## time_dim
  if (!is.character(time_dim) | length(time_dim) > 1) {
    stop("Parameter 'time_dim' must be a character string.")
  }
  if (!time_dim %in% names(dim(exp)) | !time_dim %in% names(dim(obs))) {
    stop("Parameter 'time_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## dat_dim
  if (!is.character(dat_dim)) {
    stop("Parameter 'dat_dim' must be a character vector.")
  }
  if (!all(dat_dim %in% names(dim(exp))) | !all(dat_dim %in% names(dim(obs)))) {
    stop("Parameter 'dat_dim' is not found in 'exp' or 'obs' dimension.")
  }
  ## memb
  if (!is.logical(memb) | length(memb) > 1) {
    stop("Parameter 'memb' must be one logical value.")
  }
  ## memb_dim
  if (!memb) {
    if (!is.character(memb_dim) | length(memb_dim) > 1) {
      stop("Parameter 'memb_dim' must be a character string.")
    }
    if (!memb_dim %in% names(dim(exp)) | !memb_dim %in% names(dim(obs))) {
      stop("Parameter 'memb_dim' is not found in 'exp' or 'obs' dimension.")
    }
    if (!memb_dim %in% dat_dim) {
      stop("Parameter 'memb_dim' must be one element in parameter 'dat_dim'.")
    }     
  }
  ## ncores
  if (!is.null(ncores)) {
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 |
      length(ncores) > 1) {
      stop("Parameter 'ncores' must be a positive integer.")
    }
  }
  ## exp and obs (2)
  name_exp <- sort(names(dim(exp)))
  name_obs <- sort(names(dim(obs)))
  for (i in 1:length(dat_dim)) {
    name_exp <- name_exp[-which(name_exp == dat_dim[i])]
    name_obs <- name_obs[-which(name_obs == dat_dim[i])]
  }
  if(!all(dim(exp)[name_exp] == dim(obs)[name_obs])) {
    stop(paste0("Parameter 'exp' and 'obs' must have the same length of ",
                "all dimensions except 'dat_dim'."))
  }

  ###############################
  # Sort dimension
  name_exp <- names(dim(exp))
  name_obs <- names(dim(obs))
  order_obs <- match(name_exp, name_obs)
  if (any(order_obs != sort(order_obs))) {
    obs <- Reorder(obs, order_obs)
  }

  #-----------------------------------
  # Per-paired method: If any sdate along dat_dim is NA, turn all sdate points along dat_dim into NA.
   pos <- rep(0, length(dat_dim))  # dat_dim: [dataset, member]
   for (i in 1:length(dat_dim)) {
     pos[i] <- which(names(dim(obs)) == dat_dim[i])
   }
   outrows_exp <- MeanDims(exp, pos, na.rm = FALSE) +
                  MeanDims(obs, pos, na.rm = FALSE)
   outrows_obs <- outrows_exp

   for (i in 1:length(pos)) {
     outrows_exp <- InsertDim(outrows_exp, pos[i], dim(exp)[pos[i]])
     outrows_obs <- InsertDim(outrows_obs, pos[i], dim(obs)[pos[i]])
   }
   exp_for_clim <- exp
   obs_for_clim <- obs
   exp_for_clim[which(is.na(outrows_exp))] <- NA
   obs_for_clim[which(is.na(outrows_obs))] <- NA

  #-----------------------------------

    res <- Apply(list(exp, obs, exp_for_clim, obs_for_clim),
                  target_dims = c(time_dim, dat_dim),
                  fun = .Ano_CrossValid,
                  memb_dim = memb_dim, memb = memb,
                  ncores = ncores)

  return(res)
}

.Ano_CrossValid <- function(exp, obs, exp_for_clim, obs_for_clim,
                            memb_dim = 'member', memb = TRUE, ncores = NULL) {
  # exp: [sdate, dat_dim, memb_dim]
  # obs: [sdate, dat_dim, memb_dim]
  ano_exp_list <- vector('list', length = dim(exp)[1])  #length: [sdate]
  ano_obs_list <- vector('list', length = dim(obs)[1])  

  for (tt in 1:dim(exp)[1]) {  #[sdate]
    # calculate clim
    exp_sub <- ClimProjDiags::Subset(exp_for_clim, 1, c(1:dim(exp)[1])[-tt])
    obs_sub <- ClimProjDiags::Subset(obs_for_clim, 1, c(1:dim(obs)[1])[-tt])
    clim_exp <- apply(exp_sub, c(1:length(dim(exp)))[-1], mean, na.rm = TRUE)  # average out time_dim -> [dat, memb]
    clim_obs <- apply(obs_sub, c(1:length(dim(obs)))[-1], mean, na.rm = TRUE)

    # ensemble mean
    if (!memb) {
      if (is.null(dim(clim_exp)) | length(dim(clim_exp)) == 1) {   #dim: [member]
        clim_exp <- mean(clim_exp, na.rm = TRUE)  # a number
        clim_obs <- mean(clim_obs, na.rm = TRUE)
      } else {
        pos <- which(names(dim(clim_exp)) == memb_dim)
        pos <- c(1:length(dim(clim_exp)))[-pos]
        dim_name <- names(dim(clim_exp))
        dim_exp_ori <- dim(clim_exp)
        dim_obs_ori <- dim(clim_obs)

        clim_exp <- apply(clim_exp, pos, mean, na.rm = TRUE)
        clim_obs <- apply(clim_obs, pos, mean, na.rm = TRUE)
        if (is.null(names(dim(as.array(clim_exp))))) {
          clim_exp <- as.array(clim_exp)
          clim_obs <- as.array(clim_obs)
          names(dim(clim_exp)) <- dim_name[pos]
          names(dim(clim_obs)) <- dim_name[pos]
        }
      
        # Expand it back 
        clim_exp_tmp <- array(clim_exp, dim = c(dim_exp_ori[pos], dim_exp_ori[-pos]))
        clim_obs_tmp <- array(clim_obs, dim = c(dim_obs_ori[pos], dim_obs_ori[-pos]))
        # Reorder it back to dim(clim_exp)
        tmp <- match(dim_exp_ori, dim(clim_exp_tmp))
        if (any(tmp != sort(tmp))) {
          clim_exp <- Reorder(clim_exp_tmp, tmp)
          clim_obs <- Reorder(clim_obs_tmp, tmp)
        } else {
          clim_exp <- clim_exp_tmp
          clim_obs <- clim_obs_tmp
        }
      }
    }
    # calculate ano
    ano_exp_list[[tt]] <- ClimProjDiags::Subset(exp, 1, tt, drop = 'selected') - clim_exp
    ano_obs_list[[tt]] <- ClimProjDiags::Subset(obs, 1, tt, drop = 'selected') - clim_obs
  }

  ano_exp <- array(unlist(ano_exp_list), dim = c(dim(exp)[-1], dim(exp)[1]))
  ano_exp <- Reorder(ano_exp, c(length(dim(exp)), 1:(length(dim(exp)) - 1)))
  ano_obs <- array(unlist(ano_obs_list), dim = c(dim(obs)[-1], dim(obs)[1]))
  ano_obs <- Reorder(ano_obs, c(length(dim(obs)), 1:(length(dim(obs)) - 1)))

  return(list(exp = ano_exp, obs = ano_obs))
}
