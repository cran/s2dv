#'Random walk test for skill differences
#'
#'Forecast comparison of the skill obtained with 2 forecasts (with respect to a 
#'common reference) based on Random Walks, with significance estimate at the 95%
#'confidence level, as in DelSole and Tippett (2016).
#'
#'@param skill_A A numerical array of the time series of the skill with the 
#'  forecaster A's.
#'@param skill_B A numerical array of the time series of the skill with the 
#'  forecaster B's. The dimensions should be identical as parameter 'skill_A'.
#'@param time_dim A character string indicating the name of the dimension along
#'  which the tests are computed. The default value is 'sdate'.
#'@param ncores An integer indicating the number of cores to use for parallel 
#'  computation. The default value is NULL.
#'
#'@return A list of 2:
#'\item{$score}{
#'  A numerical array with the same dimensions as the input arrays except 
#'  'time_dim'. The number of times that forecaster A has been better than 
#'  forecaster B minus the number of times that forecaster B has been better 
#'  than forecaster A (for skill positively oriented). If $score is positive 
#'  forecaster A is better than forecaster B, and if $score is negative 
#'  forecaster B is better than forecaster B.
#'}
#'\item{$signif}{
#'  A logical array with the same dimensions as the input arrays except 
#'  'time_dim'. Whether the difference is significant or not at the 5% 
#'  significance level.
#'}
#'
#'@examples
#' fcst_A <- array(c(11:50), dim = c(sdate = 10, lat = 2, lon = 2))
#' fcst_B <- array(c(21:60), dim = c(sdate = 10, lat = 2, lon = 2))
#' reference <- array(1:40, dim = c(sdate = 10, lat = 2, lon = 2))
#' skill_A <- abs(fcst_A - reference)
#' skill_B <- abs(fcst_B - reference)
#' RandomWalkTest(skill_A = skill_A, skill_B = skill_B, time_dim = 'sdate', ncores = 1)
#'
#'@import multiApply
#'@export
RandomWalkTest <- function(skill_A, skill_B, time_dim = 'sdate', ncores = NULL){
  
  ## Check inputs
  if (is.null(skill_A) | is.null(skill_B)){
    stop("Parameters 'skill_A' and 'skill_B' cannot be NULL.")
  }
  if(!is.numeric(skill_A) | !is.numeric(skill_B)){
    stop("Parameters 'skill_A' and 'skill_B' must be a numerical array.")
  }
  if (!identical(dim(skill_A),dim(skill_B))){
    stop("Parameters 'skill_A' and 'skill_B' must have the same dimensions.")
  }
  if(!is.character(time_dim)){
    stop("Parameter 'time_dim' must be a character string.")
  }
  if(!time_dim %in% names(dim(skill_A)) | !time_dim %in% names(dim(skill_B))){
    stop("Parameter 'time_dim' is not found in 'skill_A' or 'skill_B' dimensions.")
  }
  if (!is.null(ncores)){
    if (!is.numeric(ncores) | ncores %% 1 != 0 | ncores <= 0 | length(ncores) > 1){
      stop("Parameter 'ncores' must be a positive integer.")
    }
  }
  
  ## Compute the Random Walk Test
  res <- multiApply::Apply(data = list(skill_A, skill_B), 
                           target_dims = time_dim,
                           fun = .RandomWalkTest, 
                           ncores = ncores)
  return(res)
}

.RandomWalkTest <- function(skill_A, skill_B){
  
   score <- cumsum(skill_A > skill_B) - cumsum(skill_A < skill_B)
   
   ## TRUE if significant (if last value is above or below 2*sqrt(N))
   signif<- ifelse(test = (score[length(skill_A)] <  (-2)*sqrt(length(skill_A))) | (score[length(skill_A)] >  2*sqrt(length(skill_A))), 
                   yes = TRUE, no = FALSE)
   
   return(list("score"=score[length(skill_A)],"signif"=signif))
}
