context("s2dv::RMS tests")

##############################################
  # dat1
  set.seed(1)
  exp1 <- array(rnorm(120), dim = c(member = 3, sdate = 5, ftime = 2, lon = 1, lat = 4))

  set.seed(2)
  obs1 <- array(rnorm(80),  dim = c(member = 2, sdate = 5, ftime = 2, lon = 1, lat = 4))
  set.seed(2)
  na <- floor(runif(10, min = 1, max = 80))
  obs1[na] <- NA

##############################################
test_that("1. Input checks", {

  expect_error(
  RMS(c(), c()),
  "Parameter 'exp' and 'obs' cannot be NULL."
  )
  expect_error(
  RMS(c('b'), c('a')),
  "Parameter 'exp' and 'obs' must be a numeric array."
  )
  expect_error(
  RMS(c(1:10), c(2:4)),
  paste0("Parameter 'exp' and 'obs' must be at least two dimensions ",
         "containing time_dim and memb_dim.")
  )
  expect_error(
  RMS(array(1:10, dim = c(2, 5)), array(1:4, dim = c(2, 2))),
  "Parameter 'exp' and 'obs' must have dimension names."
  )
  expect_error(
  RMS(array(1:10, dim = c(a = 2, c = 5)), array(1:4, dim = c(a = 2, b = 2))),
  "Parameter 'exp' and 'obs' must have same dimension name"
  )
  expect_error(
  RMS(exp1, obs1, memb_dim = 1),
  "Parameter 'memb_dim' must be a character string."
  )
  expect_error(
  RMS(exp1, obs1, memb_dim = 'a'),
  "Parameter 'memb_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  RMS(exp1, obs1, time_dim = c('sdate', 'a')),
  "Parameter 'time_dim' must be a character string."
  )
  expect_error(
  RMS(exp1, obs1, time_dim = 'a'),
  "Parameter 'time_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  RMS(exp1, obs1, comp_dim = c('sdate', 'ftime')),
  "Parameter 'comp_dim' must be a character string."
  )
  expect_error(
  RMS(exp1, obs1, comp_dim = 'a'),
  "Parameter 'comp_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  RMS(exp1, obs1, limits = c(1,3)),
  "Paramter 'comp_dim' cannot be NULL if 'limits' is assigned."
  )
  expect_error(
  RMS(exp1, obs1, comp_dim = 'ftime', limits = c(1)),
  paste0("Parameter 'limits' must be a vector of two positive ",
                  "integers smaller than the length of paramter 'comp_dim'.")
  )
  expect_error(
  RMS(exp1, obs1, conf.lev = -1),
  "Parameter 'conf.lev' must be a numeric number between 0 and 1."
  )
  expect_error(
  RMS(exp1, obs1, conf = 1),
  "Parameter 'conf' must be one logical value."
  )
  expect_error(
  RMS(exp1, obs1, ncores = 1.5),
  "Parameter 'ncores' must be a positive integer."
  )
  expect_error(
  RMS(exp = array(1:10, dim = c(sdate = 1, member = 5, a = 1)),
       obs = array(1:4, dim = c(a = 1, sdate = 2, member = 2))),
  "Parameter 'exp' and 'obs' must have same length of all dimension expect 'memb_dim'."
  )
  expect_error(
  RMS(exp = array(1:5, dim = c(sdate = 1, member = 5, a = 1)),
       obs = array(1:2, dim = c(a = 1, sdate = 1, member = 2))),
  "The length of time_dim must be at least 2 to compute RMS."
  )



})

##############################################
test_that("2. Output checks: dat1", {

  expect_equal(
    dim(RMS(exp1, obs1)$rms),
    c(n_exp = 3, n_obs = 2, ftime = 2, lon = 1, lat = 4)
  )
  expect_equal(
    RMS(exp1, obs1)$rms[1:6],
    c(1.2815677, 2.0832803, 1.1894637, 1.3000403, 1.4053807, 0.8157563), 
    tolerance = 0.001
  )
  expect_equal(
    length(which(is.na(RMS(exp1, obs1)$conf.lower))),
    4
  )
  expect_equal(
    max(RMS(exp1, obs1)$conf.lower, na.rm = T),
    1.399509,
    tolerance = 0.001
  )
  expect_equal(
    length(which(is.na(RMS(exp1, obs1, comp_dim = 'ftime')$rms))),
    0
  )
  expect_equal(
    length(which(is.na(RMS(exp1, obs1, comp_dim = 'ftime')$conf.upper))),
    8
  )
  expect_equal(
    length(which(is.na(RMS(exp1, obs1, comp_dim = 'lat')$conf.lower))),
    36
  )
  expect_equal(
    length(which(is.na(RMS(exp1, obs1, comp_dim = 'lat', limits = c(1, 2))$conf.lower))),
    21
  )
  expect_equal(
    min(RMS(exp1, obs1, conf.lev = 0.99)$conf.upper, na.rm = TRUE),
    1.406368,
    tolerance = 0.0001
  )
  expect_equal(
    length(RMS(exp1, obs1, conf = FALSE)),
    1
  )



})

##############################################

