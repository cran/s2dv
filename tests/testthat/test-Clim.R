context("s2dv::Clim tests")

##############################################
  # dat1
  set.seed(1)
  exp1 <- array(rnorm(360), dim = c(dataset = 1, member = 3, sdate = 5, 
                                    ftime = 3, lon = 2, lat = 4))
  set.seed(2)
  obs1 <- array(rnorm(120), dim = c(dataset = 1, member = 1,
                                   ftime = 3, lon = 2, lat = 4, sdate = 5))
  # dat2
  exp2 <- exp1
  set.seed(1)
  na <- floor(runif(5, min = 1, max = 360))
  exp2[na] <- NA
  obs2 <- obs1
  set.seed(2)
  na <- floor(runif(30, min = 1, max = 120))
  obs2[na] <- NA

##############################################
test_that("1. Input checks", {

  expect_error(
  Clim(c(), c()),
  "Parameter 'exp' and 'obs' cannot be NULL."
  )
  expect_error(
  Clim(c('b'), c('a')),
  "Parameter 'exp' and 'obs' must be a numeric array."
  )
  expect_error(
  Clim(c(1:10), c(2:4)),
  paste0("Parameter 'exp' and 'obs' must be at least two dimensions ",
         "containing time_dim and dat_dim.")
  )
  expect_error(
  Clim(array(1:10, dim = c(2, 5)), array(1:4, dim = c(2, 2))),
  "Parameter 'exp' and 'obs' must have dimension names."
  )
  expect_error(
  Clim(array(1:10, dim = c(a = 2, c = 5)), array(1:4, dim = c(a = 2, b = 2))),
  "Parameter 'exp' and 'obs' must have same dimension name"
  )
  expect_error(
  Clim(exp1, obs1, method = TRUE),
  "Parameter 'method' must be one of 'clim', 'kharin' or 'NDV'."
  )
  expect_error(
  Clim(exp1, obs1, time_dim = c('sdate','ftime')),
  "Parameter 'time_dim' must be a character string."
  )
  expect_error(
  Clim(exp1, obs1, time_dim = 'asd'),
  "Parameter 'time_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  Clim(exp1, obs1, dat_dim = c(1,2)),
  "Parameter 'dat_dim' must be a character vector."
  )
  expect_error(
  Clim(exp1, obs1, dat_dim = c('member', 'dat')),
  "Parameter 'dat_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  Clim(exp1, obs1, ftime_dim = 4),
  "Parameter 'ftime_dim' must be a character string."
  )
  expect_error(
  Clim(exp1, obs1, ftime_dim = 'f'),
  "Parameter 'ftime_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  Clim(exp1, obs1, memb_dim = c('dataset', 'member')),
  "Parameter 'memb_dim' must be a character string."
  )
  expect_error(
  Clim(exp1, obs1, memb_dim = 'memb'),
  "Parameter 'memb_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  Clim(exp1, obs1, memb = 'member'),
  "Parameter 'memb' must be one logical value."
  )
  expect_error(
  Clim(exp1, obs1, na.rm = na.omit),
  "Parameter 'na.rm' must be one logical value."
  )
  expect_error(
  Clim(exp1, obs1, ncores = T),
  "Parameter 'ncores' must be a positive integer."
  )
  expect_error(
  Clim(array(1:10, dim = c(dataset = 2, member = 5, sdate = 4, ftime = 3)), 
       array(1:4, dim = c(dataset = 2, member = 2, sdate = 5, ftime = 3))),
  paste0("Parameter 'exp' and 'obs' must have same length of ",
                "all dimension expect 'dat_dim'.")
  )

})

##############################################
test_that("2. Output checks: dat1", {

  expect_equal(
    dim(Clim(exp1, obs1)$clim_exp),
    c(dataset = 1, member = 3, ftime = 3, lon = 2, lat = 4)
  )
  expect_equal(
    dim(Clim(exp1, obs1, memb = FALSE)$clim_exp),
    c(dataset = 1, ftime = 3, lon = 2, lat = 4)
  )
  expect_equal(
    dim(Clim(exp1, obs1, time_dim = 'lon')$clim_exp),
    c(dataset = 1, member = 3, sdate = 5, ftime = 3, lat = 4)
  )
  expect_equal(
    dim(Clim(exp1, obs1, method = 'kharin')$clim_exp),
    c(sdate = 5, dataset = 1, member = 3, ftime = 3, lon = 2, lat = 4)
  )
  expect_equal(
    dim(Clim(exp1, obs1, method = 'NDV')$clim_exp),
    c(sdate = 5, dataset = 1, member = 3, ftime = 3, lon = 2, lat = 4)
  )
  expect_equal(
    dim(Clim(exp1, obs1)$clim_obs),
    c(dataset = 1, member = 1, ftime = 3, lon = 2, lat = 4)
  )
  expect_equal(
    dim(Clim(exp1, obs1, method = 'kharin')$clim_obs),
    c(dataset = 1, member = 1, ftime = 3, lon = 2, lat = 4)
  )
  expect_equal(
    dim(Clim(exp1, obs1, method = 'NDV')$clim_obs),
    c(dataset = 1, member = 1, ftime = 3, lon = 2, lat = 4)
  )
  expect_equal(
    (Clim(exp1, obs1)$clim_obs)[1:5],
     c(0.14831161, -0.60462627, 0.06609153, -0.23378059, 0.50553522),
     tolerance = 0.001
  )
  expect_equal(
    (Clim(exp1, obs1, memb = FALSE)$clim_exp)[1:5],
     c(0.10084284, 0.06407350, 0.09028584, 0.17526332, 0.18004387),
     tolerance = 0.001
  )
  expect_equal(
    max(Clim(exp1, obs1)$clim_exp, na.rm = T),
    1.026186,
    tolerance = 0.001
  )
  expect_equal(
    max(Clim(exp1, obs1, method = 'kharin')$clim_exp, na.rm = T),
    2.282634,
    tolerance = 0.001
  )
  expect_equal(
    min(Clim(exp1, obs1, method = 'NDV')$clim_exp, na.rm = T),
    -4.025745,
    tolerance = 0.001
  )

})

##############################################
test_that("3. Output checks: dat2", {

  expect_equal(
    (Clim(exp2, obs2)$clim_obs)[1:5],
     c(0.23142987, -0.60462627, -0.03669491, -0.14193572, 0.61163024),
     tolerance = 0.001
  )
  expect_equal(
    (Clim(exp2, obs2)$clim_exp)[1:5],
     c(0.01054951, -0.04744191, -0.03533071, 0.14149945, 0.02359945),
     tolerance = 0.001
  )
  expect_equal(
    length(which(is.na(Clim(exp2, obs2, na.rm = FALSE)$clim_exp))),
    57
  )
  expect_equal(
    length(which(is.na(Clim(exp2, obs2, na.rm = FALSE)$clim_obs))),
    19
  )
  expect_equal(
    length(which(is.na(Clim(exp2, obs2)$clim_obs))),
    0
  )

})

##############################################

