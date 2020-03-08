context("s2dv::Trend tests")

##############################################
  # dat1
  dat1 <- array(c(-5, -7, -10:10, 12, 11, 7, 16), 
                dim = c(dat = 1, sdate = 13, ftime = 2))

  # dat2
  set.seed(10)
  dat2 <- c(1:10) + rnorm(10)

  # dat3
  set.seed(1)
  dat3 <- array(c(1:60) + rnorm(60),
                dim = c(sdate = 5, lat = 2, lon = 3, lev = 2))
  set.seed(1)
  na <- floor(runif(5, min = 1, max = 60))
  dat3[na] <- NA
  
##############################################

test_that("1. Input checks", {

  expect_error(
  Trend(c()),
  "Parameter 'data' cannot be NULL."
  )
  expect_error(
    Trend(c(NA, NA)),
    "Parameter 'data' must be a numeric array."
  )
  expect_error(
    Trend(list(a = array(rnorm(50), dim = c(dat = 5, sdate = 10)), b = c(1:4))),
    "Parameter 'data' must be a numeric array."
  )
  expect_error(
    Trend(array(1:10, dim = c(2, 5))),
    "Parameter 'data' must have dimension names."
  )
  expect_error(
    Trend(dat1, time_dim = 'a'),
    "Parameter 'time_dim' is not found in 'data' dimension."
  )
  expect_error(
    Trend(array(c(1:25), dim = c(dat = 1, date = 5, ftime = 5))),
    "Parameter 'time_dim' is not found in 'data' dimension."
  )
  expect_error(
    Trend(dat1, time_dim = 2),
    "Parameter 'time_dim' must be a character string."
  )
  expect_error(
    Trend(dat1, time_dim = c('a','sdate')),
    "Parameter 'time_dim' must be a character string."
  )
  expect_error(
    Trend(dat1, interval = 0),
    "Parameter 'interval' must be a positive number."
  )
  expect_error(
    Trend(dat1, conf = 3),
    "Parameter 'conf' must be one logical value."
  )
  expect_error(
    Trend(dat1, polydeg = 3.5),
    "Parameter 'polydeg' must be a positive integer."
  )
  expect_error(
    Trend(dat1, ncore = 3.5),
    "Parameter 'ncores' must be a positive integer."
  )
})

##############################################
test_that("2. Output checks: dat1", {

  expect_equal(
    Trend(dat1)$trend,
    array(c(-9.7692308, 0.6593407, 0.9615385, 0.7967033), 
          dim = c(stats = 2, dat = 1, ftime = 2)),
    tolerance = 0.0001
  )
  expect_equal(
    Trend(dat1)$conf.upper,
    array(c(-7.4735367, 0.9485709, 3.0167860, 1.0556402), 
          dim = c(stats = 2, dat = 1, ftime = 2)),
    tolerance = 0.0001
  )
  expect_equal(
    median(Trend(dat1)$detrended, na.rm = TRUE),
    0.1153846,
    tolerance = 0.001
  )

})

##############################################
test_that("3. Output checks: dat2", {

  expect_equal(
    Trend(dat2),
    list(trend = array(c(-0.182, 0.944), dim = c(stats = 2)), 
         conf.lower = array(c(-1.316, 0.761), dim = c(stats = 2)),
         conf.upper = array(c(0.953, 1.127), dim = c(stats = 2)),
         detrended = array(c(0.257, 0.110, -1.021, -0.193, 0.757, 0.909, 
                             -0.633, 0.267, -0.939, 0.487), dim = c(sdate = 10))),
    tolerance = 0.001
  )
  expect_equal(
    Trend(dat2, interval = 2),
    list(trend = array(c(-0.182, 0.472), dim = c(stats = 2)),
         conf.lower = array(c(-1.316, 0.381), dim = c(stats = 2)),
         conf.upper = array(c(0.953, 0.563), dim = c(stats = 2)),
         detrended = array(c(0.257, 0.110, -1.021, -0.193, 0.757, 0.909,
                             -0.633, 0.267, -0.939, 0.487), dim = c(sdate = 10))),
    tolerance = 0.001
  )
  expect_equal(
    length(Trend(dat2, conf = F)),
    2
  )
  expect_equal(
    names(Trend(dat2, conf = F)),
    c('trend', 'detrended')
  )

})

##############################################
test_that("4. Output checks: dat3", {

  expect_equal(
    median(Trend(dat3)$trend, na.rm =  TRUE),
    1.368071,
    tolerance = 0.0001
  )
  expect_equal(
    dim(Trend(dat3, polydeg = 2)$trend),
    c(stats = 3, lat = 2, lon = 3, lev = 2)
  )

})

