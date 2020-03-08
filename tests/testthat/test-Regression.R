context("s2dv::Regression tests")

##############################################
  # dat1
  set.seed(1)
  datay1 <- array(c(-39:40) + rnorm(80), 
                 dim = c(sdate = 5, ftime = 2, lon = 2, lat = 4))
  set.seed(2)
  datax1 <- array(c(1:80) + rnorm(80), 
                 dim = c(sdate = 5, ftime = 2, lon = 2, lat = 4))
  # dat2
  datay2 <- datay1
  set.seed(1)
  na <- floor(runif(5, min = 1, max = 80))
  datay2[na] <- NA
  datax2 <- datax1
  set.seed(2)
  na <- floor(runif(5, min = 1, max = 80))
  datax2[na] <- NA

  # dat3
  set.seed(1)
  datay3 <- array(c(-39:40) + rnorm(80),
                 dim = c(date = 5, ftime = 2, lon = 2, lat = 4))
  set.seed(2)
  datax3 <- array(c(1:80) + rnorm(80),
                 dim = c(date = 5, lon = 2, lat = 4,ftime = 2))


##############################################
test_that("1. Input checks", {

  expect_error(
  Regression(c(), c()),
  "Parameter 'datay' and 'datax' cannot be NULL."
  )
  expect_error(
  Regression(c('b'), c('a')),
  "Parameter 'datay' and 'datax' must be a numeric array."
  )
  expect_error(
  Regression(c(1:10), c(2:4)),
  "Parameter 'datay' and 'datax' must be at least one dimension 'time_dim'."
  )
  expect_error(
  Regression(array(1:10, dim = c(2, 5)), array(1:4, dim = c(2, 2))),
  "Parameter 'datay' and 'datax' must have dimension names."
  )
  expect_error(
  Regression(array(1:10, dim = c(a = 2, c = 5)), array(1:4, dim = c(a = 2, b = 2))),
  "Parameter 'datay' and 'datax' must have same dimension name"
  )
  expect_error(
  Regression(array(1:10, dim = c(a = 2, b = 5)), array(1:4, dim = c(a = 2, b = 2))),
  "Parameter 'datay' and 'datax' must have same length of all dimensions."
  )
  expect_error(
  Regression(datay1, datax1, time_dim = 1),
  "Parameter 'time_dim' must be a character string."
  )
  expect_error(
  Regression(datay1, datax1, time_dim = 'asd'),
  "Parameter 'time_dim' is not found in 'datay' or 'datax' dimension."
  )
  expect_error(
  Regression(datay1, datax1, na.action = TRUE),
  paste0("Parameter 'na.action' must be a function for NA values or ",
         "a numeric indicating the number of NA values allowed ",
         "before returning NA.")
  )
  expect_error(
  Regression(datay1, datax1, na.action = c(1,2)),
  paste0("Parameter 'na.action' must be a function for NA values or ",
         "a numeric indicating the number of NA values allowed ",
         "before returning NA.")
  )
  expect_error(
  Regression(datay1, datax1, formula =T),
  "Parameter 'formula' must the an object of class 'formula'."
  )
  expect_error(
  Regression(datay1, datax1, pval = 0.05),
  "Parameter 'pval' must be one logical value."
  )
  expect_error(
  Regression(datay1, datax1, conf = 0.05),
  "Parameter 'conf' must be one logical value."
  )
  expect_error(
  Regression(datay1, datax1, conf.lev = 1.5),
  "Parameter 'conf.lev' must be a numeric number between 0 and 1."
  )
  expect_error(
  Regression(datay1, datax1, ncores = T),
  "Parameter 'ncores' must be a positive integer."
  )
})


##############################################
test_that("2. Output checks: dat1", {

  expect_equal(
    dim(Regression(datay1, datax1)$regression),
    c(stats = 2, ftime = 2, lon = 2, lat = 4)
  )
  expect_equal(
    Regression(datay1, datax1)$regression[1:6],
    c(-39.0091480, 0.7290814, -39.1853129, 0.8623175, -37.4342099, 0.7844530), 
    tolerance = 0.001
  )
  expect_equal(
    length(which(is.na(Regression(datay1, datax1)$conf.low))),
    0
  )
  expect_equal(
    max(Regression(datay1, datax1)$conf.upper, na.rm = T),
    127.4267,
    tolerance = 0.001
  )
  expect_equal(
    length(Regression(datay1, datax1, conf = F)),
    3
  )
  expect_equal(
    length(Regression(datay1, datax1, pval = F)),
    4
  )
  expect_equal(
    length(Regression(datay1, datax1, pval = F, conf = F)),
    2
  )
  expect_equal(
    range(Regression(datay1, datax1, conf.lev = 0.99)$conf.low, na.rm = T),
    c(-380.888744, 0.220794),
    tolerance = 0.001
  )
  expect_equal(
    min(Regression(datay1, datax1)$p.val, na.rm = TRUE),
    0.005335,
    tolerance = 0.0001
  )
  expect_equal(
    mean(Regression(datay1, datax1, formula = y~poly(x, 2, raw = T))$p.val, na.rm = TRUE),
    0.3407307,
    tolerance = 0.0001
  )
})

##############################################
test_that("3. Output checks: dat2", {
  expect_equal(
    length(which(is.na(Regression(datay2, datax2, na.action = 0)$conf.lower))),
    14
  )
  expect_equal(
    length(which(is.na(Regression(datay2, datax2, na.action = 2)$conf.lower))),
    0
  )
  expect_equal(
    length(which(is.na(Regression(datay2, datax2, na.action = 1)$p.val))),
    2
  )
  expect_equal(
    length(which(is.na(Regression(datay2, datax2, na.action = na.pass)$p.val))),
    0
  )
  expect_equal(
    which(is.na(Regression(datay2, datax2, na.action = 1)$p.val)),
    c(3,15)
  )
  expect_equal(
    which(is.na(Regression(datay2, datax2, na.action = 1, 
                           formula = y~poly(x, 2, raw = T))$p.val)),
    c(3,15)
  )
})

##############################################
test_that("4. Output checks: dat3", {

  expect_equal(
    dim(Regression(datay3, datax3, time_dim = 'date')$regression),
    c(stats = 2, ftime = 2, lon = 2, lat = 4)
  )
  expect_equal(
    dim(Regression(datay3, datax3, time_dim = 'date')$conf.lower),
    c(stats = 2, ftime = 2, lon = 2, lat = 4)
  )
  expect_equal(
    dim(Regression(datay3, datax3, time_dim = 'date')$p.val),
    c(ftime = 2, lon = 2, lat = 4)
  )


})

##############################################

