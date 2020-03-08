context("s2dv::RMSSS tests")

##############################################
 # case 0
 set.seed(1)
 exp0 <- array(rnorm(15), dim = c(sdate = 3, member = 5))
 set.seed(2)
 obs0 <- array(rnorm(6), dim = c(sdate = 3, member = 2))

 # case 1
 set.seed(1)
 exp1 <- array(rnorm(15), dim = c(time = 3, memb = 5))
 set.seed(2)
 obs1 <- array(rnorm(6), dim = c(time = 3, memb = 2))
 
 # case 2
 set.seed(3)
 exp2 <- array(rnorm(120), dim = c(sdate = 10, dat = 1, lon = 3, lat = 2, member = 2))
 set.seed(4)
 obs2 <- array(rnorm(60), dim = c(dat = 1, sdate = 10, member = 1, lat = 2, lon = 3))

##############################################

test_that("1. Input checks", {

  expect_error(
  RMSSS(c(), c()),
  "Parameter 'exp' and 'obs' cannot be NULL."
  )
  expect_error(
  RMSSS('exp', 'obs'),
  "Parameter 'exp' and 'obs' must be a numeric array."
  )
  expect_error(
  RMSSS(c(1:10), c(2:4)),
  paste0("Parameter 'exp' and 'obs' must be at least two dimensions ",
                "containing time_dim and memb_dim.")
  )
  expect_error(
  RMSSS(array(1:10, dim = c(2, 5)), array(1:4, dim = c(2, 2))),
  "Parameter 'exp' and 'obs' must have dimension names."
  )
  expect_error(
  RMSSS(array(1:10, dim = c(a = 2, c = 5)), array(1:4, dim = c(a = 2, b = 2))),
  "Parameter 'exp' and 'obs' must have same dimension name"
  )
  expect_error(
  RMSSS(exp1, obs1, time_dim = 1),
  "Parameter 'time_dim' must be a character string."
  )
  expect_error(
  RMSSS(exp0, obs0, time_dim = 'a'),
  "Parameter 'time_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  RMSSS(exp0, obs0, memb_dim = NA),
  "Parameter 'memb_dim' must be a character string."
  )
  expect_error(
  RMSSS(exp0, obs0, memb_dim = 'memb'),
  "Parameter 'memb_dim' is not found in 'exp' or 'obs' dimension."
  )
  expect_error(
  RMSSS(exp0, obs0, pval = c(T, T)),
  "Parameter 'pval' must be one logical value."
  )
  expect_error(
  RMSSS(exp0, obs0, ncores = 1.4),
  "Parameter 'ncores' must be a positive integer."
  )
  expect_error(
  RMSSS(exp = array(1:10, dim = c(sdate = 1, member = 5, a = 1)),
        obs = array(1:4, dim = c(a = 1, sdate = 2, member = 2))),
  "Parameter 'exp' and 'obs' must have same length of all dimension expect 'memb_dim'."
  )
  expect_error(
  RMSSS(exp = array(1:10, dim = c(sdate = 1, member = 5, a = 1)),
        obs = array(1:4, dim = c(a = 1, sdate = 1, member = 2))),
  "The length of time_dim must be more than 2 to compute RMSSS."
  )
})

##############################################
test_that("1. Output checks: case 1", {

  res1_1 <- RMSSS(exp1, obs1, time_dim = 'time', memb_dim = 'memb')
  expect_equal(
    dim(res1_1$rmsss),
    c(nexp = 5, nobs = 2)
  )
  expect_equal(
    dim(res1_1$p.val),
    c(nexp = 5, nobs = 2)
  )
  expect_equal(
    mean(res1_1$rmsss),
    -0.5449538,
    tolerance = 0.00001
  )

  exp1_2 <- exp1
  exp1_2[2:4] <- NA
  obs1_2 <- obs1
  obs1_2[1:2] <- NA
  res1_2 <- RMSSS(exp1_2, obs1_2, time_dim = 'time', memb_dim = 'memb', pval = TRUE)

  expect_equal(
    length(res1_2$rmsss[which(is.na(res1_2$rmsss))]),
    7
  )
  expect_equal(
    range(res1_2$p.val, na.rm = T),
    c(0.7159769, 0.8167073),
    tolerance = 0.00001
  )

})


##############################################
test_that("2. Output checks: case 2", {

  expect_equal(
    dim(RMSSS(exp2, obs2)$rmsss),
    c(nexp = 2, nobs = 1, dat = 1, lon = 3, lat = 2) 
  )
  expect_equal(
    mean(RMSSS(exp2, obs2)$rmsss),
    -0.3912208,
    tolerance = 0.00001
  )
  expect_equal(
    range(RMSSS(exp2, obs2)$p.val),
    c(0.2627770, 0.9868412),
    tolerance = 0.00001
  )

})

##############################################


